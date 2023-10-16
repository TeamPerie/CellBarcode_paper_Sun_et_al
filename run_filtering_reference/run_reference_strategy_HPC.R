library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(CellBarcode)
library(ggrepel)
library(ROCR)
devtools::load_all("../../lib/CellBarcodeSim/")

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 
theme1 = theme0 + theme(
    axis.text.x = element_text( 
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 12, hjust = 0.5)
)

knitr::opts_chunk$set(fig.width=20, fig.height=12) 

#' # Sample info
# NOTE: input
sample_info = fread("../run_simulation_no_umi/non_umi_simu_design_matrix.tsv")
sample_info$simu_id %<>% as.character
setkey(sample_info, "simu_id")


#' #  Get True barcode
# NOTE: input
true_barcode_l = read_rds("../run_preprocess_simulation/tmp/non_umi_simulation_ref.rds")

d_true_barcode = true_barcode_l %>% rbindlist(idcol = "simu_file") 
x = d_true_barcode$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
d_true_barcode = cbind(d_true_barcode, simu_id = x)
d_true_barcode$simu_name = sample_info[d_true_barcode$simu_id, simu_name]
d_true_barcode = rename(d_true_barcode, c("seq" = "barcode_seq"))

#' # Get Simulated barcodes
# NOTE: input
bc_obj = read_rds("../run_preprocess_simulation/tmp/non_umi_simulation_bcobj.rds")
bc_names(bc_obj) %>% length

#' # Do the reference correction

random_barcode_ref = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")
hamming_barcode_ref = system.file("data", "hamming_barcodes.tsv", package = "CellBarcodeSim")
vdj_barcode_ref = system.file("data", "vdj_barcodes.tsv", package = "CellBarcodeSim")

## reference library of random barcode
x_random_ref = fread(random_barcode_ref)$seq %>% unique
## reference library of VDJ barcode
x_vdj_ref = fread(vdj_barcode_ref)[order(freq, decreasing = T)]$seq %>% unique
## reference hamming barcode
x_hamming_ref = fread(hamming_barcode_ref)$seq %>% unique
## get the clean bc object
bc_obj = bc_cure_depth(bc_obj, isUpdate = F)
# bc_obj = bc_cure_cluster(bc_obj, count_threshold = 1e7)
for (i in seq_along(true_barcode_l)) {
    #     i = 1
    simu_id_i = bc_names(bc_obj)[i] %>% str_split_fixed("_", 6) %>% extract(6)
    barcode_type_i = sample_info[simu_id == simu_id_i, barcode_type]
    if (barcode_type_i == "random") {
        x_ref_i = x_random_ref
    } else if (barcode_type_i == "hamming") {
        x_ref_i = x_hamming_ref
    } else if (barcode_type_i == "vdj") {
        x_ref_i = x_vdj_ref
    } else {
        stop("error the barcode type is not exist")
    }
    seq_trunk_i = sample_info[simu_id == simu_id_i, barcode_length]
    x_ref_i = substring(x_ref_i, 1, seq_trunk_i)
    x = bc_obj@cleanBc[[i]]
    x = x[barcode_seq %in% x_ref_i]
    bc_obj@cleanBc[[i]] = x
}

#' # Precision Recall: 
#' 
#' - Precision: TP / (TP + FP)
#' - Recall: TP / (TP + FN)
#' 
#' ROC: TPR = TP / (TP + FN), FPR = FP / (FP + TN)
#' 
#' Because there are a lot of TN, the ROC will be over optimized.

## simu name - auc
d_auc = data.table(simu_file = c(), auc = c(), batch = c())
## reads count - true barcode
d_count_true_barcode = data.table(simu_file = c(), count = c(), is_true = c())

for (i in bc_names(bc_obj)) {
    simu_file = i
    x_true = d_true_barcode[simu_file == i]
    x_true = x_true[, .(freq = sum(freq)), by = .(simu_file, barcode_seq, simu_id, simu_name)]
    x_seq = bc_2df(bc_obj[, i])

    x_merge = merge(x_true, x_seq, all=T)
    ## true barcode: 0, error: 1
    x_merge[, x_label := !is.na(freq)]
    x_merge[is.na(count), count := 0]
    x_label = x_merge$x_label %>% as.integer
    x_clone_size = x_merge$freq
    x_count = x_merge$count

    d_count_true_barcode = rbind(d_count_true_barcode, 
        data.table(simu_file = simu_file, count = x_count, is_true = x_label, clone_size = x_clone_size)
    )

    n_missing = sum(x_count == 0)
    n_ref = sum(x_label == 1)
    p_missing = n_missing / n_ref
    x_label = x_label[x_count != 0]
    x_count = x_count[x_count != 0]

    if (length(unique(x_label)) == 1) { 
        x_count = c(x_count, 0)
        x_label = c(x_label, 0)
    }

    pred = prediction(log10(x_count + 1), x_label)

    ## AUC
    y_auc = performance(pred, "auc")@y.values[[1]] %>% round(3)
    y_aucpr = performance(pred, "aucpr")@y.values[[1]] %>% round(3)
    y_aucpr = y_aucpr * (1 - p_missing)

    d_auc = rbind(d_auc, data.table(simu_file = simu_file, auc = y_auc, aucpr = y_aucpr, batch = 1))
}

## annotate d_count_true_barcode
x = d_count_true_barcode$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2)) 
d_count_true_barcode = cbind(d_count_true_barcode, simu_id = x)
d_count_true_barcode$simu_name = sample_info[d_count_true_barcode$simu_id, simu_name]

# NOTE: output
write_tsv(d_count_true_barcode, "./tmp/table_count_true_barcode_reference.tsv")

## annotate d_auc
x = d_auc$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
d_auc = cbind(d_auc, simu_id = x)
d_auc$simu_name = sample_info[d_auc$simu_id, simu_name]

# NOTE: output
write_tsv(d_auc, "./tmp/table_auc_no_umi_reference.tsv")

#' # Apply the autothreshold strategy
x_cutoff = bc_auto_cutoff(bc_obj)
bc_obj2 = bc_obj

## positive
x_highlight = true_barcode_l %>% lapply(function(x) { x$seq })
names(x_highlight)

## fix the order
i_name = names(x_cutoff)
bc_obj2 = bc_obj2[, i_name]
x_highlight = x_highlight[i_name]


#+ eval=T
d = data.table(simu_file = c(), tpr = c(), fpr = c())

for (i in seq_along(x_cutoff)) {
    simu_file = bc_names(bc_obj2)[i]
    ## Caculate the TPR and FPR
    x_detected_barcode_seq = bc_2dt(bc_obj2[, i])[count > x_cutoff[i], barcode_seq]

    x_fp = sum(!(x_detected_barcode_seq %in% x_highlight[[i]]))
    x_tn = sum(!(bc_barcodes(bc_obj2[, 1]) %in% x_highlight[[i]]))
    x_tp = sum(x_detected_barcode_seq %in% x_highlight[[i]])
    x_fn = sum(!(x_highlight[[i]] %in% x_detected_barcode_seq))

    ## TPR: TP / (TP + FN), TPR or Recall
    x_tpr = x_tp / (x_tp + x_fn)
    ## FPR: FP / (TP + TN)
    x_fpr = x_fp / (x_fp + x_tn)
    ## Precision
    x_pre = x_tp / (x_tp + x_fp)
    ## Recal
    x_rec = x_tp / (x_tp + x_fn)

    d = rbind(d, data.table(simu_file = simu_file, tpr = x_tpr, fpr = x_fpr, pre = x_pre, rec = x_rec))
}

# NOTE: output
x = d$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
d = cbind(d, simu_id = x)
d$simu_name = sample_info[d$simu_id, simu_name]

write_tsv(d, "./tmp/table_no_umi_auto_threshold_predict_rate_reference.tsv")

#' # Apply 0.0001 percentage threshold
x_cutoff = lapply(bc_obj@cleanBc, function(x) { x$count %>% sum }) %>% unlist * 0.0001
names(x_cutoff) = rownames(bc_meta(bc_obj))
bc_obj2 = bc_cure_depth(bc_obj, depth = 0, isUpdate=T)

## positive
x_highlight = true_barcode_l %>% lapply(function(x) { x$seq })
names(x_highlight)

## fix the order
i_name = names(x_cutoff)
bc_obj2 = bc_obj2[, i_name]
x_highlight = x_highlight[i_name]

x_highlight

d = data.table(simu_file = c(), tpr = c(), fpr = c())

for (i in seq_along(x_cutoff)) {
    simu_file = bc_names(bc_obj2)[i]
    ## Caculate the TPR and FPR
    x_detected_barcode_seq = bc_2dt(bc_obj2[, i])[count > x_cutoff[i], barcode_seq]

    x_fp = sum(!(x_detected_barcode_seq %in% x_highlight[[i]]))
    x_tn = sum(!(bc_barcodes(bc_obj2[, 1]) %in% x_highlight[[i]]))
    x_tp = sum(x_detected_barcode_seq %in% x_highlight[[i]])
    x_fn = sum(!(x_highlight[[i]] %in% x_detected_barcode_seq))

    ## TPR: TP / (TP + FN), TPR or Recall
    x_tpr = x_tp / (x_tp + x_fn)
    ## FPR: FP / (TP + TN)
    x_fpr = x_fp / (x_fp + x_tn)
    ## Precision
    x_pre = x_tp / (x_tp + x_fp)
    ## Recal
    x_rec = x_tp / (x_tp + x_fn)

    d = rbind(d, data.table(simu_file = simu_file, tpr = x_tpr, fpr = x_fpr, pre = x_pre, rec = x_rec))
}

# NOTE: output
x = d$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
d = cbind(d, simu_id = x)
d$simu_name = sample_info[d$simu_id, simu_name]

write_tsv(d, "./tmp/table_no_umi_p0001_threshold_predict_rate_reference.tsv")


#' # Apply 0.00001 percentage threshold
x_cutoff = lapply(bc_obj@cleanBc, function(x) { x$count %>% sum }) %>% unlist * 0.00001
names(x_cutoff) = rownames(bc_meta(bc_obj))
bc_obj2 = bc_cure_depth(bc_obj, depth = 0, isUpdate=T)

## positive
x_highlight = true_barcode_l %>% lapply(function(x) { x$seq })
names(x_highlight)

## fix the order
i_name = names(x_cutoff)
bc_obj2 = bc_obj2[, i_name]
x_highlight = x_highlight[i_name]


#+ eval=T
d = data.table(simu_file = c(), tpr = c(), fpr = c())

for (i in seq_along(x_cutoff)) {
    simu_file = bc_names(bc_obj2)[i]
    ## Caculate the TPR and FPR
    x_detected_barcode_seq = bc_2dt(bc_obj2[, i])[count > x_cutoff[i], barcode_seq]

    x_fp = sum(!(x_detected_barcode_seq %in% x_highlight[[i]]))
    x_tn = sum(!(bc_barcodes(bc_obj2[, 1]) %in% x_highlight[[i]]))
    x_tp = sum(x_detected_barcode_seq %in% x_highlight[[i]])
    x_fn = sum(!(x_highlight[[i]] %in% x_detected_barcode_seq))

    ## TPR: TP / (TP + FN), TPR or Recall
    x_tpr = x_tp / (x_tp + x_fn)
    ## FPR: FP / (TP + TN)
    x_fpr = x_fp / (x_fp + x_tn)
    ## Precision
    x_pre = x_tp / (x_tp + x_fp)
    ## Recal
    x_rec = x_tp / (x_tp + x_fn)

    d = rbind(d, data.table(simu_file = simu_file, tpr = x_tpr, fpr = x_fpr, pre = x_pre, rec = x_rec))
}

# NOTE: output
x = d$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
d = cbind(d, simu_id = x)
d$simu_name = sample_info[d$simu_id, simu_name]

write_tsv(d, "./tmp/table_no_umi_p00001_threshold_predict_rate_reference.tsv")

sessionInfo()
