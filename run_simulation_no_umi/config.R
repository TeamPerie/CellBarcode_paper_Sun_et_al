###########################################################################
#                            Workflow library                             #
###########################################################################

library(whisker)
library(magrittr)
library(stringr)
library(clustermq)
library(future)
library(promises)

plan(multicore)

## Compile template
## INPUT:
##   template_file: file with mustache {{}} template syntax
##   config: list, with 
compile = function(template_file, config) {
    template_text = readLines(template_file) %>% str_c(collapse = "\n")
    
    n_max = lapply(config, length) %>% unlist %>% max
    config = lapply(config, function(x) {rep_len(x, n_max)})

    output_file_path_v = c()

    for (i in seq_len(n_max)) {
        config_sub = list()
        config_sub = lapply(config, function(x) { x[i] })

        compiled = whisker.render(template_text, config_sub)

        output_file_name = str_c(
            config_sub$job_name,
            "_",
            basename(template_file)
            )

        output_file_path = str_c("./run/", output_file_name)

        write(compiled, file=output_file_path)

        output_file_path_v = append(output_file_path_v, output_file_path)
    }
    names(output_file_path_v) = config$job_name
    lapply(config$job_name, function(i) {
        function(exec="bash", intern=TRUE) { system(str_glue("{exec} {output_file_path_v[i]}"), intern=intern) }
    })
}

runLocal = function(compiled_script) {
    exec = "bash"
    promise_map(seq_along(compiled_script), function(i) {
        future_promise(compiled_script[[i]](exec))
    }) 
}

get_job_status = function(job_id) {                                                                                         
    job = system(str_glue("qstat -f1 {job_id}"), intern = T) %>% str_trim("both") %>% str_split_fixed(" = ", 2) %>% data.table
    job[V1 == "job_state", V2]                                                                                                
}                                                                                                                           


runPBS = function(compiled_script) {
    future_promise({
        exec = 'qsub'
        job_id = lapply(seq_along(compiled_script), function(i) {
            compiled_script[[i]](exec, intern=T) 
        })
        repeat {
            Sys.sleep(30)
            run_status = grep("C|E", vapply(job_id, get_job_status, ""))
            if (length(run_status) == length(job_id)) {
                break
            }
        }
        job_id
    })
}


dir.create("./run/", showWarnings=F, recursive=T)
dir.create("./tmp/", showWarnings=F, recursive=T)


########
#  run #
########
# TODO: Extract the reads length from fastqc result and put it in to configuration
template_file = "template_simulation_no_umi.sh"

config_run = list(
   simu_number = c(1:30)
    )
config_run$job_name = paste0(config_run$barcode_type, config_run$simu_number)

runner_run = compile(template_file, config_run)
#runPBS(runner_run)


