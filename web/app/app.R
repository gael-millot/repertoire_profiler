#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     app.R                                                           ##
##                                                                     ##
##     Web page of repertoire_profiler                                 ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################


################################ Aim

# Simple web form generated from the nextflow.config file of repertoire_profiler:
# one input per parameter of the "Parameters that must be set by the user" section, with the default value and the comment of nextflow.config as help.
# A submission writes the files in <RP_WORK_DIR>/jobs/<key>/ and the task in <RP_WORK_DIR>/todo/<key>.json.
# The web/daemon/rp_bioblend.py daemon sends the task to the repertoire_profiler Galaxy tool (galaxy.pasteur.fr -> Maestro cluster).
# The key, given to the user, is the only way to get the results back ("Check my job" tab, or the link in the e-mail).

################################ End Aim


################################ Settings (environment variables)

work_dir <- Sys.getenv("RP_WORK_DIR", "/srv/rp_work") # folder shared with the daemon
repo_dir <- Sys.getenv("RP_REPO_DIR", normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)) # repertoire_profiler clone (nextflow.config and galaxy/)
python <- Sys.getenv("RP_PYTHON", "python3")
max_upload_mb <- as.numeric(Sys.getenv("RP_MAX_UPLOAD_MB", "50"))
retention_days <- Sys.getenv("RP_RETENTION_DAYS", "15")
options(shiny.maxRequestSize = max_upload_mb * 1024^2)

config_file <- file.path(repo_dir, "nextflow.config")
fill_config <- file.path(repo_dir, "galaxy", "fill_config.py")
locked_file <- file.path(repo_dir, "galaxy", "locked_params.txt")
key_regex <- "^[A-Za-z0-9_-]{16,64}$"
mail_regex <- "^[^@[:space:]]+@[^@[:space:]]+\\.[^@[:space:]]+$"
fasta_regex <- "(\\.(fasta|fa|fas|fna|fsa|seq)|/)$" # files of the zip (or folders)

for(i0 in c("todo", "doing", "done", "error", "jobs")){
    dir.create(file.path(work_dir, i0), recursive = TRUE, showWarnings = FALSE)
}

################################ End Settings


################################ Parameters of nextflow.config

# parsed by fill_config.py, the same script that fills the config in Galaxy, so that the page and the tool always agree
spec <- jsonlite::fromJSON(paste(system2(python, c(shQuote(fill_config), "--template", shQuote(config_file), "--spec"), stdout = TRUE), collapse = "\n"), simplifyDataFrame = FALSE)
locked <- trimws(sub("#.*", "", readLines(locked_file)))
locked <- locked[locked != ""]
spec <- spec[ ! sapply(spec, function(x){x$name}) %in% locked]
sections <- unique(sapply(spec, function(x){x$section}))

param_input <- function(p){
    id <- paste0("p_", p$name)
    if(p$default %in% c("TRUE", "FALSE")){
        widget <- shiny::selectInput(id, p$name, choices = c("TRUE", "FALSE"), selected = p$default)
    }else{
        widget <- shiny::textInput(id, p$name, value = p$default, width = "100%")
    }
    shiny::div(
        class = "rp-param",
        widget,
        shiny::tags$details(shiny::tags$summary("help"), shiny::tags$pre(class = "rp-help", p$help))
    )
}

section_panel <- function(title){
    params <- Filter(function(x){x$section == title}, spec)
    shiny::tags$details(
        open = if(grepl("^Optional", title)){NULL}else{"open"}, # optional sections closed
        class = "rp-section",
        shiny::tags$summary(shiny::h4(title, style = "display: inline")),
        lapply(params, param_input)
    )
}

################################ End Parameters of nextflow.config


################################ Functions

new_key <- function(){
    # cryptographic random key: it is the only protection of the results
    bytes <- tryCatch({
        con <- file("/dev/urandom", open = "rb", raw = TRUE) # raw = TRUE: device file, not a regular file
        on.exit(close(con))
        readBin(con, "raw", 18)
    }, error = function(e){NULL})
    if(is.null(bytes) || length(bytes) != 18){
        stop("cannot read /dev/urandom to create the job key")
    }
    paste(as.character(bytes), collapse = "")
}

find_task <- function(key){
    # return list(folder, task) or NULL
    if( ! grepl(key_regex, key)){
        return(NULL)
    }
    for(i0 in c("todo", "doing", "done", "error")){
        path <- file.path(work_dir, i0, paste0(key, ".json"))
        if(file.exists(path)){
            task <- tryCatch(jsonlite::fromJSON(path), error = function(e){list(state = "unknown", message = "Status being updated, please retry.")})
            return(list(key = key, folder = i0, task = task))
        }
    }
    NULL
}

################################ End Functions


################################ UI

ui <- shiny::navbarPage(
    title = "repertoire_profiler",
    id = "tabs",
    header = shiny::tags$style(shiny::HTML("
        .rp-section {border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; margin-bottom: 10px;}
        .rp-param {margin-bottom: 6px;}
        .rp-help {white-space: pre-wrap; font-size: 12px; max-height: 250px; overflow-y: auto;}
        details > summary {cursor: pointer;}
    ")),
    shiny::tabPanel(
        "Submit",
        value = "submit",
        shiny::fluidRow(shiny::column(
            width = 10, offset = 1,
            shiny::p("Run the ", shiny::a("repertoire_profiler", href = "https://github.com/gael-millot/repertoire_profiler", target = "_blank"),
                " pipeline on the Institut Pasteur cluster. Fill the form below, which corresponds to the nextflow.config file of the pipeline, and submit.",
                " You will get a key to follow the job and download the results, which are kept ", retention_days, " days."),
            shiny::tags$details(
                open = "open",
                class = "rp-section",
                shiny::tags$summary(shiny::h4("Input files", style = "display: inline")),
                shiny::fileInput("sample_zip", paste0("Zip of the fasta files (sample_path parameter, max ", max_upload_mb, " MB)"), accept = ".zip"),
                shiny::fileInput("meta_tsv", "Metadata .tsv file (meta_path parameter, optional)", accept = c(".tsv", ".txt")),
                shiny::textInput("mail", "E-mail (optional, to be warned when the job ends)")
            ),
            lapply(sections, section_panel),
            shiny::actionButton("submit", "Submit", class = "btn-primary btn-lg"),
            shiny::br(), shiny::br()
        ))
    ),
    shiny::tabPanel(
        "Check my job",
        value = "check",
        shiny::fluidRow(shiny::column(
            width = 8, offset = 2,
            shiny::textInput("key", "Job key", width = "100%"),
            shiny::actionButton("check", "Check"),
            shiny::hr(),
            shiny::uiOutput("status")
        ))
    )
)

################################ End UI


################################ Server

server <- function(input, output, session){

    # key given in the URL (?key=...), e.g. from the e-mail
    shiny::observe({
        query <- shiny::parseQueryString(session$clientData$url_search)
        if( ! is.null(query$key)){
            shiny::updateTextInput(session, "key", value = query$key)
            shiny::updateNavbarPage(session, "tabs", selected = "check")
        }
    })

    shiny::observeEvent(input$submit, {
        fail <- function(...){
            shiny::showModal(shiny::modalDialog(title = "Submission refused", paste0(...), easyClose = TRUE))
        }
        # files
        if(is.null(input$sample_zip)){
            return(fail("Please give the zip of the fasta files."))
        }
        zip_content <- tryCatch(utils::unzip(input$sample_zip$datapath, list = TRUE)$Name, error = function(e){NULL})
        if(length(zip_content) == 0){
            return(fail("The fasta file is not a valid zip file."))
        }
        not_fasta <- zip_content[ ! grepl(fasta_regex, zip_content, ignore.case = TRUE)]
        if(length(not_fasta) > 0){
            return(fail("The zip must contain only fasta files. Not fasta: ", paste(utils::head(not_fasta, 5), collapse = ", ")))
        }
        mail <- trimws(input$mail)
        if(mail != "" && ! grepl(mail_regex, mail)){
            return(fail("Invalid e-mail address."))
        }
        # parameters: only those changed from the nextflow.config default
        values <- list()
        for(p in spec){
            value <- input[[paste0("p_", p$name)]]
            if(is.null(value)){ # input not created yet
                next
            }
            value <- trimws(value)
            if( ! identical(value, p$default)){
                values[[p$name]] <- value
            }
        }
        key <- new_key()
        job_dir <- file.path(work_dir, "jobs", key)
        dir.create(job_dir, recursive = TRUE)
        params_file <- file.path(job_dir, "params.json")
        writeLines(jsonlite::toJSON(if(length(values) == 0){structure(list(), names = character(0))}else{values}, auto_unbox = TRUE, pretty = TRUE), params_file)
        # same check as in the Galaxy tool, so that errors are given now and not after the run
        check <- suppressWarnings(system2(python, c(shQuote(fill_config), "--template", shQuote(config_file), "--locked", shQuote(locked_file), "--values", shQuote(params_file), "-o", if(.Platform$OS.type == "windows"){"NUL"}else{"/dev/null"}), stdout = TRUE, stderr = TRUE))
        if( ! is.null(attr(check, "status"))){
            unlink(job_dir, recursive = TRUE)
            return(fail(paste(check, collapse = "\n")))
        }
        file.copy(input$sample_zip$datapath, file.path(job_dir, "sample.zip"))
        if( ! is.null(input$meta_tsv)){
            file.copy(input$meta_tsv$datapath, file.path(job_dir, "metadata.tsv"))
        }
        # the task file is written last: the daemon only takes complete jobs
        task <- list(key = key, mail = mail, submitted = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), state = "queued", message = "The job is waiting to be sent to Galaxy.")
        tmp <- file.path(work_dir, "todo", paste0(key, ".json.tmp"))
        writeLines(jsonlite::toJSON(task, auto_unbox = TRUE, pretty = TRUE), tmp)
        file.rename(tmp, file.path(work_dir, "todo", paste0(key, ".json")))
        shiny::updateTextInput(session, "key", value = key)
        shiny::showModal(shiny::modalDialog(
            title = "Job submitted",
            shiny::p("Keep this key to follow your job and download the results:"),
            shiny::tags$pre(key),
            shiny::p("You can close this page and come back later with the key (Check my job tab)."),
            easyClose = TRUE
        ))
    })

    current <- shiny::eventReactive(input$check, {
        find_task(trimws(input$key))
    })

    output$status <- shiny::renderUI({
        res <- current()
        if(is.null(res)){
            return(shiny::p("Unknown key."))
        }
        items <- list(
            shiny::h4("State: ", res$task$state),
            shiny::p(res$task$message),
            shiny::p("Submitted: ", res$task$submitted)
        )
        if(res$folder == "done"){
            items <- c(items, list(shiny::downloadButton("dl_results", "Download the results (.zip)"), " ", shiny::downloadButton("dl_report", "Download the report (.html)")))
        }else if(res$folder == "error" && file.exists(file.path(work_dir, "jobs", res$key, "error.txt"))){
            items <- c(items, list(shiny::downloadButton("dl_error", "Download the error log")))
        }else if(res$folder != "error"){
            items <- c(items, list(shiny::p(shiny::em("Click Check again to refresh."))))
        }
        do.call(shiny::tagList, items)
    })

    job_file <- function(name){
        shiny::downloadHandler(
            filename = function(){paste0("repertoire_profiler_", name)},
            content = function(file){
                res <- current()
                shiny::req(res)
                file.copy(file.path(work_dir, "jobs", res$key, name), file)
            }
        )
    }
    output$dl_results <- job_file("results.zip")
    output$dl_report <- job_file("report.html")
    output$dl_error <- job_file("error.txt")
}

################################ End Server

shiny::shinyApp(ui = ui, server = server)
