suppressPackageStartupMessages(library(shiny))
options(shiny.maxRequestSize = 12 * 1024^3, shiny.fullstacktrace = FALSE)
APP_ROOT <- if (basename(getwd()) == "app") normalizePath("..", mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
setwd(APP_ROOT)
source(file.path(APP_ROOT, "R", "engine.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x
split_lines <- function(x) { z <- trimws(unlist(strsplit(x %||% "", "\\r?\\n"))); unique(z[nzchar(z)]) }
parse_tuples <- function(x, n, label) {
  lines <- split_lines(x)
  if (!length(lines)) stop("Enter at least one ", label, " tuple.")
  rows <- lapply(lines, function(line) trimws(strsplit(line, ",", fixed = TRUE)[[1]]))
  if (any(lengths(rows) != n)) stop("Each ", label, " row must contain exactly ", n, " comma-separated labels.")
  do.call(rbind, rows)
}
result_table <- function(x) {
  if (is.data.frame(x)) {
    keep <- !vapply(x, function(z) is.list(z) || is.data.frame(z), logical(1))
    if (!any(keep)) return(data.frame(note="Nested result details are preserved in the full RDS object."))
    return(x[, keep, drop=FALSE])
  }
  if (is.matrix(x)) return(as.data.frame(x))
  if (is.list(x)) {
    candidates <- x[vapply(x, function(z) is.data.frame(z) || is.matrix(z), logical(1))]
    if (length(candidates)) return(as.data.frame(candidates[[1]]))
    scalar <- x[vapply(x, function(z) is.atomic(z) && length(z) == 1L, logical(1))]
    if (length(scalar)) return(as.data.frame(scalar, stringsAsFactors = FALSE))
  }
  tryCatch(as.data.frame(x), error = function(e) data.frame(note = "Result saved as RDS; see run archive for the complete object."))
}

ui <- bslib::page_navbar(
  title = tags$span(tags$strong("ANCESTRYLAB"), tags$span("  /  RESEARCH WORKSPACE", class = "brand-sub")),
  id = "main_nav", selected = "Overview", fillable = TRUE,
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#245b8f", secondary = "#60758b"),
  header = tagList(tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML(paste0(
      ":root{--al-navy:#10243a;--al-blue:#245b8f;--al-bg:#f3f6fa;--al-line:#e2e8f0}",
      "body{background:var(--al-bg);color:#172638;font-family:Inter,Arial,sans-serif}",
      ".navbar{box-shadow:0 2px 14px #10243a12}.navbar-brand{font-size:17px;letter-spacing:.07em}",
      ".brand-sub{font-size:10px;color:#91a8bd;letter-spacing:.13em}.page-title{font-size:25px;font-weight:700;color:var(--al-navy);margin:10px 0 4px}",
      ".page-subtitle{color:#64748b;margin-bottom:22px}.al-card{background:white;border:1px solid var(--al-line);border-radius:14px;padding:20px;margin-bottom:18px;box-shadow:0 3px 14px #10243a08}",
      ".al-card h4{font-size:15px;font-weight:700;color:var(--al-navy);margin:0 0 14px}.metric{background:#fff;border:1px solid var(--al-line);border-radius:13px;padding:17px 19px;height:100%}",
      ".metric .k{font-size:11px;color:#718096;text-transform:uppercase;letter-spacing:.1em}.metric .v{font-size:25px;font-weight:700;color:var(--al-navy);margin-top:5px}.metric .d{font-size:11px;color:#8190a2;margin-top:4px}",
      ".form-label{font-weight:600;color:#33465b;font-size:13px}.btn-primary{font-weight:600}.btn{border-radius:8px}.form-control,.form-select{border-color:#d8e1eb;border-radius:8px}",
      ".nav-link{font-weight:600}.status-strip{border-left:4px solid var(--al-blue);background:#eaf2f9;padding:12px 15px;border-radius:7px;color:#27425d;font-size:13px}",
      ".small-muted{font-size:12px;color:#738196}.method-tag{display:inline-block;background:#edf3f8;color:#365a7c;border-radius:20px;padding:4px 9px;font-size:11px;font-weight:600;margin:2px}",
      ".dataTables_wrapper{font-size:13px}.tab-content{padding-top:22px}.form-text{font-size:11px;color:#7b8999}.alert{border-radius:10px}"
    )))),
    tags$div(class="container-fluid px-4 pt-2", uiOutput("global_status"))
  ),
  bslib::nav_panel("Overview", value = "Overview",
    div(class = "container-fluid px-4 py-3", div(class = "page-title", "Ancient DNA analysis workspace"), div(class = "page-subtitle", "A reproducible interface for population-genetic modelling, diagnostics and provenance."),
      fluidRow(column(3, div(class="metric",div(class="k","Active dataset"),div(class="v",textOutput("metric_dataset",inline=TRUE)),div(class="d","Genotype panel currently selected"))), column(3, div(class="metric",div(class="k","Populations"),div(class="v",textOutput("metric_pops",inline=TRUE)),div(class="d","Exact labels read from input metadata"))), column(3, div(class="metric",div(class="k","Analysis runs"),div(class="v",textOutput("metric_runs",inline=TRUE)),div(class="d","Persisted in the project results directory"))), column(3, div(class="metric",div(class="k","Engine"),div(class="v",textOutput("metric_engine",inline=TRUE)),div(class="d","ADMIXTOOLS 2 runtime status")))),
      fluidRow(column(8, div(class="al-card",h4("Workspace status"),uiOutput("workspace_status"),hr(),h4("Analysis modules"),div(span(class="method-tag","f2 / f3 / f4 / D / f4-ratio"),span(class="method-tag","FST"),span(class="method-tag","qpWave"),span(class="method-tag","qpAdm + rotation"),span(class="method-tag","qpGraph"),span(class="method-tag","FindGraph"),span(class="method-tag","GhostBuster")),p(class="small-muted mt-3","Every workflow uses explicit sample/population selections. Missing labels are reported, never silently replaced with a proxy."))), column(4,div(class="al-card",h4("Recent runs"),DT::DTOutput("recent_runs"),p(class="small-muted mt-2","Full objects, CSV tables and JSON manifests are saved per run."))))
    )
  ),
  bslib::nav_panel("Project & data", value="Project & data", div(class="container-fluid px-4 py-3",div(class="page-title","Project and dataset control"),div(class="page-subtitle","Select the exact genotype dataset, inspect its population labels and define analysis-wide settings."),
    fluidRow(column(7,div(class="al-card",h4("Genotype dataset"),textInput("prefix","Genotype file prefix",value=file.path(APP_ROOT,"data","reference"),width="100%"),p(class="small-muted","Accepted: EIGENSTRAT .geno/.snp/.ind or PLINK .bed/.bim/.fam triplet."),fileInput("eigenfiles","Upload one matching genotype triplet",multiple=TRUE,accept=c(".geno",".snp",".ind",".bed",".bim",".fam")),actionButton("inspect_data","Inspect dataset",class="btn-primary"),hr(),h4("AADR public data retrieval"),selectInput("panel","AADR panel",choices=c("1240K","HO","2M"),selected="1240K"),checkboxInput("metadata_only","Metadata only (recommended for discovery)",TRUE),actionButton("download_aadr","Download selected AADR files"),verbatimTextOutput("dataset_status"))),column(5,div(class="al-card",h4("Dataset inspection"),uiOutput("dataset_summary"),DT::DTOutput("population_table"))))),
    fluidRow(column(6,div(class="al-card",h4("Shared f2 configuration"),numericInput("block","Block size (Morgans)",.05,min=.001,step=.005),numericInput("maxmiss","Maximum missing fraction",0,min=0,max=1,step=.01),numericInput("minmaf","Minimum minor-allele frequency",0,min=0,max=.5,step=.001),numericInput("maxmaf","Maximum minor-allele frequency",.5,min=0,max=.5,step=.001),numericInput("maxmem","Maximum f2 memory (MB)",8000,min=100,step=256),checkboxInput("auto_only","Autosomes only (1–22)",TRUE),checkboxInput("transitions","Include transitions",TRUE),checkboxInput("transversions","Include transversions",TRUE),checkboxInput("apply_corr","Apply small-sample f2 correction",TRUE),selectizeInput("populations","Populations to include in f2 cache",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button"))),checkboxInput("afprod","Use allele-frequency-product mode (sparse-data sensitivity; justify in methods)",FALSE),checkboxInput("f2_adjust_pseudo","Correct for pseudohaploid genotype coding",TRUE),checkboxInput("overwrite_f2","Recompute even if a matching cache exists",FALSE),actionButton("f2","Compute / reuse f2 blocks",class="btn-primary"),div(class="small-muted mt-2","The cache is keyed to input-file hashes, population set, block size, missingness, MAF and transition filters, small-sample and pseudohaploid corrections, memory limit and ADMIXTOOLS 2 version. A single pseudohaploid target may produce missing f3 estimates when small-sample correction is enabled."))),column(6,div(class="al-card",h4("Run status and provenance"),verbatimTextOutput("status"),uiOutput("f2_status"),p(class="small-muted","A run manifest records UTC time, method, exact configuration, software version and dataset details."))))
  ),
  bslib::nav_panel("F-statistics", value="F-statistics",div(class="container-fluid px-4 py-3",div(class="page-title","F-statistics"),div(class="page-subtitle","Calculate f3, f4, D-statistics and f4-ratios with explicit data paths."),fluidRow(column(6,div(class="al-card",h4("Outgroup f3 / admixture f3"),p(class="small-muted","Tuple order: target, source 1, source 2. One tuple per line."),textAreaInput("triples","f3 population tuples",rows=7,placeholder="Target,Source_A,Source_B"),checkboxInput("f3_fixed_outgroup","These are outgroup f3 tests with the same fixed outgroup in column 3 (required with afprod)",FALSE),actionButton("run_f3","Run f3",class="btn-primary"),DT::DTOutput("f3table"))),column(6,div(class="al-card",h4("f4 / D-statistics"),p(class="small-muted","Tuple order: A, B, C, D. D-statistics use the direct genotype path (f4blockdat_from_geno), not the corrected shared f2 cache; review pseudohaploid coding and missingness carefully. f4 can reuse the f2 cache."),selectInput("fstat_kind","Statistic",choices=c("f4 statistic"="f4","D-statistic (genotype prefix required)"="D")),checkboxInput("d_allsnps","D-statistic: use all SNPs available per f4 (exploratory)",FALSE),textAreaInput("quads","f4 population tuples",rows=7,placeholder="A,B,C,D"),actionButton("run_f4","Run selected statistic",class="btn-primary"),DT::DTOutput("f4table"))))
    fluidRow(column(12,div(class="al-card",h4("F4-ratio admixture estimates"),p(class="small-muted","Five-population order: f4(1,2;3,4) / f4(1,2;5,4). This is a formal ratio estimator, not qpAdm; SNP filters are inherited from the shared f2 cache."),textAreaInput("f4ratios","Five-population tuples (one per line)",rows=4,placeholder="A,B,C,D,E"),checkboxInput("ratio_bootstrap","Use block bootstrap",FALSE),numericInput("ratio_boot_reps","Bootstrap replicates",100,min=10,max=10000,step=10),actionButton("run_f4ratio","Run f4 ratios",class="btn-primary"),DT::DTOutput("f4ratiotable"))))
  )),
  bslib::nav_panel("FST",value="FST",div(class="container-fluid px-4 py-3",div(class="page-title","Genetic differentiation"),div(class="page-subtitle","Pairwise population differentiation from the selected genotype panel."),div(class="al-card",fluidRow(column(5,selectizeInput("fst_a","Population set A",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button")))),column(5,selectizeInput("fst_b","Population set B",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button")))),column(2,br(),checkboxInput("fst_adjust_pseudo","Correct for pseudohaploid genotype coding",TRUE),checkboxInput("fst_allow_singletons","Allow populations with fewer than two individuals (not unbiased)",FALSE),checkboxInput("fst_bootstrap","Use block bootstrap instead of jackknife",FALSE),numericInput("fst_boot_reps","Bootstrap replicates (if enabled)",100,min=10,max=10000,step=10),actionButton("run_fst","Calculate FST",class="btn-primary"))),DT::DTOutput("fsttable"),p(class="small-muted mt-3","Check the upstream estimator, SNP ascertainment, missingness and uncertainty procedure before publication.")))
  ),
  bslib::nav_panel("qpWave",value="qpWave",div(class="container-fluid px-4 py-3",div(class="page-title","qpWave rank testing"),div(class="page-subtitle","Test the minimum number of ancestry streams distinguishable relative to a specified right set."),div(class="al-card",fluidRow(column(5,selectizeInput("wave_left","Left populations",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button")))),column(5,selectizeInput("wave_right","Right/reference populations",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button")))),column(2,br(),br(),numericInput("wave_fudge","Covariance diagonal fudge",0.0001,min=0,step=0.00001),checkboxInput("wave_bootstrap","Use block bootstrap",FALSE),numericInput("wave_boot_reps","Bootstrap replicates",100,min=10,max=10000,step=10),actionButton("run_qpwave","Run qpWave",class="btn-primary"))),h5("Rank / nested-model diagnostics"),DT::DTOutput("waverank"),h5("Underlying f4 statistics"),DT::DTOutput("wavef4"),verbatimTextOutput("wave_detail")))
  ),
  bslib::nav_panel("qpAdm",value="qpAdm",div(class="container-fluid px-4 py-3",div(class="page-title","qpAdm modelling"),div(class="page-subtitle","Fit an explicit target/source/right design. Treat results as conditional on the chosen populations and assumptions."),fluidRow(column(5,div(class="al-card",h4("Model specification"),selectInput("target","Target population",choices=c("Select target"="")),selectizeInput("sources","Left/source populations",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button"))),selectizeInput("right","Right/reference populations",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button"))),checkboxInput("qpadm_all_snps","Direct-genotype allsnps mode (exploratory; may inflate false positives with high missingness)",FALSE),p(class="small-muted","The direct-genotype allsnps path uses ADMIXTOOLS 2 f4blockdat_from_geno. It does not inherit the shared f2 cache's pseudohaploid correction or MAF/transition filters; prefer cached f2 for those corrections unless the allsnps trade-off is explicitly justified."),checkboxInput("qpadm_constrained","Constrain mixture coefficients",TRUE),numericInput("qpadm_fudge","Covariance diagonal fudge",0.0001,min=0,step=0.00001),checkboxInput("qpadm_fudge_twice","Apply doubled diagonal fudge (original qpAdm p-value comparability)",FALSE),checkboxInput("qpadm_return_f4","Return fitted f4 diagnostics",TRUE),checkboxInput("qpadm_bootstrap","Use block bootstrap",FALSE),numericInput("qpadm_boot_reps","Bootstrap replicates",100,min=10,max=10000,step=10),actionButton("qpadm","Run qpAdm",class="btn-primary"),hr(),h4("Model rotation"),p(class="small-muted","Evaluate alternative left/source combinations; rotation is a sensitivity analysis, not automatic model selection."),selectizeInput("rotate_left","leftright candidate populations",choices=NULL,multiple=TRUE,options=list(plugins=list("remove_button"))),actionButton("rotate","Run qpAdm rotation"))),column(7,div(class="al-card",h4("Model output"),DT::DTOutput("weights"),h5("F4 fit diagnostics"),DT::DTOutput("qpadm_f4"),h5("Rank / nested-model diagnostics"),DT::DTOutput("qpadm_rankdrop"),h5("Population-drop sensitivity"),DT::DTOutput("qpadm_popdrop"),verbatimTextOutput("qpadm_detail"),hr(),h4("Rotation output"),DT::DTOutput("rotate_table"),verbatimTextOutput("rotate_detail")))))
  ),
  bslib::nav_panel("Graphs",value="Graphs",div(class="container-fluid px-4 py-3",div(class="page-title","Admixture graph exploration"),div(class="page-subtitle","Fit a specified graph or search candidate topologies using ADMIXTOOLS 2."),fluidRow(column(6,div(class="al-card",h4("qpGraph fit"),p(class="small-muted","Upload an RDS containing a graph object in the documented ADMIXTOOLS 2 representation."),fileInput("graph_file","Graph RDS file",accept=".rds"),numericInput("starts","Random initialisations",20,min=1,max=1000),numericInput("qpgraph_seed","Random seed",1,min=0),numericInput("qpgraph_diag","Branch covariance diagonal regularization",0.0001,min=0,step=0.00001),numericInput("qpgraph_diag_f3","F3 covariance diagonal regularization",0.00001,min=0,step=0.000001),selectInput("qpgraph_lsqmode","Graph objective",choices=c("Full covariance likelihood"="FALSE","Diagonal approximation (diagnostic)"="TRUE","Unweighted residual sum (diagnostic)"="2")),checkboxInput("qpgraph_constrained","Constrain graph edge weights",TRUE),checkboxInput("qpgraph_return_fstats","Return f-statistic residuals",TRUE),checkboxInput("qpgraph_bootstrap","Use block bootstrap",FALSE),numericInput("qpgraph_boot_reps","Bootstrap replicates",100,min=10,max=10000,step=10),actionButton("run_qpgraph","Fit specified graph",class="btn-primary"),uiOutput("graph_metrics"),DT::DTOutput("graphfit"),verbatimTextOutput("graph_detail"))),column(6,div(class="al-card",h4("FindGraph search"),selectInput("outpop","Outgroup / root population",choices=c("Select outgroup"="")),numericInput("admix","Number of admixture events",1,min=0,max=8),numericInput("seconds","Search time limit (seconds)",1800,min=60),numericInput("numgraphs","Number of candidate graphs",10,min=1,max=500),numericInput("stop_gen","Search generation limit",100,min=1),numericInput("stop_gen2","Secondary generation limit",15,min=1),numericInput("max_admix","Maximum admixture events allowed",1,min=0,max=12),numericInput("stop_score","Stop at score (0 disables)",0,min=0),numericInput("reject_f4z","Reject zero-edge fits above |z|",0,min=0),checkboxInput("opt_worst_residual","Optimize worst residual instead of score",FALSE),checkboxInput("return_searchtree","Return search tree",FALSE),actionButton("find","Search graph space",class="btn-primary"),DT::DTOutput("graphs"),verbatimTextOutput("graph_search_detail")))))
  ),
  bslib::nav_panel("GhostBuster",value="GhostBuster",div(class="container-fluid px-4 py-3",div(class="page-title","GhostBuster genealogical analysis"),div(class="page-subtitle","Separate workflow based on Relate-inferred genealogical trees, genetic maps and sample/population mappings."),div(class="status-strip mb-3","Important: AADR genotype files alone are not valid GhostBuster input. Prepare Relate trees first and use a recombination map matching the genome build."),fluidRow(column(6,div(class="al-card",h4("Relate / GhostBuster inputs"),textInput("gb_dir","GhostBuster checkout directory",file.path(APP_ROOT,"external","GhostBuster")),textInput("gb_python","Python executable from GhostBuster Conda environment","python"),selectInput("gb_build","Genome build",choices=c("hg37","hg38"),selected="hg37"),textInput("gb_trees","Relate tree filename prefix (prefix + chromosome + .trees)"),textInput("gb_poplabels","Population labels file"),textInput("gb_map","Recombination map prefix/file"),textInput("gb_samples","Relate haplotype indices for --sample_id (space-separated)"),textInput("gb_chroms","Chromosomes (comma-separated)","1,2,3,4,5,6,7,8,9,10"),textInput("gb_out","Output prefix",file.path(APP_ROOT,"results","ghostbuster","run")),numericInput("gb_clusters","Number of ancestry clusters",2,min=1,max=10),numericInput("gb_repeats","EM random restarts",20,min=1,max=500),numericInput("gb_iters","EM iterations",200,min=1,max=5000),numericInput("gb_seed","Random seed",2,min=0),checkboxInput("gb_hmm","Run HMM (recommended default)",TRUE),checkboxInput("gb_pickle","Only build intermediate pickle files",FALSE),checkboxInput("gb_aggregate","Aggregate existing chromosome pickle files (omit --trees)",FALSE),actionButton("run_ghostbuster","Run GhostBuster",class="btn-primary"),hr(),h4("Downstream LD-curve dating"),selectInput("gb_ld_mode","Admixture-date model",choices=c("Run all models"="all","One pulse"="1date","Two pulses"="2date","Continuous admixture"="continuous"),selected="all"),actionButton("run_gb_ld","Run LD-curve dating"))),column(6,div(class="al-card",h4("Run log"),verbatimTextOutput("gb_log"),verbatimTextOutput("gb_ld_log"),p(class="small-muted","Record the Relate version, tree-inference parameters, sample-ID mapping, map provenance and GhostBuster commit in the project notes.")))))
  ),
  bslib::nav_panel("Run archive",value="Run archive",div(class="container-fluid px-4 py-3",div(class="page-title","Run archive and reproducibility"),div(class="page-subtitle","Browse persisted analysis runs and export result tables or manifests."),div(class="al-card",DT::DTOutput("all_runs"),downloadButton("download_archive","Download project run archive"),actionButton("refresh_runs","Refresh archive"))))
)

server <- function(input, output, session) {
  state <- shiny::reactiveValues(f2=NULL, f2_signature=NULL, status="Waiting for dataset selection.", dataset=NULL,
    labels=NULL, population_inventory=NULL, label_note=NULL, runs=data.frame(), last_result=NULL, last_run=NULL)
  results_dir <- file.path(APP_ROOT,"results")
  dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
  refresh_runs <- function() {
    dirs <- list.dirs(file.path(results_dir,"runs"),recursive=FALSE,full.names=TRUE)
    if (!length(dirs)) { state$runs <- data.frame(); return(invisible(NULL)) }
    rows <- lapply(dirs,function(d){ mf<-file.path(d,"manifest.json"); x<-tryCatch(jsonlite::read_json(mf,simplifyVector=TRUE),error=function(e)list()); data.frame(run_id=basename(d),method=x$config$method %||% "unknown",timestamp=x$timestamp_utc %||% "",status=x$status %||% "saved",path=d,stringsAsFactors=FALSE) })
    state$runs <- do.call(rbind,rows[order(vapply(rows,function(z)z$timestamp,character(1)),decreasing=TRUE)])
  }
  record <- function(method,config,result) {
    is_genealogy <- method %in% c("GhostBuster", "GhostBuster_LD_curve_dating")
    uses_f2 <- method %in% c("f3", "f4", "f4_ratio", "qpWave", "qpAdm", "qpAdm_rotate", "qpGraph", "FindGraph")
    if (identical(method,"qpAdm") && isTRUE(config$direct_genotype_allsnps)) uses_f2 <- FALSE
    ds <- list(genotype_dataset = NULL, f2_configuration = if(uses_f2) state$f2_config %||% NULL else NULL)
    if (!is_genealogy && !is.null(state$dataset)) {
      ds$genotype_dataset <- list(prefix=state$dataset$prefix,format=state$dataset$format,
        file_md5=as.list(tools::md5sum(state$dataset$files)))
    }
    saved <- persist_analysis(results_dir,method,config,result,ds)
    state$last_result <- result; state$last_run <- saved
    refresh_runs(); state$status <- paste0(method," completed and saved as ",saved$run_id)
    saved
  }
  attach_run_input <- function(saved, source_path, stored_name, original_name) {
    if (is.null(saved) || !file.exists(source_path)) stop("Could not preserve the uploaded analysis input.")
    dest <- file.path(saved$path, stored_name)
    if (!file.copy(source_path, dest, overwrite=TRUE)) stop("Could not copy the uploaded analysis input into the run archive.")
    manifest_path <- file.path(saved$path, "manifest.json")
    manifest <- jsonlite::read_json(manifest_path, simplifyVector=FALSE)
    manifest$uploaded_inputs <- list(list(original_name=original_name, stored_name=stored_name,
                                           md5=unname(tools::md5sum(dest))))
    jsonlite::write_json(manifest, manifest_path, auto_unbox=TRUE, pretty=TRUE, null="null", dataframe="rows")
    invisible(dest)
  }
  attach_external_outputs <- function(saved, output_prefix, input_files = character(), python_version = NULL) {
    if (is.null(saved)) return(invisible(NULL))
    parent <- dirname(output_prefix)
    prefix_name <- basename(output_prefix)
    files <- if (dir.exists(parent)) list.files(parent, full.names=TRUE) else character()
    files <- files[startsWith(basename(files), prefix_name) & file.exists(files) & !dir.exists(files)]
    details <- lapply(files, function(f) list(path=normalizePath(f,mustWork=TRUE), size_bytes=unname(file.info(f)$size), md5=unname(tools::md5sum(f))))
    manifest_path <- file.path(saved$path,"manifest.json")
    manifest <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
    manifest$external_outputs <- details
    input_files <- unique(input_files[file.exists(input_files) & !dir.exists(input_files)])
    manifest$external_inputs <- lapply(input_files, function(f) list(path=normalizePath(f,mustWork=TRUE),size_bytes=unname(file.info(f)$size),md5=unname(tools::md5sum(f))))
    manifest$python_version <- python_version
    jsonlite::write_json(manifest,manifest_path,auto_unbox=TRUE,pretty=TRUE,null="null",dataframe="rows")
    invisible(details)
  }
  require_current_dataset <- function() {
    if (is.null(state$dataset)) stop("Inspect a genotype dataset first.")
    info <- validate_prefix(input$prefix)
    if (!identical(info$prefix,state$dataset$prefix)) stop("The genotype prefix changed after inspection. Inspect the selected dataset again.")
    invisible(info)
  }
  require_current_f2 <- function() {
    require_current_dataset()
    if (is.null(state$f2) || is.null(state$f2_config)) stop("Compute f2 blocks first.")
    current <- list(prefix=normalizePath(input$prefix,mustWork=TRUE),
                    populations=sort(unique(split_lines(input$populations))),
                    block_size_morgans=input$block, max_missingness=input$maxmiss,
                    afprod=isTRUE(input$afprod), adjust_pseudohaploid=isTRUE(input$f2_adjust_pseudo),
                    apply_corr=isTRUE(input$apply_corr),auto_only=isTRUE(input$auto_only),
                    minmaf=input$minmaf,maxmaf=input$maxmaf,transitions=isTRUE(input$transitions),transversions=isTRUE(input$transversions),maxmem=input$maxmem)
    if (!identical(current,state$f2_config)) stop("The genotype prefix or shared f2 settings changed after computation. Recompute f2 before running this method.")
    invisible(state$f2)
  }
  output$status <- shiny::renderText(state$status)
  output$global_status <- shiny::renderUI({ tags$div(class="status-strip",state$status) })
  output$metric_dataset <- shiny::renderText(if(is.null(state$dataset))"Not selected" else basename(state$dataset$prefix))
  output$metric_pops <- shiny::renderText(if(is.null(state$labels))"—" else length(state$labels))
  output$metric_runs <- shiny::renderText(nrow(state$runs))
  output$metric_engine <- shiny::renderText(if(requireNamespace("admixtools",quietly=TRUE))paste0("ADMIXTOOLS ",utils::packageVersion("admixtools")) else "Not installed")
  output$workspace_status <- shiny::renderUI({ if(is.null(state$dataset)) return(tags$div(class="status-strip","No genotype dataset inspected. Open Project & data to upload or select a dataset.")); tags$div(class="status-strip",strong("Dataset inspected: "),state$dataset$prefix,tags$br(),state$label_note) })
  output$dataset_summary <- shiny::renderUI({ if(is.null(state$dataset)) return(tags$p(class="small-muted","Inspect a genotype prefix to view file format, sample count and population labels.")); tags$div(tags$p(tags$b("Format: "),state$dataset$format),tags$p(tags$b("Prefix: "),state$dataset$prefix),tags$p(tags$b("Groups: "),length(state$labels)),tags$p(class="small-muted",state$label_note)) })
  output$population_table <- DT::renderDT({ req(state$labels); DT::datatable(state$population_inventory,rownames=FALSE,options=list(pageLength=12,searching=TRUE,dom="ftip")) })
  output$f2_status <- shiny::renderUI({ if(is.null(state$f2)) tags$p(class="small-muted","f2 blocks not computed in this session.") else tags$div(class="status-strip","f2 blocks loaded in memory. Downstream methods can now be run.") })
  output$recent_runs <- DT::renderDT({ DT::datatable(utils::head(state$runs,5)[,intersect(c("run_id","method","status"),names(state$runs)),drop=FALSE],rownames=FALSE,options=list(dom="t",pageLength=5)) })
  output$all_runs <- DT::renderDT({ DT::datatable(state$runs[,intersect(c("run_id","method","timestamp","status"),names(state$runs)),drop=FALSE],rownames=FALSE,selection="single",options=list(pageLength=15,scrollX=TRUE)) })
  output$dataset_status <- shiny::renderText(state$dataset_status %||% "No AADR download requested.")
  output$upload_status <- shiny::renderText(state$upload_status %||% "")
  output$f3table <- DT::renderDT({ req(state$f3_result); DT::datatable(result_table(state$f3_result),options=list(pageLength=10,scrollX=TRUE)) })
  output$f4table <- DT::renderDT({ req(state$f4_result); DT::datatable(result_table(state$f4_result),options=list(pageLength=10,scrollX=TRUE)) })
  output$f4ratiotable <- DT::renderDT({ req(state$f4ratio_result); DT::datatable(result_table(state$f4ratio_result),options=list(pageLength=10,scrollX=TRUE)) })
  output$fsttable <- DT::renderDT({ req(state$fst_result); DT::datatable(result_table(state$fst_result),options=list(pageLength=15,scrollX=TRUE)) })
  output$waverank <- DT::renderDT({ req(state$wave_result); z<-state$wave_result$rankdrop %||% data.frame(); DT::datatable(result_table(z),options=list(pageLength=10,scrollX=TRUE)) })
  output$wavef4 <- DT::renderDT({ req(state$wave_result); z<-state$wave_result$f4 %||% data.frame(); DT::datatable(result_table(z),options=list(pageLength=10,scrollX=TRUE)) })
  output$wave_detail <- shiny::renderPrint({ req(state$wave_result); str(state$wave_result,max.level=2) })
  output$weights <- DT::renderDT({ req(state$qpadm_result); DT::datatable(result_table(state$qpadm_result),options=list(pageLength=12,scrollX=TRUE)) })
  output$qpadm_detail <- shiny::renderPrint({ req(state$qpadm_result); str(state$qpadm_result,max.level=3) })
  output$qpadm_f4 <- DT::renderDT({ req(state$qpadm_result); z<-state$qpadm_result$f4 %||% data.frame(); DT::datatable(result_table(z),options=list(pageLength=8,scrollX=TRUE)) })
  output$qpadm_rankdrop <- DT::renderDT({ req(state$qpadm_result); z<-state$qpadm_result$rankdrop %||% data.frame(); DT::datatable(as.data.frame(z),options=list(pageLength=8,scrollX=TRUE)) })
  output$qpadm_popdrop <- DT::renderDT({ req(state$qpadm_result); z<-state$qpadm_result$popdrop %||% data.frame(); DT::datatable(as.data.frame(z),options=list(pageLength=8,scrollX=TRUE)) })
  output$rotate_table <- DT::renderDT({ req(state$rotate_result); DT::datatable(result_table(state$rotate_result),options=list(pageLength=12,scrollX=TRUE)) })
  output$rotate_detail <- shiny::renderPrint({ req(state$rotate_result); str(state$rotate_result,max.level=2) })
  output$graph_metrics <- shiny::renderUI({ req(state$graph_result); tags$div(class="status-strip",paste0("Fit score: ",state$graph_result$score %||% "not returned","  |  Worst residual: ",state$graph_result$worst_residual %||% "not returned")) })
  output$graphfit <- DT::renderDT({ req(state$graph_result); DT::datatable(result_table(state$graph_result),options=list(pageLength=10,scrollX=TRUE)) })
  output$graph_detail <- shiny::renderPrint({ req(state$graph_result); str(state$graph_result,max.level=3) })
  output$graphs <- DT::renderDT({ req(state$find_result); DT::datatable(result_table(state$find_result),options=list(pageLength=10,scrollX=TRUE)) })
  output$graph_search_detail <- shiny::renderPrint({ req(state$find_result); str(state$find_result,max.level=2) })
  output$gb_log <- shiny::renderText(state$gb_log %||% "GhostBuster has not been run.")
  output$gb_ld_log <- shiny::renderText(state$gb_ld_log %||% "LD-curve dating has not been run.")

  shiny::observeEvent(input$inspect_data, {
    tryCatch({ info<-validate_prefix(input$prefix); pop<-read_population_labels(input$prefix); state$f2<-NULL; state$f2_config<-NULL; state$dataset<-info; state$labels<-pop$labels; state$population_inventory<-pop$inventory; state$label_note<-pop$note; shiny::updateSelectizeInput(session,"populations",choices=pop$labels,selected=intersect(isolate(input$populations),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"fst_a",choices=pop$labels,selected=intersect(isolate(input$fst_a),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"fst_b",choices=pop$labels,selected=intersect(isolate(input$fst_b),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"wave_left",choices=pop$labels,selected=intersect(isolate(input$wave_left),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"wave_right",choices=pop$labels,selected=intersect(isolate(input$wave_right),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"sources",choices=pop$labels,selected=intersect(isolate(input$sources),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"right",choices=pop$labels,selected=intersect(isolate(input$right),pop$labels),server=TRUE); shiny::updateSelectizeInput(session,"rotate_left",choices=pop$labels,selected=intersect(isolate(input$rotate_left),pop$labels),server=TRUE); shiny::updateSelectInput(session,"target",choices=c("Select target"="",pop$labels),selected=if(isolate(input$target)%in%pop$labels)isolate(input$target)else ""); shiny::updateSelectInput(session,"outpop",choices=c("Select outgroup"="",pop$labels),selected=if(isolate(input$outpop)%in%pop$labels)isolate(input$outpop)else ""); state$status<-paste("Dataset inspected:",info$format,"—",length(pop$labels),"population labels found.") },error=function(e)state$status<-paste("Dataset inspection stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$eigenfiles, {
    req(input$eigenfiles); f<-input$eigenfiles; exts<-tolower(tools::file_ext(f$name)); valid<-setequal(exts,c("geno","snp","ind"))||setequal(exts,c("bed","bim","fam"))
    if(nrow(f)!=3L||!valid){state$upload_status<-"Upload exactly one matching EIGENSTRAT triplet or one matching PLINK triplet.";return()}
    stems<-sub("\\.(geno|snp|ind|bed|bim|fam)$","",f$name,ignore.case=TRUE);if(length(unique(stems))!=1L){state$upload_status<-"All three files must share one basename.";return()}
    dest<-file.path(APP_ROOT,"data","uploads",paste0("upload_",format(Sys.time(),"%Y%m%d_%H%M%S")));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    for(i in seq_len(nrow(f)))file.copy(f$datapath[i],file.path(dest,f$name[i]),overwrite=TRUE)
    newprefix<-file.path(dest,stems[1]);shiny::updateTextInput(session,"prefix",value=newprefix);state$upload_status<-paste("Uploaded:",paste(f$name,collapse=", "));state$status<-"Upload complete. Click Inspect dataset to validate metadata."
  })
  shiny::observeEvent(input$download_aadr, {
    args<-c(file.path(APP_ROOT,"scripts","download_aadr.py"),"--panel",input$panel,"--out",file.path(APP_ROOT,"data","aadr"));if(isTRUE(input$metadata_only))args<-c(args,"--metadata-only")
    log<-tryCatch(system2("python3",args,stdout=TRUE,stderr=TRUE),error=function(e)e);state$dataset_status<-if(inherits(log,"error"))conditionMessage(log)else paste(log,collapse="\n")
  })
  shiny::observeEvent(input$f2, {
    tryCatch({ require_current_dataset(); pops<-split_lines(input$populations); withProgress(message="Estimating blockwise f2",value=0.1,{state$f2<-compute_f2(input$prefix,file.path(results_dir,"f2_cache"),pops,
          blgsize=input$block,maxmiss=input$maxmiss,afprod=input$afprod,
          adjust_pseudohaploid=input$f2_adjust_pseudo,apply_corr=input$apply_corr,
          auto_only=input$auto_only,minmaf=input$minmaf,maxmaf=input$maxmaf,
          transitions=input$transitions,transversions=input$transversions,maxmem=input$maxmem,overwrite=input$overwrite_f2);incProgress(.8)});state$status<-paste("f2 cache ready for",length(pops),"exact population labels.");state$f2_config<-list(prefix=normalizePath(input$prefix,mustWork=TRUE),populations=sort(unique(pops)),block_size_morgans=input$block,max_missingness=input$maxmiss,afprod=isTRUE(input$afprod),adjust_pseudohaploid=isTRUE(input$f2_adjust_pseudo),apply_corr=isTRUE(input$apply_corr),auto_only=isTRUE(input$auto_only),minmaf=input$minmaf,maxmaf=input$maxmaf,transitions=isTRUE(input$transitions),transversions=isTRUE(input$transversions),maxmem=input$maxmem) },error=function(e)state$status<-paste("f2 stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_f3, {
    tryCatch({require_current_f2();tuples<-parse_tuples(input$triples,3,"f3");validate_population_labels(input$prefix,as.vector(tuples));if(isTRUE(state$f2_config$afprod)){if(!isTRUE(input$f3_fixed_outgroup))stop("afprod blocks are valid for outgroup f3 only when the same outgroup is fixed in column 3.");if(length(unique(tuples[,3]))!=1L)stop("All afprod outgroup-f3 tuples must use the same third-column outgroup.")};ans<-run_f3(state$f2,tuples);state$f3_result<-ans;record("f3",list(tuples=tuples,fixed_outgroup=input$f3_fixed_outgroup),ans)},error=function(e)state$status<-paste("f3 stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_f4, {
    tryCatch({tuples<-parse_tuples(input$quads,4,"f4/D");validate_population_labels(input$prefix,as.vector(tuples));if(identical(input$fstat_kind,"D")){require_current_dataset();ans<-run_dstat(input$prefix,tuples,blgsize=input$block,auto_only=isTRUE(input$auto_only),allsnps=isTRUE(input$d_allsnps));method<-"D_statistic"}else{require_current_f2();ans<-run_f4(state$f2,tuples);method<-"f4"};state$f4_result<-ans;record(method,list(tuples=tuples,statistic=input$fstat_kind,block_size_morgans=input$block,allsnps=input$d_allsnps,auto_only=input$auto_only),ans)},error=function(e)state$status<-paste("f4/D stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_f4ratio, {
    tryCatch({require_current_f2();tuples<-parse_tuples(input$f4ratios,5,"f4-ratio");validate_population_labels(input$prefix,as.vector(tuples));ratio_boot<-if(isTRUE(input$ratio_bootstrap))as.integer(input$ratio_boot_reps)else FALSE;ans<-run_f4ratio(state$f2,tuples,boot=ratio_boot);state$f4ratio_result<-ans;record("f4_ratio",list(tuples=tuples,bootstrap=ratio_boot),ans)},error=function(e)state$status<-paste("f4-ratio stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_fst, {
    tryCatch({require_current_dataset();a<-split_lines(input$fst_a);b<-split_lines(input$fst_b);validate_population_labels(input$prefix,c(a,b));if(length(intersect(a,b)))stop("Population sets A and B must not overlap.");singletons<-state$population_inventory$population[state$population_inventory$n_individuals<2 & state$population_inventory$population%in%c(a,b)];if(length(singletons)&&!isTRUE(input$fst_allow_singletons))stop("FST is not unbiased for populations with fewer than two individuals: ",paste(singletons,collapse=", "),". Enable the explicit override only for exploratory output.");boot<-if(isTRUE(input$fst_bootstrap))as.integer(input$fst_boot_reps)else FALSE;ans<-run_fst(input$prefix,a,b,boot=boot,adjust_pseudohaploid=isTRUE(input$fst_adjust_pseudo),auto_only=isTRUE(input$auto_only),minmaf=input$minmaf,maxmaf=input$maxmaf,transitions=isTRUE(input$transitions),transversions=isTRUE(input$transversions),maxmem=input$maxmem);state$fst_result<-ans;record("FST",list(population_set_A=a,population_set_B=b,bootstrap=boot,adjust_pseudohaploid=input$fst_adjust_pseudo,allow_singletons=input$fst_allow_singletons,auto_only=input$auto_only,minmaf=input$minmaf,maxmaf=input$maxmaf,transitions=input$transitions,transversions=input$transversions,maxmem=input$maxmem),ans)},error=function(e)state$status<-paste("FST stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_qpwave, {
    tryCatch({require_current_f2();l<-split_lines(input$wave_left);r<-split_lines(input$wave_right);validate_population_labels(input$prefix,c(l,r));wave_boot<-if(isTRUE(input$wave_bootstrap))as.integer(input$wave_boot_reps)else FALSE;ans<-run_qpwave(state$f2,l,r,fudge=input$wave_fudge,boot=wave_boot);state$wave_result<-ans;record("qpWave",list(left=l,right=r,fudge=input$wave_fudge,bootstrap=wave_boot),ans)},error=function(e)state$status<-paste("qpWave stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$qpadm, {
    tryCatch({require_current_dataset();if(!isTRUE(input$qpadm_all_snps))require_current_f2();t<-trimws(input$target);s<-split_lines(input$sources);r<-split_lines(input$right);validate_population_design(t,s,r);validate_population_labels(input$prefix,c(t,s,r));qdata<-if(isTRUE(input$qpadm_all_snps)) input$prefix else state$f2;qargs<-c(list(data=qdata,left=s,right=r,target=t,constrained=isTRUE(input$qpadm_constrained),fudge=input$qpadm_fudge,fudge_twice=isTRUE(input$qpadm_fudge_twice),return_f4=isTRUE(input$qpadm_return_f4),boot=if(isTRUE(input$qpadm_bootstrap))as.integer(input$qpadm_boot_reps)else FALSE),if(isTRUE(input$qpadm_all_snps))list(allsnps=TRUE,blgsize=input$block,auto_only=isTRUE(input$auto_only))else list());ans<-do.call(run_qpadm,qargs);state$qpadm_result<-ans;record("qpAdm",list(target=t,sources=s,right=r,direct_genotype_allsnps=input$qpadm_all_snps,constrained=input$qpadm_constrained,fudge=input$qpadm_fudge,fudge_twice=input$qpadm_fudge_twice,return_f4=input$qpadm_return_f4,bootstrap=if(isTRUE(input$qpadm_bootstrap))as.integer(input$qpadm_boot_reps)else FALSE,block_size_morgans=input$block,auto_only=input$auto_only,shared_f2_adjustment_applied=!isTRUE(input$qpadm_all_snps),shared_f2_filters_applied=!isTRUE(input$qpadm_all_snps)),ans)},error=function(e)state$status<-paste("qpAdm stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$rotate, {
    tryCatch({require_current_dataset();require_current_f2();t<-trimws(input$target);l<-split_lines(input$rotate_left);r<-split_lines(input$right);validate_population_labels(input$prefix,c(t,l,r));ans<-run_qpadm_rotate(state$f2,t,l,r,full_results=TRUE);state$rotate_result<-ans;record("qpAdm_rotate",list(target=t,leftright=l,rightfix=r,full_results=TRUE),ans)},error=function(e)state$status<-paste("qpAdm rotation stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_qpgraph, {
    tryCatch({require_current_f2();req(input$graph_file);g<-readRDS(input$graph_file$datapath);graph_boot<-if(isTRUE(input$qpgraph_bootstrap))as.integer(input$qpgraph_boot_reps)else FALSE;lsqmode<-switch(input$qpgraph_lsqmode,"TRUE"=TRUE,"2"=2,FALSE);ans<-run_qpgraph(state$f2,g,numstart=input$starts,seed=as.integer(input$qpgraph_seed),diag=input$qpgraph_diag,diag_f3=input$qpgraph_diag_f3,lsqmode=lsqmode,constrained=isTRUE(input$qpgraph_constrained),return_fstats=isTRUE(input$qpgraph_return_fstats),boot=graph_boot);state$graph_result<-ans;saved<-record("qpGraph",list(graph_file=input$graph_file$name,numstart=input$starts,seed=input$qpgraph_seed,diag=input$qpgraph_diag,diag_f3=input$qpgraph_diag_f3,lsqmode=input$qpgraph_lsqmode,constrained=input$qpgraph_constrained,return_fstats=input$qpgraph_return_fstats,bootstrap=graph_boot),ans);attach_run_input(saved,input$graph_file$datapath,"input_graph.rds",input$graph_file$name)},error=function(e)state$status<-paste("qpGraph stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$admix, {
    if (!is.null(input$max_admix) && input$max_admix < input$admix) shiny::updateNumericInput(session,"max_admix",value=input$admix)
  }, ignoreInit=TRUE)
  shiny::observeEvent(input$find, {
    tryCatch({require_current_f2();validate_population_labels(input$prefix,input$outpop);ans<-run_find_graphs(state$f2,numadmix=input$admix,outpop=trimws(input$outpop),stop_gen=input$stop_gen,stop_gen2=input$stop_gen2,stop_sec=input$seconds,numgraphs=input$numgraphs,stop_score=input$stop_score,max_admix=input$max_admix,opt_worst_residual=isTRUE(input$opt_worst_residual),reject_f4z=input$reject_f4z,return_searchtree=isTRUE(input$return_searchtree));state$find_result<-ans;record("FindGraph",list(outpop=input$outpop,numadmix=input$admix,stop_gen=input$stop_gen,stop_gen2=input$stop_gen2,stop_sec=input$seconds,numgraphs=input$numgraphs,stop_score=input$stop_score,max_admix=input$max_admix,opt_worst_residual=input$opt_worst_residual,reject_f4z=input$reject_f4z,return_searchtree=input$return_searchtree),ans)},error=function(e)state$status<-paste("FindGraph stopped:",conditionMessage(e)))
  })
  shiny::observeEvent(input$run_ghostbuster, {
    tryCatch({out<-normalizePath(dirname(input$gb_out),mustWork=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE);ans<-run_ghostbuster(file.path(APP_ROOT,"scripts","run_ghostbuster.sh"),input$gb_dir,input$gb_build,input$gb_trees,input$gb_poplabels,input$gb_map,input$gb_samples,input$gb_chroms,input$gb_out,input$gb_clusters,input$gb_repeats,input$gb_iters,input$gb_seed,input$gb_hmm,input$gb_pickle,input$gb_aggregate,input$gb_python);state$gb_log<-paste(ans$log,collapse="\n");saved<-record("GhostBuster",list(genome_build=input$gb_build,trees_prefix=input$gb_trees,poplabels=input$gb_poplabels,recombination_map=input$gb_map,sample_ids=input$gb_samples,chromosomes=input$gb_chroms,output_prefix=input$gb_out,num_clusters=input$gb_clusters,n_repeats=input$gb_repeats,num_iters=input$gb_iters,seed=input$gb_seed,hmm=input$gb_hmm,only_pickle=input$gb_pickle,aggregate=input$gb_aggregate,python_bin=input$gb_python),ans);attach_external_outputs(saved,input$gb_out,ans$input_files,ans$python_version)},error=function(e){state$gb_log<-conditionMessage(e);state$status<-paste("GhostBuster stopped:",conditionMessage(e))})
  })
  shiny::observeEvent(input$run_gb_ld, {
    tryCatch({ans<-run_ghostbuster_ld_dating(input$gb_dir,input$gb_out,input$gb_build,input$gb_ld_mode,input$gb_python);state$gb_ld_log<-paste(ans$log,collapse="\n");saved<-record("GhostBuster_LD_curve_dating",list(genome_build=input$gb_build,output_prefix=input$gb_out,mode=input$gb_ld_mode,python_executable=input$gb_python),ans);attach_external_outputs(saved,input$gb_out,ans$input_files,ans$python_version)},error=function(e){state$gb_ld_log<-conditionMessage(e);state$status<-paste("LD-curve dating stopped:",conditionMessage(e))})
  })
  shiny::observeEvent(input$refresh_runs,refresh_runs())
  output$download_archive <- shiny::downloadHandler(filename=function()paste0("ancestrylab-runs-",Sys.Date(),".zip"),content=function(file){dirs<-list.dirs(file.path(results_dir,"runs"),recursive=FALSE,full.names=TRUE);if(!length(dirs))stop("No saved runs yet.");old<-getwd();on.exit(setwd(old),add=TRUE);setwd(results_dir);utils::zip(file,file.path("runs",basename(dirs)),flags="-r9X")})
  refresh_runs()
}

shiny::shinyApp(ui,server)
