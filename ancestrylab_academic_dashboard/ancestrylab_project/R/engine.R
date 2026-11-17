# AncestryLab analysis engine
# Statistical calculations are delegated to the upstream ADMIXTOOLS 2 package.
`%||%` <- function(x, y) if (is.null(x)) y else x

require_admixtools <- function() {
  if (!requireNamespace("admixtools", quietly = TRUE)) {
    stop("ADMIXTOOLS 2 is not installed. Run scripts/install.R first.", call. = FALSE)
  }
}

validate_prefix <- function(prefix) {
  if (length(prefix) != 1L || is.na(prefix) || !nzchar(trimws(prefix))) stop("Select a genotype file prefix.", call. = FALSE)
  prefix <- normalizePath(prefix, mustWork = TRUE)
  eigen_files <- paste0(prefix, c(".geno", ".snp", ".ind"))
  plink_files <- paste0(prefix, c(".bed", ".bim", ".fam"))
  if (all(file.exists(eigen_files))) return(list(prefix = prefix, format = "eigenstrat", files = eigen_files))
  if (all(file.exists(plink_files))) return(list(prefix = prefix, format = "plink", files = plink_files))
  stop("Input must be a complete EIGENSTRAT (.geno/.snp/.ind) or PLINK (.bed/.bim/.fam) triplet. No conversion or proxy substitution is performed.", call. = FALSE)
}

read_population_labels <- function(prefix) {
  info <- validate_prefix(prefix)
  if (identical(info$format, "eigenstrat")) {
    ind <- utils::read.table(paste0(info$prefix, ".ind"), header = FALSE, stringsAsFactors = FALSE,
                             comment.char = "", quote = "", fill = TRUE)
    if (ncol(ind) < 3L) stop("The .ind file must contain sample ID, sex and population columns.")
    group_values <- trimws(as.character(ind[[3L]]))
    labels <- sort(unique(group_values[!is.na(group_values)]))
    inventory <- as.data.frame(table(population = group_values), stringsAsFactors = FALSE)
    names(inventory) <- c("population", "n_individuals")
    inventory <- inventory[nzchar(inventory$population), , drop = FALSE]
    return(list(labels = labels[nzchar(labels)], inventory = inventory, format = info$format,
                note = "Population labels and individual counts read from the third column of the EIGENSTRAT .ind file."))
  }
  fam <- utils::read.table(paste0(info$prefix, ".fam"), header = FALSE, stringsAsFactors = FALSE)
  if (ncol(fam) < 2L) stop("Could not parse PLINK .fam file.")
  group_values <- as.character(fam[[1L]])
  labels <- sort(unique(group_values[!is.na(group_values)]))
  inventory <- as.data.frame(table(population = group_values), stringsAsFactors = FALSE)
  names(inventory) <- c("population", "n_individuals")
  inventory <- inventory[nzchar(inventory$population), , drop = FALSE]
  list(labels = labels[nzchar(labels)], inventory = inventory, format = info$format,
       note = "PLINK FID values and counts are shown as candidate groups. Confirm that FID encodes populations in this dataset before analysis.")
}

validate_population_labels <- function(prefix, populations) {
  populations <- unique(trimws(as.character(populations)))
  populations <- populations[!is.na(populations) & nzchar(populations)]
  if (!length(populations)) stop("No population labels were supplied.", call. = FALSE)
  available <- read_population_labels(prefix)$labels
  missing <- setdiff(populations, available)
  if (length(missing)) stop("These exact population labels are absent from the genotype metadata: ", paste(missing, collapse = ", "), ". No nearest-population substitution is allowed.", call. = FALSE)
  invisible(TRUE)
}

validate_population_design <- function(target, sources, right) {
  target <- trimws(as.character(target)); sources <- trimws(as.character(sources)); right <- trimws(as.character(right))
  vals <- c(target, sources, right)
  if (length(target) != 1L || anyNA(vals) || any(!nzchar(vals))) stop("Population names must be non-empty; target must be one label.", call. = FALSE)
  if (anyDuplicated(sources)) stop("Source populations must be unique.", call. = FALSE)
  if (anyDuplicated(right)) stop("Right/reference populations must be unique.", call. = FALSE)
  if (target %in% sources) stop("Target cannot also be a source.", call. = FALSE)
  if (target %in% right) stop("Target cannot also be a right/reference population.", call. = FALSE)
  if (length(sources) < 2L) stop("qpAdm requires at least two source populations.", call. = FALSE)
  if (length(unique(right)) < 2L) stop("Use at least two distinct right/reference populations.", call. = FALSE)
  if (length(intersect(sources, right))) stop("Do not place the same population on both the left-source and right-reference sets without an explicitly justified design.", call. = FALSE)
  invisible(TRUE)
}

validate_f2_population_labels <- function(f2_blocks, populations) {
  available <- dimnames(f2_blocks)[[1L]]
  if (is.null(available) || !length(available)) stop("The f2 object has no population dimension names; recompute it from a labelled genotype dataset.", call. = FALSE)
  requested <- unique(trimws(as.character(populations)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, available)
  if (length(missing)) stop("These exact populations are not in the current f2 cache: ", paste(missing, collapse = ", "), ". Recompute f2 with those labels selected.", call. = FALSE)
  invisible(TRUE)
}

compute_f2 <- function(prefix, cache_dir, pops, blgsize = 0.05, maxmiss = 0.0, afprod = FALSE,
                       adjust_pseudohaploid = TRUE, apply_corr = TRUE, auto_only = TRUE,
                       minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                       maxmem = 8000, overwrite = FALSE) {
  require_admixtools(); input_info <- validate_prefix(prefix)
  if (!is.numeric(blgsize) || length(blgsize) != 1L || !is.finite(blgsize) || blgsize <= 0) stop("Block size must be positive Morgans.")
  if (!is.numeric(maxmiss) || length(maxmiss) != 1L || !is.finite(maxmiss) || maxmiss < 0 || maxmiss > 1) stop("Maximum missingness must be between 0 and 1.")
  if (!is.numeric(minmaf) || !is.numeric(maxmaf) || length(minmaf) != 1L || length(maxmaf) != 1L || !is.finite(minmaf) || !is.finite(maxmaf) || minmaf < 0 || maxmaf > 0.5 || minmaf > maxmaf) stop("MAF filters must satisfy 0 <= minmaf <= maxmaf <= 0.5.")
  if (!is.numeric(maxmem) || length(maxmem) != 1L || maxmem < 100 || !is.finite(maxmem)) stop("Maximum memory must be at least 100 MB.")
  if (!isTRUE(transitions) && !isTRUE(transversions)) stop("At least one of transitions or transversions must be included.")
  validate_population_labels(prefix, pops)
  pops <- sort(unique(trimws(as.character(pops))))
  if (length(pops) < 2L) stop("Select at least two populations for f2 estimation.")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  hashes <- unname(tools::md5sum(input_info$files))
  signature <- digest::digest(list(hashes = hashes, blgsize = blgsize, maxmiss = maxmiss,
                                  pops = pops, afprod = isTRUE(afprod), adjust_pseudohaploid = adjust_pseudohaploid,
                                  apply_corr = apply_corr, auto_only = auto_only, minmaf = minmaf, maxmaf = maxmaf,
                                  transitions = transitions, transversions = transversions, maxmem = maxmem,
                                  admixtools = as.character(utils::packageVersion("admixtools"))), algo = "sha256")
  cache_file <- file.path(cache_dir, paste0("f2_blocks_", signature, ".rds"))
  if (file.exists(cache_file) && !isTRUE(overwrite)) return(readRDS(cache_file))
  args <- list(pref = input_info$prefix, blgsize = blgsize, maxmiss = maxmiss, pops = pops,
                afprod = isTRUE(afprod), adjust_pseudohaploid = adjust_pseudohaploid,
                apply_corr = apply_corr, auto_only = auto_only, minmaf = minmaf, maxmaf = maxmaf,
                transitions = transitions, transversions = transversions, maxmem = maxmem, verbose = TRUE)
  blocks <- do.call(admixtools::f2_from_geno, args)
  saveRDS(blocks, cache_file, compress = FALSE)
  manifest <- list(signature = signature, input_files = as.list(setNames(hashes, basename(input_info$files))),
    genotype_prefix = input_info$prefix, format = input_info$format, populations = pops,
    block_size_morgans = blgsize, max_missing_fraction = maxmiss, afprod = isTRUE(afprod),
    adjust_pseudohaploid = adjust_pseudohaploid, apply_corr = apply_corr, auto_only = auto_only,
    minmaf = minmaf, maxmaf = maxmaf, transitions = transitions, transversions = transversions,
    maxmem_mb = maxmem, admixtools_version = as.character(utils::packageVersion("admixtools")),
    R = R.version.string, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(manifest, file.path(cache_dir, paste0("f2_blocks_", signature, ".manifest.json")),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")
  blocks
}

run_f3 <- function(f2_blocks, triples) {
  require_admixtools(); triples <- as.matrix(triples)
  if (ncol(triples) != 3L) stop("Each f3 tuple must have exactly three populations: target, source1, source2.")
  validate_f2_population_labels(f2_blocks, as.vector(triples))
  do.call(rbind, lapply(seq_len(nrow(triples)), function(i) {
    z <- trimws(as.character(triples[i, ])); as.data.frame(admixtools::qp3pop(f2_blocks, pop1=z[1], pop2=z[2], pop3=z[3]))
  }))
}

run_f4 <- function(f2_blocks, quadruples) {
  require_admixtools(); quadruples <- as.matrix(quadruples)
  if (ncol(quadruples) != 4L) stop("Each f4 tuple must have exactly four populations: A, B, C, D.")
  validate_f2_population_labels(f2_blocks, as.vector(quadruples))
  as.data.frame(admixtools::qpdstat(f2_blocks, pop1=quadruples, f4mode=TRUE))
}

run_f4ratio <- function(f2_blocks, quintets, boot = FALSE) {
  require_admixtools(); quintets <- as.matrix(quintets)
  if (ncol(quintets) != 5L) stop("Each f4-ratio tuple must contain exactly five populations.")
  validate_f2_population_labels(f2_blocks, as.vector(quintets))
  # The SNP filters are already frozen in the cached f2 array; qpf4ratio
  # cannot retroactively change chromosome/MAF filters on that array.
  admixtools::qpf4ratio(f2_blocks, pops = quintets, boot = boot)
}

run_dstat <- function(prefix, quadruples, blgsize = 0.05, auto_only = TRUE,
                      allsnps = FALSE, poly_only = FALSE) {
  require_admixtools(); validate_prefix(prefix); quadruples <- as.matrix(quadruples)
  if (ncol(quadruples) != 4L) stop("Each D-statistic tuple must have exactly four populations: A, B, C, D.")
  validate_population_labels(prefix, as.vector(quadruples))
  # Supply the complete tuple matrix in one call. With allsnps=FALSE this
  # lets qpdstat use the same SNP set across all requested D-statistics.
  # f4blockdat_from_geno has a narrower filter set than f2_from_geno.
  as.data.frame(admixtools::qpdstat(prefix, pop1=quadruples, f4mode=FALSE,
    blgsize=blgsize, auto_only=auto_only, allsnps=allsnps, poly_only=poly_only))
}

run_qpwave <- function(f2_blocks, left, right, ...) {
  require_admixtools()
  if (length(left) < 2L || length(right) < 2L) stop("qpWave needs at least two left and two right populations.")
  if (length(intersect(left, right))) stop("qpWave left and right sets must not overlap.")
  validate_f2_population_labels(f2_blocks, c(left, right))
  admixtools::qpwave(f2_blocks, left = left, right = right, ...)
}

run_qpadm <- function(data, target, sources, right, ...) {
  require_admixtools(); validate_population_design(target, sources, right)
  if (is.array(data) && length(dim(data)) == 3L) validate_f2_population_labels(data, c(target, sources, right))
  # With a cached f2 array, qpadm uses those precomputed blocks. With a
  # genotype prefix, supported ... arguments are forwarded to f4blockdat_from_geno
  # (for example allsnps=TRUE). The caller must choose one data path explicitly.
  admixtools::qpadm(data, left = sources, right = right, target = target, ...)
}

run_qpadm_rotate <- function(f2_blocks, target, leftright, rightfix, full_results = TRUE, ...) {
  require_admixtools()
  if (length(leftright) < 2L || length(rightfix) < 2L) stop("qpAdm rotation requires at least two leftright candidates and two fixed right populations.")
  if (target %in% leftright || target %in% rightfix) stop("The target must be separate from leftright and rightfix populations.")
  validate_f2_population_labels(f2_blocks, c(target, leftright, rightfix))
  admixtools::qpadm_rotate(f2_blocks, leftright = leftright, target = target,
                           rightfix = rightfix, full_results = full_results, ...)
}

run_qpgraph <- function(f2_blocks, graph, ...) {
  require_admixtools()
  if (is.null(graph)) stop("Provide an ADMIXTOOLS 2 graph object in RDS format.")
  leaves <- tryCatch(admixtools::get_leafnames(graph), error = function(e) character())
  if (length(leaves)) validate_f2_population_labels(f2_blocks, leaves)
  admixtools::qpgraph(f2_blocks, graph, ...)
}

run_find_graphs <- function(f2_blocks, numadmix = 1L, outpop, stop_gen = 100L,
                            stop_gen2 = 15L, stop_sec = 3600, numgraphs = 10L,
                            stop_score = 0, max_admix = numadmix,
                            opt_worst_residual = FALSE, reject_f4z = 0,
                            return_searchtree = FALSE, ...) {
  require_admixtools()
  if (length(outpop) != 1L || !nzchar(outpop)) stop("Specify exactly one outgroup/root population.")
  validate_f2_population_labels(f2_blocks, outpop)
  if (numadmix < 0 || stop_gen < 1 || stop_gen2 < 1 || numgraphs < 1 || stop_sec < 1 || max_admix < numadmix || stop_score < 0 || reject_f4z < 0) stop("Graph search limits are invalid; max_admix must be at least numadmix and score/z thresholds cannot be negative.")
  admixtools::find_graphs(f2_blocks, numadmix = as.integer(numadmix), outpop = outpop,
                          stop_gen = as.integer(stop_gen), stop_gen2 = as.integer(stop_gen2),
                          stop_sec = stop_sec, numgraphs = as.integer(numgraphs),
                          stop_score = stop_score, max_admix = as.integer(max_admix),
                          opt_worst_residual = opt_worst_residual, reject_f4z = reject_f4z,
                          return_searchtree = return_searchtree, ...)
}

run_fst <- function(prefix, pop1, pop2, boot = FALSE, adjust_pseudohaploid = TRUE,
                    auto_only = TRUE, minmaf = 0, maxmaf = 0.5,
                    transitions = TRUE, transversions = TRUE, maxmem = 8000, ...) {
  require_admixtools(); info <- validate_prefix(prefix)
  validate_population_labels(prefix, c(pop1, pop2))
  if (!isTRUE(transitions) && !isTRUE(transversions)) stop("At least one of transitions or transversions must be included.")
  # fst() uses the Hudson estimator and block-jackknife by default; boot may
  # be TRUE or an integer number of bootstrap resamples. Genotype filters are
  # forwarded to f2_from_geno by ADMIXTOOLS 2.
  admixtools::fst(info$prefix, pop1 = pop1, pop2 = pop2, boot = boot,
    adjust_pseudohaploid = adjust_pseudohaploid, auto_only = auto_only,
    minmaf = minmaf, maxmaf = maxmaf, transitions = transitions,
    transversions = transversions, maxmem = maxmem, ...)
}

write_run_manifest <- function(path, config, dataset = list(), status = "completed", warnings = character()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  payload <- list(run_id = basename(dirname(path)), timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                  status = status, config = config, dataset = dataset, warnings = warnings,
                  R = R.version.string,
                  admixtools_version = if (requireNamespace("admixtools", quietly = TRUE)) as.character(utils::packageVersion("admixtools")) else NA_character_,
                  platform = R.version$platform)
  jsonlite::write_json(payload, path, auto_unbox = TRUE, pretty = TRUE, null = "null", dataframe = "rows")
  invisible(payload)
}

persist_analysis <- function(results_dir, method, config, result, dataset = list()) {
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", gsub("[^A-Za-z0-9]+", "-", method), "_", substr(digest::digest(list(config, Sys.time()), algo="sha1"), 1, 7))
  run_dir <- file.path(results_dir, "runs", run_id); dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(run_dir, "result.rds"))
  tab <- tryCatch(as.data.frame(result), error = function(e) NULL)
  if (!is.null(tab) && nrow(tab)) tryCatch(utils::write.csv(tab, file.path(run_dir, "result.csv"), row.names = FALSE), error = function(e) NULL)
  write_run_manifest(file.path(run_dir, "manifest.json"), c(list(method = method), config), dataset)
  list(run_id = run_id, path = run_dir, table = tab)
}

run_ghostbuster <- function(script, ghostbuster_dir, genome_build, trees_prefix, poplabels,
                            recombination_map, sample_ids, chromosomes, output_prefix,
                            num_clusters = 2L, n_repeats = 20L, num_iters = 200L,
                            seed = 2L, hmm = TRUE, only_pickle = FALSE,
                            aggregate = FALSE, python_bin = "python") {
  if (!file.exists(script) || !file.exists(poplabels)) stop("GhostBuster launcher and population-label file must exist.")
  if (!dir.exists(ghostbuster_dir) || !file.exists(file.path(ghostbuster_dir, "ghost_buster.py"))) stop("Select a GhostBuster checkout containing ghost_buster.py.")
  if (!genome_build %in% c("hg37", "hg38")) stop("Genome build must be hg37 or hg38.")
  if (!nzchar(sample_ids) || !nzchar(chromosomes) || !nzchar(python_bin) || !nzchar(output_prefix)) stop("Haplotype indices, chromosomes, Python executable and output prefix are required.")
  ids <- unlist(strsplit(trimws(sample_ids), "[[:space:],]+"))
  chr <- unlist(strsplit(gsub("[[:space:]]", "", chromosomes), ",", fixed = TRUE))
  if (!length(ids) || any(!grepl("^[0-9]+$", ids)) || anyDuplicated(ids)) stop("Relate --sample_id values must be unique non-negative integer haplotype indices separated by spaces.")
  if (!length(chr) || any(!grepl("^[0-9]+$", chr)) || anyDuplicated(chr)) stop("Chromosomes must be unique integer values separated by commas.")
  if (!isTRUE(aggregate)) {
    missing_trees <- chr[!file.exists(paste0(trees_prefix, chr, ".trees"))]
    if (length(missing_trees)) stop("Relate .trees files missing for chromosomes: ", paste(missing_trees, collapse = ", "))
  }
  if (!file.exists(recombination_map) && !dir.exists(recombination_map)) {
    map_exists <- vapply(chr, function(cc) {
      candidates <- c(paste0(recombination_map, cc), paste0(recombination_map, cc, ".txt"),
                      paste0(recombination_map, "chr", cc), paste0(recombination_map, "chr", cc, ".txt"))
      any(file.exists(candidates))
    }, logical(1))
    missing_maps <- chr[!map_exists]
    if (length(missing_maps)) stop("Recombination map files missing for chromosomes: ", paste(missing_maps, collapse = ", "))
  }
  if (num_clusters < 1 || n_repeats < 1 || num_iters < 1 || seed < 0) stop("GhostBuster EM settings are outside valid ranges.")
  args <- c(script, ghostbuster_dir, genome_build, trees_prefix, poplabels, recombination_map,
            sample_ids, chromosomes, output_prefix, as.integer(num_clusters), as.integer(n_repeats),
            as.integer(num_iters), as.integer(seed), as.character(hmm), as.character(only_pickle),
            as.character(aggregate), python_bin)
  log <- system2("bash", args, stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status") %||% 0L
  if (status != 0L) stop(paste(c("GhostBuster exited with status", status, log), collapse = "\n"))
  tree_files <- if (isTRUE(aggregate)) character() else paste0(trees_prefix, chr, ".trees")
  map_files <- if (file.exists(recombination_map)) recombination_map else if (dir.exists(recombination_map)) {
    list.files(recombination_map, full.names = TRUE)
  } else unlist(lapply(chr, function(cc) {
    candidates <- c(paste0(recombination_map, cc), paste0(recombination_map, cc, ".txt"),
                    paste0(recombination_map, "chr", cc), paste0(recombination_map, "chr", cc, ".txt"))
    candidates[file.exists(candidates)][1]
  }))
  input_files <- unique(c(poplabels, tree_files, map_files, file.path(ghostbuster_dir, "ghost_buster.py")))
  input_files <- input_files[!is.na(input_files) & file.exists(input_files) & !dir.exists(input_files)]
  python_version <- tryCatch(system2(python_bin, "--version", stdout = TRUE, stderr = TRUE), error = function(e) as.character(e$message))
  list(log = log, output_prefix = output_prefix, input_files = input_files,
       python_version = python_version)
}

run_ghostbuster_ld_dating <- function(ghostbuster_dir, output_prefix, genome_build, mode = "all", python_bin = "python") {
  script <- file.path(ghostbuster_dir, "ld_curve_dating.py")
  if (!file.exists(script)) stop("ld_curve_dating.py was not found in the GhostBuster checkout.")
  if (!nzchar(output_prefix) || !nzchar(genome_build) || !nzchar(python_bin)) stop("Output prefix, genome build and Python executable are required.")
  if (!genome_build %in% c("hg37", "hg38")) stop("Genome build must be hg37 or hg38.")
  if (!mode %in% c("all", "1date", "2date", "continuous")) stop("Unknown LD-curve dating mode.")
  args <- c(script, "--output", output_prefix, "--genome_build", genome_build)
  if (!identical(mode, "all")) args <- c(args, "--mode", mode)
  log <- system2(python_bin, args, stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status") %||% 0L
  if (status != 0L) stop(paste(c("LD-curve dating exited with status", status, log), collapse = "\n"))
  list(log = log, output_prefix = output_prefix, mode = mode,
       input_files = c(script),
       python_version = tryCatch(system2(python_bin, "--version", stdout = TRUE, stderr = TRUE), error = function(e) as.character(e$message)))
}
