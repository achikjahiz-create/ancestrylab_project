source("R/engine.R")

# Design validation: catch structurally invalid qpAdm models before computation.
stopifnot(isTRUE(tryCatch({ validate_population_design("T", c("A", "B"), c("R1", "R2")); TRUE }, error=function(e) FALSE)))
stopifnot(inherits(tryCatch(validate_population_design("T", c("A", "A"), c("R1", "R2")), error=function(e) e), "error"))
stopifnot(inherits(tryCatch(validate_population_design("T", c("A"), c("R1", "R2")), error=function(e) e), "error"))
stopifnot(inherits(tryCatch(validate_population_design("T", c("A", "B"), c("A", "R2")), error=function(e) e), "error"))
stopifnot(inherits(tryCatch(validate_population_design("T", c("A", "B"), c("T", "R2")), error=function(e) e), "error"))

# Exact population-label validation: a nearby or misspelled label must never pass.
tmp <- tempfile("ancestrylab_test_"); dir.create(tmp)
prefix <- file.path(tmp, "synthetic")
file.create(paste0(prefix, ".geno")); file.create(paste0(prefix, ".snp"))
writeLines(c("sample_1 U Ancient_A", "sample_2 U Ancient_B", "sample_3 U Modern_C"), paste0(prefix, ".ind"))
stopifnot(setequal(read_population_labels(prefix)$labels, c("Ancient_A", "Ancient_B", "Modern_C")))
stopifnot(isTRUE(tryCatch({ validate_population_labels(prefix, c("Ancient_A", "Modern_C")); TRUE }, error=function(e) FALSE)))
stopifnot(inherits(tryCatch(validate_population_labels(prefix, "Ancient_A_proxy"), error=function(e) e), "error"))
stopifnot(inherits(tryCatch(validate_population_labels(prefix, "Ancient_A "), error=function(e) e), "error") == FALSE) # whitespace is normalized
unlink(tmp, recursive=TRUE)

# The current f2 cache must contain the exact labels requested by a model.
f2_test <- array(0, dim=c(3,3,1), dimnames=list(c("Ancient_A","Ancient_B","Modern_C"),
                                                    c("Ancient_A","Ancient_B","Modern_C"), "block1"))
stopifnot(isTRUE(tryCatch({ validate_f2_population_labels(f2_test,c("Ancient_A","Modern_C")); TRUE }, error=function(e) FALSE)))
stopifnot(inherits(tryCatch(validate_f2_population_labels(f2_test,"NotInCache"), error=function(e) e), "error"))

unlink(tmp, recursive=TRUE)

cat("AncestryLab engine smoke tests passed.\n")
