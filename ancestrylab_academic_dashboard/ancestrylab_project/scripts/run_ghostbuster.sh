#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 7 ]]; then
  echo "Usage: $0 GHOSTBUSTER_DIR GENOME_BUILD TREES_PREFIX POPLABELS RECOMB_MAP SAMPLE_IDS CHROMOSOMES [OUTPUT_PREFIX] [NUM_CLUSTERS] [N_REPEATS] [NUM_ITERS] [SEED] [HMM] [ONLY_PICKLE] [AGGREGATE] [PYTHON_BIN]" >&2
  exit 2
fi
GB_DIR=$1; BUILD=$2; TREES=$3; POPLABELS=$4; REC=$5; SAMPLE_IDS=$6; CHROMS=$7
OUT=${8:-results/ghostbuster/run}; CLUSTERS=${9:-2}; REPEATS=${10:-20}; ITERS=${11:-200}; SEED=${12:-2}
HMM=${13:-true}; ONLY_PICKLE=${14:-false}; AGGREGATE=${15:-false}; PYTHON_BIN=${16:-python}
HMM=${HMM,,}; ONLY_PICKLE=${ONLY_PICKLE,,}; AGGREGATE=${AGGREGATE,,}
[[ -f "$POPLABELS" ]] || { echo "Population-label file not found: $POPLABELS" >&2; exit 3; }
if [[ ! -f "$REC" && ! -d "$REC" ]]; then
  IFS=',' read -r -a MAP_CHRS <<< "$CHROMS"
  for chr in "${MAP_CHRS[@]}"; do
    [[ -f "${REC}${chr}" || -f "${REC}${chr}.txt" || -f "${REC}chr${chr}" || -f "${REC}chr${chr}.txt" ]] || { echo "Recombination map file not found for chromosome $chr (prefix: $REC)" >&2; exit 3; }
  done
fi
[[ -d "$GB_DIR" ]] || { echo "GhostBuster source directory not found: $GB_DIR" >&2; exit 3; }
[[ -f "$GB_DIR/ghost_buster.py" ]] || { echo "Not a GhostBuster checkout: $GB_DIR" >&2; exit 3; }
if [[ "$AGGREGATE" != "true" ]]; then
  IFS=',' read -r -a CHR_ARRAY <<< "$CHROMS"
  for chr in "${CHR_ARRAY[@]}"; do
    [[ -f "${TREES}${chr}.trees" ]] || { echo "Relate tree file not found: ${TREES}${chr}.trees" >&2; exit 3; }
  done
fi
[[ "$CLUSTERS" =~ ^[1-9][0-9]*$ ]] || { echo "NUM_CLUSTERS must be a positive integer" >&2; exit 4; }
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || { echo "N_REPEATS must be a positive integer" >&2; exit 4; }
[[ "$ITERS" =~ ^[1-9][0-9]*$ ]] || { echo "NUM_ITERS must be a positive integer" >&2; exit 4; }
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "SEED must be a non-negative integer" >&2; exit 4; }
mkdir -p "$(dirname "$OUT")"
read -r -a IDS <<< "$SAMPLE_IDS"
CMD=("$PYTHON_BIN" "$GB_DIR/ghost_buster.py" --genome_build "$BUILD" \
  --poplabels "$POPLABELS" --rec "$REC" --sample_id "${IDS[@]}" \
  --chrs "$CHROMS" --output "$OUT" --num_clusters "$CLUSTERS" \
  --n_repeats "$REPEATS" --num_iters "$ITERS" --seed "$SEED")
if [[ "$AGGREGATE" != "true" ]]; then CMD+=(--trees "$TREES"); fi
if [[ "$HMM" == "false" ]]; then CMD+=(--hmm False); fi
if [[ "$ONLY_PICKLE" == "true" ]]; then CMD+=(--only_make_pickle_files); fi
printf 'Executing:'; printf ' %q' "${CMD[@]}"; printf '\n'
"${CMD[@]}"
