#!/bin/bash
# Build chunk runlists and submit map + dependent merge jobs for each sigma.
#   usage: submit_all.sh [-m macro[,macro...]] [chunk_gb=20] [sigma ...]
#     -m      subset of lad_tracking_eff,lad_hodo_eff,lad_hodo_dist,proton_tof_plot
#             (default: the first three)
#     sigma   default: 5sigma_new 10sigma_new (the runlist tag: all_C3_runlist_SHMS_13p5_submit_<tag>.dat)
#   e.g.  submit_all.sh -m lad_hodo_eff,lad_hodo_dist 20 10sigma
while getopts "m:" opt; do
  case $opt in
    m) export LADTOF_MACROS=$OPTARG ;;
    *) echo "usage: $0 [-m macro[,macro...]] [chunk_gb] [sigma ...]" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))
export LADTOF_MACROS=${LADTOF_MACROS:-}
source "$(dirname "$0")/config.sh" || exit 1
chunk_gb=${1:-20}
[ $# -gt 0 ] && shift
if [ $# -gt 0 ]; then sigmas=("$@"); else sigmas=(5sigma_new 10sigma_new); fi

# lad_tracking_eff needs ~18 GB with 3 chi2 cuts and >25 GB with 4; the hodo macros peak at ~2.5 GB and
# proton_tof_plot at ~1 GB, so runs without lad_tracking_eff ask for much less
# memory and start sooner.
mem=4G
[[ " ${MACROS[*]} " == *" lad_tracking_eff "* ]] && mem=${LADTOF_MEM:-40G}
echo "macros: ${MACROS[*]}  (map --mem=$mem)"

mkdir -p "$WORK_BASE/logs"
cd "$LADTOF_DIR" || exit 1
setup_root
# Build the macro libraries once here so concurrent jobs never race ACLiC.
for m in "${MACROS[@]}" slurm/fix_sig; do
  root -l -b -q -e ".L ${m}.C+" >/dev/null 2>&1 && [ -f "${m}_C.so" ] || { echo "compile of $m failed"; exit 1; }
done

for sig in "${sigmas[@]}"; do
  # A fresh submission recomputes the selected macros: drop their old caches
  # (map_job.sh skips macros whose chunk cache exists, which is only meant for
  # resubmitting failed tasks of the SAME submission).
  for m in "${MACROS[@]}"; do rm -f "$WORK_BASE/$sig/cache/${m}"_chunk_*.root "$WORK_BASE/$sig/cache/${m}_merged.root"; done
  n=$(bash "$LADTOF_DIR/slurm/make_chunks.sh" "$sig" "$chunk_gb") || exit 1
  # --export=ALL passes LADTOF_DIR, WORK_BASE and LADTOF_MACROS to the jobs.
  map=$(sbatch --parsable --mem=$mem --array=0-$((n - 1)) --output="$WORK_BASE/logs/map_%x_%A_%a.log" \
    --export=ALL,SIGMA=$sig "$LADTOF_DIR/slurm/map_job.sh") || exit 1
  merge=$(sbatch --parsable --dependency=afterok:$map --output="$WORK_BASE/logs/merge_%x_%j.log" \
    --export=ALL,SIGMA=$sig "$LADTOF_DIR/slurm/merge_job.sh") || exit 1
  echo "$sig: $n chunks -> map job $map, merge job $merge"
done
