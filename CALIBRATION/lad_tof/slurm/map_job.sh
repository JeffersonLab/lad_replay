#!/bin/bash
# Split ("map") step: one array task per chunk runlist. Runs the selected macros
# (LADTOF_MACROS, default all three; see config.sh) on the chunk with a histogram
# cache; the caches are what the merge step sums.
# The per-chunk plot outputs are by-products and can be deleted.
# Normally submitted by submit_all.sh, which exports LADTOF_DIR/WORK_BASE and sets
# the log path (--output). By hand, from CALIBRATION/lad_tof:
#   source slurm/config.sh
#   sbatch --array=0-<N-1> --output=$WORK_BASE/logs/map_%x_%A_%a.log --export=ALL,SIGMA=<sigma> slurm/map_job.sh
# Tasks skip a macro whose cache already exists, so failed tasks can simply be resubmitted.
#SBATCH --job-name=ladtof_map
#SBATCH --partition=submit
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00

# Slurm runs a spooled copy of this script, so config.sh is found via LADTOF_DIR
# (exported at submission), not via this script's own path.
source "${LADTOF_DIR:?export LADTOF_DIR (source slurm/config.sh) before sbatch}/slurm/config.sh" || exit 1
sig=${SIGMA:?set SIGMA via --export}
i=$(printf %03d "$SLURM_ARRAY_TASK_ID")
W=$WORK_BASE/$sig
list=$W/chunks/chunk_$i.dat
# Slurm may allocate more CPUs than requested (12-20 seen). RDataFrame keeps one copy
# of every histogram per thread, so lad_tracking_eff's memory grows with the thread
# count (~18 GB at 12 threads with 4 chi2 cuts; >40 GB at 20 threads with 5). Cap it;
# the macros keep only ~4-8 cores busy anyway. Override with LADTOF_THREADS.
nthreads=${SLURM_CPUS_PER_TASK:-4}
max_threads=${LADTOF_THREADS:-8}
[ "$nthreads" -gt "$max_threads" ] && nthreads=$max_threads
mkdir -p "$W/cache" "$W/chunk_out"

cd "$LADTOF_DIR" || exit 1
setup_root
echo "host=$(hostname) sigma=$sig chunk=$i threads=$nthreads files=$(grep -c . "$list") start=$(date)"
rc=0
for m in "${MACROS[@]}"; do
  cache=$W/cache/${m}_chunk_$i.root
  if [ -s "$cache" ]; then echo "[$m] cache exists, skipping"; continue; fi
  tmp=$W/cache/tmp_${m}_chunk_$i.root   # renamed on success so a crash never leaves a partial cache
  rm -f "$tmp"
  /usr/bin/time -f "TIME $m wall=%e s maxrss=%M kB cpu=%P" \
    root -l -b -q "${m}.C+(\"$list\",\"$W/chunk_out/${m}_chunk_$i.root\",$nthreads,\"$tmp\")" 2>&1 |
    grep -av 'no dictionary for class' | tr '\r' '\n' | grep -avE '^\s*$|^\|=|%\]'
  if [ -s "$tmp" ]; then # the macros write the cache only after a completed event loop
    mv "$tmp" "$cache"
  else
    echo "[$m] ERROR: no cache written"; rc=1
  fi
done
echo "end=$(date) rc=$rc"
exit $rc
