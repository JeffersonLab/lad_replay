#!/bin/bash
# Merge ("reduce") step for one sigma: sum the per-chunk caches with hadd, rewrite
# the merged cache's signature to the FULL runlist (fix_sig.C), then run each
# macro on the full runlist with that cache -> cache HIT, the event loop is
# skipped and only fits/plots are made. Full runlist: ~5 GB peak, ~1 min per macro.
# Normally submitted by submit_all.sh. Can also run directly on the login node:
#   SIGMA=<sigma> bash slurm/merge_job.sh
#SBATCH --job-name=ladtof_merge
#SBATCH --partition=submit
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00

# In a batch job LADTOF_DIR comes from the submitting shell; when run directly,
# fall back to this script's own location.
source "${LADTOF_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/config.sh" || exit 1
sig=${SIGMA:?set SIGMA}
W=$WORK_BASE/$sig
full=$(runlist_for "$sig")
nchunks=$(ls "$W"/chunks/chunk_*.dat | wc -l)

cd "$LADTOF_DIR" || exit 1
setup_root
rc=0
for m in "${MACROS[@]}"; do
  echo "===== $m ($sig)"
  caches=("$W"/cache/${m}_chunk_*.root)
  if [ ${#caches[@]} -ne "$nchunks" ]; then
    echo "ERROR: ${#caches[@]} caches for $nchunks chunks; resubmit the missing map tasks first"; rc=1; continue
  fi
  merged=$W/cache/${m}_merged.root
  hadd -f -k "$merged" "${caches[@]}" >/dev/null || { echo "ERROR: hadd failed"; rc=1; continue; }
  root -l -b -q "$LADTOF_DIR/slurm/fix_sig.C+(\"$merged\",\"$full\")" 2>&1 | grep fix_sig
  out=${OUTBASE[$m]}_${sig}.root
  mkdir -p "$(dirname "$out")"
  log=$(root -l -b -q "${m}.C+(\"$full\",\"$out\",1,\"$merged\")" 2>&1 | grep -av 'no dictionary for class')
  echo "$log" | grep -aE 'cache (HIT|MISS)|rror'
  if echo "$log" | grep -q 'cache HIT'; then echo "wrote $LADTOF_DIR/$out"; else echo "ERROR: expected a cache HIT"; rc=1; fi
done
exit $rc
