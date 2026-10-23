#!/bin/bash
# Merge a runlist from the chunk caches of other runlists ("pieces") instead of mapping it again. Used when two
# runlists share most files: map the shared part once, map each remainder, then combine. The target runlist
# (all_C3_runlist_SHMS_13p5_submit_<target>.dat) must be exactly the union of the pieces' files.
#   SIGMA=<target> PIECES="<piece> <piece> ..." sbatch --export=ALL ... slurm/combine_merge.sh
# Links every piece's chunk caches into $WORK_BASE/<target>/cache (renumbered) with matching dummy chunk lists,
# then runs merge_job.sh for <target>: hadd, fix_sig to the target runlist, and a cache-HIT run of each macro.
#SBATCH --job-name=ladtof_combine
#SBATCH --partition=submit
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00

source "${LADTOF_DIR:?export LADTOF_DIR (source slurm/config.sh)}/slurm/config.sh" || exit 1
tgt=${SIGMA:?set SIGMA to the target runlist tag}
pieces=(${PIECES:?set PIECES})
T=$WORK_BASE/$tgt
# Only the selected macros' caches are replaced: other macros' merged caches in the target stay usable for redraws.
rm -rf "$T/chunks"
mkdir -p "$T/chunks" "$T/cache"
for m in "${MACROS[@]}"; do rm -f "$T/cache/${m}"_chunk_*.root "$T/cache/${m}_merged.root"; done
n=0
for p in "${pieces[@]}"; do
  for c in "$WORK_BASE/$p"/chunks/chunk_*.dat; do
    i=$(basename "$c" .dat | sed 's/chunk_//')
    k=$(printf %03d $n)
    cp "$c" "$T/chunks/chunk_$k.dat"
    for m in "${MACROS[@]}"; do
      src=$WORK_BASE/$p/cache/${m}_chunk_$i.root
      [ -s "$src" ] || { echo "ERROR: missing $src"; exit 1; }
      ln -s "$src" "$T/cache/${m}_chunk_$k.root"
    done
    n=$((n + 1))
  done
done
echo "$tgt: $n chunks from ${pieces[*]}"
# the target runlist must be the union of the pieces (same files, any order)
diff <(cat "$T"/chunks/chunk_*.dat | sort) <(grep -vE '^\s*(#|$)' "$(runlist_for "$tgt")" | sed -E 's/^\s+|\s+$//g' | sort) >/dev/null ||
  { echo "ERROR: pieces do not cover $(runlist_for "$tgt") exactly"; exit 1; }
bash "$LADTOF_DIR/slurm/merge_job.sh"
