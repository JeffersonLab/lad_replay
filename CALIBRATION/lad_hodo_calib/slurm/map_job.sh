#!/bin/bash
# Slurm array task: run one map macro on one chunk list. Submitted by submit_map.sh.
#SBATCH --partition=submit
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=04:00:00
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1
W=${CAL_WORK:?}
i=$(printf %03d "$SLURM_ARRAY_TASK_ID")
cd "${CAL_DIR:?}" || exit 1
echo "host=$(hostname) chunk=$i files=$(grep -c . $W/chunks/chunk_$i.dat) start=$(date)"
root -l -b -q "${CAL_MACRO:?}+(\"$W/chunks/chunk_$i.dat\",\"$W/out/tmp_chunk_$i.root\",4${CAL_ARGS:+,$CAL_ARGS})" 2>&1 | grep -av 'no dictionary'
[ -s "$W/out/tmp_chunk_$i.root" ] && mv "$W/out/tmp_chunk_$i.root" "$W/out/chunk_$i.root"
echo "end=$(date)"
