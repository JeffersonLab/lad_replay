#!/bin/bash
# Slurm array task for lad_photon_timing.C: one task per chunk list.
# Submit with slurm/submit_photon_timing.sh (it makes the chunks and the merge job).
#SBATCH --job-name=lad_photon_timing
#SBATCH --partition=submit
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=02:00:00
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1
W=${PT_WORK:?}
i=$(printf %03d "$SLURM_ARRAY_TASK_ID")
cd "${LADTOF_DIR:?}" || exit 1
echo "host=$(hostname) chunk=$i files=$(grep -c . $W/chunks/chunk_$i.dat) start=$(date)"
root -l -b -q "lad_photon_timing.C+(\"$W/chunks/chunk_$i.dat\",\"$W/out/tmp_chunk_$i.root\",4,${PT_GOODVTX:-false})" 2>&1 | grep -av 'no dictionary'
[ -s "$W/out/tmp_chunk_$i.root" ] && mv "$W/out/tmp_chunk_$i.root" "$W/out/chunk_$i.root"
echo "end=$(date)"
