#!/bin/bash
# Run lad_photon_timing.C over a set of replay files on Slurm, then hadd.
#   usage: slurm/submit_photon_timing.sh <runlist|dir> <out.root> [files_per_chunk=10]
# Must be run from CALIBRATION/lad_tof. With PT_GOODVTX=true in the environment,
# events without a reaction point inside |z| < 20 cm are dropped. Compiles the macro once first so the
# array tasks don't race on ACLiC.
in=$1; out=$2; nper=${3:-10}
[ -n "$in" ] && [ -n "$out" ] || { echo "usage: $0 <runlist|dir> <out.root> [files_per_chunk]" >&2; exit 1; }
export LADTOF_DIR=$(pwd)
tag=$(basename "$out" .root)
export PT_WORK=/ceph/submit/data/user/${USER:0:1}/$USER/hallc/lad_tof_scratch/photon_timing/$tag
mkdir -p $PT_WORK/chunks $PT_WORK/out $PT_WORK/logs
rm -f $PT_WORK/chunks/*.dat $PT_WORK/out/*.root
if [ -d "$in" ]; then ls "$in"/*.root; else grep -vE '^\s*(#|$)' "$in"; fi |
  awk -v n=$nper -v d=$PT_WORK/chunks '{printf "%s\n",$0 > sprintf("%s/chunk_%03d.dat",d,int((NR-1)/n))}'
nchunk=$(ls $PT_WORK/chunks | wc -l)
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1
root -l -b -q -e '.L lad_photon_timing.C+' >/dev/null 2>&1 || { echo "compile failed" >&2; exit 1; }
jid=$(sbatch --parsable --array=0-$((nchunk-1)) --output=$PT_WORK/logs/map_%A_%a.log --export=ALL slurm/photon_timing_job.sh)
mid=$(sbatch --parsable --dependency=afterany:$jid --job-name=lad_photon_merge --partition=submit --mem=4G --time=01:00:00 \
  --output=$PT_WORK/logs/merge_%j.log --wrap="source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1; ls $PT_WORK/out/chunk_*.root | wc -l; hadd -f $LADTOF_DIR/$out $PT_WORK/out/chunk_*.root")
echo "chunks=$nchunk map=$jid merge=$mid work=$PT_WORK"
