#!/bin/bash
# Run a map macro over replay files on Slurm (one task per chunk), then hadd.
#   usage: slurm/submit_map.sh <macro.C> <runlist|dir> <out.root> [files_per_chunk=5] [extra macro args]
# Run from the directory that holds the macro. The macro signature must be (list, out, nthreads[, extra...]).
# extra args are passed verbatim, e.g. '"P"'. Compiles once first so the array tasks don't race on ACLiC.
macro=$1; in=$2; out=$3; nper=${4:-5}; extra=$5
[ -n "$macro" ] && [ -n "$in" ] && [ -n "$out" ] || { echo "usage: $0 <macro.C> <runlist|dir> <out.root> [files_per_chunk] [args]" >&2; exit 1; }
export CAL_DIR=$(pwd) CAL_MACRO=$macro CAL_ARGS=$extra
tag=$(basename "$out" .root)
export CAL_WORK=/ceph/submit/data/user/${USER:0:1}/$USER/hallc/lad_calib_2026/work/$tag
mkdir -p $CAL_WORK/chunks $CAL_WORK/out $CAL_WORK/logs
rm -f $CAL_WORK/chunks/*.dat $CAL_WORK/out/*.root
if [ -d "$in" ]; then ls "$in"/*.root; else grep -vE '^\s*(#|$)' "$in"; fi |
  awk -v n=$nper -v d=$CAL_WORK/chunks '{printf "%s\n",$0 > sprintf("%s/chunk_%03d.dat",d,int((NR-1)/n))}'
nchunk=$(ls $CAL_WORK/chunks | wc -l)
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1
root -l -b -q -e ".L $macro+" 2>&1 | grep -E 'error|Error' && { echo "compile failed" >&2; exit 1; }
jid=$(sbatch --parsable --job-name=${tag} --array=0-$((nchunk-1)) --output=$CAL_WORK/logs/map_%A_%a.log --export=ALL $(dirname $(readlink -f $0))/map_job.sh)
mid=$(sbatch --parsable --dependency=afterany:$jid --job-name=${tag}_merge --partition=submit --mem=4G --time=01:00:00 \
  --output=$CAL_WORK/logs/merge_%j.log --wrap="source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh >/dev/null 2>&1; ls $CAL_WORK/out/chunk_*.root | wc -l; hadd -f $CAL_DIR/$out $CAL_WORK/out/chunk_*.root")
echo "chunks=$nchunk map=$jid merge=$mid work=$CAL_WORK"
