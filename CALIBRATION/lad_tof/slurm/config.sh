#!/bin/bash
# Shared settings for the subMIT Slurm split/merge workflow (see ../README_subMIT_slurm.md).
# Sourced by every script in this folder. Nothing here is user-specific: paths
# follow this file's location and $USER, and both can be overridden from the
# environment (submit_all.sh exports them to the batch jobs).

# This lad_replay checkout's CALIBRATION/lad_tof (the folder above slurm/).
LADTOF_DIR=${LADTOF_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
RUNLIST_DIR=$LADTOF_DIR/../files/run-lists
# Scratch area for chunk runlists, per-chunk caches/outputs and logs: the
# user's own Ceph area (never /home or /work).
WORK_BASE=${WORK_BASE:-/ceph/submit/data/user/${USER:0:1}/$USER/hallc/lad_tof_scratch/production}
ROOT_SETUP=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
export LADTOF_DIR WORK_BASE

# Macro name -> final output file (relative to $LADTOF_DIR); "_<sigma>.root" is appended.
ALL_MACROS=(lad_tracking_eff lad_hodo_eff lad_hodo_dist proton_tof_plot lad_gem_resid)
declare -A OUTBASE=(
  [lad_tracking_eff]=files/tracking_eff/tracking_eff_C3_SHMS_13p5_v5_PH
  [lad_hodo_eff]=files/hodo_eff/hodo_eff_C3_SHMS_13p5_v4_PH
  [lad_hodo_dist]=files/hodo_dist/hodo_dist_C3_SHMS_13p5_v4_PH
  [proton_tof_plot]=files/proton_tof/proton_tof_C3_SHMS_13p5_PH
  [lad_gem_resid]=files/gem_resid/gem_resid_C3_SHMS_13p5_P
)
# Run when no selection is given. proton_tof_plot is opt-in (its canvases are
# also produced by lad_tracking_eff), and so is lad_gem_resid (GEM residual study).
DEFAULT_MACROS=(lad_tracking_eff lad_hodo_eff lad_hodo_dist)

# Macros to run: LADTOF_MACROS (comma- or space-separated) selects a subset;
# unset = DEFAULT_MACROS. submit_all.sh -m sets it and exports it to the jobs.
MACROS=(${LADTOF_MACROS//,/ })
[ ${#MACROS[@]} -eq 0 ] && MACROS=("${DEFAULT_MACROS[@]}")
for m in "${MACROS[@]}"; do
  # return (not exit) when sourced interactively, so a typo doesn't close the user's shell
  [ -n "${OUTBASE[$m]}" ] || { echo "unknown macro '$m' (choose from: ${ALL_MACROS[*]})" >&2; return 1 2>/dev/null || exit 1; }
done

# Full runlist for a sigma tag (5sigma / 10sigma).
runlist_for() { echo "$RUNLIST_DIR/all_C3_runlist_SHMS_13p5_submit_$1.dat"; }

setup_root() { source "$ROOT_SETUP" >/dev/null 2>&1; }
