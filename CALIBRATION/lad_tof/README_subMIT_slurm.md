# Running lad_tracking_eff / lad_hodo_eff / lad_hodo_dist / proton_tof_plot on subMIT (Slurm)

> **Status (October 2026).** The current replays are the ones made with the 1D-cluster tracking
> road, in `ROOTfiles/LAD_COIN/PRODUCTION/{5sigma_new,10sigma_new}/`, with runlists
> `all_C3_runlist_SHMS_13p5_submit_{5sigma_new,10sigma_new}.dat` (214 files each). The older
> `5sigma/` and `10sigma/` replays below have been deleted. Wherever this README says `5sigma` or
> `10sigma` (runlist names, `$WORK_BASE/<sigma>`, output names), use `5sigma_new` / `10sigma_new`;
> `slurm/submit_all.sh` now uses those by default. The workflow files and runlists described below
> as "not on GitHub" are now in the repository.

This README explains how `lad_tracking_eff.C`, `lad_hodo_eff.C`, `lad_hodo_dist.C` and
`proton_tof_plot.C` are run on the subMIT cluster as Slurm batch jobs. It covers the 5σ and 10σ
LAD_COIN replays (GEM strip zero suppression at 5σ and 10σ) in:

```
/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/5sigma
/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/10sigma
```

Those ROOT files are readable by every subMIT user, so **anyone can reproduce the outputs without
replaying anything**. The **"Setting up as a new user"** section below lists everything to clone,
copy and link. The original three macros are unchanged; everything else is a new runlist, a
wrapper script in `slurm/`, or `proton_tof_plot.C`.

## Summary

- **Runlists.** The ifarm runlist `all_C3_runlist_SHMS_13p5.dat` points at `/volatile/...`.
  Two copies with the paths rewritten to the subMIT 5σ/10σ directories were added:
  `../files/run-lists/all_C3_runlist_SHMS_13p5_submit_{5sigma,10sigma}.dat`. Each has 215 files,
  793 GB (5σ) and 544 GB (10σ).
- **Split, merge, plot.** Each runlist is split into chunks of about 20 GB. One Slurm job per chunk
  runs the macros with their built-in *histogram cache* (4th macro argument). The per-chunk caches
  are then summed with `hadd`. The merged cache is relabelled for the full runlist, and each macro is
  run once on the full runlist with that cache. The macro sees a cache HIT, skips the event loop,
  and only fits and draws.
- **Time.** A full run of all three default macros on both sigmas is 75 jobs of ~5 min each, plus
  queue time. That took ~40 min in total on 2026-09-28.

## Setting up as a new user

This works for any subMIT user. The example user is `aviarm`, but every command uses `$USER`. It
assumes a subMIT account with Slurm access and a Ceph area (sections 1 and 2.0 of
`/home/submit/ehingerl/LAD_submit_getting_started.md`).

**You do not need hcana, LADlib or a replay.** The macros only need ROOT, taken from CVMFS. The
data are Lucas's ROOT files above, read in place.

### What you need, at a glance

| What | Where it comes from | Why |
|---|---|---|
| `lad_replay`, branch `lad-1Dcluster-tracking` of `LucasEhinger/lad_replay` | `git clone` | `lad_tracking_eff.C`, `lad_hodo_eff.C`, `lad_hodo_dist.C` and the base runlists are **only on this branch**, not on JeffersonLab `master` |
| `CALIBRATION/lad_tof/slurm/` (5 `.sh` + 2 `.C`) | copy from Lucas's checkout | batch workflow; not on GitHub |
| `CALIBRATION/lad_tof/proton_tof_plot.C` | copy | standalone proton-tof canvas; not on GitHub |
| `CALIBRATION/lad_tof/lad_gem_resid.C` | copy (optional) | GEM residual / hit-efficiency study; not on GitHub |
| `CALIBRATION/lad_tof/README_subMIT_slurm.md` | copy | this file |
| `CALIBRATION/files/run-lists/all_C3_runlist_SHMS_13p5_submit_{5sigma,10sigma}.dat` | copy | subMIT paths of the 5σ/10σ ROOT files; not on GitHub |
| ROOT 6.30.02 | CVMFS LCG_105 view | the scripts source it themselves |
| Symlinks | **none needed** | see "Soft links" below |

### 1. Get the code

The scripts must go in `/work`, not `/home` (10 GB quota) or `/ceph`. Use a fresh clone if you
don't have `lad_replay` yet:

```bash
mkdir -p /work/submit/$USER/hallc/software
cd /work/submit/$USER/hallc/software
git clone --branch lad-1Dcluster-tracking https://github.com/LucasEhinger/lad_replay.git
```

If you already have a `lad_replay` checkout (e.g. from the getting-started guide), add Lucas's fork
and switch to the branch instead. Commit or stash your own changes first:

```bash
cd /work/submit/$USER/hallc/software/lad_replay
git remote add lucas https://github.com/LucasEhinger/lad_replay.git
git fetch lucas
git switch -c lad-1Dcluster-tracking lucas/lad-1Dcluster-tracking
```

Submodules (`UTIL_OL`) are only needed for replaying, not for this.

### 2. Copy the files that aren't on GitHub

These live only in Lucas's checkout, which is readable by everyone:

```bash
SRC=/work/submit/ehingerl/hallc/software/lad_replay
DST=/work/submit/$USER/hallc/software/lad_replay

cp $SRC/CALIBRATION/lad_tof/README_subMIT_slurm.md \
   $SRC/CALIBRATION/lad_tof/proton_tof_plot.C       $DST/CALIBRATION/lad_tof/
mkdir -p $DST/CALIBRATION/lad_tof/slurm
cp $SRC/CALIBRATION/lad_tof/slurm/*.sh \
   $SRC/CALIBRATION/lad_tof/slurm/*.C               $DST/CALIBRATION/lad_tof/slurm/
cp $SRC/CALIBRATION/files/run-lists/all_C3_runlist_SHMS_13p5_submit_{5,10}sigma.dat \
                                                    $DST/CALIBRATION/files/run-lists/
```

Don't copy these:

- `*_C.so`, `*_C.d`, `*_ACLiC_dict_rdict.pcm`: compiled macros that ROOT rebuilds in your
  checkout.
- `CALIBRATION/lad_tof/files/*/*.root`: Lucas's outputs. You make your own.

After copying, `git status` in `$DST` shows exactly these new files as untracked. Nothing tracked
changes.

### 3. Soft links

**No new symlinks are needed.** In detail:

- **Input data.** The two `_submit_*sigma.dat` runlists contain absolute paths to Lucas's ROOT
  files. They're read in place and never go through a `ROOTfiles` link. If those files ever move
  (e.g. to `/ceph/submit/data/group/lad`), regenerate the runlists with step 1 of the full
  procedure below, using the new directory.
- **Scratch.** Chunk runlists, caches and logs go to
  `$WORK_BASE = /ceph/submit/data/user/<first letter of $USER>/$USER/hallc/lad_tof_scratch/production`.
  `slurm/config.sh` works that out from `$USER`, and the scripts create it. It never touches Lucas's
  area.
- **Outputs.** The final ROOT files (~150 MB for everything) go to your own checkout's
  `CALIBRATION/lad_tof/files/<macro>/`. That's in `/work`, which is fine at this size.
- **The `ROOTfiles` / `REPORT_OUTPUT` / `raw` links** from the getting-started guide (section 2.6)
  are for **replaying** and aren't used here. It's fine to have them or not. Don't point
  `CALIBRATION/lad_tof/files` or your `ROOTfiles` at Lucas's directories, because you can't write
  there.

### 4. Check the setup (login node, seconds)

```bash
cd /work/submit/$USER/hallc/software/lad_replay/CALIBRATION/lad_tof
source slurm/config.sh
echo $LADTOF_DIR          # -> /work/submit/<you>/hallc/software/lad_replay/CALIBRATION/lad_tof
echo $WORK_BASE           # -> /ceph/submit/data/user/<x>/<you>/hallc/lad_tof_scratch/production
head -1 ../files/run-lists/all_C3_runlist_SHMS_13p5_submit_10sigma.dat | xargs ls -l   # readable?
sacctmgr -n show assoc user=$USER format=user,account      # you have a Slurm account
```

### 5. Small test run first (recommended, ~5 min)

This runs 2 tiny files through the full split → merge → plot chain, with the two light macros and
4 GB jobs:

```bash
B=/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/10sigma
printf '%s\n' $B/LAD_COIN_22569_0_0_-1.root $B/LAD_COIN_23470_1_1_-1.root \
  > ../files/run-lists/all_C3_runlist_SHMS_13p5_submit_test.dat
bash slurm/submit_all.sh -m lad_hodo_eff,proton_tof_plot 0 test    # chunk size 0 = one file per chunk
squeue -u $USER                                                     # wait until empty
grep -h 'wrote\|ERROR' $WORK_BASE/logs/merge_*.log
ls files/hodo_eff/*_test.root files/proton_tof/*_test.root
```

A "sigma" argument is just the `<tag>` in `all_C3_runlist_SHMS_13p5_submit_<tag>.dat`, so any
runlist named that way can be run. Afterwards, delete the test runlist, the `*_test.root` outputs
and `$WORK_BASE/test`.

### 6. Full run

Continue with "Step-by-step" below, from step 2. The result is your own equivalent of Lucas's
output files. They should be identical: the histograms don't depend on the chunking, and the final
fits run single-threaded in the merge step. To check, compare against Lucas's copy (a few seconds;
it prints `differing(>1e-9 rel)=0` if they're identical):

```bash
L=/work/submit/ehingerl/hallc/software/lad_replay/CALIBRATION/lad_tof
root -l -b -q 'slurm/compare.C+("files/hodo_eff/hodo_eff_C3_SHMS_13p5_v4_PH_10sigma.root","'$L'/files/hodo_eff/hodo_eff_C3_SHMS_13p5_v4_PH_10sigma.root")'
```

A difference means the macros, the runlists or the data differ between the two checkouts. For
example, Lucas may have edited a macro since you copied it.

## LAD time-of-flight convention

The ToF offset is set per spectrometer in `PARAM/LAD/LADKINE/lladkine.param`:
`lglobal_time_offset_shms = -1728.0` and `lglobal_time_offset_hms = -1752.3` ns. With these, photons
arrive at tof = L/c (tof − L/c = 0) for both P and H. Older replays used one shared −1710 ns, which
put the photon peak at 18 ns (P) and 42 ns (H).

The macros don't use `goodhit_hit_tof_*` as stored. `lad_tof_offset.h` recomputes it for each hit as
hittime − t_vertex + offset, using the offsets above (RF-corrected ToF is shifted by the same amount).
Old and new replays therefore give identical histograms. All windows are written in this
convention: they are the former SHMS windows shifted by −18 ns, and they now apply to H as well.

| region (tof − L/c, ns) | now | before (old −1710 offset, SHMS) |
|---|---|---|
| histogram range | [−168, 157] | [−150, 175] |
| proton peak | [12, 32] | [30, 50] |
| in-time sidebands | [−43, 12] ∪ [32, 107] | [−25, 30] ∪ [50, 125] |
| out-of-time sidebands | [−168, −118] ∪ [107, 157] | [−150, −100] ∪ [125, 175] |

Vertex z: `t_vertex` now uses the vertex z in the electron path length (z·cosθ_e shorter), and the
RF-corrected ToF takes its z/c from the same z. Events whose reaction point is outside |z| < 20 cm
(`lvertex_zmax`) get no LAD ToF (`lvertex_bad_z_reject = 1`). This covers ~9% of SHMS LAD hits, and
rejecting them narrows the photon peak. The z used is written as `X.ladkin.z_tof`. For replays made
before this change, `lad_tof_offset.h` applies the same correction itself, including which RF bucket
is picked, so old and new replays still agree.

If the offsets in `lladkine.param` change, update `LAD_TOF_OFFSET_SHMS/HMS` in `lad_tof_offset.h` to
match. Do the same for `VERTEX_ZMAX` / `REJECT_BAD_Z` if `lvertex_zmax` / `lvertex_bad_z_reject` change. The offsets come from the photon peak: `lad_photon_timing.C` (fill, via
`slurm/submit_photon_timing.sh`) and `lad_photon_timing_fit.py` (fit and suggested offsets). The cache
signatures include the offsets, so caches made under the old convention are rebuilt automatically.

## Why it's done this way

These numbers were measured on subMIT (ROOT 6.30.02, LCG_105):

| | peak memory | time |
|---|---|---|
| `lad_tracking_eff`, event loop | ~17–18 GB with 3 χ² cuts, >25 GB with 4, at any thread count (map jobs request 40G; override with `LADTOF_MEM`) | fixed ~2–3 min, plus ~15–30 s per GB |
| `lad_hodo_eff` / `lad_hodo_dist`, event loop | ~1.8 GB | ~35–55 s per 7 GB |
| `proton_tof_plot`, event loop | ~0.85 GB | ~40–65 s per 7 GB |
| any macro, plotting from a cache | < 1 GB (4 files) – ~5 GB (full 215-file runlist) | ~15–60 s |

- **Not on the login node.** Each user on a login node (`submit04` etc.) is capped at **12 GB of
  RAM and 20 CPUs in total**, and that includes the VS Code server. `lad_tracking_eff` needs ~17 GB
  even with 1 thread. Running it there OOMs the whole session and disconnects VS Code.
- **Many small jobs, not a few large ones.** The macros scale poorly with threads. A 16-thread run
  averaged only ~4 busy cores. Chunks of ~20 GB with ~4 threads use the cluster far better.
- **Skimming doesn't help.** A lossless skim, keeping only the branches the macros read and events
  with a vertex or proton hit, keeps 87% of events and ~3000 of 4556 branches. It took longer
  (5 min/GB, single-threaded) than the analysis itself. The macros are limited by CPU, not I/O.

## Validation

Validation was done on 4 files (6.9 GB, 10σ). The files are in
`/ceph/submit/data/user/e/ehingerl/hallc/lad_tof_scratch/validate/`.

- **Splitting is exact.** The merged cache of 2 chunks of 2 files matches the cache of a direct
  4-file run bin for bin, for all three macros (7948 / 638 / 560 histograms, max difference 0).
  The plots made from the merged cache are identical to plots made from the direct run's own cache.
- **Caveat: the three original macros' fitted plots depend on the thread count.** With implicit MT
  on, ROOT parallelizes the fit sums, and the changed summation order moves the tof background fits.
  Identical histograms fitted with N threads (a direct run) and with 1 thread (the merge step)
  therefore differ in the background-subtracted ratio panels (`*_proton_track_ratio_*`, one
  `apvsummary` panel, `*_ratio_*_cut*`). That's a few percent in 236 of 55k (tracking_eff) and 168
  of 66k (hodo_*) drawn histograms. All other plots are identical. This was confirmed by redrawing
  the same cache with and without MT. The macros are unchanged here, so the batch outputs match a
  1-thread run. Adding `ROOT::DisableImplicitMT();` after the event loop in a macro would make its
  plots independent of the thread count. `proton_tof_plot.C` does this, and its direct and merged
  outputs are identical.
- **Portability.** The workflow was tested end to end from a fresh `git clone` of the branch plus
  only the files copied in "Setting up as a new user", with all paths derived automatically.

## Step-by-step: how to reproduce

Paths are set in `slurm/config.sh`. They follow the checkout location and `$USER` automatically.
To put scratch elsewhere, `export WORK_BASE=...` before running the scripts.

### 0. Environment (login node)

```bash
cd /work/submit/$USER/hallc/software/lad_replay/CALIBRATION/lad_tof
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
source slurm/config.sh        # sets LADTOF_DIR and WORK_BASE in this shell
```

Don't run the macros themselves on the login node (see above). Compiling them, merging and
plotting from caches are all fine there.

### 1. Make the subMIT runlists (already done; copied in "Setting up")

This is how they were made. Rerun it with a different `B` if the ROOT files move:

```bash
cd ../files/run-lists
B=/ceph/submit/data/user/e/ehingerl/hallc/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION
for s in 5sigma 10sigma; do
  sed "s#^/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION#$B/$s#" \
      all_C3_runlist_SHMS_13p5.dat > all_C3_runlist_SHMS_13p5_submit_$s.dat
done
cd ../../lad_tof
```

### 2. Submit everything (one command)

```bash
bash slurm/submit_all.sh 20            # 20 = chunk size in GB; defaults to both 5sigma and 10sigma
# only one sigma:           bash slurm/submit_all.sh 20 10sigma
# only some of the macros:  bash slurm/submit_all.sh -m lad_hodo_eff,lad_hodo_dist 20
#                           bash slurm/submit_all.sh -m lad_tracking_eff 20 5sigma
# only the proton-tof plot: bash slurm/submit_all.sh -m proton_tof_plot 20
```

`-m` takes a comma-separated subset of
`lad_tracking_eff,lad_hodo_eff,lad_hodo_dist,proton_tof_plot`. The default is the first three;
`proton_tof_plot` is opt-in, because `lad_tracking_eff` already draws the same canvases. Macros you
leave out are not run and their existing caches and outputs are left alone.

This command does the following:

1. Compiles the selected macros and `slurm/fix_sig.C` once, so parallel jobs don't race ACLiC.
2. Deletes the selected macros' old caches for each sigma, so a fresh submission never reuses
   stale results.
3. Runs `slurm/make_chunks.sh <sigma> 20`, which writes
   `$WORK_BASE/<sigma>/chunks/chunk_NNN.dat` (46 chunks for 5σ, 29 for 10σ).
4. Submits `slurm/map_job.sh` as an array, one task per chunk. Each task gets 4 CPUs (Slurm
   usually allocates 12) and a 4 h limit. It gets 24 GB if `lad_tracking_eff` is selected, and
   4 GB otherwise, since the hodo macros peak at ~2.5 GB.
5. Submits `slurm/merge_job.sh` with `--dependency=afterok` on the array, so it starts only after
   every chunk has succeeded.

Logs go to `$WORK_BASE/logs/`. The jobs get `LADTOF_DIR`, `WORK_BASE` and the macro selection
(`LADTOF_MACROS`, set by `-m`) from the submitting shell via `--export=ALL`. The equivalent manual
commands for one sigma are:

```bash
source slurm/config.sh                                # exports LADTOF_DIR and WORK_BASE
export LADTOF_MACROS=lad_hodo_eff,lad_hodo_dist       # or leave unset for the default three
mkdir -p $WORK_BASE/logs
n=$(bash slurm/make_chunks.sh 10sigma 20)
map=$(sbatch --parsable --mem=4G --array=0-$((n-1)) --output=$WORK_BASE/logs/map_%x_%A_%a.log \
      --export=ALL,SIGMA=10sigma slurm/map_job.sh)
sbatch --dependency=afterok:$map --output=$WORK_BASE/logs/merge_%x_%j.log \
      --export=ALL,SIGMA=10sigma slurm/merge_job.sh
```

If you submit by hand, delete the selected macros' old `$WORK_BASE/<sigma>/cache/<macro>_chunk_*.root`
first. `map_job.sh` skips a macro whose chunk cache already exists (see step 4).

### 3. Monitor

```bash
squeue -u $USER                                             # pending/running jobs
ls $WORK_BASE/10sigma/cache/*_chunk_*.root | wc -l          # finished caches (one per macro per chunk)
tail $WORK_BASE/logs/map_ladtof_map_<jobid>_<task>.log      # one chunk's log (TIME lines give wall/maxrss)
sacct -j <jobid> -o JobID,State,Elapsed,MaxRSS              # accounting
```

### 4. If a chunk task fails

The merge job won't start: its dependency is `afterok`, so it stays pending (`DependencyNeverSatisfied`).
Cancel it, resubmit only the failed task(s), then resubmit the merge. Tasks skip any macro whose
cache already exists. If the original run used `-m`, export the same selection first.

```bash
source slurm/config.sh
export LADTOF_MACROS=<same selection as the original run, or leave unset for the default three>
scancel <merge_jobid>
map=$(sbatch --parsable --array=<failed task ids, e.g. 7,12> --output=$WORK_BASE/logs/map_%x_%A_%a.log \
      --export=ALL,SIGMA=10sigma slurm/map_job.sh)
sbatch --dependency=afterok:$map --output=$WORK_BASE/logs/merge_%x_%j.log \
      --export=ALL,SIGMA=10sigma slurm/merge_job.sh
```

For OOM failures, raise `--mem` (e.g. `sbatch --mem=32G ...`).

### 5. Merge by hand (optional)

The merge fits on the login node. For the full runlist it peaks at ~5 GB and takes ~3 min for all
three macros, which is within the 12 GB per-user cap if not much else is running:

```bash
SIGMA=10sigma bash slurm/merge_job.sh
LADTOF_MACROS=lad_hodo_dist SIGMA=10sigma bash slurm/merge_job.sh     # example of only one macro
```

For each macro, the merge step does the following:

```bash
hadd -f -k $WORK_BASE/10sigma/cache/lad_tracking_eff_merged.root $WORK_BASE/10sigma/cache/lad_tracking_eff_chunk_*.root
root -l -b -q 'slurm/fix_sig.C+("<merged cache>","../files/run-lists/all_C3_runlist_SHMS_13p5_submit_10sigma.dat")'
root -l -b -q 'lad_tracking_eff.C+("../files/run-lists/all_C3_runlist_SHMS_13p5_submit_10sigma.dat","files/tracking_eff/tracking_eff_C3_SHMS_13p5_v5_PH_10sigma.root",1,"<merged cache>")'
```

The log must say `cache HIT`. A MISS means the signature didn't match, and the macro would start a
full single-job event loop.

About `fix_sig.C`: each macro stores a signature string in its cache (binning, cuts, tracking
variants, and `;runlist=<std::hash of the runlist>`). The chunk caches carry their chunk's hash.
`fix_sig.C` replaces only the `;runlist=` field with the hash of the full runlist, computed exactly
as the macros compute it. The merged cache is then accepted for the full runlist.

### 6. Outputs

Final outputs, in your checkout's `CALIBRATION/lad_tof/` (these are the macros' default names with
`_<sigma>` appended):

```
files/tracking_eff/tracking_eff_C3_SHMS_13p5_v5_PH_{5sigma,10sigma}.root
files/hodo_eff/hodo_eff_C3_SHMS_13p5_v4_PH_{5sigma,10sigma}.root
files/hodo_dist/hodo_dist_C3_SHMS_13p5_v4_PH_{5sigma,10sigma}.root
files/proton_tof/proton_tof_C3_SHMS_13p5_PH_{5sigma,10sigma}.root      # only when run with -m proton_tof_plot
```

Merged caches are in `$WORK_BASE/<sigma>/cache/<macro>_merged.root`. To change only plotting code,
rerun step 5's last command against the merged cache. It takes about a minute per macro, with no
event loop.

**Only the `_c_proton_tof` canvas:** `proton_tof_plot.C` (in this folder) makes
`<spec>_c_proton_tof` on its own. It reads the raw replay ROOT files from a runlist and takes the
same arguments as the other macros (runlist, output, threads, cache), so it runs standalone or in
the batch workflow (`-m proton_tof_plot`, alone or with others). It fills only the 62 histograms
this canvas needs, with the same selection as `lad_tracking_eff.C`: event vertex, `isProton_1` hits
on planes 001/101, and paddle-centre path-length correction. It then draws the same four pads and
fits. It covers the standard and 1D tracking variants, and leaves out `x`, `noTrackVertex` and
`noTrackVertex_x`.

It's light: ~0.85 GB and ~40–65 s for 6.9 GB of input with 4 threads, against ~17 GB for
`lad_tracking_eff`. A few runs work on the login node. For the full dataset (~550–800 GB), use
batch (`-m proton_tof_plot`), because it's roughly 1–2 h in a single process.

```bash
root -l -b -q 'proton_tof_plot.C+("../files/run-lists/<runlist>.dat","proton_tof.root",4)'
root -l -b -q 'proton_tof_plot.C+("<runlist>","proton_tof.root",4,"","png_dir")'         # also one PNG per canvas
# redraw only (seconds), from this macro's cache OR an existing lad_tracking_eff cache:
root -l -b -q -e '.L proton_tof_plot.C+' -e 'proton_tof_draw("'$WORK_BASE'/10sigma/cache/lad_tracking_eff_merged.root","proton_tof_10sigma.root","png_10sigma")'
```

It writes 60 canvases (P/H × 3 χ² cuts × 10 variants) under
`<spec>/proton_id/chi2cut_<val>/<variant>/`. It was validated with `slurm/compare.C` on the 4-file
set in three ways. First, its 62 histograms are identical to `lad_tracking_eff`'s. Second, the
split-and-merged cache is identical to a direct run. Third, direct (4 threads) and merged outputs
are identical canvas for canvas. It always fits single-threaded (see the caveat under
Validation).

Scratch you can delete: `$WORK_BASE/<sigma>/chunk_out/` holds per-chunk plot files, which are
by-products. `$WORK_BASE/<sigma>/cache/*_chunk_*.root` can go too once the merge is done.


## Files

| file | purpose |
|---|---|
| `slurm/config.sh` | paths (derived from the checkout location and `$USER`; override with `LADTOF_DIR`/`WORK_BASE`), macro → output names, macro selection (`LADTOF_MACROS`) |
| `slurm/make_chunks.sh` | split a runlist into ~N GB chunk runlists |
| `slurm/map_job.sh` | Slurm array task: the selected macros on one chunk, writing caches |
| `slurm/merge_job.sh` | hadd caches, fix signature, plot from merged cache (selected macros) |
| `slurm/submit_all.sh` | compile, clear old caches, chunk, submit map array + dependent merge per sigma (`-m` to select macros) |
| `proton_tof_plot.C` | standalone `<spec>_c_proton_tof` from raw ROOT files (runlist) or from a cache; batch-compatible (no x / noTrackVertex variants) |
| `lad_gem_resid.C` | P GEM residuals (every cluster vs the vertex→hodoscope line, and layer 0 vs vertex→layer 1) and a χ²-free GEM hit efficiency per layer/axis; opt-in, run with `bash slurm/submit_all.sh -m lad_gem_resid 20` (4G, ~30 s per chunk); output `files/gem_resid/gem_resid_C3_SHMS_13p5_P_<sigma>.root` |
| `slurm/fix_sig.C` | rewrite a merged cache's runlist hash |
| `slurm/compare.C` | bin-by-bin comparison of two ROOT files (validation) |
