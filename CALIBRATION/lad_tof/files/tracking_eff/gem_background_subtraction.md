# Background Subtraction for the GEM Hit-Position Plots

Tracking-efficiency macro (lad_tracking_eff.C) — the two derived GEM panels IT − OOT and peak − IT − OOT.

## The model

The pad-2 (proton+track) corrected-tof fit (trapgaus) splits the tof spectrum into three additive components:

- Flat F — constant in tof; the random / accidental coincidence background.
- Trapezoid T(t) — fixed corners at t = −75, −25, 50, 125 ns, plateau height (p1 − p0) above the flat line; the in-time, non-signal background.
- Gaussian G(t) — centered at 41 ns; the signal peak.

The GEM (x, y) histograms are filled separately for four tof windows: OOT (out-of-time sidebands [−150,−100] ∪ [125,175]), IT (in-time sidebands [−25,30] ∪ [50,125]), peak ([30,50]), and "all". Call the corresponding position histograms H_OOT, H_IT, H_PK.

Core assumption: each physics component has one (x, y) shape that does not depend on the tof window — only its amount changes between windows, and the fit gives those amounts. So the windows are subtracted against one another with scale factors matched to those amounts.

## Scale factor 1 — flat background (exact, from widths)

The flat term is uniform in tof, so its accidental counts in any window are proportional to that window's total tof width. OOT is pure flat (the trapezoid is zero outside [−75,125] and the gaussian is far away), so H_OOT is a clean spatial template of the accidental background. To predict the flat contribution inside another window R, scale H_OOT by the width ratio:

```
s_flat^R = w_R / w_OOT            (exact — tof widths only)

w_OOT = 50 + 50 = 100 ns
w_IT  = 55 + 75 = 130 ns   ->   s_flat^IT = 130 / 100 = 1.30
w_PK  = 20      =  20 ns    ->   s_flat^PK =  20 / 100 = 0.20
```

This is the exact scale factor — no fit values enter, only the region widths. In the code these are sflat_it and sflat_pk.

## Scale factor 2 — trapezoid background (from integrals / areas)

The IT window contains flat + trapezoid and, by construction, no peak ([30,50] is excluded from IT). So the difference H_IT − s_flat^IT · H_OOT is a spatial template of the trapezoid component alone, holding an amount equal to the trapezoid's area over the IT window. To predict the trapezoid amount inside the peak window, scale that template by the ratio of the fitted trapezoid's areas in the two windows:

```
s_trap = ∫_PK T(t) dt  /  ∫_IT T(t) dt
```

The trapezoid amplitude cancels in this ratio, so only the fixed corner geometry matters (the code integrates the fitted TF1 over each region's intervals — itrPK, itrIT). Numerically:

```
∫_PK T = A · 20                         (peak is entirely on the plateau)
∫_IT T = A · 55 + A · 37.5 = A · 92.5    ([-25,30] plateau + [50,125] falling edge)
s_trap = 20 / 92.5 = 0.2162
```

## The two subtracted panels

```
IT − OOT        :  H_IT − s_flat^IT · H_OOT
                =  trapezoid (in-time background) spatial distribution

peak − IT − OOT :  H_PK − s_flat^PK · H_OOT − s_trap · ( H_IT − s_flat^IT · H_OOT )
                =  gaussian (signal) spatial distribution
```

In words: IT − OOT removes the flat / accidental floor from the in-time window, leaving where the in-time background strikes the GEMs. peak − IT − OOT removes the flat floor and the in-time trapezoid from the peak window, leaving the spatial distribution of the true signal protons. Each removed term is scaled to the amount the fit assigns to the target window — the region width for the flat piece, the trapezoid area for the in-time piece. In the macro these two derived histograms are built by mkITmOOT and mkPeakBS.

## Caveat

The trapezoid template comes from the raw IT histogram, whose inner edges (30 and 50 ns) sit ~1.8–2.2σ from the 41 ns peak, so a small gaussian tail leaks into IT. peak − IT − OOT therefore slightly over-subtracts the signal tails. The effect is small; making it exact would require modelling the gaussian leakage into the IT window rather than treating IT as pure flat + trapezoid.
