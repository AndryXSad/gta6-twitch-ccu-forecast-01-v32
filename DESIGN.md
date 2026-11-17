# v32 — design record

## What v32 is

v32 is a MEMBER-REPAIR round, not a combine-tuning round. The v31
combine cell already beat v30 on the cached member rows (practice
0.117532, LOO 0.208604/0.210005) — but v31's program did not reproduce
those rows (direct grade 0.122963, worse than v30), because its
Y5/Z0/Y1 members were mis-built. This round fixes the member
construction and applies the v31 combine cell unchanged. The result
reproduces the cached-row scores exactly in direct grading.

## The member gap, and what fixed it

v31's diagnosis (BUILD-v31.md, "The member gap") established:

- The clone member reproduced the cached rows almost exactly (sig max
  diff 0.00014; mu max diff 0.019, only at p00) — so the OLS view and
  the shape-sig formula were right in v31's code.
- Y5/Z0/Y1 did not reproduce: practice mu max diff 0.214 / 0.206 /
  0.190; sig max diff 0.185 / 0.171 / 0.171.
- Ruled out by v31's 18 diagnostic scripts: per-member OLS features,
  OLS window, OLS recentering, OLS library filter, shape template
  (median/mean), (c, s) fit variants, per-hour OLS shape view, mix
  weights (grids + least-squares, 2-view blends), sd conventions
  (ddof 0/1), EA filter variants, size clamp, library composition.
- Sharpest fact: the y5 sig **cannot be fit by any fixed weights** on
  the validated view sigs (least-squares max residual 0.42 over 943
  rows) — the y5/z0/y1 shape views differ STRUCTURALLY, not just in
  weights. And the y5 mu error pattern was consistent with the shape
  amplitude (s) being too small for big targets and too big for p00.

The missing structural detail, recovered this round from the harness
that generated the cached rows (v22's work/h22.py, which v23 verified
reproduces v21's rows byte-identically): **the v4-family spread blend
is a sqrt-variance blend** — `sig = sqrt(sum(w * sig^2) / sum(w))` —
not a linear weighted mean of the view spreads. v31's mix function
used a linear mean, so no fixed weights could ever fit the y5 spread.
The second missing piece is that the shape views differ per member in
several small choices at once:

| choice | v4 family (y5/z0/y1) | clone |
| --- | --- | --- |
| shape pool | windowed (+-2.5 log units) | whole library |
| fit_hours / min_fits | 10 / 6 | 20 / 5 |
| residual sd | population (ddof=0) | sample (ddof=1) |
| OLS pool | window counted first | hour filter counted first |
| OLS sd | population | sample |
| spread blend | sqrt-variance | 0.81 x linear mean |

Plus the per-member EA filters and shape features: y5 = tier filter +
[lnV1, g12]; z0 = naive filter + [lnV1, g12]; y1 = naive filter +
[lnV1, lnPREM]; clone = hyphen filter + [lnV1, lnPREM, g12]. And the
OLS fallback: when the v4 OLS view is missing at an hour, the fallback
is (lnV1, 0.5) at the OLS weight; the clone has no such fallback.

This file implements each of those choices per member. Validation
(work/h32.py): all four members reproduce the cached rows to
**0.000000** max |diff| on practice (92 rows) AND on LOO (851 rows),
mu and sig. The combine then scores 0.117532 / 0.208604 / 0.210005 /
median-only 0.182041 from the rows, and the direct grade reproduces
0.117532 exactly (per-case 0.1070/0.1366/0.0778/0.1488).

## The packaged cell (v31, unchanged)

Shift (0.60 + 0.20*size), late base 0.45; size shift -0.04; early ramp
+0.045, mid ramp -0.14, late ramp +0.05; mid median 0.6*Y5 + 0.4*clone;
early median Y1; late median median(Y5, Z0, Y1); k_dis 0.85; spread
0.975*(0.45*sigY5 + 0.35*sigZ0 + 0.20*sigClone); late spread ramp
0.975 -> 0.84 REPLACES the 0.975 base factor in the late band;
W5 spread *= 1 + 1.25*size*sd; W1 sig_lo *= 1 + 0.05*size; asym
(0.90,1.10)/(1.10,0.80)/(0.95,0.85); cap 0.315; disagreement sd is the
population sd (ddof=0) of the four member medians.

## Results (direct-graded, exact RULES scorer, "Average viewers")

- Practice **0.117532** (p00 0.1070 / p01 0.1366 / p02 0.0778 /
  p03 0.1488) vs v30 0.118710 — **-0.001178**.
- LOO-large **0.208604** vs 0.209744 — **-0.001140**.
- LOO-mega **0.210005** vs 0.211400 — **-0.001395**.
- Median-only practice **0.182041** (v30: 0.184501) — the eighth round
  in a row the median-only loss moves. Median-only LOO-large 0.316268
  (v30: 0.317142), LOO-mega 0.316741 (v30: 0.317647).

Two direct grades were used (this candidate, before and after the
three authorship edits — identical scores both times). The harness
sweep numbers matched the direct grade to six decimals, because the
harness scores through the same `combine`.

## Honest caveats

1. **Selection bias.** Every parameter was chosen on the same four
   practice launches and the same 37-launch LOO that all previous
   rounds selected on. The margins (-0.001178 practice, -0.0011 LOO)
   are thin, and part of them is selection noise. The practice board
   is only a proxy: placements are held until finalize, and the final
   is graded on GTA VI's real release day.
2. **The practice board is flooded.** As of 17 Nov ~21:08 WIB, 198
   entries score below 0.122315 — the top is dominated by replay /
   memorization entries (67 entries at exactly 1e-06). v32 at 0.117532
   would sit around rank #199 on the practice board, nowhere near the
   top 20, let alone the top 10 the owner now requires for submission.
   The final rewards honest generalization (v32's LOO), not memorized
   practice curves — but whether prizes follow the practice rank or the
   final rank is decided by the platform's finalize step, and that is
   not in this repo's control.
3. **Extrapolation.** The final launch (lnV1 ~ 14+) sits beyond the
   library's largest launch (13.68). The size shift (-0.04*size) and
   the size-indexed disagreement shift (0.60 + 0.20*size) both
   extrapolate hardest exactly there, and no gate in this repo can
   validate them at that size. The shift parameters in the packaged
   cell are the most extreme ever used (0.60 base, 0.20 slope, mid
   ramp -0.14) — the frontier had not flattened in v31's search, but
   the parameters were getting extreme, so the search stopped there.
4. **0.100 is practically unreachable honestly.** The median is the
   bottleneck (best median-only practice ever is now 0.182041), the
   per-hour OLS residual sd is 0.35–0.45 log, and the practice answers
   are 2–2.5× rougher than the comparables. 0.100 would need a near-
   perfect smooth model. v32 gains -0.001178, the biggest single-round
   practice gain since v25 — but it comes from the v31 combine cell,
   not from a new median idea; the member repair only unlocked it.
5. **Not submitted.** The owner's rule (17 Nov ~21:21 WIB) is to submit
   only a candidate that reaches the top 10. v32 does not. It is a
   verified fallback, better than v30 on all three gates.
