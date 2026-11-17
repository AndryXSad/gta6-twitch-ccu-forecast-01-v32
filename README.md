# GTA VI launch-day Twitch CCU forecaster — v32

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v32/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v32/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

A make.fun `gta6-twitch-ccu-forecast-01` entry. It forecasts, for each of
hours 2–24 of Grand Theft Auto VI's release day (2026-11-19 UTC), five
quantiles (5%, 25%, 50%, 75%, 95%) of the concurrent Twitch viewers in the
"Grand Theft Auto VI" category, from the first hour's print and the
launch-day curves of past big game launches shipped inside each case.

## How it works

`main.py` (Python 3.13, standard library only, no network) builds four
member predictors in log space from the 37 comparable past launches:

1. **Y5** — v4-family blend, tier EA filter (drop tier == ea), shape
   regressions on `[lnV1, g12]`. `g12` is the "unlock momentum" feature:
   lnV1 minus log viewers 12 hours before release day.
2. **Z0** — v4-family blend, naive EA filter (drop launches whose unlock
   note contains "early access"), shape regressions on `[lnV1, g12]`.
3. **Y1** — v4-family blend, naive EA filter, shape regressions on
   `[lnV1, lnPREM]`.
4. **Clone** — v9-clone baseline: whole-library shape views, sample-sd
   spreads, weighted-mean blend (0.6 OLS / 0.2 console-shape / 0.2
   all-shape) with spread 0.81× the weighted mean of the view spreads.

A v4-family member mixes three views per hour (0.5 OLS, 0.2
console-only shape, 0.3 all-launch shape): the median is the weighted
mean, the spread is the **sqrt-variance blend**
`sqrt(sum(w * sig^2) / sum(w))`. That last detail is the whole point of
this round: v31 mixed the spreads with a linear mean and its members
missed the cached rows by up to ~0.21 in log-median — no fixed linear
weights can fit a sqrt-variance spread. The shape views also differ per
member (fit_hours 10 vs 20, min_fits 6 vs 5, windowed vs whole-library
pools, population vs sample sd, OLS pool counted before vs after the
hour filter). This file implements each choice per member, and all four
members reproduce the cached rows to 0.000000 on practice AND on LOO.

The **log-median is chosen per hour band** — the median is the bottleneck:

| band | hours | median source | ramp |
| ---- | ----- | ------------- | ---- |
| early | 2–6 | Y1 log-median | 0.00 → +0.045 |
| mid | 7–14 | 0.6·Y5 + 0.4·clone | 0.00 → −0.14 |
| late | 15–24 | median of Y5 / Z0 / Y1 log-medians | 0.00 → +0.05 |

On top of that, a **size-dependent shift**: above lnV1 11.5 the median is
pulled down by 0.04 per log unit of size, in every band, and a **deeper,
size-indexed disagreement shift** `mu += (0.60 + 0.20 * size) * sd` of the
four member log-medians (late band base 0.45). When the members disagree
the median sits systematically low — shifting it up buys practice AND
LOO at once.

The **log-spread** is 0.975× the 3-way weighted mean of the member spreads
(45% Y5 + 35% Z0 + 20% clone), widened by a **disagreement factor**
`1 + 0.85 * sd(member medians)`, shaped per band with separate
lower/upper scales: early (0.90, 1.10), mid (1.10, 0.80), late
(0.95, 0.85), narrowed across the late band (spread scale 0.975 → 0.84
from hour 15 to hour 24 — this REPLACES the 0.975 base factor in the
late band), and finally **capped at 0.315**. Two size levers from v31:
**W1** grows only the lower spread with size (`slo *= 1 + 0.05 * size`)
and **W5** lets size amplify the disagreement effect on both sides
(`*= 1 + 1.25 * size * sd`). Quantiles are lognormal.

## Results

Direct-graded with the exact RULES scorer on the four practice cases
("Average viewers" answers): **practice pinball loss 0.117532**
(p00 0.1070 / p01 0.1366 / p02 0.0778 / p03 0.1488),
LOO-large **0.208604**, LOO-mega **0.210005**. That beats v30
(0.118710 / 0.209744 / 0.211400) on practice by −0.001178, with
LOO-large −0.001140 and LOO-mega −0.001395. Median-only practice is
**0.182041** (LOO-large 0.316268) — the median moved again, the
eighth round in a row. See `DESIGN.md` for the full reasoning record
and the honest caveats. Not submitted: the owner's rule is to submit
only a candidate that reaches the top 10.

## Run it

```sh
MAKEFUN_INPUT_DIR=<dir with *.in> MAKEFUN_OUTPUT_DIR=<dir> python3 main.py
```

Each `*.in` case produces a matching `*.out`:
`{"hours": [{"hour": h, "q": [q05, q25, q50, q75, q95]}]}` for h = 2..24.
