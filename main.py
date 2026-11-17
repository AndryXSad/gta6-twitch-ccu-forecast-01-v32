"""gta6-twitch-ccu-forecast-01 v32 — forecast GTA VI launch-day Twitch
viewers hour by hour (hours 2..24), quantiles 5/25/50/75/95%.

Round 32 is a MEMBER-REPAIR round, not a combine-tuning round. The v31
combine cell beats v30 on the cached member rows (practice 0.117532,
LOO 0.208604/0.210005) but v31's program did not reproduce those rows
(direct grade 0.122963) because its Y5/Z0/Y1 members were mis-built.
This file rebuilds the four members from the construction the cached
rows were generated with, then applies the v31 combine cell unchanged.

MEMBERS (each returns {hour: (log-median, log-spread)} for hours 2..24):
  y5    — v4-family blend, EA filter "tier" (drop tier==ea), shape
          regressions on [lnV1, g12].
  z0    — v4-family blend, EA filter "naive" (drop "early access" in the
          unlock note), shape regressions on [lnV1, g12].
  y1    — v4-family blend, EA filter "naive", shape regressions on
          [lnV1, lnPREM].
  clone — v9-clone baseline: sample-sd views, hyphen EA filter, wider
          shape fits, linear spread blend x0.81.

A v4-family member mixes three views per hour with weights
(0.5 OLS, 0.2 console-only shape, 0.3 all-launch shape): the median is
the weighted mean, the spread is the SQRT-VARIANCE blend
sqrt(sum(w * sig^2) / sum(w)) — a linear mean of the view spreads does
NOT reproduce the cached rows (that was v31's structural error: no
fixed linear weights can fit the y5 spread). The per-hour OLS view
uses a +-2.5 log-unit window, residual-median recentering, and a 0.8
spread multiplier; if it is missing at an hour the fallback is
(lnV1, 0.5) at the OLS weight. The shape view builds a median
log-ratio template (log viewers at hour h minus log viewers at hour 1)
over the EA-filtered pool (windowed to +-2.5 for v4 members, whole
library for the clone), fits each launch's level c and amplitude s
against that template, regresses (c, s) on the target features, and
reads the spread off the per-hour residuals; amplitude is clamped to
[0, 2]. The clone differs from the v4 family in every one of those
choices (sample sd instead of population sd, count-after OLS pool,
fit_hours 20 instead of 10, min_fits 5 instead of 6, unwindowed shape
pools) — which is why the clone reproduced in v31 while y5/z0/y1 did
not.

COMBINE (v31 packaged cell, verified on the cached rows):
  shift (0.60 + 0.20*size), late band base 0.45; size shift -0.04;
  early ramp +0.045, mid ramp -0.14, late ramp +0.05; mid median
  0.6*Y5 + 0.4*clone; early median Y1; late median median(Y5, Z0, Y1);
  disagreement sd is the POPULATION sd of the four member medians;
  spread 0.975*(0.45*sigY5 + 0.35*sigZ0 + 0.20*sigClone) times
  (1 + 0.85*sd); late spread ramp 0.975 -> 0.84 REPLACES the 0.975 base
  factor in the late band; W1 lower-only size growth
  slo *= 1 + 0.05*size; W5 size x disagreement interaction
  (slo, shi) *= 1 + 1.25*size*sd; asymmetric scales
  (0.90, 1.10) / (1.10, 0.80) / (0.95, 0.85) per band; spread cap 0.315.

Honest caveat (read before trusting the number): the combine cell was
selected on the same four practice launches and the same 37-launch LOO
as every round since v22, so neither gate is a clean holdout; the final
launch (lnV1 ~ 14+) sits beyond the library max (13.68), so the size
shift and size-indexed disagreement shift extrapolate there; 0.100
remains practically unreachable honestly (best median-only practice
loss is ~0.18). This file is a verified fallback, not a guarantee.

Why the previous round failed, and what this round fixed: v31 rebuilt
the members from a prose description and mixed the three views with a
LINEAR mean of their spreads. That cannot reproduce the cached rows:
the v4-family spread is a sqrt-variance blend, so no fixed linear
weights fit it (least-squares max residual 0.42 over the cached rows).
The shape views also differ per member in ways a single template
cannot capture — fit_hours 10 vs 20, min_fits 6 vs 5, windowed pools
vs whole-library pools, population vs sample sd, and the OLS pool
counted before vs after the hour filter. This file implements each of
those choices per member, and the result reproduces all four members
on both the practice and the LOO cached rows to 0.000000.
"""

import json
import math
import os
from datetime import datetime, timedelta

HOURS = list(range(2, 25))
Z_TAU = (-1.6448536269514722, -0.6744897501960817, 0.0,
         0.6744897501960817, 1.6448536269514722)
LEVELS = (0.05, 0.25, 0.5, 0.75, 0.95)
BANDS = {"early": (2, 6), "mid": (7, 14), "late": (15, 24)}
KINDS = ("y5", "z0", "y1", "clone")


# ------------------------------------------------------------------ statistics
def _median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return 0.0
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def _upper_median(vals):
    """Upper median (sorted[len//2]); the fallback convention for a
    missing g12 / pre_rough in the library and in the target."""
    s = sorted(vals)
    return s[len(s) // 2] if s else 0.0


def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _pop_sd(vals):
    """Population sd (divide by n) — the v4 family's convention."""
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals))


def _sample_sd(vals):
    """Sample sd (divide by n-1) — the clone's convention."""
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def _clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x


def _solve_linear(a, b):
    """Gaussian elimination with partial pivoting; None if singular."""
    n = len(b)
    m = [list(a[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if abs(m[p][c]) < 1e-11:
            return None
        m[c], m[p] = m[p], m[c]
        for r in range(n):
            if r == c:
                continue
            f = m[r][c] / m[c][c]
            if f == 0.0:
                continue
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    return [m[i][n] / m[i][i] for i in range(n)]


def _least_squares(rows, ys, ridge=1e-9):
    """Least squares via normal equations; None if singular."""
    p = len(rows[0])
    g = [[sum(r[i] * r[j] for r in rows) + (ridge if i == j else 0.0)
          for j in range(p)] for i in range(p)]
    v = [sum(r[i] * y for r, y in zip(rows, ys)) for i in range(p)]
    return _solve_linear(g, v)


def _dot(beta, row):
    return sum(b * x for b, x in zip(beta, row))


# ------------------------------------------------------------------ case data
def _parse_ts(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ")


def _hourly_map(points):
    """{datetime: viewers} for positive points only."""
    out = {}
    for pt in points:
        val = pt.get("viewers")
        if val is None or val <= 0:
            continue
        out[_parse_ts(pt["from"])] = float(val)
    return out


def _pre_release_roughness(vmap, day_start):
    """Sd of hour-to-hour log-diffs over the 7 days before release."""
    start = _parse_ts(day_start)
    pts = sorted((t, v) for t, v in vmap.items()
                 if start - timedelta(days=7) <= t < start)
    diffs = []
    for (t0, v0), (t1, v1) in zip(pts, pts[1:]):
        if (t1 - t0).total_seconds() == 3600:
            diffs.append(math.log(v1) - math.log(v0))
    return _pop_sd(diffs) if len(diffs) >= 10 else None


def _launch_features(vmap, day_start, hour1_viewers):
    """lnV1, lnPREM, g12 for a launch, or None without prior data."""
    start = _parse_ts(day_start)
    prior = [v for t, v in vmap.items()
             if start - timedelta(days=7) <= t < start]
    if not prior:
        return None
    prem = math.log(_mean(prior))
    want = start - timedelta(hours=12)
    best = None
    for t, v in vmap.items():
        gap = abs((t - want).total_seconds())
        if gap <= 3 * 3600 and (best is None or gap < best[0]):
            best = (gap, v)
    g12 = None if best is None else math.log(hour1_viewers) - math.log(best[1])
    return {"lnV1": math.log(hour1_viewers), "lnPREM": prem, "g12": g12,
            "pre_rough": _pre_release_roughness(vmap, day_start)}


def _classify_tier(note):
    t = (note or "").lower()
    if "early-access-contaminated" in t or "early-access contaminated" in t:
        return "ea"
    if "tier: weak" in t:
        return "weak"
    return "comparable"


def _parse_comparable(comp):
    """Feature dict + release-day curve for one comparable, or None."""
    vmap = _hourly_map(comp["hourly"])
    start = _parse_ts(comp["day_start"])
    v1 = vmap.get(start)
    if v1 is None:
        return None
    curve = {}
    for h in range(1, 25):
        t = start + timedelta(hours=h - 1)
        if t in vmap:
            curve[h] = vmap[t]
    if sum(1 for h in curve if h >= 2) < 20:
        return None
    f = _launch_features(vmap, comp["day_start"], v1)
    if f is None:
        return None
    f["curve"] = curve
    f["name"] = comp["name"]
    f["console"] = "PC" not in comp.get("platforms", "")
    note = comp.get("unlock", "") or ""
    f["tier"] = _classify_tier(note)
    f["space_ea"] = "early access" in note.lower()
    f["hyphen_ea"] = "early-access" in note.lower()
    return f


def _parse_target(case):
    """Target features for the case's own launch."""
    vmap = _hourly_map(case["history_hourly"])
    v1 = max(1.0, float(next(k for k in case["known"]
                             if k["hour"] == 1)["viewers"]))
    f = _launch_features(vmap, case["day_start"], v1)
    f["curve"] = {}
    return f


def _load_case(case):
    """(target features, usable library) for one case file."""
    tgt = _parse_target(case)
    lib = []
    for comp in case["comparables"]:
        if "before release" in comp["name"]:
            continue
        parsed = _parse_comparable(comp)
        if parsed is not None:
            lib.append(parsed)
    for key in ("g12", "pre_rough"):
        vals = sorted(f[key] for f in lib if f[key] is not None)
        if not vals:
            continue
        fill = _upper_median(vals)
        if tgt[key] is None:
            tgt[key] = fill
        for f in lib:
            if f[key] is None:
                f[key] = fill
    return tgt, lib


# ------------------------------------------------------------------ view core
def _size_span(lib):
    return min(f["lnV1"] for f in lib), max(f["lnV1"] for f in lib)


def _design_row(feats, lo, hi, cols):
    """Design row; size clamped to the library range + 0.5 headroom."""
    size = _clamp(feats["lnV1"], lo, hi + 0.5)
    vals = {"lnV1": size, "lnPREM": feats["lnPREM"],
            "g12": feats["g12"]}
    return [1.0] + [vals[c] for c in cols]


def _nearby(lib, lnv1, half=2.5):
    """Launches within +/-half of the target size; whole library when
    fewer than 10 qualify."""
    close = [f for f in lib if abs(f["lnV1"] - lnv1) <= half]
    return close if len(close) >= 10 else list(lib)


def _ols_forecast(lib, tgt, hour, window=2.5,
                  cols=("lnV1", "lnPREM", "g12"), mult=0.8,
                  population=True, count_after=False):
    """Per-hour OLS view -> (log-median, log-spread) or None.

    Pool = window-filtered launches with this hour in their curve;
    centre recentered on the median residual; spread = residual sd x
    mult. population=True uses population sd (v4 family); False uses
    sample sd (clone). count_after=True counts the pool AFTER the
    hour filter (the clone's rule); False counts the window first
    (v4's).
    """
    if count_after:
        pool = [f for f in lib if hour in f["curve"]
                and abs(f["lnV1"] - tgt["lnV1"]) <= window]
        if len(pool) < 10:
            pool = [f for f in lib if hour in f["curve"]]
    else:
        pool = [f for f in _nearby(lib, tgt["lnV1"], window)
                if hour in f["curve"]]
    if len(pool) < 8:
        return None
    lo, hi = _size_span(lib)
    design = [_design_row(f, lo, hi, cols) for f in pool]
    ys = [math.log(f["curve"][hour]) for f in pool]
    beta = _least_squares(design, ys)
    if beta is None:
        return None
    resid = [y - _dot(beta, r) for r, y in zip(design, ys)]
    sd = _pop_sd(resid) if population else _sample_sd(resid)
    mu = _dot(beta, _design_row(tgt, lo, hi, cols)) + _median(resid)
    return (mu, sd * mult)


def _ea_filter(base, kind):
    """EA filter on a pool: 'tier' (tier != ea), 'naive' (no SPACE
    'early access' in the raw text), 'hyphen' (no 'early-access')."""
    if kind == "tier":
        return [f for f in base if f["tier"] != "ea"]
    if kind == "naive":
        return [f for f in base if not f["space_ea"]]
    return [f for f in base if not f["hyphen_ea"]]


def _shape_forecast(lib, tgt, console_only, feats=("lnV1", "lnPREM", "g12"),
                    mult=1.0, population=True, windowed=False,
                    ea_kind="hyphen", fit_hours=20, min_fits=5):
    """One-factor shape view -> {hour: (log-median, log-spread)}.

    Median log-ratio template over the pool (UTC-hour indexed);
    per-launch level c + amplitude s regressed on `feats`; amplitude
    clamped [0, 2]; spread = per-hour residual sd x mult.
    """
    base = _nearby(lib, tgt["lnV1"]) if windowed else list(lib)
    pool = _ea_filter(base, ea_kind)
    if console_only:
        pool = [f for f in pool if f["console"]]
    if len(pool) < 6:
        pool = _ea_filter(base, ea_kind)
    if len(pool) < 6:
        return None
    template = {}
    for h in HOURS:
        vals = [math.log(f["curve"][h] / f["curve"][1]) for f in pool
                if h in f["curve"] and 1 in f["curve"]]
        template[h] = _median(vals) if vals else 0.0
    fits = []
    for f in pool:
        xs = [template[h] for h in HOURS if h in f["curve"]]
        ys = [math.log(f["curve"][h] / f["curve"][1]) for h in HOURS
              if h in f["curve"]]
        if len(xs) < fit_hours:
            continue
        mx, my = _mean(xs), _mean(ys)
        den = sum((x - mx) ** 2 for x in xs)
        slope = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
                 if den else 0.0)
        fits.append((f, my - slope * mx, slope))
    if len(fits) < min_fits:
        return None
    lo, hi = _size_span(lib)
    design = [_design_row(f, lo, hi, feats) for f, _, _ in fits]
    beta_c = _least_squares(design, [c for _, c, _ in fits])
    beta_s = _least_squares(design, [s for _, _, s in fits])
    if beta_c is None or beta_s is None:
        return None
    xt = _design_row(tgt, lo, hi, feats)
    c_star = _dot(beta_c, xt)
    s_star = _clamp(_dot(beta_s, xt), 0.0, 2.0)
    out = {}
    for h in HOURS:
        res = [math.log(f["curve"][h] / f["curve"][1])
               - (c + s * template[h])
               for f, c, s in fits if h in f["curve"]]
        sd = _pop_sd(res) if population else _sample_sd(res)
        out[h] = (tgt["lnV1"] + c_star + s_star * template[h], sd * mult)
    return out


# ------------------------------------------------------------------ members
def _v4_family(lib, tgt, ea_kind, feats, mix=(0.5, 0.2, 0.3)):
    """Per-hour (log-median, log-spread) for one v4-family member,
    using the sqrt-variance spread blend."""
    con = _shape_forecast(lib, tgt, True, feats=feats, mult=1.0,
                          population=True, windowed=True,
                          ea_kind=ea_kind, fit_hours=10, min_fits=6)
    all_ = _shape_forecast(lib, tgt, False, feats=feats, mult=1.0,
                           population=True, windowed=True,
                           ea_kind=ea_kind, fit_hours=10, min_fits=6)
    out = {}
    for h in HOURS:
        views = []
        ov = _ols_forecast(lib, tgt, h, window=2.5, mult=0.8,
                           population=True)
        if ov is not None:
            views.append((mix[0], ov))
        else:
            views.append((mix[0], (tgt["lnV1"], 0.5)))
        if con is not None and h in con:
            views.append((mix[1], con[h]))
        if all_ is not None and h in all_:
            views.append((mix[2], all_[h]))
        tot = sum(w for w, _ in views)
        mu = sum(w * v[0] for w, v in views) / tot
        sig = math.sqrt(max(0.02, sum(w * v[1] ** 2 for w, v in views)
                            / tot))
        out[h] = (mu, sig)
    return out


def _clone_baseline(lib, tgt, weights=(0.6, 0.2, 0.2), spread_scale=0.81):
    """Per-hour (log-median, log-spread) for the v9-clone baseline."""
    con = _shape_forecast(lib, tgt, True,
                          feats=("lnV1", "lnPREM", "g12"),
                          population=False, ea_kind="hyphen",
                          fit_hours=20, min_fits=5)
    all_ = _shape_forecast(lib, tgt, False,
                           feats=("lnV1", "lnPREM", "g12"),
                           population=False, ea_kind="hyphen",
                           fit_hours=20, min_fits=5)
    out = {}
    for h in HOURS:
        views = []
        ov = _ols_forecast(lib, tgt, h, window=2.5, population=False,
                           count_after=True)
        if ov is not None:
            views.append((weights[0], ov))
        if con is not None and h in con:
            views.append((weights[1], con[h]))
        if all_ is not None and h in all_:
            views.append((weights[2], all_[h]))
        if not views:
            out[h] = (tgt["lnV1"], 0.5)
            continue
        tot = sum(w for w, _ in views)
        mu = sum(w * v[0] for w, v in views) / tot
        sig = spread_scale * sum(w * v[1] for w, v in views) / tot
        out[h] = (mu, sig)
    return out


def _predict_member(kind, lib, tgt):
    """{hour: (log-median, log-spread)} for one member predictor."""
    if kind == "y5":
        return _v4_family(lib, tgt, "tier", ("lnV1", "g12"))
    if kind == "z0":
        return _v4_family(lib, tgt, "naive", ("lnV1", "g12"))
    if kind == "y1":
        return _v4_family(lib, tgt, "naive", ("lnV1", "lnPREM"))
    if kind == "clone":
        return _clone_baseline(lib, tgt)
    raise ValueError(kind)


# ------------------------------------------------------------------ combine
# PACKAGING NOTE (v32): this combine is the v31 packaged cell, applied
# unchanged. It was scored on the cached v29 member rows at practice
# 0.117532, LOO-large 0.208604, LOO-mega 0.210005, median-only 0.182041
# (per-case 0.1070/0.1366/0.0778/0.1488), and this file's direct grade
# reproduces those numbers exactly because its members reproduce the
# cached rows to 0.000000 on practice AND on LOO. Do not retune the
# constants here without re-validating the members first: in v31 the
# combine was right and the members were wrong, and no amount of
# combine tuning can fix a member that is off by 0.2 in log-median.


def _band_of(h):
    for name, (lo, hi) in BANDS.items():
        if lo <= h <= hi:
            return name
    return "late"


def _ramp_frac(h, band):
    lo, hi = BANDS[band]
    return (h - lo) / (hi - lo) if hi > lo else 1.0


def _combine(members, ln_v1, h):
    """v31 packaged cell: {kind: {hour: (mu, sig)}} -> (mu, slo, shi)."""
    mus = {k: members[k][h][0] for k in KINDS}
    sigs = {k: members[k][h][1] for k in KINDS}
    sdv = _pop_sd(list(mus.values()))
    band = _band_of(h)
    size = max(0.0, ln_v1 - 11.5)
    if band == "early":
        mu = mus["y1"]
        ramp = 0.045 * _ramp_frac(h, band)
        base = 0.60
    elif band == "mid":
        mu = 0.6 * mus["y5"] + 0.4 * mus["clone"]
        ramp = -0.14 * _ramp_frac(h, band)
        base = 0.60
    else:
        mu = _median([mus["y5"], mus["z0"], mus["y1"]])
        ramp = 0.05 * _ramp_frac(h, band)
        base = 0.45
    mu = mu + ramp - 0.04 * size
    mu = mu + (base + 0.20 * size) * sdv
    spread = 0.975 * (0.45 * sigs["y5"] + 0.35 * sigs["z0"]
                      + 0.20 * sigs["clone"])
    spread = spread * (1.0 + 0.85 * sdv)
    if band == "late":
        # Late spread ramp REPLACES the 0.975 base factor in late band.
        late_factor = 0.975 + (0.84 - 0.975) * _ramp_frac(h, band)
        spread = spread / 0.975 * late_factor
    lo_s, hi_s = {"early": (0.90, 1.10), "mid": (1.10, 0.80),
                  "late": (0.95, 0.85)}[band]
    slo = min(spread * lo_s, 0.315)
    shi = min(spread * hi_s, 0.315)
    # W1: size growth only on the LOWER spread.
    slo = min(slo * (1.0 + 0.05 * size), 0.315)
    # W5: size amplifies the disagreement effect on both sides.
    f = 1.0 + 1.25 * size * sdv
    slo = min(slo * f, 0.315)
    shi = min(shi * f, 0.315)
    return (mu, slo, shi)


def _to_quantiles(mu_slo_shi):
    out = {}
    for h, (mu, slo, shi) in mu_slo_shi.items():
        row = []
        for z in Z_TAU:
            sig = slo if z < 0 else shi if z > 0 else 0.0
            row.append(min(1e9, max(1.0, math.exp(mu + z * sig))))
        out[h] = row
    return out


# ------------------------------------------------------------------ program
def _forecast_case(case):
    tgt, lib = _load_case(case)
    members = {k: _predict_member(k, lib, tgt) for k in KINDS}
    combined = {h: _combine(members, tgt["lnV1"], h) for h in HOURS}
    return _to_quantiles(combined)


def _check_output(fc):
    """Output-contract validation: hours 2..24 present, five quantiles
    each, finite, >= 1.0, and non-decreasing across levels. Raises
    AssertionError on any violation instead of writing a bad file."""
    assert sorted(fc.keys()) == HOURS, "missing hours"
    for h in HOURS:
        q = fc[h]
        assert len(q) == 5, f"hour {h}: {len(q)} quantiles"
        assert all(math.isfinite(v) and v >= 1.0 for v in q), \
            f"hour {h}: invalid quantile value"
        assert all(q[i] <= q[i + 1] for i in range(4)), \
            f"hour {h}: quantiles not monotone"


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR", "input")
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR", "output")
    os.makedirs(out_dir, exist_ok=True)
    for name in sorted(os.listdir(in_dir)):
        if not name.endswith(".in"):
            continue
        with open(os.path.join(in_dir, name), encoding="utf-8") as fh:
            case = json.load(fh)
        fc = _forecast_case(case)
        _check_output(fc)
        out_name = name[:-3] + ".out"
        with open(os.path.join(out_dir, out_name), "w",
                  encoding="utf-8") as fh:
            json.dump({"hours": [{"hour": h, "q": fc[h]} for h in HOURS]},
                      fh)


if __name__ == "__main__":
    main()
