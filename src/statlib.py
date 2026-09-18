"""Small statistics toolkit: OLS with HC3 errors, binned/spline controls,
binned mutual information with a permutation-calibrated floor, bootstrap CIs.

Deliberately explicit rather than pulling in statsmodels/sklearn: every number
reported in RESULTS.md should be traceable to a few lines of arithmetic.
"""

from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------ regression

def ols(X, y):
    """Least squares with HC3 (heteroskedasticity-robust) standard errors."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    # HC3: leverage-corrected sandwich
    hat = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    w = (resid / np.clip(1.0 - hat, 1e-12, None)) ** 2
    cov = XtX_inv @ (X.T * w) @ X @ XtX_inv
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, beta / se, np.nan)
    return dict(beta=beta, se=se, t=tstat, r2=float(r2), resid=resid, n=n, k=k)


def normal_two_sided_p(t):
    """Two-sided normal p-value (n is large everywhere we use this)."""
    from math import erfc, sqrt
    return np.array([erfc(abs(v) / sqrt(2.0)) for v in np.atleast_1d(t)])


def design(*cols, intercept=True):
    """Stack columns / column blocks into a design matrix (2D blocks allowed)."""
    blocks = []
    n = None
    for c in cols:
        c = np.asarray(c, float)
        c = c.reshape(len(c), -1)
        n = len(c) if n is None else n
        assert len(c) == n, "design blocks must share their first dimension"
        blocks.append(c)
    if intercept:
        blocks.insert(0, np.ones((n, 1)))
    return np.hstack(blocks)


def quantile_bins(x, nbins):
    """Bin indices from quantiles of x (ties collapsed)."""
    x = np.asarray(x, float)
    edges = np.unique(np.quantile(x[np.isfinite(x)], np.linspace(0, 1, nbins + 1)))
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    return idx, edges


def bin_dummies(x, nbins):
    """Indicator design matrix for quantile bins of x (drop first for rank)."""
    idx, edges = quantile_bins(x, nbins)
    nb = len(edges) - 1
    D = np.zeros((len(x), nb - 1))
    for b in range(1, nb):
        D[idx == b, b - 1] = 1.0
    return D, idx, edges


def group_mean_residual(y, idx):
    """Residual of y after removing per-bin means -- exactly the residual of an
    OLS on bin dummies, computed in O(n) instead of O(n k^2)."""
    y = np.asarray(y, float)
    nb = idx.max() + 1
    cnt = np.bincount(idx, minlength=nb).astype(float)
    sm = np.bincount(idx, weights=y, minlength=nb)
    mean = np.where(cnt > 0, sm / np.maximum(cnt, 1), 0.0)
    return y - mean[idx]


def binned_partial(x, y, zcond, nbins):
    """Partial correlation and exact delta-R^2 for adding x to a bin model in z.

    Returns corr(res_x, res_y | z bins), the R^2 of the bin model for y, and the
    increment in R^2 from adding x on top of it.
    """
    idx, _ = quantile_bins(zcond, nbins)
    rx = group_mean_residual(x, idx)
    ry = group_mean_residual(y, idx)
    y = np.asarray(y, float)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_within = float(np.sum(ry ** 2))
    denom = float(np.sqrt(np.sum(rx ** 2) * ss_within))
    rho = float(np.sum(rx * ry) / denom) if denom > 0 else np.nan
    r2_base = 1.0 - ss_within / ss_tot if ss_tot > 0 else np.nan
    delta = (rho ** 2) * ss_within / ss_tot if ss_tot > 0 else np.nan
    n_eff = len(y)
    t = rho * np.sqrt(max(n_eff - nbins - 1, 1) / max(1.0 - rho ** 2, 1e-30))
    return dict(n_bins=int(nbins), partial_corr=rho, r2_bins_only=r2_base,
                delta_r2=float(delta), t=float(t),
                p=float(normal_two_sided_p([t])[0]))


def r2_of(y, yhat):
    y = np.asarray(y, float)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1.0 - np.sum((y - yhat) ** 2) / ss_tot) if ss_tot > 0 else np.nan


def nested_r2(y, base_cols, extra_cols):
    """R^2 of base model, of base+extra, and the increment."""
    Xb = design(*base_cols)
    Xf = design(*base_cols, *extra_cols)
    fb, ff = ols(Xb, y), ols(Xf, y)
    return dict(r2_base=fb["r2"], r2_full=ff["r2"], delta_r2=ff["r2"] - fb["r2"],
                beta_extra=ff["beta"][Xb.shape[1]:].tolist(),
                se_extra=ff["se"][Xb.shape[1]:].tolist(),
                t_extra=ff["t"][Xb.shape[1]:].tolist(),
                p_extra=normal_two_sided_p(ff["t"][Xb.shape[1]:]).tolist())


def bootstrap_ci(stat_fn, n, n_boot=400, seed=0, alpha=0.05):
    rng = np.random.default_rng(seed)
    vals = np.array([stat_fn(rng.integers(0, n, n)) for _ in range(n_boot)])
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2], axis=0)
    return float(lo), float(hi)


# ------------------------------------------------------- mutual information

def _entropy_from_counts(counts):
    tot = counts.sum()
    p = counts[counts > 0] / tot
    return float(-(p * np.log(p)).sum())


def mutual_info_binned(x, y, nb=16):
    """Plug-in MI in nats from quantile bins (bias assessed separately)."""
    ix, _ = quantile_bins(x, nb)
    iy, _ = quantile_bins(y, nb)
    joint = np.histogram2d(ix, iy, bins=[np.arange(nb + 1) - 0.5] * 2)[0]
    return (_entropy_from_counts(joint.sum(axis=1)) +
            _entropy_from_counts(joint.sum(axis=0)) - _entropy_from_counts(joint))


def conditional_mutual_info(x, y, zcond, nb=12, nz=8):
    """I(X;Y|Z) in nats, quantile-binned in all three variables."""
    iz, _ = quantile_bins(zcond, nz)
    total, n = 0.0, len(x)
    for b in range(iz.max() + 1):
        m = iz == b
        if m.sum() < 50:
            continue
        total += (m.sum() / n) * mutual_info_binned(x[m], y[m], nb=nb)
    return float(total)


def cmi_with_null(x, y, zcond, nb=12, nz=8, n_shuffle=12, seed=0):
    """I(X;Y|Z) together with its within-stratum-shuffled null (the bias floor)."""
    rng = np.random.default_rng(seed)
    obs = conditional_mutual_info(x, y, zcond, nb=nb, nz=nz)
    iz, _ = quantile_bins(zcond, nz)
    nulls = []
    for _ in range(n_shuffle):
        ysh = y.copy()
        for b in range(iz.max() + 1):
            m = np.flatnonzero(iz == b)
            if len(m) > 1:
                ysh[m] = y[rng.permutation(m)]
        nulls.append(conditional_mutual_info(x, ysh, zcond, nb=nb, nz=nz))
    nulls = np.array(nulls)
    return dict(cmi=obs, null_mean=float(nulls.mean()), null_std=float(nulls.std()),
                excess=float(obs - nulls.mean()),
                z_score=float((obs - nulls.mean()) / nulls.std()) if nulls.std() > 0
                else float("inf"))
