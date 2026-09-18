"""Ensemble generators shared by src/ensembles.py and src/homogeneity_check.py.

Kept in its own module so that importing a generator does not execute an analysis
script.  All generators take an explicit seed; see EXPERIMENT.md for the values used.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import svcore as sv  # noqa: E402,F401  (used by generators)


def gaussian_solenoidal_gradients(n_grid=64, n_real=4, k0=4.0, seed=1):
    """Velocity gradients of a divergence-free isotropic Gaussian random field."""
    rng = np.random.default_rng(seed)
    N = n_grid
    k1 = np.fft.fftfreq(N, d=1.0 / N)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    K2s = np.where(K2 == 0, 1.0, K2)
    kmag = np.sqrt(K2)
    # E(k) ~ k^4 exp(-2 k^2 / k0^2)  ->  amplitude ~ sqrt(E/k^2)
    amp = np.where(K2 > 0, kmag ** 2 * np.exp(-K2 / k0 ** 2), 0.0)
    grads = []
    for _ in range(n_real):
        uh = []
        for _c in range(3):
            w = rng.normal(size=(N, N, N))
            uh.append(np.fft.fftn(w) * amp)
        uh = np.stack(uh)                                     # (3, N, N, N)
        kdotu = KX * uh[0] + KY * uh[1] + KZ * uh[2]
        uh = uh - np.stack([KX, KY, KZ]) * kdotu / K2s        # project solenoidal
        G = np.empty((N, N, N, 3, 3))
        Kv = (KX, KY, KZ)
        for i in range(3):
            for j in range(3):
                G[..., i, j] = np.fft.ifftn(1j * Kv[j] * uh[i]).real
        grads.append(G.reshape(-1, 3, 3))
    G = np.concatenate(grads)
    # enforce tracelessness of S up to round-off (spectral projection is exact
    # to machine precision; this only removes accumulated FFT round-off)
    tr = np.trace(G, axis1=-2, axis2=-1)
    G = G - np.eye(3) * (tr / 3.0)[:, None, None]
    return G


def random_traceless_symmetric(n, rng):
    M = rng.normal(size=(n, 3, 3))
    S = 0.5 * (M + np.swapaxes(M, -1, -2))
    tr = np.trace(S, axis1=-2, axis2=-1)
    return S - np.eye(3) * (tr / 3.0)[:, None, None]


def isotropic_dirs(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def ensemble_naive(n, seed=2):
    """E2: independent random S and omega, magnitudes lognormal."""
    rng = np.random.default_rng(seed)
    S = random_traceless_symmetric(n, rng)
    S = S / np.linalg.norm(S, axis=(-2, -1))[:, None, None]
    S = S * rng.lognormal(0.0, 1.0, size=n)[:, None, None]
    w = isotropic_dirs(n, rng) * rng.lognormal(0.0, 1.0, size=n)[:, None]
    return S, w


def ensemble_zeta_designed(n, seed=3):
    """E3: zeta ~ U(-1,1), strain shape random, orientation isotropic."""
    rng = np.random.default_rng(seed)
    S = random_traceless_symmetric(n, rng)
    S = S / np.linalg.norm(S, axis=(-2, -1))[:, None, None]      # ||S||_F = 1
    zeta = rng.uniform(-0.995, 0.995, size=n)
    nW2 = (1.0 + zeta) / (1.0 - zeta)                            # since ||S||^2 = 1
    w = isotropic_dirs(n, rng) * np.sqrt(2.0 * nW2)[:, None]
    return S, w, zeta


def ensemble_aligned(n, seed=4, which=2, kappa=4.0):
    """E4: omega tilted towards eigenvector `which` (0=min,1=mid,2=max) of S."""
    rng = np.random.default_rng(seed)
    S = random_traceless_symmetric(n, rng)
    S = S / np.linalg.norm(S, axis=(-2, -1))[:, None, None]
    lam, vec = np.linalg.eigh(S)                                 # ascending
    axis = vec[..., which]
    v = isotropic_dirs(n, rng)
    d = axis * kappa + v                                         # tilt then renormalize
    d = d / np.linalg.norm(d, axis=-1, keepdims=True)
    zeta = rng.uniform(-0.995, 0.995, size=n)
    nW2 = (1.0 + zeta) / (1.0 - zeta)
    return S, d * np.sqrt(2.0 * nW2)[:, None]


