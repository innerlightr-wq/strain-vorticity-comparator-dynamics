"""Core coordinate algebra for the strain--vorticity / Thales-allocation audit.

Conventions (fixed once, used everywhere):

    (grad u)_ij = d u_i / d x_j
    S = (grad u + grad u^T)/2          symmetric, traceless iff div u = 0
    W = (grad u - grad u^T)/2          antisymmetric  (called Omega in the papers)
    omega_i = eps_ijk d_j u_k          vorticity, so W_ij = -1/2 eps_ijk omega_k
    |omega|^2 = 2 ||W||_F^2            (exact, see DERIVATIONS.md)

Coordinates:

    Qtot = ||grad u||_F^2 = ||S||_F^2 + ||W||_F^2      (Pythagoras)
    a    = ||S||_F^2 / Qtot      strain allocation
    b    = ||W||_F^2 / Qtot      rotation allocation      a + b = 1
    zeta = (||W||^2 - ||S||^2)/(||W||^2 + ||S||^2) = b - a
    h    = sqrt(a b)             Thales altitude
    L    = (b - a)/2 = zeta/2    Thales lateral displacement
    D    = 1/2 - h               Thales coherence deficit
    A    = omega.S omega / (||S||_F |omega|^2)           alignment coordinate
    P    = omega.S omega                                  enstrophy production

No SI units anywhere; every quantity below is either dimensionless or scales as a
power of inverse time.
"""

from __future__ import annotations

import numpy as np

SQRT_2_3 = np.sqrt(2.0 / 3.0)          # sharp bound on |A|
SHARP_P_CONST = 4.0 * np.sqrt(2.0) / 9.0   # sharp bound on P / ||grad u||_F^3


# ---------------------------------------------------------------- tensor algebra

def sym(G):
    """Symmetric part of a stack of 3x3 matrices (..., 3, 3)."""
    return 0.5 * (G + np.swapaxes(G, -1, -2))


def skew(G):
    """Antisymmetric part of a stack of 3x3 matrices."""
    return 0.5 * (G - np.swapaxes(G, -1, -2))


def frob2(M):
    """Squared Frobenius norm, summed over the last two axes."""
    return np.sum(M * M, axis=(-1, -2))


def vorticity_from_skew(W):
    """omega from W_ij = -1/2 eps_ijk omega_k, i.e. omega = (W_32-W_23, ...)."""
    W = np.asarray(W)
    return np.stack([W[..., 2, 1] - W[..., 1, 2],
                     W[..., 0, 2] - W[..., 2, 0],
                     W[..., 1, 0] - W[..., 0, 1]], axis=-1)


def skew_from_vorticity(w):
    """W_ij = -1/2 eps_ijk omega_k."""
    w = np.asarray(w)
    z = np.zeros(w.shape[:-1])
    return 0.5 * np.stack([
        np.stack([z, w[..., 2], -w[..., 1]], axis=-1),
        np.stack([-w[..., 2], z, w[..., 0]], axis=-1),
        np.stack([w[..., 1], -w[..., 0], z], axis=-1),
    ], axis=-2)


# ---------------------------------------------------------------- coordinates

def coordinates(S, w, eps=0.0):
    """All audit coordinates for stacks of strain tensors S and vorticities w.

    Returns a dict of arrays. Entries that are undefined (||S||=0 or |w|=0) are
    returned as NaN rather than silently regularized; `eps` is accepted only to
    make the degenerate-denominator convention explicit and defaults to zero.
    """
    S = np.asarray(S, dtype=float)
    w = np.asarray(w, dtype=float)

    nS2 = frob2(S)                       # ||S||_F^2
    w2 = np.sum(w * w, axis=-1)          # |omega|^2
    nW2 = 0.5 * w2                       # ||W||_F^2 = |omega|^2 / 2
    Qtot = nS2 + nW2                     # ||grad u||_F^2

    nS = np.sqrt(nS2)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(Qtot > eps, nS2 / Qtot, np.nan)
        b = np.where(Qtot > eps, nW2 / Qtot, np.nan)
        zeta = np.where(Qtot > eps, (nW2 - nS2) / Qtot, np.nan)

        P = np.einsum("...i,...ij,...j->...", w, S, w)
        denom = nS * w2
        A = np.where(denom > eps, P / denom, np.nan)

    h = np.sqrt(a * b)
    L = 0.5 * (b - a)
    D = 0.5 - h
    return dict(nS2=nS2, nW2=nW2, w2=w2, Qtot=Qtot, a=a, b=b, zeta=zeta,
                h=h, L=L, D=D, A=A, P=P,
                p_norm=np.where(Qtot > 0, P / Qtot ** 1.5, np.nan))


def coordinates_from_gradient(G):
    """Audit coordinates directly from a stack of velocity-gradient tensors."""
    S = sym(G)
    w = vorticity_from_skew(skew(G))
    return coordinates(S, w)


# ------------------------------------------------- deterministic zeta relations

def a_of_zeta(zeta):
    return 0.5 * (1.0 - np.asarray(zeta))


def b_of_zeta(zeta):
    return 0.5 * (1.0 + np.asarray(zeta))


def h_of_zeta(zeta):
    """Thales altitude as a function of zeta: h = (1/2) sqrt(1 - zeta^2)."""
    return 0.5 * np.sqrt(1.0 - np.asarray(zeta) ** 2)


def L_of_zeta(zeta):
    return 0.5 * np.asarray(zeta)


def D_of_zeta(zeta):
    return 0.5 * (1.0 - np.sqrt(1.0 - np.asarray(zeta) ** 2))


def allocation_factor(zeta):
    """g(zeta) = 2 b sqrt(a) = (1 + zeta) sqrt((1 - zeta)/2).

    The exact allocation factor in P = ||grad u||_F^3 * g(zeta) * A.
    Maximal at zeta = 1/3 with g = 4/(3 sqrt 3).
    """
    z = np.asarray(zeta, dtype=float)
    return (1.0 + z) * np.sqrt(np.clip((1.0 - z) / 2.0, 0.0, None))


def production_envelope(zeta):
    """Sharp bound on |P| / ||grad u||_F^3 at fixed zeta."""
    return SQRT_2_3 * allocation_factor(zeta)


# ------------------------------------------------- strain-state invariants

def strain_eigen(S):
    """Eigenvalues (ascending) and eigenvectors of a stack of symmetric S."""
    lam, vec = np.linalg.eigh(S)
    return lam, vec


def lund_rogers_s(S):
    """Normalized third strain invariant s = -3 sqrt(6) det S / ||S||_F^3.

    Lund & Rogers (1994) convention, verified against explicit states:
      s = -1  <->  lambda ~ (2,-1,-1): one extensional, two compressive axes
      s = +1  <->  lambda ~ (1,1,-2):  two extensional, one compressive axis
                   (the state isotropic turbulence prefers, s* ~ +0.5)
    Undefined (NaN) for S = 0.
    """
    S = np.asarray(S, dtype=float)
    det = np.linalg.det(S)
    n3 = frob2(S) ** 1.5
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n3 > 0, -3.0 * np.sqrt(6.0) * det / n3, np.nan)


def alignment_cosines(S, w):
    """cos^2 of the angle between omega and each strain eigenvector (ascending)."""
    lam, vec = strain_eigen(S)
    w = np.asarray(w, dtype=float)
    wn = w / np.linalg.norm(w, axis=-1, keepdims=True)
    c = np.einsum("...i,...ij->...j", wn, vec)
    return lam, c ** 2
