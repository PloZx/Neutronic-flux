#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualisation 3D du flux ϕ
- Fenêtre 1 : cas stationnaire (3D)
- Fenêtre 2 : évolution temporelle (3D)
Protocole d'entrée (stdin) :
  1) une ligne avec n valeurs (le vecteur stationnaire), sans en-tête
  2) répétitions de blocs :
       "# step K t=<valeur>"
       "<n valeurs séparées par des espaces>"
"""

import argparse
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

# ---------- Même logique de géométrie que le code C ----------
def build_map_like_c(m: int, alpha: float, L:float):
    nx = m - 2
    h = (L / (m - 1)) * alpha
    map_idx = -np.ones((nx, nx), dtype=int)
    k = 0
    iL_ix, jT_iy = None, -1
    for iy in range(nx):
        y = (iy + 1) * h
        for ix in range(nx):
            x = (ix + 1) * h
            in_notch = (x >= (1.5/2) * L * alpha) and (x <= L * alpha) and (y >= 0.0) and (y <= (0.75/2) * L * alpha)
            if in_notch:
                if iL_ix is None or ix < iL_ix:
                    iL_ix = ix
                if iy > jT_iy:
                    jT_iy = iy
            else:
                map_idx[iy, ix] = k
                k += 1
    iL_full = (iL_ix + 1) if iL_ix is not None else nx
    jT_full = (jT_iy + 1) if jT_iy >= 0 else 1
    return map_idx, k, iL_full, jT_full


def scatter_on_full_grid(phi_vec: np.ndarray, m: int, map_idx: np.ndarray) -> np.ndarray:
    Z = np.zeros((m, m), dtype=float)
    it = np.nditer(map_idx, flags=['multi_index'])
    while not it.finished:
        idx = int(it[0])
        iy, ix = it.multi_index
        if idx >= 0:
            Z[iy + 1, ix + 1] = phi_vec[idx]
        it.iternext()
    return Z

# Prépare X,Y et “découpe” en deux blocs surface (comme la stationnaire)
def prepare_surfaces_geometry(m, alpha, iL_full, jT_full, L):
    X = np.linspace(0.0, 2.0/2.0 * L * alpha, m)
    Y = np.linspace(0.0, 2.0/2.0 * L * alpha, m)
    Xm, Ym = np.meshgrid(X, Y)
    # Bloc gauche (incluant frontière verticale interne)
    sl_left = np.s_[:, :iL_full + 1]
    # Bloc haut-droit (incluant frontière horizontale)
    sl_top  = np.s_[jT_full:, iL_full:]
    return Xm, Ym, sl_left, sl_top

# ---------- Fenêtre 1 : 3D stationnaire ----------
def make_stationary_plot_3d(m, alpha, phi_vec, map_idx, iL_full, jT_full, L):
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    phi = np.abs(phi_vec.copy())
    mx = phi.max() if phi.size else 0.0
    if mx > 0:
        phi /= mx

    Z = scatter_on_full_grid(phi, m, map_idx)
    Xm, Ym, sl_left, sl_top = prepare_surfaces_geometry(m, alpha, iL_full, jT_full, L)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    s1 = ax.plot_surface(Xm[sl_left], Ym[sl_left], Z[sl_left],
                         linewidth=0, antialiased=False, shade=False,
                         cmap="viridis", vmin=0.0, vmax=1.0)
    s2 = ax.plot_surface(Xm[sl_top],  Ym[sl_top],  Z[sl_top],
                         linewidth=0, antialiased=False, shade=False,
                         cmap="viridis", vmin=0.0, vmax=1.0)

    ax.plot3D([(1.5 / 2) * L * alpha, (1.5 / 2) * L * alpha], [0.0, (0.75 / 2) * L * alpha], [0.0, 0.0],
              color="black", linewidth=2)
    ax.plot3D([(1.5 / 2) * L * alpha, (2.0 / 2.0) * L * alpha], [(0.75 / 2) * L * alpha, (0.75 / 2) * L * alpha],
              [0.0, 0.0],
              color="black", linewidth=2)

    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("ϕ normalisé (|ϕ|/max)")
    ax.set_title(f"Flux stationnaire — mode fondamental (m={m})")
    ax.set_zlim(0.0, 1.0)
    ax.set_box_aspect((2/2 * L * alpha, 2/2 * L * alpha, 1.0))

    fig.colorbar(s2, shrink=0.82, aspect=22).set_label("ϕ normalisé (|ϕ|/max)")

    plt.show(block=False)  # on laisse la fenêtre ouverte
    return fig, ax, Xm, Ym, (sl_left, sl_top)

# ---------- Fenêtre 2 : 3D temporelle ----------
def live_temporal_plot_3d(m, alpha, map_idx, geom, L,pending_header=None):
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    Xm, Ym, (sl_left, sl_top) = geom
    n_expected = (map_idx >= 0).sum()

    plt.ion()
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    zmax_t = 1e-12 # Echelle qui bougera avec le temps
    def draw_frame(phi_vec, step, t):
        nonlocal zmax_t
        phi = np.abs(phi_vec)
        Z = scatter_on_full_grid(phi, m, map_idx)
        zcur = float(Z.max()) if Z.size else 0.0
        if zcur > zmax_t * 1.02:  # augmente l'échelle si on dépasse de 2%
            zmax_t = zcur

        ax.clear()
        s1 = ax.plot_surface(Xm[sl_left], Ym[sl_left], Z[sl_left],
                             linewidth=0, antialiased=False, shade=False,
                             cmap="viridis", vmin=0.0, vmax=max(zmax_t, 1e-12))
        s2 = ax.plot_surface(Xm[sl_top],  Ym[sl_top],  Z[sl_top],
                             linewidth=0, antialiased=False, shade=False,
                             cmap="viridis", vmin=0.0, vmax=max(zmax_t, 1e-12))

        ax.plot3D([(1.5 / 2) * L * alpha, (1.5 / 2) * L * alpha], [0.0, (0.75 / 2) * L * alpha], [0.0, 0.0],
                  color="black", linewidth=2)
        ax.plot3D([(1.5 / 2) * L * alpha, (2.0 / 2.0) * L * alpha], [(0.75 / 2) * L * alpha, (0.75 / 2) * L * alpha],
                  [0.0, 0.0],
                  color="black", linewidth=2)

        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("ϕ non-normalisé [m\u207b\u00b2.s\u207b\u00b9]")
        ax.set_title(f"Évolution temporelle — t={t:.2e}  (step {step})")
        ax.set_zlim(0.0, max(zmax_t, 1e-12))
        ax.set_box_aspect((2/2 * L * alpha, 2/2 * L * alpha, max(zmax_t, 1e-12)))

        # pour éviter d'empiler les colorbars à chaque frame, on n'en met pas ici
        plt.pause(0.03)

    # Si on a déjà lu un header, traite-le d'abord
    if pending_header is not None:
        parts = pending_header.strip().split()
        step0 = int(parts[2]) if len(parts) >= 3 else 0
        t0 = 0.0
        if len(parts) >= 4 and '=' in parts[3]:
            try: t0 = float(parts[3].split('=')[1])
            except Exception: t0 = 0.0
        vec_line = sys.stdin.readline()
        phi0 = np.fromstring(vec_line, sep=' ')
        if phi0.size == n_expected:
            draw_frame(phi0, step0, t0)

    # Boucle temps réel
    for line in sys.stdin:
        if not line:
            break
        if not line.startswith("# step"):
            # ignorer toute autre ligne (au cas où)
            continue

        parts = line.strip().split()
        try:
            step = int(parts[2])
        except Exception:
            step = 0
        t = 0.0
        if len(parts) >= 4 and '=' in parts[3]:
            try: t = float(parts[3].split('=')[1])
            except Exception: t = 0.0

        vec_line = sys.stdin.readline()
        if not vec_line:
            break
        phi = np.fromstring(vec_line, sep=' ')
        if phi.size != n_expected:
            print(f"⚠️  taille reçue {phi.size} ≠ n attendu {n_expected}", file=sys.stderr)
            continue

        draw_frame(phi, step, t)

    plt.ioff()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualisation 3D (stationnaire + temporelle) de ϕ.")
    parser.add_argument("--L", type=float, required=True)
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--alpha", type=float, default=math.sqrt(5.85 / 2.0))
    args = parser.parse_args()

    m = args.m
    L = args.L
    alpha = args.alpha

    # Géométrie et n attendu
    map_idx, n, iL_full, jT_full = build_map_like_c(m, alpha, L)

    # 1) Lire le vecteur stationnaire (n valeurs), avant tout header "# step"
    nums = []
    pending_header = None
    for line in sys.stdin:
        if line.startswith("# step"):
            pending_header = line
            break
        if not line.strip():
            continue
        arr = np.fromstring(line, sep=' ')
        if arr.size > 0:
            nums.append(arr)
            if sum(x.size for x in nums) >= n:
                break

    if sum(x.size for x in nums) >= n:
        phi_stationary = np.concatenate(nums)[:n]
        fig1, ax1, Xm1, Ym1, slices = make_stationary_plot_3d(
            m, alpha, phi_stationary, map_idx, iL_full, jT_full,L)
    else:
        # aucun stationnaire fourni, on affichera juste la fenêtre temporelle
        slices = prepare_surfaces_geometry(m, alpha, iL_full, jT_full, L)[:2], (np.s_[:, :iL_full+1], np.s_[jT_full:, iL_full:])

    # 2) Lancer la fenêtre 3D temporelle
    # Prépare la “géométrie” (Xm, Ym et les deux slices)
    Xm, Ym, sl_left, sl_top = prepare_surfaces_geometry(m, alpha, iL_full, jT_full, L)
    live_temporal_plot_3d(m, alpha, map_idx, (Xm, Ym, (sl_left, sl_top)),L, pending_header=pending_header)

if __name__ == "__main__":
    main()
