"""Side-by-side GIF for lane neg-forward with the web's Length Fit overlay.

python render.py OUT.gif TITLE label1=run1.json label2=run2.json
Faces coloured like web/src/curvature.ts 'Length Fit' (ball K of the 3D shape minus ball K of the flow
metric, over sensitivity 0.4 |rho|, clipped to [-1, 1], blue-white-red; grey = cone 1-rings); deep
folds (dihedral < -0.15) as thick black edges.  Per-panel label: t, state, len rms/max, K err, deep.
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

SENS = 0.4


def kabsch(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (X - cx) @ R.T + cy


def deep_edges(X, F):
    inc = {}
    for fi, t in enumerate(F):
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            inc.setdefault((min(a, b), max(a, b)), []).append(fi)
    n = np.cross(X[F[:, 1]] - X[F[:, 0]], X[F[:, 2]] - X[F[:, 0]])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-300
    vol = np.einsum("ij,ij->i", X[F[:, 0]], np.cross(X[F[:, 1]], X[F[:, 2]])).sum()
    s = 1.0 if vol > 0 else -1.0
    out = []
    for (a, b), fl in inc.items():
        if len(fl) != 2:
            continue
        f1, f2 = fl
        t = F[f1]
        if not any(t[k] == a and t[(k + 1) % 3] == b for k in range(3)):
            f1, f2 = f2, f1
        d = s * np.arctan2(np.dot(np.cross(n[f1], n[f2]), X[b] - X[a]), np.dot(n[f1], n[f2]))
        if d < -0.15:
            out.append((X[a], X[b]))
    return out


def load(p):
    d = json.load(open(p))
    F = np.array(d["faces"])
    fr = []
    for f in d["frames"]:
        X = np.array(f["x"]).reshape(-1, 3)
        kv = np.array([np.nan if v is None else v for v in f["kv"]]) if f["kv"] else np.full(len(X), np.nan)
        fr.append(dict(f, X=X, kv=kv))
    return d, F, fr


def main():
    out, title = sys.argv[1], sys.argv[2]
    runs = []
    for a in sys.argv[3:]:
        lab, p = a.split("=", 1)
        runs.append((lab,) + load(p))
    n = max(len(r[3]) for r in runs)
    # align every frame to the previous one of the same run (display only)
    for r in runs:
        for i in range(1, len(r[3])):
            r[3][i]["X"] = kabsch(r[3][i]["X"], r[3][i - 1]["X"])
    lim = max(np.abs(f["X"]).max() for r in runs for f in r[3])
    cmap = plt.get_cmap("RdBu_r")
    fig = plt.figure(figsize=(6 * len(runs), 6.6))
    axs = [fig.add_subplot(1, len(runs), i + 1, projection="3d") for i in range(len(runs))]

    def upd(i):
        for ax, (lab, d, F, frs) in zip(axs, runs):
            ax.clear()
            f = frs[min(i, len(frs) - 1)]
            X = f["X"]
            vt = np.clip(f["kv"] / SENS, -1, 1)
            ft = np.nanmean(np.stack([vt[F[:, 0]], vt[F[:, 1]], vt[F[:, 2]]]), axis=0) if np.isfinite(vt).any() else np.full(len(F), np.nan)
            cols = [(0.6, 0.6, 0.6, 1) if not np.isfinite(v) else cmap(0.5 + 0.5 * v) for v in ft]
            ax.add_collection3d(Poly3DCollection(X[F], facecolors=cols, edgecolors=(0, 0, 0, 0.08), linewidths=0.2))
            de = deep_edges(X, F)
            if de:
                ax.add_collection3d(Line3DCollection(de, colors="k", linewidths=2.5))
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
            ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
            ax.view_init(elev=25, azim=35)
            ax.set_title(
                f"{lab}\nt={f['t']:.3g}  {f['state']}  deep={f['n_deep']}\n"
                f"len rms {100*f['len_rel_rms']:.2f}% max {100*f['len_rel_max']:.2f}%  K err {100*f['k_err']:+.0f}%",
                fontsize=10)
        fig.suptitle(title + "  (Length Fit overlay: red/blue = 3D curvature above/below flow target; black = deep folds)", fontsize=9)
        return []

    anim = animation.FuncAnimation(fig, upd, frames=n, interval=700)
    anim.save(out, writer=animation.PillowWriter(fps=1.4))
    print("wrote", out)


if __name__ == "__main__":
    main()
