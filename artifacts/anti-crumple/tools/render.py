"""Side-by-side GIF + worst-frame PNG strip for anti-crumple runs.

python render.py OUT_PREFIX TITLE label1=run1.json label2=run2.json [...]
Faces shaded by a fixed light; deep folds (dihedral < -0.15) drawn as red edges; inward-dent
vertices (> 1% bbox diag inside the convex hull) not drawn separately (see table).
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


def kabsch(X, Y):
    """rigid motion of X best matching Y"""
    cx, cy = X.mean(0), Y.mean(0)
    H = (X - cx).T @ (Y - cy)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return (X - cx) @ R.T + cy


def dihedral_edges(X, F):
    inc = {}
    for fi, t in enumerate(F):
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            inc.setdefault((min(a, b), max(a, b)), []).append(fi)
    n = np.cross(X[F[:, 1]] - X[F[:, 0]], X[F[:, 2]] - X[F[:, 0]])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    vol = np.einsum("ij,ij->i", X[F[:, 0]], np.cross(X[F[:, 1]], X[F[:, 2]])).sum()
    s = 1.0 if vol > 0 else -1.0
    deep = []
    for (a, b), fl in inc.items():
        if len(fl) != 2:
            continue
        f1, f2 = fl
        t = F[f1]
        trav = any(t[k] == a and t[(k + 1) % 3] == b for k in range(3))
        if not trav:
            f1, f2 = f2, f1
        d = s * np.arctan2(np.dot(np.cross(n[f1], n[f2]), X[b] - X[a]), np.dot(n[f1], n[f2]))
        if d < -0.15:
            deep.append((a, b))
    return deep, n


def load(p):
    d = json.load(open(p))
    F0 = np.array(d["faces"])
    fr = []
    for f in d["frames"]:
        X = np.array(f["x"]).reshape(-1, 3)
        F = np.array(f["faces"]) if f["faces"] is not None else F0
        fr.append(dict(f, X=X, F=F))
    return d, fr


LIGHT = np.array([0.4, -0.5, 0.75])
LIGHT /= np.linalg.norm(LIGHT)


def draw(ax, X, F, title, lim):
    ax.cla()
    deep, n = dihedral_edges(X, F)
    sh = 0.35 + 0.6 * np.clip(n @ LIGHT, 0, 1)
    cols = np.stack([0.55 * sh, 0.7 * sh, 0.95 * sh, np.ones_like(sh)], 1)
    pc = Poly3DCollection(X[F], facecolors=cols, edgecolors=(0, 0, 0, 0.08), linewidths=0.2)
    ax.add_collection3d(pc)
    if deep:
        ax.add_collection3d(Line3DCollection([X[[a, b]] for a, b in deep], colors="red", linewidths=2.0))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=-60)
    ax.set_axis_off()
    ax.set_title(title, fontsize=9)


def main():
    out, title = sys.argv[1], sys.argv[2]
    runs = []
    for a in sys.argv[3:]:
        lab, p = a.split("=", 1)
        d, fr = load(p)
        runs.append((lab, d, fr))
    # common frame: align every run's frames to the first run's final frame (same mesh)
    ref = runs[0][2][-1]["X"]
    ref = ref - ref.mean(0)
    for _, _, fr in runs:
        R = kabsch(fr[-1]["X"], ref)
        # chain: same rigid motion for the whole run (frames are already chained by Kabsch)
        cx = fr[-1]["X"].mean(0)
        H = (fr[-1]["X"] - cx).T @ ref
        U, S, Vt = np.linalg.svd(H)
        dd = np.sign(np.linalg.det(Vt.T @ U.T))
        Rm = Vt.T @ np.diag([1, 1, dd]) @ U.T
        for f in fr:
            f["X"] = (f["X"] - f["X"].mean(0)) @ Rm.T
    lim = max(np.abs(f["X"]).max() for _, _, fr in runs for f in fr) * 0.8
    nfr = len(runs[0][2])
    fig = plt.figure(figsize=(4.2 * len(runs), 4.6))
    axes = [fig.add_subplot(1, len(runs), i + 1, projection="3d") for i in range(len(runs))]
    sup = fig.suptitle("", fontsize=10)

    def lab(lab_, f):
        return f"{lab_}  t={f['t']:.4f}\ndeep {f['n_deep']}  dent {f['dent_depth']:.3f}  rms {100*f['len_rel_rms']:.2f}%  jump {f['jump']:.3f}"

    def upd(i):
        for ax, (lb, d, fr) in zip(axes, runs):
            draw(ax, fr[i]["X"], fr[i]["F"], lab(lb, fr[i]), lim)
        sup.set_text(f"{title} — frame {i}/{nfr-1}")
        return []

    anim = animation.FuncAnimation(fig, upd, frames=nfr, interval=600)
    anim.save(out + ".gif", writer=animation.PillowWriter(fps=2), dpi=70)
    plt.close(fig)
    # worst frames by the first run: deep folds, then dent, then jump
    fr0 = runs[0][2]
    score = [(f["n_deep"], f["dent_depth"], f["jump"]) for f in fr0]
    worst = sorted(range(1, nfr), key=lambda i: score[i], reverse=True)[:4]
    worst.sort()
    fig = plt.figure(figsize=(3.6 * len(worst), 3.8 * len(runs)))
    for r, (lb, d, fr) in enumerate(runs):
        for c, i in enumerate(worst):
            ax = fig.add_subplot(len(runs), len(worst), r * len(worst) + c + 1, projection="3d")
            draw(ax, fr[i]["X"], fr[i]["F"], lab(lb, fr[i]), lim)
    fig.suptitle(f"{title} — worst frames of {runs[0][0]}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out + "_worst.png", dpi=80)
    plt.close(fig)
    print(out + ".gif", out + "_worst.png")


if __name__ == "__main__":
    main()
