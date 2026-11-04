"""Render GIFs (3 seeds side by side) and frame-0 crease PNGs from true_polytope_compare JSON.
python artifacts/true-polytope/render.py artifacts/true-polytope/data artifacts/true-polytope"""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL import Image

data, out = sys.argv[1], sys.argv[2]
SEEDS = [("true_polytope", "True Polytope"), ("polytope", "Polytope"), ("gluing", "Cone Gluing")]
LIGHT = np.array([0.4, 0.3, 0.85]); LIGHT /= np.linalg.norm(LIGHT)


def draw(ax, X, F, cones, crease=None, title=""):
    X = X - X.mean(0)
    s = np.abs(X).max()
    X = X / s
    T = X[F]
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-300
    sh = 0.35 + 0.65 * np.clip(n @ LIGHT, 0, 1)
    col = np.stack([0.55 * sh, 0.7 * sh, 0.95 * sh, np.ones_like(sh)], 1)
    pc = Poly3DCollection(T, facecolors=col, edgecolors=(0, 0, 0, 0.12), linewidths=0.2)
    ax.add_collection3d(pc)
    if crease is not None and len(crease):
        ax.add_collection3d(Line3DCollection(X[np.array(crease)], colors="crimson", linewidths=1.2))
    ax.scatter(*X[cones].T, color="gold", s=25, depthshade=False, edgecolors="k")
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=22, azim=35); ax.set_axis_off()
    ax.set_title(title, fontsize=8)


def frame_title(name, fr):
    return f"{name}  k={fr['k']} t={fr['t']:.3g}\nlen_rms={fr['len_rel_rms']:.1e} deep={fr['n_deep']} reflex={fr['n_reflex']}"


configs = sorted({f.rsplit("_", 1)[0].replace("_true", "") for f in os.listdir(data)})
configs = sorted({f[: f.index("_" + s)] for f in os.listdir(data) for s, _ in SEEDS if f.endswith("_" + s + ".json") and not (s == "polytope" and f.endswith("true_polytope.json"))})
paths = []
for c in configs:
    runs = {}
    for s, _ in SEEDS:
        p = os.path.join(data, f"{c}_{s}.json")
        if os.path.exists(p):
            runs[s] = json.load(open(p))
    nfr = max(len(r["frames"]) for r in runs.values())
    imgs = []
    for k in range(nfr):
        fig = plt.figure(figsize=(10.5, 3.9), dpi=80)
        for i, (s, lab) in enumerate(SEEDS):
            ax = fig.add_subplot(1, 3, i + 1, projection="3d")
            if s not in runs:
                ax.set_axis_off(); ax.set_title(lab + ": failed", fontsize=8); continue
            r = runs[s]
            fr = r["frames"][min(k, len(r["frames"]) - 1)]
            F = np.array(fr["faces"] if fr["faces"] is not None else r["faces"])
            draw(ax, np.array(fr["vertices"]), F, r["cones"], r["crease"] if (k == 0 and s == "true_polytope") else None, frame_title(lab, fr))
        fig.suptitle(c, fontsize=10)
        fig.canvas.draw()
        imgs.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3]))
        if k == 0:
            png = os.path.join(out, f"{c}_frame0.png")
            fig.savefig(png); paths.append(png)
        plt.close(fig)
    gif = os.path.join(out, f"{c}.gif")
    imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=[1200] + [350] * (len(imgs) - 2) + [1500], loop=0)
    paths.append(gif)
print("\n".join(os.path.abspath(p) for p in paths))
