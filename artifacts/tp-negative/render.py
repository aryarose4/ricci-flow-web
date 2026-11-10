"""Frame PNGs + GIF from crumple-harness JSON (rho<0 True Polytope).
python artifacts/tp-negative/render.py artifacts/tp-negative NAME"""
import json, sys, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
d, name = sys.argv[1], sys.argv[2]
J = json.load(open(os.path.join(d, name + ".json")))
F = np.array(J["faces"]); cones = J["cones"]
L = np.array([0.4, 0.3, 0.85]); L /= np.linalg.norm(L)
def img(fr):
    X = np.array(fr["x"]).reshape(-1, 3); X = X - X.mean(0); X /= np.abs(X).max()
    T = X[F]; n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-300
    sh = 0.3 + 0.7 * np.abs(n @ L)
    fig = plt.figure(figsize=(4.5, 4.5), dpi=100); ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(T, facecolors=np.stack([.55*sh, .7*sh, .95*sh, np.ones_like(sh)], 1), edgecolors=(0, 0, 0, .15), linewidths=.2))
    ax.scatter(*X[cones].T, color="gold", s=25, depthshade=False, edgecolors="k")
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1); ax.set_axis_off(); ax.view_init(25, 35)
    ax.set_title(f"{name} k={fr['k']} t={fr['t']:.3g} {fr['state']}\nlen rms {fr['len_rel_rms']:.1e}", fontsize=8)
    fig.canvas.draw(); im = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3]); plt.close(fig); return im
ims = [img(f) for f in J["frames"]]
for k in [0, 5, 10, len(ims) - 1]:
    ims[k].save(os.path.join(d, f"{name}_k{k}.png"))
ims[0].save(os.path.join(d, name + ".gif"), save_all=True, append_images=ims[1:], duration=350, loop=0)
