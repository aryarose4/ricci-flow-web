"""Per-frame table: python table.py run.json [...]"""
import json
import sys

for p in sys.argv[1:]:
    d = json.load(open(p))
    print(f"## {p}  wall {d['wall']:.1f}s  hidden {d['n_hidden']}  quality {d['quality']!r}")
    print(f"notice: {d['notice']}")
    print("| k | t | state | start | deep | reflex | dent | len rms % | len max % | K err % | K rms % |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for f in d["frames"]:
        print(f"| {f['k']} | {f['t']:.3g} | {f['state']} | {f['start']} | {f['n_deep']} | {f['n_reflex']} | {f['dent_depth']:.3f} | "
              f"{100*f['len_rel_rms']:.2f} | {100*f['len_rel_max']:.2f} | {100*f['k_err']:+.1f} | {100*f['k_err_rms']:.1f} |")
    print()
