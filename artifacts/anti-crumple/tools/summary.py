"""Markdown summary per case: python summary.py case_label run.json [run.json ...] (frame 0 excluded from maxima)"""
import json, sys
def row(p, wall=None):
    d = json.load(open(p)); fr = d['frames'][1:]
    js = [f['jump'] for f in fr if f['jump'] >= 0]
    nd = [f['n_deep'] for f in fr]
    fail = [f['t'] for f in d['frames'] if f['state'] == 'failed']
    return '| %s | %d | %d | %.3f | %.3f | %.3f / %.2f | %.3f / %.4f | %s | %s | %s |' % (
        d['mode'], max(nd), sum(1 for x in nd if x), max(f['dent_depth'] for f in fr), sum(f['h_neg_frac'] for f in fr) / len(fr),
        100 * max(f['len_rel_rms'] for f in fr), 100 * max(f['len_rel_max'] for f in fr), max(js), sum(js) / len(js),
        ('%.1f' % wall) if wall else '%.1f*' % d['wall'], '%.3f' % fail[0] if fail else '-', d['quality'])
print('| mode | max deep | frames w/ deep | max dent | mean H<0 frac | len rms / max % (max) | jump max / mean | wall s | failed_t | quality |')
print('|---|---|---|---|---|---|---|---|---|---|')
for a in sys.argv[1:]:
    p, _, w = a.partition('@')
    print(row(p, float(w) if w else None))
