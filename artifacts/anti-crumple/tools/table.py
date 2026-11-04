import json,sys
for p in sys.argv[1:]:
    d=json.load(open(p))
    print(p, d['mode'], 'wall %.1f'%d['wall'], 'hidden',d['n_hidden'], d['quality'], d.get('notice',''))
    print(' k        t   rms%   max%  refl deep  mindih  dent  dentfr  Hneg   jump state')
    for f in d['frames']:
        print('%2d %8.4f %6.3f %6.3f %5d %4d %7.3f %5.3f %6.3f %5.3f %6.4f %s'%(f['k'],f['t'],100*f['len_rel_rms'],100*f['len_rel_max'],f['n_reflex'],f['n_deep'],f['min_dih'],f['dent_depth'],f['dent_frac'],f['h_neg_frac'],f['jump'],f['state']))
