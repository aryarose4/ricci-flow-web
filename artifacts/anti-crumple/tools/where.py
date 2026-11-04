import json,sys,numpy as np
sys.path.insert(0,'tools')
from render import dihedral_edges
d=json.load(open(sys.argv[1])); F0=np.array(d['faces']); cones=set(d['cones'])
nb={}
for t in F0:
    for k in range(3): nb.setdefault(t[k],set()).update(t)
ring1=set().union(*[nb[c] for c in cones])
for k in map(int,sys.argv[2:]):
    f=d['frames'][k]; X=np.array(f['x']).reshape(-1,3); F=np.array(f['faces']) if f['faces'] else F0
    deep,_=dihedral_edges(X,F)
    print(k, 'deep',len(deep),'cone-incident',sum(a in cones or b in cones for a,b in deep),'in cone 1-ring',sum(a in ring1 and b in ring1 for a,b in deep))
