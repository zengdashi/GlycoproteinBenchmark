"""The explicit four-atom registry from the frozen reference compiler."""
from collections import defaultdict

def pkey(position, atom):
    return f"protein/P1/{position}/{atom}"

def integer(value):
    try:
        return int(value)
    except (ValueError, TypeError):
        return None

def registered_torsions(bonds, bond_kinds, expected_meta, nodes):
    torsions, angles = [], []
    graph=defaultdict(set)
    for a,b,o in bonds:graph[a].add(b);graph[b].add(a)
    def elem(k):return expected_meta.get(k,{}).get('element')
    for bi,(a,b,o)in enumerate(bonds):
     if bond_kinds[bi]not in ['attachment','glycosidic']:continue
     ak,bk=a.split('/'),b.split('/');ends=[a,b];don=[k for k in ends if k.startswith('glycan/')and k.split('/')[-1]in ['C1','C2']]
     record={'bond_index':bi,'kind':bond_kinds[bi],'bond_atoms':[a,b],'phi':None,'psi':None,'omega':None,'status':'unsupported_definition'}
     # A donor carbon has the glycan ring oxygen adjacent and a non-ring exocyclic partner.
     donors=[k for k in don if any(n in nodes[k.split('/')[2]]['ring_atom_keys']and elem(n)=='O'for n in graph[k])]
     if len(donors)==1:
      dc=donors[0];accept=b if dc==a else a;no=nodes[dc.split('/')[2]];ros=[k for k in graph[dc]if elem(k)=='O'and k in no['ring_atom_keys']];ac=[k for k in graph[accept]if k!=dc and elem(k)=='C']
      if len(ros)==1 and len(ac)==1:
       ca=ac[0];record['phi']=[ros[0],dc,accept,ca];an=accept.split('/')[-1]
       if accept.startswith('protein/'):
        prev=pkey(int(accept.split('/')[2]),'CB'if an=='ND2'else'CA');record['psi']=[dc,accept,ca,prev]if prev in graph[ca]else None
       else:
        cname=ca.split('/')[-1];number=integer(cname[1:])if cname.startswith('C')else None
        ns=[k for k in graph[ca]if k!=accept and elem(k)=='C'];preferred=[k for k in ns if number is not None and k.split('/')[-1]=='C'+str(number-1)]
        prev=preferred[0]if len(preferred)==1 else ns[0]if len(ns)==1 else None
        if prev:record['psi']=[dc,accept,ca,prev]
        if ca not in nodes[ca.split('/')[2]]['ring_atom_keys']and prev:
         ringos=[k for k in graph[prev]if elem(k)=='O'and k in nodes[prev.split('/')[2]]['ring_atom_keys']]
         if len(ringos)==1:record['omega']=[accept,ca,prev,ringos[0]]
       record['status']='registered'if record['phi']and record['psi']else'partial_definition'
     torsions.append(record)
     if record['phi']:
      quad=record['phi'];angles.append({'kind':record['kind'],'atoms':[quad[1],quad[2],quad[3]],'source':'adjacent atom graph around explicit linking heteroatom'})
    return torsions
