"""CPU nonbonded/contact/surface diagnostics on explicit canonical input graphs.
No coordinate repair, docking, minimization or inference is performed.
"""
from __future__ import annotations
from collections import defaultdict,deque,Counter,OrderedDict
from itertools import combinations
import hashlib,json,math,time
import numpy as np
try:
 from scipy.spatial import cKDTree
 KD_BACKEND='scipy_cKDTree'
except ImportError:
 KD_BACKEND='numpy_brute_fallback'
 class cKDTree:
  def __init__(self,data):self.data=np.asarray(data,float)
  def query_ball_tree(self,other,r):
   return [np.flatnonzero(np.sum((other.data-v)**2,axis=1)<=r*r).tolist()for v in self.data]

VERSION='glyco_physics_v1_20261113'
CONFIG={'contact_thresholds_A':[4.,4.5,5.],'contact_union_max_A':5.,'graph_neighbor_exclusion_max_bonds':3,'vdw_radii_A':{'C':1.7,'N':1.55,'O':1.52,'S':1.8,'P':1.8,'F':1.47,'CL':1.75,'BR':1.85,'I':1.98},'severe_overlap_A':.6,'sensitivity_overlap_A':.4,'cross_reference_radius_A':15.,'cross_error_thresholds_A':[.5,1.,2.,4.],'SASA_probe_A':1.4,'SASA_algorithm':'LeeRichards','SASA_n_slices':20,'SASA_threads':1,'SASA_min_exposed_area_A2':1.}
_CACHE=OrderedDict()
def _gid(k):
 p=k.split('/');return p[1]if len(p)>1 and p[0]=='glycan'else None

def _node(k):
 p=k.split('/');return p[2]if len(p)>2 and p[0]=='glycan'else None

def _pos(k):
 p=k.split('/');return int(p[2])if len(p)>2 and p[0]=='protein'else None

def _valid(x):
 try:return len(x)==3 and all(math.isfinite(float(v))for v in x)
 except Exception:return False

def _edge(a,b):return(a,b)if a<b else(b,a)
def _elements(ref,pred=None):
 out={k:str(v['element']).upper()for k,v in ref.get('expected_atom_meta',{}).items()}
 if pred:
  for k,v in pred.get('atom_meta',{}).items():
   if k not in out and v.get('element'):out[k]=str(v['element']).upper()
 return out

def _excluded_pairs(ref):
 graph=defaultdict(set)
 for a,b,*_ in ref['bonds']:graph[a].add(b);graph[b].add(a)
 excluded=set()
 for key in sorted(k for k in ref['expected_atom_meta']if k.startswith('glycan/')):
  seen={key};todo=deque([(key,0)])
  while todo:
   a,d=todo.popleft()
   if d==3:continue
   for b in graph[a]:
    if b not in seen:seen.add(b);todo.append((b,d+1));excluded.add(_edge(key,b))
 return excluded

def _near_gp_pairs(atoms,pkeys,gkeys,excluded,radius):
 pp=[k for k in pkeys if k in atoms and _valid(atoms[k])];gg=[k for k in gkeys if k in atoms and _valid(atoms[k])]
 if not pp or not gg:return set()
 pt=cKDTree(np.asarray([atoms[k]for k in pp],float));gt=cKDTree(np.asarray([atoms[k]for k in gg],float));out=set()
 for gi,neighbors in enumerate(gt.query_ball_tree(pt,radius)):
  g=gg[gi]
  for pi in neighbors:
   p=pp[pi]
   if _edge(g,p)not in excluded:out.add((_pos(p),_node(g),_gid(g)))
 return out

def _minimum(atoms,aa,bb,excluded=None):
 aa=[k for k in aa if k in atoms and _valid(atoms[k])];bb=[k for k in bb if k in atoms and _valid(atoms[k])]
 if not aa or not bb:return None,None
 xa=np.asarray([atoms[k]for k in aa]);xb=np.asarray([atoms[k]for k in bb]);sq=np.sum((xa[:,None,:]-xb[None,:,:])**2,axis=2)
 if excluded:
  for i,a in enumerate(aa):
   for j,b in enumerate(bb):
    if _edge(a,b)in excluded:sq[i,j]=np.inf
 idx=np.unravel_index(int(np.argmin(sq)),sq.shape);d=float(np.sqrt(sq[idx]))
 return(None,None)if not math.isfinite(d)else(d,[aa[idx[0]],bb[idx[1]]])

def _prepare(ref):
 token=(ref['case_id'],ref.get('reference_source',{}).get('task_id'),ref.get('reference_source',{}).get('reference_only_sha256'))
 if token in _CACHE:_CACHE.move_to_end(token);return _CACHE[token]
 ra=ref['atoms'];pkeys=sorted(k for k in ra if k.startswith('protein/'));gkeys=sorted(k for k in ra if k.startswith('glycan/'));excluded=_excluded_pairs(ref);res=defaultdict(list);node=defaultdict(list)
 for k in pkeys:res[_pos(k)].append(k)
 for k in gkeys:node[_node(k)].append(k)
 near=_near_gp_pairs(ra,pkeys,gkeys,excluded,5.)
 minima={pair:_minimum(ra,node[pair[1]],res[pair[0]],excluded)for pair in near}
 cross=[]
 if pkeys and gkeys:
  pp=np.asarray([ra[k]for k in pkeys]);gg=np.asarray([ra[k]for k in gkeys]);pt=cKDTree(pp);gt=cKDTree(gg)
  for gi,neighbors in enumerate(gt.query_ball_tree(pt,15.)):
   g=gkeys[gi]
   for pi in sorted(neighbors):
    p=pkeys[pi]
    if _edge(g,p)not in excluded:cross.append((g,p,float(np.linalg.norm(gg[gi]-pp[pi]))))
 cross.sort(key=lambda z:(z[0],z[1]));x={'protein_keys':pkeys,'glycan_keys':gkeys,'protein_residues':dict(res),'node_keys':dict(node),'excluded':excluded,'reference_contact_minima':minima,'cross_pairs':cross,'reference_sasa':None};_CACHE[token]=x
 while len(_CACHE)>8:_CACHE.popitem(last=False)
 return x

def _prf(rows,t):
 ref={(r['protein_residue_position'],r['node_id'])for r in rows if r['reference_min_distance_A']is not None and r['reference_min_distance_A']<=t};pred={(r['protein_residue_position'],r['node_id'])for r in rows if r['prediction_min_distance_A']is not None and r['prediction_min_distance_A']<=t};tp=len(ref&pred);fp=len(pred-ref);fn=len(ref-pred)
 return {'reference_contacts':len(ref),'prediction_contacts':len(pred),'true_positive':tp,'false_positive':fp,'false_negative':fn,'precision':tp/len(pred)if pred else None,'recall':tp/len(ref)if ref else None,'F1':2*tp/(2*tp+fp+fn)if ref else None,'status':'ok'if ref else'not_applicable_no_reference_nonbonded_contacts','no_reference_contact_prediction_also_empty':not pred if not ref else None}

def contacts(ref,pred,cache):
 ra,pa=ref['atoms'],pred['atoms'];union=set(cache['reference_contact_minima'])|_near_gp_pairs(pa,cache['protein_keys'],cache['glycan_keys'],cache['excluded'],5.);rows=[]
 for pos,nid,gid in sorted(union):
  rd,rnear=cache['reference_contact_minima'].get((pos,nid,gid),(None,None))
  if rd is None:rd,rnear=_minimum(ra,cache['node_keys'][nid],cache['protein_residues'][pos],cache['excluded'])
  pd,pnear=_minimum(pa,cache['node_keys'][nid],cache['protein_residues'][pos],cache['excluded'])
  rows.append({'protein_residue_position':pos,'node_id':nid,'glycan_id':gid,'reference_min_distance_A':rd,'prediction_min_distance_A':pd,'reference_nearest_atom_keys':rnear,'prediction_nearest_atom_keys':pnear,'excluded_covalent_graph_distance_leq':3})
 sites={}
 for gid,s in ref['sites'].items():
  rr=[r for r in rows if r['glycan_id']==gid];nonanchor=[r for r in rr if r['protein_residue_position']!=int(s['position'])]
  sites[gid]={'by_threshold_A':{str(t):_prf(rr,t)for t in CONFIG['contact_thresholds_A']},'exclude_attachment_residue_by_threshold_A':{str(t):_prf(nonanchor,t)for t in CONFIG['contact_thresholds_A']},'union_rows':len(rr)}
 return {'status':'ok','score_mask':'reference-observed common protein and glycan heavy atoms','distance_union_cutoff_A':5.,'rows':rows,'sites':sites,'prediction_missing_reference_mask_atoms':[k for k in cache['protein_keys']+cache['glycan_keys']if k not in pa],'graph_neighbor_exclusion_max_bonds':3}

def clashes(ref,pred,cache):
 pa=pred['atoms'];els=_elements(ref,pred);allowed=set(ref['expected_atom_meta'])|{f"protein/P1/{len(ref['protein_sequence'])}/OXT"};keys=sorted(k for k in pa if k in allowed and _valid(pa[k]));gkeys=[k for k in keys if k.startswith('glycan/')];unsupported=[k for k in keys if els.get(k)not in CONFIG['vdw_radii_A']]
 if unsupported:return {'status':'unsupported_element_radii','atom_keys':unsupported}
 if not keys or not gkeys:return {'status':'not_applicable_empty_glycan'}
 x=np.asarray([pa[k]for k in keys]);gx=np.asarray([pa[k]for k in gkeys]);tree=cKDTree(x);gt=cKDTree(gx);seen=set();rows=[];maxdistance=2*max(CONFIG['vdw_radii_A'].values())
 for gi,neighbors in enumerate(gt.query_ball_tree(tree,maxdistance)):
  a=gkeys[gi]
  for j in neighbors:
   b=keys[j];pair=_edge(a,b)
   if a==b or pair in seen or pair in cache['excluded']:continue
   seen.add(pair);d=float(np.linalg.norm(gx[gi]-x[j]));r=CONFIG['vdw_radii_A'][els[a]]+CONFIG['vdw_radii_A'][els[b]];ov=r-d
   if ov<=0:continue
   gb=_gid(b);kind='glycan_protein'if gb is None else'within_glycan'if _gid(a)==gb else'between_glycans';gids=sorted({_gid(a)}|({gb}if gb else set()))
   rows.append({'atom_1':a,'atom_2':b,'distance_A':d,'vdw_radii_sum_A':r,'overlap_A':ov,'kind':kind,'glycan_ids':gids,'severe_gt_0_6_A':ov>.6,'sensitivity_gt_0_4_A':ov>.4})
 sites={gid:{'positive_overlap_pairs':0,'severe_gt_0_6_A':0,'sensitivity_gt_0_4_A':0,'maximum_overlap_A':0.}for gid in ref['sites']};by=defaultdict(Counter)
 for r in rows:
  by[r['kind']]['positive_overlap_pairs']+=1;by[r['kind']]['severe_gt_0_6_A']+=r['severe_gt_0_6_A'];by[r['kind']]['sensitivity_gt_0_4_A']+=r['sensitivity_gt_0_4_A']
  for gid in r['glycan_ids']:
   sites[gid]['positive_overlap_pairs']+=1;sites[gid]['severe_gt_0_6_A']+=r['severe_gt_0_6_A'];sites[gid]['sensitivity_gt_0_4_A']+=r['sensitivity_gt_0_4_A'];sites[gid]['maximum_overlap_A']=max(sites[gid]['maximum_overlap_A'],r['overlap_A'])
 severe=sum(r['severe_gt_0_6_A']for r in rows)
 return {'status':'ok','scope':'all predicted expected heavy atoms involving >=1 glycan; protein-protein clashes not evaluated','graph_exclusion_max_bonds':3,'radii_A':CONFIG['vdw_radii_A'],'pairs':rows,'sites':sites,'by_kind':{k:dict(v)for k,v in by.items()},'severe_overlap_pair_count':severe,'severe_pairs_per_1000_glycan_atoms':severe*1000/len(gkeys),'severe_pairs_per_1000_total_heavy_atoms':severe*1000/len(keys),'glycan_atom_denominator':len(gkeys),'total_atom_denominator':len(keys)}

def cross_distance(ref,pred,cache):
 pa=pred['atoms'];groups=defaultdict(list)
 for g,p,rd in cache['cross_pairs']:
  err=None if g not in pa or p not in pa or not _valid(pa[g])or not _valid(pa[p])else abs(float(np.linalg.norm(np.asarray(pa[g])-np.asarray(pa[p])))-rd)
  groups[_node(g)].append((g,p,rd,err))
 nodes={};sites={}
 for nid,node in ref['nodes'].items():
  pairs=groups[nid];errs=[r[3]for r in pairs];counts=[sum(e is not None and e<t for e in errs)for t in CONFIG['cross_error_thresholds_A']];n=len(errs);score=sum(counts)/(4*n)if n else None;sig=hashlib.sha256(json.dumps([[g,p,rd]for g,p,rd,_ in pairs],separators=(',',':')).encode()).hexdigest()
  nodes[nid]={'glycan_id':node['glycan_id'],'status':'ok'if n else'not_applicable_no_reference_pairs','pair_count':n,'missing_prediction_pair_count':sum(e is None for e in errs),'absolute_distance_errors_A':errs,'preserved_counts':counts,'score':score,'reference_pair_order_sha256':sig}
 for gid,s in ref['sites'].items():
  ns=[nodes[n]for n in s['node_ids']if nodes[n]['score']is not None];sites[gid]={'node_macro_score':float(np.mean([n['score']for n in ns]))if ns else None,'nodes_scored':len(ns),'nodes_total':len(s['node_ids']),'pair_micro_score':sum(sum(n['preserved_counts'])for n in ns)/(4*sum(n['pair_count']for n in ns))if ns else None,'pair_count':sum(n['pair_count']for n in ns)}
 return {'status':'ok','name':'custom_GP_distance_preservation_not_standard_lDDT','reference_radius_A':15.,'thresholds_A':CONFIG['cross_error_thresholds_A'],'graph_exclusion_max_bonds':3,'nodes':nodes,'sites':sites,'reference_pair_order':'lexical glycan_atom then protein_atom; reconstructed from immutable reference and cutoff/exclusion registry','missing_prediction_pairs_count_as_not_preserved':True}

def multisite(ref,pred,cache):
 ra,pa=ref['atoms'],pred['atoms'];rows=[]
 for ga,gb in combinations(sorted(ref['sites']),2):
  sa,sb=ref['sites'][ga],ref['sites'][gb];aa=sa['attachment'];ba=sb['attachment'];ka=f"glycan/{ga}/{aa['glycan_node_id']}/{aa['glycan_atom']}";kb=f"glycan/{gb}/{ba['glycan_node_id']}/{ba['glycan_atom']}";terminal_a=[k for n in sa['node_ids']if not ref['nodes'][n]['children']for k in cache['node_keys'][n]];terminal_b=[k for n in sb['node_ids']if not ref['nodes'][n]['children']for k in cache['node_keys'][n]]
  rd,_=_minimum(ra,[ka],[kb]);pd,_=_minimum(pa,[ka],[kb]);rt,rta=_minimum(ra,terminal_a,terminal_b);pt,pta=_minimum(pa,terminal_a,terminal_b);rows.append({'glycan_1':ga,'glycan_2':gb,'root_reference_distance_A':rd,'root_prediction_distance_A':pd,'root_absolute_error_A':abs(pd-rd)if pd is not None and rd is not None else None,'terminal_reference_min_distance_A':rt,'terminal_prediction_min_distance_A':pt,'terminal_absolute_error_A':abs(pt-rt)if pt is not None and rt is not None else None,'reference_nearest_terminal_atoms':rta,'prediction_nearest_terminal_atoms':pta})
 node_contacts=[]
 for ga,gb in combinations(sorted(ref['sites']),2):
  for na in ref['sites'][ga]['node_ids']:
   for nb in ref['sites'][gb]['node_ids']:
    rd,_=_minimum(ra,cache['node_keys'][na],cache['node_keys'][nb],cache['excluded']);pd,_=_minimum(pa,cache['node_keys'][na],cache['node_keys'][nb],cache['excluded'])
    if(rd is not None and rd<=5.)or(pd is not None and pd<=5.):node_contacts.append({'glycan_1':ga,'glycan_2':gb,'node_1':na,'node_2':nb,'reference_min_distance_A':rd,'prediction_min_distance_A':pd})
 metrics={}
 for t in CONFIG['contact_thresholds_A']:
  rs={(v['node_1'],v['node_2'])for v in node_contacts if v['reference_min_distance_A']is not None and v['reference_min_distance_A']<=t};ps={(v['node_1'],v['node_2'])for v in node_contacts if v['prediction_min_distance_A']is not None and v['prediction_min_distance_A']<=t};tp=len(rs&ps);fp=len(ps-rs);fn=len(rs-ps);metrics[str(t)]={'reference_contacts':len(rs),'prediction_contacts':len(ps),'TP':tp,'FP':fp,'FN':fn,'precision':tp/len(ps)if ps else None,'recall':tp/len(rs)if rs else None,'F1':2*tp/(2*tp+fp+fn)if rs else None}
 return {'status':'ok'if rows else'not_applicable_single_site','site_pairs':rows,'glycan_node_contact_union':node_contacts,'glycan_node_contact_metrics_by_threshold_A':metrics}

def _sasa(atoms,pkeys,gkeys,elements):
 import freesasa
 keys=pkeys+gkeys;missing=[k for k in keys if k not in atoms or not _valid(atoms[k])]
 if missing:return {'status':'model_or_reference_missing_mask_atoms','missing_atoms':missing}
 bad=[k for k in keys if elements.get(k)not in CONFIG['vdw_radii_A']]
 if bad:return {'status':'unsupported_element_radii','atom_keys':bad}
 xyz=np.asarray([atoms[k]for k in keys],float);unique=np.unique(xyz,axis=0)
 if len(unique)!=len(xyz):return {'status':'not_assessed_duplicate_coordinates','duplicate_coordinate_count':len(xyz)-len(unique)}
 params=freesasa.Parameters({'algorithm':freesasa.LeeRichards,'probe-radius':1.4,'n-slices':20,'n-threads':1});radii=[CONFIG['vdw_radii_A'][elements[k]]for k in keys];full=freesasa.calcCoord(xyz.reshape(-1).tolist(),radii,params);plain=freesasa.calcCoord(xyz[:len(pkeys)].reshape(-1).tolist(),radii[:len(pkeys)],params);res=defaultdict(lambda:[0.,0.])
 for i,k in enumerate(pkeys):res[_pos(k)][0]+=float(plain.atomArea(i));res[_pos(k)][1]+=float(full.atomArea(i))
 out={}
 for pos,(p,pg)in sorted(res.items()):
  delta=p-pg;out[str(pos)]={'protein_only_fixed_coordinates_A2':p,'protein_with_given_glycans_A2':pg,'delta_A2':delta,'shielding_fraction_raw':delta/p if p>=1. else None,'exposed_for_relative_score':p>=1.}
 return {'status':'ok','protein_atom_count':len(pkeys),'glycan_atom_count':len(gkeys),'per_residue':out,'protein_SASA_P_A2':sum(v[0]for v in res.values()),'protein_SASA_P_plus_G_A2':sum(v[1]for v in res.values())}

def surface(ref,pred,cache):
 try:import freesasa
 except ImportError:return {'status':'unsupported_missing_freesasa','reference':None,'prediction':None}
 if cache['reference_sasa']is None:cache['reference_sasa']=_sasa(ref['atoms'],cache['protein_keys'],cache['glycan_keys'],_elements(ref))
 rr=cache['reference_sasa'];pp=_sasa(pred['atoms'],cache['protein_keys'],cache['glycan_keys'],_elements(ref,pred));diff=[]
 if rr['status']=='ok'and pp['status']=='ok':
  for pos,r in rr['per_residue'].items():
   p=pp['per_residue'][pos];a=r['shielding_fraction_raw'];b=p['shielding_fraction_raw'];diff.append({'position':int(pos),'reference_relative':a,'prediction_relative':b,'relative_absolute_error':abs(a-b)if a is not None and b is not None else None,'reference_delta_A2':r['delta_A2'],'prediction_delta_A2':p['delta_A2'],'delta_absolute_error_A2':abs(r['delta_A2']-p['delta_A2'])})
 vals=[v['relative_absolute_error']for v in diff if v['relative_absolute_error']is not None]
 return {'status':'ok'if rr['status']=='ok'and pp['status']=='ok'else'not_assessed_subitem','mask':'fixed reference observed common protein+given glycan heavy atoms; no protein-only prediction','probe_radius_A':1.4,'radii_A':CONFIG['vdw_radii_A'],'algorithm':'FreeSASA LeeRichards 20 slices 1 thread','minimum_exposed_protein_SASA_A2':1.,'reference':rr,'prediction':pp,'comparison':{'relative_MAE':float(np.mean(vals))if vals else None,'relative_residues_scored':len(vals),'relative_residues_total':len(cache['protein_residues']),'per_residue':diff},'scope':'static resolved-fragment surface geometry, not natural complete glycan shield or functional evidence'}

def evaluate_physics(ref,pred):
 start=time.monotonic();out={'metric_version':VERSION,'config':CONFIG,'coordinates_modified':False,'neighbor_search_backend':KD_BACKEND};cache=_prepare(ref)
 for name,fn in [('contacts',contacts),('clashes',clashes),('cross_distance',cross_distance),('multisite',multisite),('surface',surface)]:
  t=time.monotonic()
  try:out[name]=fn(ref,pred,cache)
  except Exception as e:out[name]={'status':'error','error':type(e).__name__+': '+str(e)}
  out[name]['elapsed_seconds']=time.monotonic()-t
 out['elapsed_seconds']=time.monotonic()-start;out['status']='complete_with_subitem_errors'if any(v.get('status')=='error'for v in [out[n]for n in ['contacts','clashes','cross_distance','multisite','surface']])else'complete'
 return out
