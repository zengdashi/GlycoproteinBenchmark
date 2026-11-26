"""Fixed-graph glycoprotein geometry, CPU only, no coordinate repair or selection.
Proper-rotation row convention: transformed = X @ rotation + translation.
"""
from __future__ import annotations
import math
from collections import defaultdict, deque
from itertools import combinations
import numpy as np
METRIC_VERSION='glyco_geometry_v1_20261113'
EPS=1e-10
CONFIG={'protein_fit':'fixed_observed_CA_proper_Kabsch_no_outlier_removal','site_cutoff_A':6.0,'fallback_sequence_window':4,'lddt_radius_A':15.0,'lddt_thresholds_A':[.5,1.,2.,4.],'lddt_chemistry_penalty':False,'ring_classes':'diagnostic_poles_not_IUPAC'}

def _xyz(atoms,k):
 v=atoms.get(k)
 if isinstance(v,dict):v=v.get('xyz',v.get('coordinate',v.get('coordinates')))
 try:a=np.asarray(v,dtype=float)
 except (ValueError,TypeError):return None
 return a if a.shape==(3,) and np.isfinite(a).all() else None

def _pos(k):
 try:return int(k.split('/')[2]) if k.startswith('protein/P1/') else None
 except (ValueError,IndexError):return None

def _name(k):return k.rsplit('/',1)[-1]
def _node(k):
 p=k.split('/');return p[2] if len(p)==4 and p[0]=='glycan' else None

def _sort(k):
 p=_pos(k);return (0,p,_name(k)) if p is not None else (1,k,'')

def _points(a,keys):
 v=[_xyz(a,k)for k in keys]
 return None if any(x is None for x in v) else np.asarray(v,dtype=float).reshape((-1,3))

def _good(t):return isinstance(t,dict)and t.get('status')=='ok'

def proper_kabsch(x,y):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
 if x.ndim!=2 or x.shape[1:]!=(3,)or x.shape!=y.shape:return {'status':'shape_error'}
 if len(x)<3:return {'status':'insufficient_fit_points','fit_count':len(x)}
 if not np.isfinite(x).all()or not np.isfinite(y).all():return {'status':'nonfinite_fit_points','fit_count':len(x)}
 xc=x-x.mean(0);yc=y-y.mean(0);rx=int(np.linalg.matrix_rank(xc,tol=1e-8));ry=int(np.linalg.matrix_rank(yc,tol=1e-8))
 if min(rx,ry)<2:return {'status':'degenerate_fit_points','fit_count':len(x),'rank_pred':rx,'rank_ref':ry}
 u,sv,vt=np.linalg.svd(xc.T@yc);d=np.eye(3);d[2,2]=-1. if np.linalg.det(u@vt)<0 else 1.
 R=u@d@vt;t=y.mean(0)-x.mean(0)@R;e=x@R+t-y;sse=float(np.sum(e*e))
 return {'status':'ok','fit_count':len(x),'rotation':R.tolist(),'translation':t.tolist(),'convention':'row_vector_x@R+t','determinant':float(np.linalg.det(R)),'fit_sse_A2':sse,'fit_rmsd_A':math.sqrt(max(0.,sse/len(x))),'rank_pred':rx,'rank_ref':ry,'singular_values':sv.tolist()}

def fit_keys(ra,pa,keys):
 keys=list(dict.fromkeys(keys));mr=[k for k in keys if _xyz(ra,k)is None];mp=[k for k in keys if _xyz(pa,k)is None]
 if mr:return {'status':'reference_missing_or_nonfinite','fit_atom_keys':keys,'missing_reference_keys':mr,'fit_count':len(keys)}
 if mp:return {'status':'model_missing_or_nonfinite','fit_atom_keys':keys,'missing_prediction_keys':mp,'fit_count':len(keys)}
 z=proper_kabsch(_points(pa,keys),_points(ra,keys));z['fit_atom_keys']=keys;return z

def score_transformed(ra,pa,keys,transform,details=False):
 keys=list(dict.fromkeys(keys));z={'atom_count':len(keys),'sse_A2':None,'rmsd_A':None}
 if not keys:z['status']='not_applicable_empty_mask';return z
 mr=[k for k in keys if _xyz(ra,k)is None];mp=[k for k in keys if _xyz(pa,k)is None]
 if mr or mp:z.update(status='reference_missing_or_nonfinite'if mr else'model_missing_or_nonfinite',missing_reference_keys=mr,missing_prediction_keys=mp);return z
 if not _good(transform):z.update(status='fit_unavailable',fit_status=transform.get('status')if transform else None);return z
 e=_points(pa,keys)@np.asarray(transform['rotation'])+np.asarray(transform['translation'])-_points(ra,keys);sq=np.sum(e*e,axis=1)
 z.update(status='ok',sse_A2=float(sq.sum()),rmsd_A=float(np.sqrt(sq.mean())))
 if details:z['atom_squared_errors_A2']={k:float(q)for k,q in zip(keys,sq)}
 return z

def circular_error_degrees(a,b):return None if a is None or b is None else float(abs((float(a)-float(b)+180)%360-180))

def angle_degrees(p):
 if p is None or len(p)!=3:return None
 a,b,c=np.asarray(p,dtype=float);u=a-b;v=c-b;den=float(np.linalg.norm(u)*np.linalg.norm(v))
 return None if den<EPS else float(np.degrees(np.arccos(np.clip(np.dot(u,v)/den,-1.,1.))))

def dihedral_degrees(p):
 if p is None or len(p)!=4:return None
 p0,p1,p2,p3=np.asarray(p,dtype=float);b0=-(p1-p0);b1=p2-p1;b2=p3-p2;n=float(np.linalg.norm(b1))
 if n<EPS:return None
 b1=b1/n;v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
 if min(float(np.linalg.norm(v)),float(np.linalg.norm(w)))<EPS:return None
 return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))

def geometry_record(ra,pa,keys,kind,definition=None):
 keys=list(keys);r=_points(ra,keys);p=_points(pa,keys)
 if kind=='bond':rv=None if r is None else float(np.linalg.norm(r[0]-r[1]));pv=None if p is None else float(np.linalg.norm(p[0]-p[1]))
 else:
  fn=angle_degrees if kind=='angle'else dihedral_degrees;rv=None if r is None else fn(r);pv=None if p is None else fn(p)
 err=circular_error_degrees(pv,rv)if kind=='torsion'else(None if pv is None or rv is None else abs(pv-rv))
 return {'kind':kind,'atom_keys':keys,'definition':definition,'reference':rv,'prediction':pv,'absolute_error':err,'unit':'A'if kind=='bond'else'degree','status':'ok'if rv is not None and pv is not None else('reference_missing_or_degenerate'if rv is None else'model_missing_or_degenerate')}

def cremer_pople(points):
 """CP Fourier coordinates; normal=R_sin cross R_cos, given zero-based ring order.
 Specific IUPAC labels are deliberately not inferred from a chair pole alone.
 """
 p=np.asarray(points,dtype=float);n=len(p)
 if n not in(5,6):return {'status':'unsupported_ring_size','ring_size':n}
 if p.shape!=(n,3)or not np.isfinite(p).all():return {'status':'nonfinite_or_shape_error','ring_size':n}
 p=p-p.mean(0);j=np.arange(n);rs=np.sum(p*np.sin(2*np.pi*j/n)[:,None],axis=0);rc=np.sum(p*np.cos(2*np.pi*j/n)[:,None],axis=0);normal=np.cross(rs,rc);norm=float(np.linalg.norm(normal))
 if norm<EPS:return {'status':'degenerate_ring_plane','ring_size':n}
 normal/=norm;z=p@normal;qc=float(np.sqrt(2/n)*np.sum(z*np.cos(4*np.pi*j/n)));qs=float(-np.sqrt(2/n)*np.sum(z*np.sin(4*np.pi*j/n)));q2=float(np.hypot(qc,qs));phi=float(np.degrees(np.arctan2(qs,qc))%360)if q2>1e-8 else None
 out={'status':'ok','ring_size':n,'q2_A':q2,'phi2_degree':phi,'normal_convention':'R_sin_cross_R_cos; zero_based_input_ring_order','normal_vector':normal.tolist(),'mean_plane_heights_A':z.tolist(),'iupac_conformer':None}
 if n==6:
  q3=float(np.sum(z*((-1.)**j))/np.sqrt(6));q=float(np.hypot(q2,q3));theta=float(np.degrees(np.arctan2(q2,q3)))if q>1e-8 else None
  cl='near_planar_diagnostic'if q<.05 else('chair_pole_0_diagnostic'if theta<=15 else('chair_pole_180_diagnostic'if theta>=165 else'nonchair_region_diagnostic'))
  out.update(q3_A=q3,Q_A=q,theta_degree=theta,diagnostic_class=cl)
 else:out.update(Q_A=q2,diagnostic_class='near_planar_diagnostic'if q2<.05 else'five_member_pseudorotation')
 return out

def ca_lddt(ra,pa,keys,radius=15.,thresholds=(.5,1.,2.,4.)):
 keys=list(keys);n=len(keys)
 if not n:return {'status':'not_applicable_empty_CA_mask','score':None,'pair_count':0,'per_residue':[]}
 y=_points(ra,keys)
 if y is None:return {'status':'reference_missing_or_nonfinite','score':None,'pair_count':0,'per_residue':[]}
 vals=[_xyz(pa,k)for k in keys];good=np.asarray([p is not None for p in vals]);x=np.asarray([p if p is not None else[np.nan]*3 for p in vals]);counts=np.zeros((n,len(thresholds)),dtype=np.int64);den=np.zeros(n,dtype=np.int64);gc=np.zeros(len(thresholds),dtype=np.int64);total=0
 for start in range(0,n,128):
  stop=min(n,start+128);dref=np.linalg.norm(y[start:stop,None,:]-y[None,:,:],axis=2);ii,jj=np.nonzero((dref<radius)&(np.arange(start,stop)[:,None]<np.arange(n)[None,:]))
  if not len(ii):continue
  i=ii+start;dr=dref[ii,jj];dp=np.linalg.norm(x[i]-x[jj],axis=1);valid=good[i]&good[jj];dif=np.abs(dp-dr);den+=np.bincount(np.concatenate([i,jj]),minlength=n)
  for ti,t in enumerate(thresholds):
   passed=valid&(dif<t);np.add.at(counts[:,ti],i[passed],1);np.add.at(counts[:,ti],jj[passed],1);gc[ti]+=int(passed.sum())
  total+=len(i)
 scores=np.divide(counts.sum(1),den*len(thresholds),out=np.full(n,np.nan),where=den>0)
 return {'status':'ok'if total else'not_applicable_no_reference_pairs','score':float(gc.sum()/(total*len(thresholds)))if total else None,'pair_count':total,'thresholds_A':list(thresholds),'radius_A':radius,'preserved_pair_counts':gc.tolist(),'sequence_separation':0,'reference_inclusion_operator':'distance<radius','preservation_operator':'absolute_error<threshold','chemical_quality_penalty':False,'missing_model_CA_count':int((~good).sum()),'per_residue':[{'position':_pos(k),'score':float(scores[i])if np.isfinite(scores[i])else None,'reference_pair_count':int(den[i]),'preserved_counts':counts[i].tolist()}for i,k in enumerate(keys)]}

def _site_data(ref,case):
 if ref.get('sites'):return ref['sites']
 return {g['glycan_id']:{'root_node_id':g['root_node_id'],'node_ids':[n['node_id']for n in g['nodes']],'position':g['attachment']['protein_residue_index'],'link_type':g.get('link_type')}for g in case.get('glycans',[])}

def _node_data(ref,case):
 out=dict(ref.get('nodes')or{})
 for g in case.get('glycans',[]):
  adj=defaultdict(list)
  for b in g.get('bonds',[]):adj[b['node_1']].append(b['node_2']);adj[b['node_2']].append(b['node_1'])
  depths={g['root_node_id']:0};children=defaultdict(list);queue=deque([g['root_node_id']])
  while queue:
   a=queue.popleft()
   for b in adj[a]:
    if b not in depths:depths[b]=depths[a]+1;children[a].append(b);queue.append(b)
  for n in g['nodes']:
   z=dict(out.get(n['node_id'],{}));z.setdefault('glycan_id',g['glycan_id']);z.setdefault('ccd_id',n.get('ccd_id'));z.setdefault('depth',depths.get(n['node_id']));z.setdefault('children',children[n['node_id']]);out[n['node_id']]=z
 return out

def _ring_keys(node):
 k=node.get('ring_atom_keys',[])
 return list(k[0])if k and isinstance(k[0],list)else list(k)

def _graph_bonds(ref,case):
 rows=[]
 for b in ref.get('bonds',[]):
  if isinstance(b,(list,tuple))and len(b)>=2:rows.append((b[0],b[1],b[2]if len(b)>2 else'SING'))
  elif isinstance(b,dict):
   a=b.get('atom_1',b.get('a'));z=b.get('atom_2',b.get('b'))
   if a and z:rows.append((a,z,b.get('order',b.get('bond_order','SING'))))
 for g in case.get('glycans',[]):
  gid=g['glycan_id'];at=g['attachment'];rows.append((f"protein/P1/{at['protein_residue_index']}/{at['protein_atom']}",f"glycan/{gid}/{at['glycan_node_id']}/{at['glycan_atom']}",at.get('bond_order','SING')))
  for b in g.get('bonds',[]):rows.append((f"glycan/{gid}/{b['node_1']}/{b['atom_1']}",f"glycan/{gid}/{b['node_2']}/{b['atom_2']}",b.get('bond_order','SING')))
 unique=[];seen=set()
 for a,b,o in rows:
  k=tuple(sorted((a,b)))
  if k not in seen:unique.append((a,b,o));seen.add(k)
 return unique

def torsion_registry(ref,case,nodes=None,bonds=None):
 """Register explicit C1/C2 donor O-link and ASN/SER/THR root torsions.

 Root psi uses Cdonor-Xacceptor-CG/CB-CB/CA; this explicit convention
 must not be pooled with a different Privateer convention without conversion.
 """
 if ref.get('torsion_registry'):
  result=[]
  for row in ref['torsion_registry']:
   ba=row.get('bond_atoms',row.get('bond_atom_keys',[]));gkeys=[k for k in ba if k.startswith('glycan/')]
   gid=gkeys[0].split('/')[1] if gkeys else None
   kind='protein_glycan' if any(k.startswith('protein/')for k in ba)else'glycan_glycan'
   status=row.get('status','unsupported')
   r={'glycan_id':gid,'kind':kind,'bond_atom_keys':ba,'torsions':[{'name':n,'atom_keys':row[n]}for n in ['phi','psi','omega']if row.get(n)],'template_status':'registered'if status=='registered'else'unsupported','template':'reference_compiler_explicit_four_atom_registry','reference_registry_status':status}
   if status!='registered':r['unsupported_reason']=row.get('reason',row.get('unsupported_reason',status))
   result.append(r)
  return result
 nodes=_node_data(ref,case)if nodes is None else nodes;bonds=_graph_bonds(ref,case)if bonds is None else bonds
 adj=defaultdict(set)
 for a,b,_ in bonds:adj[a].add(b);adj[b].add(a)
 out=[]
 def nk(gid,nid,name):return f'glycan/{gid}/{nid}/{name}'
 for g in case.get('glycans',[]):
  gid=g['glycan_id'];at=g['attachment'];edges=[{'kind':'protein_glycan','a':f"protein/P1/{at['protein_residue_index']}/{at['protein_atom']}",'b':nk(gid,at['glycan_node_id'],at['glycan_atom'])}]
  edges += [{'kind':'glycan_glycan','a':nk(gid,b['node_1'],b['atom_1']),'b':nk(gid,b['node_2'],b['atom_2'])}for b in g.get('bonds',[])]
  for e in edges:
   a,b=e['a'],e['b'];donor=a if _node(a)and _name(a)in('C1','C2')else b if _node(b)and _name(b)in('C1','C2')else None;acceptor=b if donor==a else a if donor==b else None
   rec={'glycan_id':gid,'kind':e['kind'],'bond_atom_keys':[a,b],'torsions':[],'template_status':'unsupported'}
   if donor is None:rec['unsupported_reason']='no_registered_C1_or_C2_donor';out.append(rec);continue
   ring=_ring_keys(nodes.get(_node(donor),{}));ro=[k for k in ring if _name(k).startswith('O')and k in adj[donor]]
   if not ro and donor in ring:
    i=ring.index(donor);ro=[k for k in[ring[(i-1)%len(ring)],ring[(i+1)%len(ring)]]if _name(k).startswith('O')]
   if len(ro)!=1:rec['unsupported_reason']='donor_ring_oxygen_not_uniquely_registered';out.append(rec);continue
   if acceptor.startswith('protein/'):
    pos=_pos(acceptor);name=_name(acceptor)
    if name=='ND2':ac,ap=f'protein/P1/{pos}/CG',f'protein/P1/{pos}/CB'
    elif name in('OG','OG1'):ac,ap=f'protein/P1/{pos}/CB',f'protein/P1/{pos}/CA'
    else:rec['unsupported_reason']='unregistered_protein_attachment_atom';out.append(rec);continue
   else:
    an=_name(acceptor)
    if not(an.startswith('O')and an[1:].isdigit()):rec['unsupported_reason']='unregistered_acceptor_atom_name';out.append(rec);continue
    index=int(an[1:]);anid=_node(acceptor);ac=nk(gid,anid,f'C{index}');ap=nk(gid,anid,f'C{index-1}')if index>1 else None
    if ap is None:rec['unsupported_reason']='acceptor_C1_psi_not_registered';out.append(rec);continue
   rec.update(template_status='registered',donor_atom_key=donor,template='ringO-donorC-acceptorX-acceptorC; psi preceding acceptor carbon; explicit_root_variant')
   rec['torsions']=[{'name':'phi','atom_keys':[ro[0],donor,acceptor,ac]},{'name':'psi','atom_keys':[donor,acceptor,ac,ap]}]
   if not acceptor.startswith('protein/')and _name(acceptor)=='O6':
    omega=[acceptor,ac,ap,nk(gid,_node(acceptor),'O5')]
    if omega[3]in adj[omega[2]]:rec['torsions'].append({'name':'omega','atom_keys':omega})
   out.append(rec)
 return out

def _root_sidechain(ra,pa,case,gid,pos,by_res,gfit,sfit):
 ccd=next((g['attachment'].get('protein_ccd_id')for g in case.get('glycans',[])if g['glycan_id']==gid),None)
 defs={'ASN':[('chi1',['N','CA','CB','CG']),('chi2_OD1',['CA','CB','CG','OD1'])],'SER':[('chi1',['N','CA','CB','OG'])],'THR':[('chi1_OG1',['N','CA','CB','OG1'])]}
 sk=[k for k in by_res.get(pos,[])if _name(k)not in('N','CA','C','O','OXT')]
 return {'root_ccd_id':ccd,'chi_angles':[dict(name=name,**geometry_record(ra,pa,[f'protein/P1/{pos}/{an}'for an in atoms],'torsion','explicit_root_sidechain'))for name,atoms in defs.get(ccd,[])],'sidechain_protein_frame':score_transformed(ra,pa,sk,gfit),'sidechain_site_frame':score_transformed(ra,pa,sk,sfit)}

def centroid_axis_diagnostic(ra,pa,keys,transform):
 """Centroid and unsigned PCA-axis diagnostic; near-degenerate axis is NA."""
 z=score_transformed(ra,pa,keys,transform)
 if z['status']!='ok':return {'status':z['status'],'centroid_error_A':None,'axis_error_degree':None}
 y=_points(ra,keys);x=_points(pa,keys)@np.asarray(transform['rotation'])+np.asarray(transform['translation'])
 out={'status':'ok','centroid_error_A':float(np.linalg.norm(x.mean(0)-y.mean(0))),'axis_error_degree':None,'axis_sign':'unsigned_acute_angle','axis_stability_threshold_relative':1e-6}
 if len(keys)<3:out['axis_status']='insufficient_points';return out
 vals=[];axes=[];gaps=[]
 for v in[x,y]:
  c=v-v.mean(0);e,u=np.linalg.eigh(c.T@c/len(c));order=np.argsort(e)[::-1];e=e[order];u=u[:,order];gap=float((e[0]-e[1])/e[0])if e[0]>1e-12 else 0.;vals.append(e.tolist());axes.append(u[:,0]);gaps.append(gap)
 out.update(eigenvalues_pred_ref_A2=vals,relative_top_eigenvalue_gaps_pred_ref=gaps)
 if min(gaps)<=1e-6:out['axis_status']='near_degenerate_axis';return out
 out['axis_status']='ok';out['axis_error_degree']=float(np.degrees(np.arccos(np.clip(abs(np.dot(axes[0],axes[1])),-1.,1.))))
 return out

def evaluate_geometry(ref,pred,case):
 """Evaluate continuous fixed-mask geometry without modifying inputs.

 Missing model atoms remain missing, not removed from primary score masks.
 TM-score is intentionally delegated to a separate fixed-alignment adapter.
 """
 ra=ref.get('atoms',{});pa=pred.get('atoms',{});sites=_site_data(ref,case);nodes=_node_data(ref,case)
 ca=sorted([k for k in ra if k.startswith('protein/P1/')and k.endswith('/CA')],key=_sort)
 gfit=fit_keys(ra,pa,ca);pg=score_transformed(ra,pa,ca,gfit)
 protein={'ca_global':pg,'global_transform':gfit,'ca_lddt':ca_lddt(ra,pa,ca),'per_residue':[],'sites':{},'tm_score':{'status':'external_fixed_sequence_adapter_required','score':None}}
 by_res=defaultdict(list)
 for k in sorted([k for k in ra if k.startswith('protein/P1/')],key=_sort):by_res[_pos(k)].append(k)
 for pos,ks in sorted(by_res.items()):
  cak=f'protein/P1/{pos}/CA';cac=score_transformed(ra,pa,[cak],gfit)if cak in ra else{'status':'reference_CA_unobserved','rmsd_A':None}
  protein['per_residue'].append({'position':pos,'CA_displacement_A':cac.get('rmsd_A'),'CA_status':cac['status'],'common_heavy':score_transformed(ra,pa,ks,gfit)})
 out={'schema':METRIC_VERSION,'case_id':ref.get('case_id',case.get('case_id')),'definitions':dict(CONFIG),'protein':protein,'sites':{},'nodes':{},'linkages':[],'angles':[],'tm_score_included':False,'coordinates_modified':False}
 for gid,st in sites.items():
  nids=list(st.get('node_ids')or[k for k,v in nodes.items()if v.get('glycan_id')==gid]);gkeys=[k for k in ra if k.startswith(f'glycan/{gid}/')];rid=st['root_node_id'];rkeys=_ring_keys(nodes.get(rid,{}));pos=int(st.get('position',st.get('protein_residue_index')))
  wca=[k for k in ca if abs(_pos(k)-pos)<=4];wfit=fit_keys(ra,pa,wca);pocket={pos};gy=_points(ra,gkeys)
  if gy is not None and len(gy):
   for rp,ks in by_res.items():
    coords=_points(ra,ks)
    if coords is not None and len(coords)and np.any(np.sum((coords[:,None]-gy[None,:])**2,axis=2)<=36.):pocket.add(rp)
  pca=[k for k in ca if _pos(k)in pocket];rp=_points(ra,pca);rank=int(np.linalg.matrix_rank(rp-rp.mean(0),tol=1e-8))if rp is not None and len(rp)else 0;fallback=len(pca)<3 or rank<2
  if fallback:pca=sorted(set(pca)|set(wca),key=_sort)
  sfit=fit_keys(ra,pa,pca);sfit['reference_defined_site_residue_positions']=sorted(pocket);sfit['sequence_window_fallback_used']=fallback
  fits={'protein':gfit,'site':sfit,'root':fit_keys(ra,pa,rkeys),'self':fit_keys(ra,pa,gkeys)}
  frames={name:dict(score_transformed(ra,pa,gkeys,tr),transform=tr)for name,tr in fits.items()}
  site={'glycan_id':gid,'position':pos,'root_node_id':rid,'link_type':st.get('link_type'),'node_ids':nids,'frames':frames,'reference_atom_count':len(gkeys),'max_depth':max([nodes.get(n,{}).get('depth',0)or 0 for n in nids],default=0),'branch_node_count':sum(len(nodes.get(n,{}).get('children',[]))>=2 for n in nids)}
  site['centroid_axis_protein_frame']=centroid_axis_diagnostic(ra,pa,gkeys,gfit)
  rootinfo={'position':pos,'window_positions':[_pos(k)for k in wca],'window_global_CA':score_transformed(ra,pa,wca,gfit),'window_self_CA':score_transformed(ra,pa,wca,wfit),'window_transform':wfit,'spatial_neighborhood_CA_global':score_transformed(ra,pa,pca,gfit)}
  rootinfo.update(_root_sidechain(ra,pa,case,gid,pos,by_res,gfit,sfit));protein['sites'][gid]=rootinfo
  for nid in nids:
   nd=nodes.get(nid,{});nk=[k for k in gkeys if _node(k)==nid];rk=_ring_keys(nd);nt=fit_keys(ra,pa,nk);rt=fit_keys(ra,pa,rk)
   nr={'node_id':nid,'glycan_id':gid,'ccd_id':nd.get('ccd_id'),'depth':nd.get('depth'),'children':nd.get('children',[]),'is_terminal':not nd.get('children'),'is_branch_point':len(nd.get('children',[]))>=2,'frame_errors':{name:score_transformed(ra,pa,nk,tr)for name,tr in fits.items()},'monosaccharide_self':score_transformed(ra,pa,nk,nt),'ring_self':score_transformed(ra,pa,rk,rt),'ring_atom_keys':rk,'ring_protein_frame':score_transformed(ra,pa,rk,gfit)}
   if rk:
    ry,px=_points(ra,rk),_points(pa,rk);nr['reference_pucker']=cremer_pople(ry)if ry is not None else{'status':'reference_missing'};nr['prediction_pucker']=cremer_pople(px)if px is not None else{'status':'model_missing'};rr,pp=nr['reference_pucker'],nr['prediction_pucker']
    if rr.get('status')==pp.get('status')=='ok':
     nr['pucker_comparison']={'Q_abs_error_A':abs(rr['Q_A']-pp['Q_A']),'theta_abs_error_degree':abs(rr['theta_degree']-pp['theta_degree'])if rr.get('theta_degree')is not None and pp.get('theta_degree')is not None else None,'phi2_circular_error_degree':circular_error_degrees(rr.get('phi2_degree'),pp.get('phi2_degree')),'diagnostic_class_agreement':rr['diagnostic_class']==pp['diagnostic_class']}
   if _good(gfit)and nk and _points(ra,nk)is not None and _points(pa,nk)is not None:
    xx=_points(pa,nk)@np.asarray(gfit['rotation'])+np.asarray(gfit['translation']);nr['centroid_protein_frame_error_A']=float(np.linalg.norm(xx.mean(0)-_points(ra,nk).mean(0)))
   out['nodes'][nid]=nr
  site['depth_shells']=[]
  for depth in sorted({nodes.get(n,{}).get('depth')for n in nids if nodes.get(n,{}).get('depth')is not None}):
   for cumulative in(False,True):
    selected=[n for n in nids if nodes.get(n,{}).get('depth')is not None and(nodes[n]['depth']<=depth if cumulative else nodes[n]['depth']==depth)];ks=[k for k in gkeys if _node(k)in selected]
    site['depth_shells'].append({'depth':depth,'cumulative':cumulative,'node_ids':selected,'frame_errors':{name:score_transformed(ra,pa,ks,tr)for name,tr in fits.items()}})
  out['sites'][gid]=site
 bonds=_graph_bonds(ref,case);adj=defaultdict(set)
 for a,b,o in bonds:adj[a].add(b);adj[b].add(a)
 for center in sorted(adj):
  for a,b in combinations(sorted(adj[center]),2):
   if any(k.startswith('glycan/')for k in(a,center,b)):out['angles'].append(geometry_record(ra,pa,[a,center,b],'angle','reacted_graph_neighbor_triplet'))
 for rec in torsion_registry(ref,case,nodes,bonds):
  rec['bond_geometry']=geometry_record(ra,pa,rec['bond_atom_keys'],'bond','provided_covalent_graph');rec['torsions']=[dict(name=t['name'],**geometry_record(ra,pa,t['atom_keys'],'torsion',rec.get('template')))for t in rec['torsions']]
  if rec['kind']=='protein_glycan':
   a,b=rec['bond_atom_keys'];pk=a if a.startswith('protein/')else b;gk=b if pk==a else a;pos=_pos(pk);name=_name(pk);an='CG'if name=='ND2'else'CB'if name in('OG','OG1')else None
   rec['root_attachment_angle']=geometry_record(ra,pa,[f'protein/P1/{pos}/{an}',pk,gk],'angle','N_CG-ND2-C1_or_O_CB-OGx-C1')if an else{'status':'unsupported_attachment_atom'}
  out['linkages'].append(rec)
 pv=[v['frames']['protein'].get('rmsd_A')for v in out['sites'].values()];sv=[v['frames']['site'].get('rmsd_A')for v in out['sites'].values()]
 out['summary']={'protein_CA_RMSD_A':pg.get('rmsd_A'),'protein_CA_lDDT':protein['ca_lddt'].get('score'),'site_count':len(sites),'glycan_protein_frame_macro_RMSD_A':float(np.mean(pv))if pv and all(v is not None for v in pv)else None,'glycan_site_frame_macro_RMSD_A':float(np.mean(sv))if sv and all(v is not None for v in sv)else None,'protein_frame_sites_scored':sum(v is not None for v in pv),'site_frame_sites_scored':sum(v is not None for v in sv),'linkage_count':len(out['linkages']),'registered_linkage_count':sum(x['template_status']=='registered'for x in out['linkages'])}
 return out
