"""Scoring functions extracted from the benchmark evaluation scripts.
Atom identifiers and coordinate policies are documented in ../README.md.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import gzip, hashlib, json, math, copy
import numpy as np
VERSION='glyco_eval_common_v1_20261113'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def load(p):
 op=gzip.open if str(p).endswith('.gz')else open
 with op(p,'rt',encoding='utf-8')as f:return json.load(f)

def inventory(ref,pred):
 expected=set(ref['expected_atom_meta']);observed=set(pred['atoms']);pkeys={k for k in expected if k.startswith('protein/')};gkeys={k for k in expected if k.startswith('glycan/')};pos=len(ref['protein_sequence']);optional={f'protein/P1/{pos}/OXT'}
 missing=pkeys-observed;mg=gkeys-observed;extra=observed-expected-optional;nonfinite=[k for k,v in pred['atoms'].items()if len(v)!=3 or not all(math.isfinite(float(z))for z in v)];ca_missing=[f'protein/P1/{i}/CA'for i in range(1,pos+1)if f'protein/P1/{i}/CA'not in observed];dupes=pred.get('duplicate_keys',[]);em=[]
 for k in expected&observed:
  want=ref['expected_atom_meta'][k]['element'].upper();got=pred.get('atom_meta',{}).get(k,{}).get('element')
  if got and got.upper()!=want:em.append({'atom':k,'expected':want,'observed':got})
 return {'version':VERSION,'status':'pass'if not(missing or mg or extra or nonfinite or ca_missing or dupes or em)else'fail','expected_protein_heavy_atoms':len(pkeys),'expected_glycan_heavy_atoms':len(gkeys),'observed_protein_heavy_atoms':sum(k.startswith('protein/')for k in observed),'observed_glycan_heavy_atoms':sum(k.startswith('glycan/')for k in observed),'missing_protein_atoms':sorted(missing),'missing_glycan_atoms':sorted(mg),'unexpected_atoms':sorted(extra),'optional_terminal_OXT_observed':sorted(optional&observed),'duplicate_keys':dupes,'nonfinite_atoms':nonfinite,'missing_CA_positions':[int(k.split('/')[2])for k in ca_missing],'element_mismatches':em,'global_hard_fail':bool(nonfinite or ca_missing or dupes),'independent_sequence_CCD_identity_checked_here':False,'sequence_identity_basis':'frozen original v1 protein_sequence report retained separately','mask':'complete standard protein and reacted glycan identity, not reference-only observed mask'}

def _gid(k):
 v=k.split('/')if isinstance(k,str)else[]
 return v[1]if len(v)>1 and v[0]=='glycan'else None

def local_v1(ref,report,inv):
 report_path=report if isinstance(report,(str,Path))else None
 if report_path is not None:report=load(report_path)
 out={};seq=report.get('protein_sequence',{});global_fail=inv.get('global_hard_fail',False)or seq.get('status')=='fail'or report.get('finite_coordinates')is False
 centers=report.get('stereochemistry',{}).get('centers',[]);gb=report.get('glycosidic_geometry',{}).get('bonds',[]);wb=report.get('within_monosaccharide_geometry',{}).get('bond_distances',[]);gi=report.get('glycan_heavy_atom_identity',{});unlocalized_mapping_errors=gi.get('mapping_errors',[])
 if unlocalized_mapping_errors:global_fail=True
 for gid,site in ref['sites'].items():
  pos=int(site['position']);prefix=f'protein/P1/{pos}/';belongs=lambda k:_gid(k)==gid or(isinstance(k,str)and k.startswith(prefix))
  cs=[c for c in centers if belongs(c.get('atom'))];bonds=[b for b in gb+wb if belongs(b.get('atom_1'))or belongs(b.get('atom_2'))];badcent=[c for c in cs if c.get('pass')is False];badbonds=[b for b in bonds if b.get('broad_anomaly')is True or b.get('within_fixed_routine_range')is False or b.get('distance_A')is None]
  miss=[k for k in inv.get('missing_glycan_atoms',[])+inv.get('missing_protein_atoms',[])if belongs(k)];extra=[k for k in inv.get('unexpected_atoms',[])if belongs(k)];forbidden=[k for k in gi.get('forbidden_leaving_heavy_atoms_observed',[])if belongs(k)]
  unknown=not cs or not bonds or report.get('stereochemistry',{}).get('status')in [None,'unassessed','error']
  failed=bool(global_fail or badcent or badbonds or miss or extra or forbidden)
  out[gid]={'status':'fail'if failed else'unknown'if unknown else'pass','profile':'legacy_v1_local_components_only','global_hard_fail':global_fail,'centers_checked':len(cs),'stereo_failures':badcent,'bonds_checked':len(bonds),'bond_failures':badbonds,'missing_atoms':miss,'extra_atoms':extra,'forbidden_leaving_atoms':forbidden,'unknown_required_legacy_coverage':unknown,'centers':cs,'bonds':bonds,'new_angles_pucker_clashes_not_in_this_profile':True}
 return {'version':VERSION,'profile':'legacy_v1_local_components_only','legacy_candidate_status':report.get('status'),'global_sequence_status':seq.get('status'),'global_hard_fail':global_fail,'sites':out,'input_report_SHA256':sha(report_path)if report_path is not None else None,'not_whole_molecule_extended_validity':True}

@lru_cache(maxsize=16)
def _parent_centers(sequence):
 from rdkit import Chem
 mol=Chem.MolFromFASTA(sequence);Chem.AssignStereochemistry(mol,cleanIt=True,force=True);centers=[]
 for atom in mol.GetAtoms():
  if not atom.HasProp('_CIPCode'):continue
  info=atom.GetPDBResidueInfo();pos=info.GetResidueNumber();name=info.GetName().strip();neigh=[]
  for n in atom.GetNeighbors():
   ni=n.GetPDBResidueInfo();neigh.append(f'protein/P1/{ni.GetResidueNumber()}/{ni.GetName().strip()}')
  centers.append({'atom':f'protein/P1/{pos}/{name}','expected_parent_CIP':atom.GetProp('_CIPCode'),'expected_chiral_tag':int(atom.GetChiralTag()),'neighbor_keys_in_RDKit_order':neigh})
 # Coordinate determinant-to-neighbor-order tag calibration on a complete 4-atom fragment.
 rw=Chem.RWMol();c=Chem.Atom('C');c.SetNumExplicitHs(1);c.SetNoImplicit(True);rw.AddAtom(c)
 for sym in ['N','C','O']:
  j=rw.AddAtom(Chem.Atom(sym));rw.AddBond(0,j,Chem.BondType.SINGLE)
 toy=rw.GetMol();toy.UpdatePropertyCache(strict=False);conf=Chem.Conformer(4);conf.Set3D(True)
 for i,v in enumerate([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]):conf.SetAtomPosition(i,v)
 toy.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(toy,replaceExistingTags=True);positive_tag=int(toy.GetAtomWithIdx(0).GetChiralTag())
 return centers,positive_tag

def protein_stereochemistry(ref,pred):
 try:centers,positive_tag=_parent_centers(ref['protein_sequence'])
 except Exception as e:return {'status':'unsupported','error':str(e),'graph_basis':'standard unmodified parent protein MolFromFASTA'}
 rows=[]
 for c in centers:
  atom=c['atom'];ns=c['neighbor_keys_in_RDKit_order'];keys=[atom]+ns;missing=[k for k in keys if k not in pred['atoms']]
  r={**c,'observed_parent_CIP':None,'signed_normalized_triple_product':None,'pass':None}
  if missing:r.update(status='model_missing_neighbor_atoms',missing_atoms=missing)
  elif len(ns)!=3:r.update(status='unsupported_neighbor_count',neighbor_count=len(ns))
  else:
   p=np.asarray([pred['atoms'][k]for k in keys],float);v=p[1:]-p[0];den=float(np.prod(np.linalg.norm(v,axis=1)));val=float(np.linalg.det(v)/den)if den>1e-12 else 0.
   r['signed_normalized_triple_product']=val
   if not np.isfinite(p).all():r.update(status='nonfinite',pass_value=False)
   elif den<=1e-12 or abs(val)<1e-3:r.update(status='near_coplanar_or_collapsed',pass_value=False)
   else:
    tag=positive_tag if val>0 else(1 if positive_tag==2 else 2);same=tag==c['expected_chiral_tag'];actual=c['expected_parent_CIP']if same else('S'if c['expected_parent_CIP']=='R'else'R');r.update(status='pass'if same else'inverted',observed_parent_CIP=actual,pass_value=same)
   r['pass']=r.pop('pass_value',None)
  rows.append(r)
 fail=[r for r in rows if r['pass']is False];unknown=[r for r in rows if r['pass']is None]
 return {'version':VERSION,'status':'fail'if fail else'unknown'if unknown else'pass','graph_basis':'standard unmodified parent protein MolFromFASTA; neighbor-order parity; NOT full glycosylated product CIP ranking','calibrated_positive_determinant_RDKit_tag':positive_tag,'near_coplanar_normalized_volume_threshold':1e-3,'centers_total':len(rows),'centers_evaluable':len(rows)-len(unknown),'inverted_count':sum(r['status']=='inverted'for r in rows),'degenerate_count':sum(r['status']=='near_coplanar_or_collapsed'for r in rows),'failure_count':len(fail),'unknown_count':len(unknown),'centers':rows,'coordinates_imputed_or_modified':False}

def convert(v1,seq):
 v=copy.deepcopy(v1);changes=[]
 for r in v['centers']:
  if r['atom'].endswith('/CB')and seq[int(r['atom'].split('/')[2])-1]=='I':
   old=r['status'];r['neighbor_keys_in_RDKit_order']=[k[:-3]+('CG2'if k.endswith('CG1')else'CG1')if k.endswith(('/CG1','/CG2'))else k for k in r['neighbor_keys_in_RDKit_order']]
   if r.get('signed_normalized_triple_product')is not None:r['signed_normalized_triple_product']*=-1
   if old in['pass','inverted']:
    r['status']='pass'if old=='inverted'else'inverted';r['pass']=not r['pass'];r['observed_parent_CIP']=r['expected_parent_CIP']if r['pass']else('R'if r['expected_parent_CIP']=='S'else'S')
   r['v1_ILE_name_parity_corrected']=True;changes.append(dict(atom=r['atom'],v1_status=old,v2_status=r['status']))
 fail=[r for r in v['centers']if r.get('pass')is False];unknown=[r for r in v['centers']if r.get('pass')is None]
 v.update(version="parent_protein_stereo_CCD_ILE_name_v2_20261113",status='fail'if fail else'unknown'if unknown else'pass',graph_basis='Standard parent MolFromFASTA with ILE CG1/CG2 name binding corrected to official CCD CG1-CD1, CB=S; not full glycoprotein product CIP',inverted_count=sum(r['status']=='inverted'for r in v['centers']),degenerate_count=sum(r['status']=='near_coplanar_or_collapsed'for r in v['centers']),failure_count=len(fail),unknown_count=len(unknown),ILE_name_correction=changes,old_status=v1['status'],coordinates_imputed_or_modified=False)
 return v

def _three_value_and(states):
 if any(s=='fail'for s in states):return 'fail'
 return 'pass'if states and all(s in('pass','not_applicable_known_no_centers')for s in states)else'unknown'

def _root_parent_stereo(ref,ps,position):
 prefix=f'protein/P1/{position}/';expected={x['atom']for x in ref.get('protein_standard_parent_stereocenters',[])if x.get('atom','').startswith(prefix)}
 rows=[x for x in ps.get('centers',[])if x.get('atom','').startswith(prefix)];lookup={x.get('atom'):x for x in rows}
 if not expected:
  # Standard, nonchiral anchor identity must be explicit to count as inapplicable.
  seq=ref.get('protein_sequence','');known_nonchiral=1<=position<=len(seq)and seq[position-1]=='G'
  return {'status':'not_applicable_known_no_centers'if known_nonchiral else'unknown','centers_expected':0,'centers_seen':len(rows),'failures':sum(x.get('pass')is False for x in rows),'unknown':sum(x.get('pass')is None for x in rows)}
 failed=sum(lookup.get(k,{}).get('pass')is False for k in expected);unknown=sum(k not in lookup or lookup[k].get('pass')is None for k in expected)
 return {'status':'fail'if failed else'unknown'if unknown else'pass','centers_expected':len(expected),'centers_seen':sum(k in lookup for k in expected),'failures':failed,'unknown':unknown}

def _finite_number(x):
 return isinstance(x,(int,float,np.integer,np.floating)) and not isinstance(x,(bool,np.bool_)) and math.isfinite(float(x))

def _joint_bit(local_status,error,tau,model_hard_failure=False):
 if model_hard_failure or local_status=='fail':return 0
 if _finite_number(error):
  if float(error)>tau:return 0
  return 1 if local_status=='pass'else None
 return None

