"""Fixed-correspondence TM-score, truth-free candidate diversity, and scalar views.

No structure prediction, relaxation, candidate repair, or quality-based selection.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from itertools import combinations
from functools import lru_cache
import hashlib,json,math,os,re,subprocess,tempfile
import numpy as np
from .geometry_metrics import fit_keys,score_transformed,geometry_record,_points,_pos,_name

SUPPLEMENTAL_VERSION='glyco_supplemental_v1_20261113'
AA3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

def _sha_bytes(b):return hashlib.sha256(b).hexdigest()
@lru_cache(maxsize=8)
def _executable_sha(path,mtime_ns,size):
 return _sha_bytes(Path(path).read_bytes())

def _ca_mask(ref):
 return sorted([k for k in ref.get('atoms',{})if k.startswith('protein/P1/')and k.endswith('/CA')],key=lambda k:_pos(k))

def _pdb_ca(atoms,keys,sequence):
 lines=[];rounding=[]
 for serial,key in enumerate(keys,1):
  pos=_pos(key)
  if pos is None or not 1<=pos<=len(sequence)or sequence[pos-1]not in AA3:raise ValueError('Reference CA position lacks a registered standard amino acid')
  xyz=_points(atoms,[key])
  if xyz is None:raise ValueError('Missing or nonfinite fixed-mask CA atom')
  x,y,z=xyz[0]
  if max(abs(x),abs(y),abs(z))>=9999:raise ValueError('PDB interface coordinate exceeds registered width')
  strings=[f'{q:8.3f}'for q in(x,y,z)]
  if any(len(q)!=8 for q in strings):raise ValueError('PDB interface coordinate exceeds fixed numeric field width')
  readback=np.array([float(q)for q in strings]);rounding.append(float(np.max(np.abs(readback-xyz[0]))))
  lines.append(f'ATOM  {serial:5d}  CA  {AA3[sequence[pos-1]]:>3s} A{serial:4d}    {strings[0]}{strings[1]}{strings[2]}  1.00  0.00           C  ')
 return '\n'.join(lines+['TER','END'])+'\n',max(rounding,default=0.)

def evaluate_tm(ref,pred,workdir=None,executable=None,timeout_seconds=60):
 """Run TM-align with immutable fixed CA mask and explicit matching alignment.

 Both CA-only interface PDBs have identical observed-mask sequence/order. The
 -I fasta fixes correspondence; -L normalizes to reference evaluable CA count.
 A fixed aligned-length check rejects any apparent correspondence dropping.
 """
 exe=Path(executable or os.environ.get('GLYCO_TMALIGN','/usr/bin/TMalign'))
 keys=_ca_mask(ref);n=len(keys);sequence=ref.get('protein_sequence','')
 out={'schema':SUPPLEMENTAL_VERSION,'status':'not_run','score':None,'tm_score':None,'fixed_matched_CA_count':n,'normalization_CA_count':n,'reference_CA_positions':[_pos(k)for k in keys],'fixed_correspondence':True,'free_structure_alignment_allowed':False,'coordinate_interface':'CA-only_PDB_3_decimals; original coordinates unchanged','coordinates_modified':False}
 if n<3:out['status']='insufficient_reference_CA';return out
 if _points(ref.get('atoms',{}),keys)is None:out['status']='reference_missing_or_nonfinite';return out
 missing=[k for k in keys if _points(pred.get('atoms',{}),[k])is None]
 if missing:out.update(status='model_missing_or_nonfinite',missing_prediction_keys=missing);return out
 if not exe.is_file():out.update(status='tool_unavailable',executable=str(exe));return out
 try:
  seq=''.join(sequence[_pos(k)-1]for k in keys)
  if any(a not in AA3 for a in seq):raise ValueError('Unregistered residue in fixed-mask sequence')
  ptext,pround=_pdb_ca(pred['atoms'],keys,sequence);rtext,rround=_pdb_ca(ref['atoms'],keys,sequence)
 except (ValueError,IndexError,TypeError)as e:out.update(status='interface_input_error',error=str(e));return out
 align=f'>prediction_fixed_CA_mask\n{seq}\n>reference_fixed_CA_mask\n{seq}\n'
 out.update(interface_prediction_pdb_sha256=_sha_bytes(ptext.encode()),interface_reference_pdb_sha256=_sha_bytes(rtext.encode()),fixed_alignment_sha256=_sha_bytes(align.encode()),max_per_coordinate_rounding_A=max(pround,rround))
 parent=Path(workdir)if workdir is not None else None
 if parent is not None:parent.mkdir(parents=True,exist_ok=True)
 st=exe.stat();out.update(executable=str(exe),executable_sha256=_executable_sha(str(exe.resolve()),st.st_mtime_ns,st.st_size))
 with tempfile.TemporaryDirectory(prefix='tm_fixed_',dir=str(parent)if parent else None)as tmp:
  d=Path(tmp);pp=d/'prediction_CA.pdb';rr=d/'reference_CA.pdb';aa=d/'fixed_correspondence.fasta';pp.write_text(ptext);rr.write_text(rtext);aa.write_text(align)
  cmd=[str(exe),str(pp),str(rr),'-I',str(aa),'-L',str(n)]
  out['argv_template']=[str(exe),'prediction_CA.pdb','reference_CA.pdb','-I','fixed_correspondence.fasta','-L',str(n)]
  try:
   proc=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout_seconds,env={**os.environ,'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
  except subprocess.TimeoutExpired as e:out.update(status='tool_timeout',timeout_seconds=timeout_seconds,stdout_excerpt=str(e.stdout or'')[-2000:]);return out
  except OSError as e:out.update(status='tool_execution_error',error=type(e).__name__+': '+str(e));return out
  out.update(exit_code=proc.returncode,raw_stdout=proc.stdout,raw_stderr=proc.stderr,temporary_interface_files_deleted=True)
 if proc.returncode!=0:out['status']='tool_nonzero_exit';return out
 aligned=re.search(r'Aligned\s+length\s*=\s*(\d+)',proc.stdout,re.I)
 out['reported_aligned_length']=int(aligned.group(1))if aligned else None
 if aligned is None or int(aligned.group(1))!=n:out['status']='fixed_alignment_count_not_confirmed';return out
 candidates=[]
 for line in proc.stdout.splitlines():
  m=re.search(r'TM-score\s*=\s*([0-9.eE+\-]+)',line)
  if m:candidates.append((float(m.group(1)),line))
 out['raw_tm_score_lines']=[line for _,line in candidates]
 chosen=[x for x in candidates if 'user' in x[1].lower() and any(v in x[1].lower()for v in ['specified','assign','normal'])]
 if not chosen:
  # Both interface chains have exactly n residues, so chain-2 normalization is
  # identical to requested reference length if this old writer omits a third line.
  chosen=[x for x in candidates if re.search(r'(Chain[_ ]?2|Structure[_ ]?2)',x[1],re.I)]
  origin='chain_2_length_equals_explicit_reference_mask_length'
 else:origin='user_specified_L'
 if len(chosen)!=1:out['status']='normalization_line_not_uniquely_identified';return out
 value,line=chosen[0]
 if not math.isfinite(value)or not -1e-6<=value<=1.000001:out['status']='invalid_tm_value';return out
 out.update(status='ok',score=value,tm_score=value,normalization_line=line,normalization_origin=origin)
 return out

def _compact_transform(t,definition):
 out={k:v for k,v in t.items()if k!='fit_atom_keys'}
 if 'fit_atom_keys'in t:out['fit_atom_keys_sha256']=_sha_bytes('\n'.join(t['fit_atom_keys']).encode());out['fit_mask_count']=len(t['fit_atom_keys'])
 out['fit_mask_definition']=definition;return out

def _input_sites(ref):return ref.get('sites',{})

def _input_glycan_keys(ref,views):
 meta=ref.get('expected_atom_meta')or{}
 keys=[k for k in meta if k.startswith('glycan/')]
 if keys:return sorted(keys),'known_reacted_input_expected_atom_meta'
 keys=[]
 for nid,nd in ref.get('nodes',{}).items():
  for atom in nd.get('expected_reacted_heavy_atom_names',[]):keys.append(f"glycan/{nd['glycan_id']}/{nid}/{atom}")
 if keys:return sorted(set(keys)),'known_reacted_input_node_atom_names'
 return sorted({k for v in views for k in v.get('atoms',{})if k.startswith('glycan/')}),'union_prediction_glycan_keys; expected_input_mask_not_supplied'

def _precomputed_contact_jaccard(a,b):
 ca=a.get('precomputed_contacts');cb=b.get('precomputed_contacts')
 if not isinstance(ca,dict)or not isinstance(cb,dict):return {'status':'precomputed_contact_sets_unavailable','jaccard':None}
 da=ca.get('definition_hash',ca.get('definition_id'));db=cb.get('definition_hash',cb.get('definition_id'))
 if not da or da!=db:return {'status':'contact_definitions_missing_or_incomparable','jaccard':None}
 try:
  sa={tuple(map(str,x))for x in ca['pairs']};sb={tuple(map(str,x))for x in cb['pairs']}
 except (KeyError,TypeError):return {'status':'invalid_precomputed_contact_schema','jaccard':None}
 union=sa|sb;inter=sa&sb
 return {'status':'ok'if union else'not_applicable_both_sets_empty','jaccard':len(inter)/len(union)if union else None,'intersection_n':len(inter),'union_n':len(union),'candidate_i_contacts_n':len(sa),'candidate_j_contacts_n':len(sb),'definition_hash':da}

def evaluate_candidate_pairs(ref,candidate_views):
 """Within-case diversity using only input metadata and prediction coordinates.

 This function deliberately never reads ref['atoms']. Input windows use known
 attachment positions +/-4; global anchor uses every input protein CA position.
 Contact Jaccard is computed only from compatible caller-supplied contact sets.
 """
 if isinstance(candidate_views,dict):views=[dict(v,candidate_index=v.get('candidate_index',k))for k,v in candidate_views.items()]
 else:views=list(candidate_views)
 sites=_input_sites(ref);length=len(ref.get('protein_sequence',''));global_keys=[f'protein/P1/{p}/CA'for p in range(1,length+1)];glykeys,maskorigin=_input_glycan_keys(ref,views)
 ids=[v.get('candidate_index',i)for i,v in enumerate(views)];out={'schema':SUPPLEMENTAL_VERSION,'case_id':ref.get('case_id'),'candidate_count':len(views),'expected_five_candidate_count':5,'all_five_views_present':len(views)==5,'candidate_ids':ids,'expected_pairs_for_five':10,'reference_coordinates_used':False,'anchors':'all_input_CA_and_input_site_plusminus4_CA; never_reference_6A','glycan_scoring_mask_origin':maskorigin,'pairs':[],'pair_count':0}
 if len(set(map(str,ids)))!=len(ids):out['status']='duplicate_candidate_ids';return out
 registry=ref.get('torsion_registry',[])
 for i,j in combinations(range(len(views)),2):
  a=views[i].get('atoms',{});b=views[j].get('atoms',{});gf=fit_keys(b,a,global_keys)
  pair={'candidate_i':ids[i],'candidate_j':ids[j],'global_prediction_to_prediction_transform':_compact_transform(gf,'all_input_positions_1_to_L_CA'),'sites':{},'torsion_differences':[],'contact_jaccard':_precomputed_contact_jaccard(views[i],views[j])}
  for gid,site in sites.items():
   pos=int(site.get('position',site.get('protein_residue_index')));window=[f'protein/P1/{p}/CA'for p in range(max(1,pos-4),min(length,pos+4)+1)];keys=[k for k in glykeys if k.startswith(f'glycan/{gid}/')];wf=fit_keys(b,a,window);sf=fit_keys(b,a,keys)
   pair['sites'][gid]={'global_protein_frame':score_transformed(b,a,keys,gf),'input_window_frame':score_transformed(b,a,keys,wf),'self_shape':score_transformed(b,a,keys,sf),'input_window_transform':_compact_transform(wf,'input_attachment_position_plusminus4_CA'),'self_shape_transform':_compact_transform(sf,'known_input_glycan_heavy_atoms')}
  for link in registry:
   for name in ['phi','psi','omega']:
    if link.get(name):
     rec=geometry_record(b,a,link[name],'torsion','same_input_chemical_registry_candidate_i_vs_j');rec.update(name=name,bond_atoms=link.get('bond_atoms'),source='prediction_i_vs_prediction_j');pair['torsion_differences'].append(rec)
  for field in ['global_protein_frame','input_window_frame','self_shape']:
   vals=[x[field]['rmsd_A']for x in pair['sites'].values()];pair[field+'_macro_RMSD_A']=float(np.mean(vals))if vals and all(v is not None for v in vals)else None;pair[field+'_sites_scored']=sum(v is not None for v in vals)
  angles=[t['absolute_error']for t in pair['torsion_differences']];valid=[v for v in angles if v is not None]
  pair['registered_torsion_count']=len(angles);pair['torsions_scored']=len(valid);pair['torsion_mean_absolute_circular_difference_degree']=float(np.mean(valid))if angles and len(valid)==len(angles)else None
  pair['status']='ok'if pair['global_protein_frame_macro_RMSD_A']is not None else'partial_or_missing_geometry'
  out['pairs'].append(pair)
 out['pair_count']=len(out['pairs']);out['status']='ok'if len(views)==5 and all(x['status']=='ok'for x in out['pairs'])else'incomplete_or_partial_pairs'
 return out

