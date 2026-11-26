"""Model-output atom naming from the frozen prediction validators."""
import collections
from pathlib import Path
def pkey(pos,atom):return ('protein','P1',str(pos),atom)

def gkey(gid,node,atom):return ('glycan',gid,node,atom)

def missing(s):return s in ('','.', '?',None)

def choose(*v):return next((x for x in v if not missing(x)),None)

def cif_dictionary(path):
 try:
  import gemmi
  b=gemmi.cif.read_file(str(path)).sole_block();out={}
  for category in ('_atom_site.','_struct_conn.'):
   for key,values in b.get_mmcif_category(category).items():out[category+key]=['?' if v is None else '.' if v is False else v for v in values]
  return out
 except ImportError:
  from Bio.PDB.MMCIF2Dict import MMCIF2Dict
  return MMCIF2Dict(str(path))

def read_structure(path,model=1):
 path=Path(path);atoms=[];connections=[];counts=collections.Counter()
 if path.suffix.lower() in ('.cif','.mmcif'):
  d=cif_dictionary(path);n=len(d.get('_atom_site.id',[]))
  def col(k,i,default='?'):return (d[k][i] if k in d else default)
  models=[]
  for i in range(n):
   mid=col('_atom_site.pdbx_PDB_model_num',i,'1')
   if mid not in models:models.append(mid)
   if len(models)<model or mid!=models[model-1]:continue
   alt=col('_atom_site.label_alt_id',i,'.')
   if alt not in ('.','?','A','1',''):continue
   row={'index':len(atoms),'cif_row':i,'serial':str(col('_atom_site.id',i)),'record':col('_atom_site.group_PDB',i),'name':col('_atom_site.label_atom_id',i),'resname':col('_atom_site.label_comp_id',i),'chain':col('_atom_site.label_asym_id',i),'resid':col('_atom_site.label_seq_id',i),'auth_chain':col('_atom_site.auth_asym_id',i),'auth_resid':col('_atom_site.auth_seq_id',i),'auth_name':col('_atom_site.auth_atom_id',i),'element':str(col('_atom_site.type_symbol',i)).upper(),'xyz':[float(col('_atom_site.Cartn_'+axis,i)) for axis in 'xyz']}
   atoms.append(row)
  nc=len(d.get('_struct_conn.id',[]))
  for i in range(nc):
   def c(k,default='?'):return (d['_struct_conn.'+k][i] if '_struct_conn.'+k in d else default)
   typ=c('conn_type_id')
   if typ.lower() not in ('covale','disulf','modres'):continue
   ends=[]
   for j in (1,2):
    t='ptnr'+str(j)+'_';ends.append({'chain':c(t+'label_asym_id'),'resid':c(t+'label_seq_id'),'name':c(t+'label_atom_id'),'auth_chain':c(t+'auth_asym_id'),'auth_resid':c(t+'auth_seq_id'),'auth_name':c(t+'label_atom_id')})
   connections.append({'source':'mmCIF_struct_conn','endpoints':ends,'id':c('id')})
  return atoms,connections,{'format':'mmCIF','models_seen':models,'selected_model_1based':model,'struct_conn_covalent_rows':nc}
 # RFAA stock PDB preserves atomized-token order but duplicates names/resids.
 model_id=0;explicit_models=False;selected=True;models=[]
 for line in path.read_text().splitlines():
  if line.startswith('MODEL'):
   explicit_models=True;model_id+=1;selected=model_id==model;models.append(model_id);continue
  if line.startswith('ENDMDL'):selected=False;continue
  if line.startswith(('ATOM  ','HETATM')) and selected:
   alt=line[16:17]
   if alt not in (' ','A','1'):continue
   name=line[12:16].strip();element=line[76:78].strip().upper() if len(line)>=78 else ''
   if not element:element=''.join(c for c in name if c.isalpha())[0:1].upper()
   atoms.append({'index':len(atoms),'serial':line[6:11].strip(),'record':line[:6].strip(),'name':name,'resname':line[17:20].strip(),'chain':line[21:22].strip(),'resid':line[22:26].strip(),'auth_chain':line[21:22].strip(),'auth_resid':line[22:26].strip(),'auth_name':name,'element':element,'xyz':[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]})
  if line.startswith('CONECT'):
   ids=[line[i:i+5].strip() for i in range(6,len(line),5) if line[i:i+5].strip()]
   if ids:connections.extend({'source':'PDB_CONECT','serials':[ids[0],v]} for v in ids[1:])
 return atoms,connections,{'format':'PDB','models_seen':models or [1],'selected_model_1based':model,'CONECT_rows':len(connections)}

def aliases(atom):
 result=[]
 for ch,re,na in [(atom.get('chain'),atom.get('resid'),atom.get('name')),(atom.get('auth_chain'),atom.get('auth_resid'),choose(atom.get('auth_name'),atom.get('name')))]:
  if not any(missing(x) for x in (ch,re,na)):result.append((str(ch),str(re),str(na)))
 return set(result)

def standard_mapping(method,case,mapping,atoms):
 native={};custom_res=set();errors=[];proxy_positions=set();mode='named native residues'
 if method=='rfaa':
  mode='RFAA ordinary residue names plus atomized-AA token/element-order mapping'
  rows=sorted(mapping['atomized_atoms'],key=lambda r:r.get('native_token_index_0based',r.get('native_token_index_0based_expected')))
  het=[a for a in atoms if a['record']=='HETATM' and a['element'] not in ('H','D')]
  if len(het)!=len(rows):errors.append({'code':'RFAA_atomized_emission_count_mismatch','expected':len(rows),'observed':len(het)})
  for a,r in zip(het,rows):
   if a['element']!=r['native_output_atom_name'].upper():errors.append({'code':'RFAA_atomized_order_or_element_mismatch','serial':a['serial'],'expected':r['native_output_atom_name'],'observed':a['element']})
   if r['kind']=='protein_atom':k=pkey(r['protein_residue_index_1based'],r['atom_name']);proxy_positions.add(int(r['protein_residue_index_1based']))
   else:k=gkey(r['glycan_id'],r['glycan_node_id'],r['atom_name'])
   native[a['index']]=k
  ordinary=mapping['ordinary_protein_positions_1based_by_native_token']
  for a in atoms:
   if a['record']!='ATOM' or a['element'] in ('H','D'):continue
   try:i=int(a['resid'])-1
   except ValueError:errors.append({'code':'RFAA_invalid_ordinary_residue_index','serial':a['serial']});continue
   if not 0<=i<len(ordinary):errors.append({'code':'RFAA_ordinary_token_out_of_bounds','serial':a['serial']});continue
   native[a['index']]=pkey(ordinary[i],a['name'])
  return native,errors,proxy_positions,mode
 look={};protein_chains=set()
 if method=='boltz2' and mapping.get('schema','').startswith('boltz2-glycoAA'):
  protein_chains.add(mapping['protein_rule']['native_chain'])
  for r in mapping['atoms']:
   pos=int(r['native_protein_position_1based']);custom_res.add((r['native_chain'],str(pos)));proxy_positions.add(pos)
   k=pkey(r['canonical_position_1based'],r['canonical_atom_name']) if r['canonical_kind']=='protein' else gkey(r['glycan_id'],r['node_id'],r['canonical_atom_name'])
   look[(r['native_chain'],str(pos),r['native_atom_name'])]=k
 else:
  protein_chains.update(p.get('af3_chain_id') for p in mapping.get('proteins',[]))
  protein_chains.update((mapping.get('protein_chain_map') or {}).values())
  for g in mapping.get('glycans',[]):
   chain=g.get('af3_chain_id',g.get('native_chain'))
   for n in g['nodes']:
    pos=n.get('af3_residue_index',n.get('draft_preorder_position_1based'))
    if pos is None:pos=int(n['native_monosaccharide_index_0based'])+1
    source=next(x for gly in case['glycans'] if gly['glycan_id']==g['glycan_id'] for x in gly['nodes'] if x['node_id']==n['node_id'])
    for name in source['expected_reacted_heavy_atom_names']+source['leaving_atom_names']:look[(chain,str(pos),name)]=gkey(g['glycan_id'],n['node_id'],name)
 for a in atoms:
  if a['element'] in ('H','D'):continue
  found={look[v] for v in aliases(a) if v in look}
  if len(found)>1:errors.append({'code':'ambiguous_native_atom_identity','serial':a['serial']});continue
  if found:native[a['index']]=next(iter(found));continue
  locs={(ch,re) for ch,re,na in aliases(a)}
  if locs & custom_res:errors.append({'code':'unexpected_atom_in_Boltz_custom_residue','serial':a['serial'],'native_name':a['name']});continue
  protein_alias=next(((ch,re,na) for ch,re,na in aliases(a) if ch in protein_chains and re.isdigit()),None)
  if protein_alias:native[a['index']]=pkey(protein_alias[1],protein_alias[2])
  else:errors.append({'code':'unmapped_output_heavy_atom','serial':a['serial'],'native_name':a['name'],'chain':a['chain'],'resid':a['resid']})
 return native,errors,proxy_positions,mode

AA_HEAVY={'ASN':{'N','CA','C','O','CB','CG','OD1','ND2'},'SER':{'N','CA','C','O','CB','OG'},'THR':{'N','CA','C','O','CB','OG1','CG2'}}
def canonical_mapping(method,case,mapping,atoms):
 if mapping.get('schema')!='boltz-parent-atomized-map-1':return standard_mapping(method,case,mapping,atoms)
 if method!='boltz2' or mapping['case_id']!=case['case_id']:raise ValueError('Parent atom map method/case mismatch')
 if mapping['chemical_input_sha256']!=case['chemical_input_sha256']:raise ValueError('Parent atom map chemical hash mismatch')
 pchain=mapping['protein_rule']['native_chain'];look={};expectgly={};observedgly=set();observedp=set();expectedp=set();glychains={}
 for g in case['glycans']:
  at=g['attachment'];expectedp.update((at['protein_residue_index'],n) for n in AA_HEAVY[at['protein_ccd_id']])
  for n in g['nodes']:
   for a in n['expected_reacted_heavy_atom_names']:expectgly[(g['glycan_id'],n['node_id'],a)]=n['ccd_id']
 for r in mapping['atoms']:
  ch=r['native_chain'];pos=int(r['native_residue_position_1based']);name=r['native_atom_name'];alias=(ch,str(pos),name)
  if alias in look:raise ValueError('Duplicate native atom identity')
  if not name or len(name)>4:raise ValueError('Invalid native atom name')
  if r['canonical_kind']=='protein':
   cp=int(r['canonical_position_1based']);ca=r['canonical_atom_name']
   if ch!=pchain or pos!=cp or name!=ca or (cp,ca) not in expectedp:raise ValueError('Unexpected anchor protein atom')
   key=pkey(cp,ca);observedp.add((cp,ca))
  elif r['canonical_kind']=='glycan':
   gid=r['glycan_id'];ck=(gid,r['node_id'],r['canonical_atom_name'])
   if ck not in expectgly or r['canonical_ccd']!=expectgly[ck]:raise ValueError('Unexpected canonical glycan atom')
   if ch==pchain or pos!=1:raise ValueError('Glycan must be separate NONPOLYMER chain at residue1')
   if gid in glychains and glychains[gid]!=ch:raise ValueError('Glycan spans multiple ligand chains')
   glychains[gid]=ch;key=gkey(*ck);observedgly.add(ck)
  else:raise ValueError('Unknown map atom kind')
  look[alias]=key
 if observedp!=expectedp or observedgly!=set(expectgly) or len(set(glychains.values()))!=len(case['glycans']):raise ValueError('Map does not exactly cover all anchors and glycans')
 native={};errors=[]
 for a in atoms:
  if a['element'] in ('H','D'):continue
  possibilities=aliases(a);found={look[v] for v in possibilities if v in look}
  if len(found)>1:errors.append({'code':'ambiguous_parent_atom_identity','serial':a['serial']});continue
  if found:native[a['index']]=next(iter(found));continue
  protein={(int(re),na) for ch,re,na in possibilities if ch==pchain and re.isdigit() and 1<=int(re)<=len(case['protein']['sequence'])}
  if len(protein)==1:
   pos,name=next(iter(protein));native[a['index']]=pkey(pos,name)
  else:errors.append({'code':'unmapped_or_ambiguous_output_heavy_atom','serial':a['serial'],'native_name':a['name'],'chain':a['chain'],'resid':a['resid']})
 return native,errors,set(),'canonical protein plus explicit per-atom multi-site N/O glycan map'
