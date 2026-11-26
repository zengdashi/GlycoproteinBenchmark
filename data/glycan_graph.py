"""Validate atom identities, glycan connectivity and frozen chemical fingerprints."""
from __future__ import annotations
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'glyco-connectivity-1.0'
AA = dict(zip('ACDEFGHIKLMNPQRSTVWY', [
    'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
    'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']))
ANCHORS = {'ASN': ('N', 'ND2', 'N-Glycosylation'),
           'SER': ('S', 'OG', 'O-Glycosylation'),
           'THR': ('T', 'OG1', 'O-Glycosylation'),
           'TYR': ('Y', 'OH', 'O-Glycosylation')}
FORBIDDEN_KEYS = {'xyz', 'coordinates', 'Cartn_x', 'Cartn_y', 'Cartn_z',
                  'pdb_id', 'auth_seq_id', 'auth_asym_id', 'label_asym_id',
                  'protein_auth_seq_id', 'mmcif_path', 'source_path'}

class AnnotationError(ValueError):
    """An annotation is inconsistent or cannot be safely converted."""

def require(ok, message):
    if not ok:
        raise AnnotationError(message)

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def component_digest(component):
    return digest({k:component[k] for k in ('ccd_id','name','type','atoms','bonds')})

def _no_truth_or_source_fields(value, path='document'):
    if isinstance(value,dict):
        for key,v in value.items():
            require(key not in FORBIDDEN_KEYS, f'{path}: source/coordinate field {key} is not permitted')
            _no_truth_or_source_fields(v,path+'.'+key)
    elif isinstance(value,list):
        for i,v in enumerate(value):
            _no_truth_or_source_fields(v,f'{path}[{i}]')

def _string(value,label):
    require(isinstance(value,str) and bool(value),label+' must be a nonempty string')
    return value

def _unique_strings(value,label):
    require(isinstance(value,list),label+' must be a list')
    require(all(isinstance(x,str) and x for x in value),label+' contains invalid strings')
    require(len(value)==len(set(value)),label+' contains duplicates')
    return set(value)

def _one_based(value,length,label):
    require(type(value) is int and 1 <= value <= length,label+' must be a 1-based sequence index')
    return value

def _component(component,key):
    require(component.get('ccd_id')==key,'CCD dictionary key/identity mismatch: '+key)
    require(re.fullmatch('[A-Z0-9]+',key) is not None,'Invalid CCD ID '+key)
    require(component.get('chemical_definition_sha256')==component_digest(component),'Chemical definition SHA256 mismatch for '+key)
    atoms={};adj=defaultdict(set)
    require(isinstance(component.get('atoms'),list) and component['atoms'],'Empty CCD atom dictionary '+key)
    for atom in component['atoms']:
        name=_string(atom.get('atom_name'),key+' atom_name')
        require(name not in atoms,'Duplicate CCD atom '+key+':'+name)
        require(atom.get('element') in {'H','D','C','N','O','S','P','F','CL','BR','I','SE'},'Unsupported/missing CCD element '+key+':'+name)
        require(type(atom.get('formal_charge')) is int,'Missing CCD formal charge '+key+':'+name)
        require(type(atom.get('leaving_atom')) is bool,'Missing CCD leaving flag '+key+':'+name)
        require(atom.get('stereo_config') in {'R','S','N'},'Invalid/missing CCD atom stereochemistry '+key+':'+name)
        require(type(atom.get('aromatic')) is bool,'Missing CCD aromatic flag '+key+':'+name)
        atoms[name]=atom
    seen=set()
    require(isinstance(component.get('bonds'),list),'Missing CCD bonds '+key)
    for bond in component['bonds']:
        a,b=bond.get('atom_1'),bond.get('atom_2')
        require(a in atoms and b in atoms and a!=b,'Invalid CCD bond endpoint '+key)
        pair=frozenset((a,b));require(pair not in seen,'Duplicate CCD bond '+key);seen.add(pair)
        require(bond.get('bond_order') in {'SING','DOUB','TRIP','AROM'},'Unsupported CCD bond order '+key)
        require(bond.get('stereo_config') in {'N','E','Z'},'Invalid/missing CCD bond stereochemistry '+key)
        require(type(bond.get('aromatic')) is bool,'Missing CCD bond aromatic flag '+key)
        adj[a].add(b);adj[b].add(a)
    return {'component':component,'atoms':atoms,'adj':adj}

def _chemical_key(case):
    """Match the original frozen input-key definition, excluding source IDs."""
    glycans=[]
    for g in case['glycans']:
        nodes={n['node_id']:n for n in g['nodes']};adj=defaultdict(list)
        for b in g['bonds']:
            adj[b['node_1']].append((b['node_2'],b['atom_1'],b['atom_2'],b['bond_order']))
            adj[b['node_2']].append((b['node_1'],b['atom_2'],b['atom_1'],b['bond_order']))
        def walk(n,parent=None):
            kids=[[a,b,o,walk(v,n)] for v,a,b,o in adj[n] if v!=parent]
            kids.sort(key=lambda x:json.dumps(x,sort_keys=True))
            return [nodes[n]['ccd_id'],kids]
        a=g['attachment']
        glycans.append([a['protein_residue_index'],a['protein_ccd_id'],a['protein_atom'],a['glycan_atom'],a['bond_order'],walk(a['glycan_node_id'])])
    glycans.sort(key=lambda x:json.dumps(x,sort_keys=True))
    return {'sequence':case['protein']['sequence'],
            'modifications':sorted(case['protein']['modifications'],key=lambda m:(m['position'],m['ccd_id'])),
            'glycans':glycans,'additional_covalent_bonds':[]}

def chemical_input_digest(case):
    return digest(_chemical_key(case))

def _validate_case(case,chemistry):
    cid=_string(case.get('case_id'),'case_id')
    require(re.fullmatch(r'GLYCO_[0-9a-f]{16}',cid) is not None,'Invalid anonymous case_id '+cid)
    protein=case.get('protein');require(isinstance(protein,dict),cid+': missing protein')
    require(protein.get('id')=='P1',cid+': only anonymous single-chain P1 is supported')
    seq=protein.get('sequence');require(isinstance(seq,str) and seq and set(seq)<=set(AA),cid+': sequence must contain only canonical uppercase 20AA')
    query=case.get('msa_query');require(isinstance(query,dict),cid+': missing msa_query')
    seqsha=hashlib.sha256(seq.encode()).hexdigest()
    require(query.get('sequence_sha256')==seqsha,cid+': MSA query sequence SHA mismatch')
    require(query.get('sequence_id')=='GPSEQ_'+seqsha[:16],cid+': MSA query ID does not match sequence')
    mods={};require(isinstance(protein.get('modifications'),list),cid+': missing protein modifications list')
    for mod in protein['modifications']:
        pos=_one_based(mod.get('position'),len(seq),cid+': protein modification')
        require(pos not in mods,cid+': duplicate protein modification position')
        require(mod.get('ccd_id') in chemistry,cid+': missing modified-protein CCD definition')
        mods[pos]=mod['ccd_id']
    glycans=case.get('glycans');require(isinstance(glycans,list) and glycans,cid+': no glycan components')
    allnodes=set();gids=set();attachsites=set();counts=Counter()
    for g in glycans:
        gid=_string(g.get('glycan_id'),cid+': glycan_id')
        require(gid not in gids,cid+': duplicate glycan_id');gids.add(gid)
        local={};nodes=g.get('nodes');require(isinstance(nodes,list) and nodes,cid+': empty glycan '+gid)
        for n in nodes:
            nid=_string(n.get('node_id'),cid+': node_id')
            require(nid not in allnodes,cid+': duplicate glycan node ID');allnodes.add(nid);local[nid]=n
            require(n.get('ccd_id') in chemistry,cid+': missing glycan CCD definition '+str(n.get('ccd_id')))
            _unique_strings(n.get('leaving_atom_names'),cid+': leaving_atom_names')
            _unique_strings(n.get('expected_reacted_heavy_atom_names'),cid+': reacted heavy atoms')
            _unique_strings(n.get('observed_atom_names'),cid+': observed atom names')
        root=g.get('root_node_id');require(root in local,cid+': invalid glycan root')
        adj=defaultdict(set);external=defaultdict(set);edges=set();bonds=g.get('bonds')
        require(isinstance(bonds,list),cid+': missing internal glycan bonds')
        def endpoint(nid,atom):
            require(nid in local,cid+': internal bond refers to another/unknown glycan node')
            _string(atom,cid+': bond atom');n=local[nid];ccd=chemistry[n['ccd_id']]
            require(atom in ccd['atoms'],cid+': named atom absent from CCD '+n['ccd_id']+':'+atom)
            require(atom not in n['leaving_atom_names'],cid+': bound atom is marked for deletion')
            require(atom in n['observed_atom_names'],cid+': covalent endpoint not observed')
            require(atom not in external[nid],cid+': one atom has multiple external covalent edges')
            external[nid].add(atom)
        for b in bonds:
            n1,n2=b.get('node_1'),b.get('node_2')
            require(n1!=n2,cid+': inter-residue self bond')
            require(b.get('bond_order')=='SING',cid+': unsupported inter-residue bond order')
            pair=frozenset(((n1,b.get('atom_1')),(n2,b.get('atom_2'))))
            require(pair not in edges,cid+': duplicate internal covalent edge');edges.add(pair)
            endpoint(n1,b.get('atom_1'));endpoint(n2,b.get('atom_2'))
            adj[n1].add(n2);adj[n2].add(n1)
        visited=set();todo=[root]
        while todo:
            n=todo.pop()
            if n not in visited:visited.add(n);todo.extend(adj[n]-visited)
        require(visited==set(local) and len(bonds)==len(local)-1,cid+': glycan must be a connected acyclic graph')
        at=g.get('attachment');require(isinstance(at,dict),cid+': missing protein attachment')
        require(at.get('protein_chain_id')=='P1',cid+': attachment protein is not P1')
        require(at.get('glycan_node_id')==root,cid+': root differs from attached glycan node')
        require(at.get('bond_order')=='SING',cid+': unsupported protein-glycan bond order')
        pos=_one_based(at.get('protein_residue_index'),len(seq),cid+': attachment position')
        expected=mods.get(pos,AA[seq[pos-1]])
        require(at.get('protein_ccd_id')==expected,cid+': attachment residue CCD differs from sequence/PTM')
        if expected in ANCHORS:
            aa,atom,linktype=ANCHORS[expected]
            require(seq[pos-1]==aa and at.get('protein_atom')==atom,cid+': wrong protein attachment atom')
            require(g.get('link_type')==linktype,cid+': glycosylation type/anchor mismatch')
        else:
            require(expected in mods.values(),cid+': unreviewed non-N/O protein attachment')
            require(at.get('protein_atom') in chemistry[expected]['atoms'],cid+': modified protein attachment atom absent from CCD')
            require(g.get('link_type') in {'N-Glycosylation','O-Glycosylation'},cid+': unsupported glycosylation class')
        site=(pos,at.get('protein_atom'))
        require(site not in attachsites,cid+': multiple glycans use the same protein atom');attachsites.add(site)
        endpoint(root,at.get('glycan_atom'))
        for nid,n in local.items():
            cc=chemistry[n['ccd_id']];atoms=cc['atoms'];gone=set(n['leaving_atom_names'])
            require(gone<=set(atoms),cid+': removed atom is absent from CCD')
            required=set()
            for bound in external[nid]:
                candidates={a for a in cc['adj'][bound] if atoms[a]['leaving_atom'] and atoms[a]['element'] not in {'H','D'}}
                require(len(candidates)<=1,cid+': ambiguous heavy leaving group')
                for old in candidates:
                    require(old not in external[nid],cid+': leaving atom is another bond endpoint')
                    require(sum(atoms[x]['element'] not in {'H','D'} for x in cc['adj'][old])==1,cid+': nonterminal heavy leaving group')
                    required.add(old);required.update(x for x in cc['adj'][old] if atoms[x]['element'] in {'H','D'})
            require(gone==required,cid+': per-node leaving-atom inventory differs from explicit reaction')
            expectedheavy={a for a,x in atoms.items() if x['element'] not in {'H','D'}}-gone
            require(set(n['expected_reacted_heavy_atom_names'])==expectedheavy,cid+': reacted heavy-atom inventory mismatch')
            require(set(n['observed_atom_names'])<=set(atoms),cid+': observed atom absent from CCD')
            require(not (gone&set(n['observed_atom_names'])),cid+': experimentally observed leaving atom')
        counts['glycans']+=1;counts['sugar_nodes']+=len(nodes);counts['internal_bonds']+=len(bonds);counts[g['link_type']]+=1
    fingerprint=chemical_input_digest(case)
    require(case.get('chemical_input_sha256')==fingerprint,cid+': chemical input SHA256 mismatch')
    require(cid=='GLYCO_'+fingerprint[:16],cid+': case ID does not match chemical input')
    return counts

def validate(document):
    require(isinstance(document,dict) and document.get('schema_version')==SCHEMA,'Unsupported annotation schema')
    _no_truth_or_source_fields(document)
    components=document.get('ccd_components');require(isinstance(components,dict) and components,'Missing CCD definitions')
    companion=document.get('ccd_companion')
    if companion is not None:
        require(isinstance(companion,dict),'Invalid CCD companion descriptor')
        filename=_string(companion.get('filename'),'CCD companion filename')
        require(Path(filename).name==filename,'CCD companion filename must be a local basename')
        require(isinstance(companion.get('sha256'),str) and re.fullmatch('[0-9a-f]{64}',companion['sha256']) is not None,'Invalid CCD companion SHA256')
    chemistry={k:_component(v,k) for k,v in components.items()}
    cases=document.get('cases');require(isinstance(cases,list) and cases,'Missing cases')
    ids=set();queries=set();counts=Counter()
    for case in cases:
        cid=case.get('case_id');require(cid not in ids,'Duplicate case ID');ids.add(cid)
        counts.update(_validate_case(case,chemistry));queries.add(case['msa_query']['sequence_id'])
    datasets=document.get('datasets');require(isinstance(datasets,dict) and datasets,'Missing dataset memberships')
    ds={}
    for name,value in datasets.items():
        require(isinstance(value,dict),'Invalid dataset membership '+name)
        members=_unique_strings(value.get('case_ids'),'Dataset '+name)
        require(members<=ids,'Dataset references unknown case IDs '+name)
        ds[name]={'cases':len(members),'queries':len({c['msa_query']['sequence_id'] for c in cases if c['case_id'] in members})}
    return {'status':'pass','schema_version':SCHEMA,'cases':len(cases),'unique_queries':len(queries),
            'CCD_definitions':len(components),'datasets':ds,**dict(counts),
            'source_ids_and_target_coordinates_present':False,'native_featurisation_tested':False,'model_inference_run':False}

def resolve_ccd_companion(document,annotations_path,explicit=None):
    """Validate a companion before use; JSON-only graph validation needs no file."""
    descriptor=document.get('ccd_companion')
    require(isinstance(descriptor,dict),'AF3 export requires a declared, hash-pinned CCD companion')
    path=Path(explicit) if explicit is not None else Path(annotations_path).parent/descriptor['filename']
    require(path.is_file(),'AF3 CCD companion file is missing')
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    require(actual==descriptor['sha256'],'AF3 CCD companion SHA256 mismatch')
    return path
