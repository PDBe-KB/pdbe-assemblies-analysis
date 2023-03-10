import re
import csv
from scripts.constants import UNIPROT_PATTERN, RFAM_PATTERN, RIBOSOMES_RFAM_MAPPING
from collections import Counter

SYMMETRY_MAPPING = {}
with open("data/symmetry_reference.csv") as f:
  csv_reader = csv.reader(f, delimiter=",")
  next(csv_reader)
  for line in csv_reader:
    SYMMETRY_MAPPING[line[0]] = line[1]

def check_for_protein(assembly_string):
  if re.search(UNIPROT_PATTERN, assembly_string):
    return True
  elif "Protein_" in assembly_string:
    return True
  elif "antibody_" in assembly_string:
    return True
  else:
    return False

def check_for_rna(assembly_string):
  if re.search(RFAM_PATTERN, assembly_string):
    return True
  elif "RNA_" in assembly_string:
    return True
  else:
    return False

def check_for_dna(assembly_string):
  if "DNA" in assembly_string:
    return True
  else:
    return False

def validate_uniprot(assembly_string, query_type):
  """
  Returns whether any or all components in an assembly are mapped
  to UniProt
  """  
  assembly_components = assembly_string.split(",")
  uniprot_matches = [bool(re.search(UNIPROT_PATTERN, component)) for component in assembly_components]
  return query_type(uniprot_matches)

def assembly_composition(valid_protein, valid_RNA, valid_DNA):
    """
    Returns the assembly polymer composition
    """
    composition = set()
    
    if valid_protein:
        composition.add("protein")
        
    if valid_RNA:
        composition.add("RNA")
        
    if valid_DNA:
        composition.add("DNA")
        
    composition = list(composition)
    
    return ",".join(composition)

def get_asymmetrical_assemblies(assembly, sym_operator):
  asymmetrical_assemblies = []
  assembly = assembly.split(",")
  sym_operator = sym_operator.split(",")
  for assembly, sym_operator in zip(assembly, sym_operator):
    if sym_operator == "no-sym":
      asymmetrical_assemblies.append(assembly)
  return ", ".join(asymmetrical_assemblies)

def get_extended_symmetry_operators(assembly, sym_operator):
  assemblies = []
  assembly = assembly.split(",")
  sym_operator = sym_operator.split(",")
  for assembly, sym_operator in zip(assembly, sym_operator):
    assemblies.append(f"{assembly}|{sym_operator}")
  return ", ".join(assemblies)

def most_frequent_symmetry(symmetry_list):
  symmetry_list = symmetry_list.split(",")
  occurence_count = Counter(symmetry_list)
  return occurence_count.most_common(1)[0][0]

def get_unique_symmetries(symmetry_list):
  result = set(symmetry_list.split(","))
  return ", ".join(list(result))

def get_assembly_type(assembly_string):
    """
    Returns the assembly type whether
    it's a homomeric, heteromeric or
    monomeric
    """
    assembly_components = assembly_string.split(",")
    if len(assembly_components) > 1:
        assembly_type = "heteromeric"
    else:
        *_, stoic = assembly_string.split("_")
        if stoic == str(1):
          assembly_type = "monomeric"
        else:
          assembly_type = "homomeric"
    return assembly_type    

def group_UniProt_accessions(data):
  grouped_UNP_accessions = {}
  for elem in data:
    accession, stoichiometry = elem.split("_")
    grouped_UNP_accessions.setdefault(accession, []).append(stoichiometry)
  return grouped_UNP_accessions

def group_subassemblies(data):
  grouped_subassemblies = {}
  grouped_superassemblies = {}
  for elem in data:
    grouped_subassemblies.setdefault(elem[0], []).append(elem[1])
    grouped_superassemblies.setdefault(elem[1], []).append(elem[0])
  return grouped_subassemblies, grouped_superassemblies

def get_sym_op(data):
  assemblies_list = data.split(",")
  sym_op_list = [SYMMETRY_MAPPING.get(assembly, "no-sym") for assembly in assemblies_list]
  return ",".join(sym_op_list)

def validate_consistent_symmetry(data):
  symmetry_list = data.split(",")
  if len(set(symmetry_list)) == 1:
    return True
  else:
    return False

def count_unique_symmetry(data):
  symmetry_list = data.split(",")
  return len(set(symmetry_list))

def count_assemblies(data):
  assemblies_list = data.split(",")
  return len(assemblies_list)

def count_unique_pdb(data):
  assemblies_list = data.split(",")
  pdb_list = [assembly.split("_")[0] for assembly in assemblies_list]
  return len(set(pdb_list))

def check_complete_ribosome(assembly_string):
  for ribosome_name, ribosome_components in RIBOSOMES_RFAM_MAPPING.items():
    if all(component in assembly_string for component in ribosome_components):
      return ribosome_name
  return None

def get_sym_variant(ref_symmetry, assemblies):
  assemblies = assemblies.split(",")
  variants = []
  for assembly in assemblies:
    _, sym = assembly.split("|")
    if sym != ref_symmetry:
      variants.append(assembly)
  return ", ".join(variants)

def dict_compare(d1, d2):
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    modified = {o : (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    same = set(o for o in shared_keys if d1[o] == d2[o])
    return added, removed, modified, same




