UNIPROT_PATTERN = r"([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?"
RFAM_PATTERN = r"RF\d{5}"

RIBOSOMES_RFAM_MAPPING = {
    "Bacterial complete ribosome": ["RF00177", "RF02541"],
    "Eukaryotic complete ribosome": ["RF01960", "RF02543"],
    "Archeal complete ribosome": ["RF01989", "RF02540"],
}