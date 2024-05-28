import argparse
from Bio.PDB import *


scriptparser = argparse.ArgumentParser()
scriptparser.add_argument("--filename",
                    help="structure file",
                    type=str, required=True)

scriptparser.add_argument("--file_format",
                    help="can be pdb or mmcif, if not specified it is parsed from file",
                    type=str, required=False)

args = scriptparser.parse_args()

file_format = args.file_format
if file_format=="pdb":
    parser = PDBParser()
elif file_format=="mmcif":
    parser = MMCIFParser()
else:
    file_format = args.filename.split(".")[-1]
    if file_format == "pdb":
        parser = PDBParser()
    elif file_format == "mmcif":
        parser = MMCIFParser()
    elif file_format=="cif":
        parser = MMCIFParser()
    else:
        print("file must be .pdb or .cif")
        exit()


class VeryHighSelect(Select):
    def accept_atom(self, atom):
        if atom.get_bfactor() >90:
            return 1
        else:
            return 0

class HighSelect(Select):
    def accept_atom(self, atom):
        if atom.get_bfactor() >70:
            return 1
        else:
            return 0

class LowSelect(Select):
    def accept_atom(self, atom):
        if atom.get_bfactor() >50:
            return 1
        else:
            return 0

jobname = args.filename.split(".")[0]
structure = parser.get_structure('structure_id', args.filename)


io = PDBIO()

# Write out the selected atoms to a new PDB file
file_name_very_confident = jobname + "_very_confident.pdb"
file_name_confident = jobname + "_confident.pdb"
file_name_low_confidence = jobname + "_low_confidence.pdb"

io.set_structure(structure)
io.save(file_name_very_confident, VeryHighSelect())
io.save(file_name_confident, HighSelect())
io.save(file_name_low_confidence, LowSelect())