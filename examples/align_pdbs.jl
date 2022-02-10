using IterativeClosestPoint
using BioStructures

ICP = IterativeClosestPoint

pdb_path = "C:\\Users\\Tom\\Downloads"
receptor1 = "Vmn1r78"
receptor2 = "Rhodopsin"

struc1 = read(joinpath(pdb_path,receptor1*".pdb"),PDB)
struc2 = read(joinpath(pdb_path,receptor2*".pdb"),PDB)

coords1 = hcat( [ struc1["A"].residues["$i"].atoms[" CA "].coords for i=1:length(struc1["A"].residues)]...)
coords2 = hcat( [ struc2["A"].residues["$i"].atoms[" CA "].coords for i=1:length(struc2["A"].residues)]...)

alignment = ICP.tm_align(coords1, coords2, )
tmscore = ICP.tm_score(coords1, coords2, alignment)