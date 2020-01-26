from Entries_calc_Pit import fasta_to_entry_pit
from Species_calc_Pit import *

fasta_to_entry_pit("uniprot_sprot.fasta", "Entries sprot.csv")  # insert sprot.fasta input file and the output file name
fasta_to_entry_pit("uniprot_trembl.fasta", "Entries trembl.csv")  # insert trembl.fasta input file and the output file name
entry_to_species()  # insert entries input file and the output file name
species_pit()  # insert species file and the output file name
