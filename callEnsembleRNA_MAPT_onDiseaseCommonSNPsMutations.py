import Functions_EnsembleRNA as frna
import pandas as pd


# Go through all the mutations and SNPs and get number of structures for vector 0 and vector 1
# Read in the mutation and common SNPs around MAPT hairpin file
start_Coord=46010343
WT_seq="AACGTCCAGTCCAAGTGTGGCTCAAAGGATAATATCAAACACGTCCCGGGAGGCGGCAGTGTGAGTACCTTCACACGTCCCATGCGCCGTGCTGTGGCTTGAATTATTAGGAAGTGGTGTGAGTGCGTACACTTGCGAGACACTGCATAGAATAAATCCTTCTTGGGCTCTCAGGATCTGGCTGCGACCTCTGGGTGAATG"
muts_snps = pd.read_csv("../processed_data/MutsAndSNPsAroundTauHairpin_ForLaederachPrimerSet.bed",sep="\t",header=None)

with open("../tmp/Structures_forMuts.csv","w") as sw:
    for i in range(muts_snps.shape[0]):
        # Write the coordinates of mutation in the first line
        positionMut=muts_snps.iloc[i,1]-(start_Coord-1)
        mut = muts_snps.iloc[i,5]
        name_Mut = muts_snps.iloc[i,3]
        print(name_Mut)
        nums = frna.getPercentageChangesInEnsemble(WT_seq,mut,positionMut,"../tmp/BuildEnsembleFromMutationsSNPs/MAPThairpin_LaederachPrimerSet_SuboptimalStructuresWithMutsAndSNPs_deltaEnergy-2_Top500_UniqueDBs.db",[50,80])
        toWrite = [str(name_Mut),str(nums[0]),str(nums[1]),str(nums[2]),str(nums[3])]
        sw.write(",".join(toWrite))
        sw.write("\n")