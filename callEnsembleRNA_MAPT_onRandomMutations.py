import Functions_EnsembleRNA as frna
import pandas as pd

# Generate random mutations

WT_seq="AACGTCCAGTCCAAGTGTGGCTCAAAGGATAATATCAAACACGTCCCGGGAGGCGGCAGTGTGAGTACCTTCACACGTCCCATGCGCCGTGCTGTGGCTTGAATTATTAGGAAGTGGTGTGAGTGCGTACACTTGCGAGACACTGCATAGAATAAATCCTTCTTGGGCTCTCAGGATCTGGCTGCGACCTCTGGGTGAATG"

# Open a file to write 
with open("../tmp/Structures_forRandomMuts.csv","w") as rw:
# Go through every position
    for pos in range(len(WT_seq)):
        # Look for WT base at position
        WT_base = WT_seq[pos]
        # Get possible mutations at position
        MUT_bases = [i for i in ["A","G","T","C"] if i!=WT_base]
        print(MUT_bases)
        for mut in MUT_bases:
            name_MUT = "Position_"+str(pos)+"_"+WT_base+"-"+mut
            print(name_MUT)
            nums = frna.getPercentageChangesInEnsemble(WT_seq,mut,pos,"../tmp/BuildEnsembleFromMutationsSNPs/MAPThairpin_LaederachPrimerSet_SuboptimalStructuresWithMutsAndSNPs_deltaEnergy-2_Top500_UniqueDBs.db",[50,80])
            toWrite = [str(name_MUT),str(nums[0]),str(nums[1]),str(nums[2]),str(nums[3])]
            rw.write(",".join(toWrite))
            rw.write("\n")