import Functions_EnsembleRNA as frna
import pandas as pd
import sys

withDMS=sys.argv[1]

if withDMS=="True":
    with open("../tmp/MAPT_hairpin_CoordsForLaederachPrimer_IncludeDMS_1000Iterations.txt","w") as fdms:
        for i in range(1000):
            numClusters = frna.getNumberPerClusterEnsemble("../data/MAPT_hairpin_CoordsForLaederachPrimer.fa","../tmp/BuildEnsembleFromMutationsSNPs/MAPThairpin_LaederachPrimerSet_SuboptimalStructuresWithMutsAndSNPs_deltaEnergy-2_Top500_UniqueDBs.db",[55,75],"../tmp/MAPT_amplicon_HCC1500_Take2/Pipeline_MAPT-Exon10Plus300basesintoIntron10.shape")
            fdms.write(str(numClusters[0])+"\t"+str(numClusters[1]) +"\n")
else:
    with open("../tmp/MAPT_hairpin_CoordsForLaederachPrimer_NoDMS_1000Iterations.txt","w") as f:
        for i in range(1000):
            numClusters = frna.getNumberPerClusterEnsemble("../data/MAPT_hairpin_CoordsForLaederachPrimer.fa","../tmp/BuildEnsembleFromMutationsSNPs/MAPThairpin_LaederachPrimerSet_SuboptimalStructuresWithMutsAndSNPs_deltaEnergy-2_Top500_UniqueDBs.db",[55,75])
            f.write(str(numClusters[0])+"\t"+str(numClusters[1]) +"\n")