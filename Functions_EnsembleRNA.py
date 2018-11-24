import subprocess
import pandas as pd

def getNumberPerClusterEnsemble(fastaseqfile,projectedDBfile,rangeToCluster,optionalDMSFile=""):
    with open(fastaseqfile) as fw:
        lines = [line.strip() for line in fw]
    seq = lines[1]
    print(seq)
    print(len(seq))
    # Get 1000 suboptimal structures the fastasequence
    if optionalDMSFile!="":
        dbfileToOpen="seq_1000SuboptimalStructures_WithDMS.txt"
        with open("../tmp/"+dbfileToOpen,"w") as gws:
            subprocess.call(["RNAsubopt","--shape",optionalDMSFile,"-p","1000", "-i",fastaseqfile], stdout=gws)
    else:
        dbfileToOpen="seq_1000SuboptimalStructures_WithoutDMS.txt"
        with open("../tmp/"+dbfileToOpen,"w") as gw:
            subprocess.call(["RNAsubopt","-p","1000", "-i",fastaseqfile], stdout=gw)
    # Read in 1000 suboptimal structures 
    with open("../tmp/"+dbfileToOpen) as f:
        wholefile = f.read().strip()    
    # Split by sequences based on new line
    splitByNewLine = wholefile.split("\n")
    print(len(splitByNewLine))
    # Only grab lines that have . or ( or ) as the start character
    splitByNewLine_OnlyDBs = [i for i in splitByNewLine if i[0] in [".","(",")"]]
    print(len(splitByNewLine_OnlyDBs))
    # Only grab the first number of characters in the WT_seq 
    onlyDBs = [i[0:len(seq)] for i in splitByNewLine_OnlyDBs]
    print(len(onlyDBs))
    dbfileToWrite = dbfileToOpen.split(".")[0]+".db"
    # Write the suboptimal DB structures into a file
    with open("../tmp/"+dbfileToWrite,"w") as fw:
        #fw.write("\n".join(DBs_unique))
        fw.write("\n".join(onlyDBs))
        fw.write("\n")
    # Call ensemblerna on previously created db file
    if optionalDMSFile!="":
        subprocess.call(["ensemblerna", "-d","../tmp/"+dbfileToWrite,"-md",projectedDBfile,"-r",str(rangeToCluster[0]),str(rangeToCluster[1]),fastaseqfile,"../tmp/wtDMSdata_EnsembleRNA_output"])
        name_seq = fastaseqfile.split("/")[2].split(".")[0]
        seq_DB = pd.read_csv("../tmp/wtDMSdata_EnsembleRNA_output/"+name_seq+".csv",sep=",",header=0)
    else:
        subprocess.call(["ensemblerna", "-d","../tmp/"+dbfileToWrite,"-md",projectedDBfile,"-r",str(rangeToCluster[0]),str(rangeToCluster[1]),fastaseqfile,"../tmp/wtNoDMSdata_EnsembleRNA_output"])
        name_seq = fastaseqfile.split("/")[2].split(".")[0]
        seq_DB = pd.read_csv("../tmp/wtNoDMSdata_EnsembleRNA_output/"+name_seq+".csv",sep=",",header=0)
    # Read in seq db file
    Num_map_vector0 = seq_DB[seq_DB["vectorization"]=='[0]']["frequency"].values[0]
    Num_map_vector1 = seq_DB[seq_DB["vectorization"]=='[1]']["frequency"].values[0]
    return(Num_map_vector0,Num_map_vector1)

def getPercentageChangesInEnsemble(WTseq,mutation,mutposition,projectedDBfile,rangeToCluster):
    # First get Mut Seq
    mutSeq = WTseq[0:mutposition]+mutation+WTseq[mutposition+1:]
    # Write MutSeq into fastafile
    with open("../tmp/mutationSeq.fa","w") as fw:
        fw.write("> Mutation_Position"+str(mutposition)+"_"+WTseq[mutposition]+"->"+mutation)
        fw.write("\n")
        fw.write(mutSeq)
        fw.write("\n")
        # Get 1000 suboptimal structures for the mutation
    with open("../tmp/SeqWithMutation_1000SuboptimalStructures.txt","w") as gw:
        subprocess.call(["RNAsubopt", "-p","1000", "-i","../tmp/mutationSeq.fa"], stdout=gw)
    # Read in 1000 suboptimal structures 
    with open("../tmp/SeqWithMutation_1000SuboptimalStructures.txt") as f:
        wholefile = f.read().strip()    
    # Split by sequences based on new line
    splitByNewLine = wholefile.split("\n")
    print(len(splitByNewLine))
    # Only grab lines that have . or ( or ) as the start character
    splitByNewLine_OnlyDBs = [i for i in splitByNewLine if i[0] in [".","(",")"]]
    print(len(splitByNewLine_OnlyDBs))
    # Only grab the first number of characters in the WT_seq 
    onlyDBs = [i[0:len(WTseq)] for i in splitByNewLine_OnlyDBs]
    print(len(onlyDBs))
    # Only grab unique DBs
    #DBs_unique = list(set(onlyDBs))
    #print(len(DBs_unique))
    # Write the suboptimal DB structures into a file
    with open("../tmp/SeqWithMutation_1000SuboptimalStructures_UniqueDBs.db","w") as fw:
        #fw.write("\n".join(DBs_unique))
        fw.write("\n".join(onlyDBs))
        fw.write("\n")
    # Call ensemblerna on this mutation db file
    subprocess.call(["ensemblerna", "-d","../tmp/SeqWithMutation_1000SuboptimalStructures_UniqueDBs.db","-md",projectedDBfile,"-r",str(rangeToCluster[0]),str(rangeToCluster[1]),"../tmp/mutationSeq.fa","../tmp/mutation_EnsembleRNA_output"])
    # Read in map db csv
    map_DB = pd.read_csv("../tmp/mutation_EnsembleRNA_output/mutationSeq_map.csv",sep=",",header=0)
    # Read in mut csv
    mut_DB = pd.read_csv("../tmp/mutation_EnsembleRNA_output/mutationSeq.csv",sep=",",header=0)
    # Return the numbers 
    Num_map_vector0 = map_DB[map_DB["vectorization"]=='[0]']["frequency"].values[0]
    Num_map_vector1 = map_DB[map_DB["vectorization"]=='[1]']["frequency"].values[0]
    Num_mut_vector0 = mut_DB[mut_DB["vectorization"]=='[0]']["frequency"].values[0]
    Num_mut_vector1 = mut_DB[mut_DB["vectorization"]=='[1]']["frequency"].values[0]
    return(Num_map_vector0,Num_map_vector1,Num_mut_vector0,Num_mut_vector1)