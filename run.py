import ncbi_genome_download as ngd
import os, re, gzip
from ete3 import NCBITaxa, Tree
import  sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from progressbar import ProgressBar
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

'''
SUEDE: Strain-level UniquE Dna probE finder
@author: Yilei Fu
@Email: yf20@rice.edu
'''
ncbi = NCBITaxa()                                                               # Build NCBI Taxon Database
pbar = ProgressBar()
# workpath = os.path.join("../" + "NCBITaxa/")
# controlpath = os.path.join("./example_control")
# ntpath =  /home/Users/yf20/ncbi_database/newnt/nt

def parseArgs(argv):
    """Function for parsing Arguments
    
    Args:
        argv: The arguments sent into this program
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser(description = "SUEDE: Strain-level UniquE \
        Dna probE finder  !!!IMPORTANT!!!: Do not use remote with blastn if you are \
        using blast 2.9.0. This is a bug in blast -remote, will be fixed in \
        blast 2.10.0. Blast nt database download instructions: \
        https://www.ncbi.nlm.nih.gov/books/NBK537770/ ")
    parser.add_argument("-w", "--work_dir",  type=str, 
    help="Directory of work directory, default: ../NCBITaxa/", 
    default="../NCBITaxa/")
    parser.add_argument("-c", "--config_dir", type=str, required=True,
    help = "Directory of configuration file")
    parser.add_argument("-d", "--database_path", type=str,
    help = "Directory of NCBI nt database, defalult: remote (!!!Do not use \
         remote with blastn 2.9.0, there is a bug in that version, \
            will be fixed in 2.10.0!!!) ")
    parser.add_argument("-t", "--tree_help", 
    help = "Use parsnp to generate a tree to evaluate positive control? \
        defalut: False", action="store_true")
    parser.add_argument("-p", "--processes", type = int, 
    help = "Threads number, default: 70", default=70)
    parser.add_argument("-m", "--MUMS_only", 
    help = "Only get MUMs, not blast against database? \
    defalut: False", action="store_true")

    args = parser.parse_args(argv)
    return args


def getTaxid(namelist):
    """Get Taxon id from NCBI taxon databse
    Use try except to return wrongly input name, and get accession ids fron the 
    name
    Args:
        namelist: the list of target names
    Returns:
        return the list of taxon ids from the list name
    """    
    accessid = []
    for i in namelist:
        name2taxid = ncbi.get_name_translator([i])
        if name2taxid == {}:
            print("Wrong Taxon name: " + i + "!")
            exit()
            return
        else:
            accessid.append(name2taxid)
    return accessid

def ungz_all_fasta(top_dir):
    """ungip all fasta.gz files
    ncbi database have default gz download filetype, recursively ungzip
    Args:
        top_dir: the directory of downloaded files
    Returns:
        nonthing 
    """   
    if os.path.isdir(top_dir):
        if os.listdir(top_dir) == []:
            return
        for i in os.listdir(top_dir):
            if i.split(".")[-2:] == ["fna", "gz"]:
                os.system("gunzip " + top_dir + "/" + i)
                continue
            ungz_all_fasta(top_dir  + "/" + i)
    else:
        return


def rm_not_fasta(current_dir):
    """delete all non-fasta files in the directory

    Args:
        current_dir: the directory
    Returns:
        nonthing 
    """   
    remove_list = []
    for i in os.listdir(current_dir):
        if i[-4:] != ".fna":
            remove_list.append(i)
    for i in remove_list:
        os.system("rm -r " + current_dir + i)


def download_db(taxnamelist, group):
    """download the NCBI fasta files from name
    use ncbi-gneome-download to download sequences
    Args:
        taxnamelist: the list of taxon names
        group: the name of taxon groups
    Returns:
        nonthing 
    """   
    for i in getTaxid(taxnamelist):
        Taxon = list(i.keys())[0]
        Taxonid = str(list(i.values())[0][0])
        outdir = workpath + "_".join(Taxon.split(" ")) + "/"
        try:
            os.mkdir(outdir)
        except FileExistsError:
            print("Path exists: "+ outdir)
        print("#############################################################")
        print("Downloading complete sequence in fasta from NCBI database...\n" + 
            "Taxon: " + Taxon + "\n" + 
            "Taxon id: " + Taxonid + "\n" + 
            "Directory: " + outdir)
        print("Executing: " + "ncbi-genome-download -t " + Taxonid +  \
            " -F fasta -l complete "  +" -o " + outdir + " " + \
                group)
        os.system("ncbi-genome-download -t " + Taxonid +  \
            " -F fasta -l complete "  +" -o " + outdir + " " + \
                group)
        print("...Taxon " + Taxon + " downloaded complete!")
        print("Unzip and re-organizing...")
        ungz_all_fasta(outdir)
        for i in os.walk(outdir):
            for j in i[2]:
                if j[-4:] == ".fna":
                    os.system("cp " + i[0]+"/"+j + " " + outdir)
        rm_not_fasta(outdir)



def getinfo(config):
    """get taxon informations form the inpuut configuration file

    Args:
        config: the input configuration file
    Returns:
        positive: the list of positive taxon names
        positive_paths: the path of positive taxons
        negative: the list of negative taxon names
        negative_paths: the path of negative taxons
        whole: combined list
        group: taxon group
    """   
    output = [[]]
    for x in config:
        output[-1].append(x)
        if x == '\n':
            output.append([])
    output[0].remove('\n')
    positive = []
    negative = []
    positive_paths = []
    negative_paths = []
    whole = []
    group = output[0][0][:-1]
    for i in output[1][:-1]:
        if i[-1] == '\n':
            positive.append(i[:-1])
            positive_paths.append(workpath + "_".join(i[:-1].split(" ")) + "/")
            whole.append(i[:-1])
        else:
            positive.append(i)
            positive_paths.append(workpath + "_".join(i.split(" ")) + "/")
            whole.append(i)
    for i in output[2]:
        if i[-1] == '\n':
            negative.append(i[:-1])
            negative_paths.append(workpath + "_".join(i[:-1].split(" ")) + "/")
            whole.append(i[:-1])
        else:
            negative.append(i)
            negative_paths.append(workpath + "_".join(i.split(" ")) + "/")
            whole.append(i)
    for i in positive:
        if i in negative:
            print("Positive and Negative control sets cannot have\
 same element! Exiting")
            exit
    return(positive, positive_paths, negative, negative_paths, whole, group)



def detectdataexist(path):
    """detect if the data exist

    Args:
        path: the directory of input file
    Returns:
        True/False
    """   
    if "positive" in os.listdir(path) and "negative" in os.listdir(path):
        return True
    else:
        return False


def generate_read_list(reads, MUML):
    """generate the list of reads with names

    Args:
        reads: reads
    Returns:
        name_list: names
        read_list: reads
    """   
    name_list = []
    read_list = []
    tempstr = ""
    for i in range(len(reads)):
        if "#SequenceCount" in reads[i]:
            sc = int(reads[i].split()[1])
        if reads[i][0] == "#" or reads[i][0] == "=":
            continue
        if reads[i][0] == ">":
            name_list.append(reads[i][:-1])
            # print(reads[i])
        elif i == len(reads) - 2:
            tempstr = tempstr+reads[i][:-1]
            read_list.append([i for i in re.split("A|C|T|G|-",tempstr) if i] )
            tempstr = ""
        elif reads[i+1][0] == ">":
            tempstr = tempstr+reads[i][:-1]
            read_list.append([i for i in re.split("A|C|T|G|-",tempstr) if i] )
            tempstr = ""
        elif reads[i+1][0] == "=":
            tempstr = tempstr+reads[i][:-1]
            read_list.append([i for i in re.split("A|C|T|G|-",tempstr) if i] )
            tempstr = ""
        else:
            tempstr = tempstr+reads[i][:-1]
    final_reads = []
    for i in range(len(read_list)):
        if i % sc == 0:
            final_reads+=read_list[i]
    clusters = {}
    for i in range(len(final_reads)):
        if len(final_reads[i]) > MUML:
            clusters.update({"cluster"+str(i):final_reads[i]})
    return clusters

def percenttostrain(filepath, blastresultfile, clusters):
    """from file to get a dictionary that contains the alignment score of
    strains VS MUMs

    Args:
        filepath: path to files
        blastresultfile: path to the blast results
        clusters: MUMs clusters(parsnp result)
    Returns:
        percent2strain: a dictionary stors the percentage of one MUM against 
        strains
        strain: a list of strains
    """   
    percent2strain = []
    strain = []
    for i in blastresultfile:
        with open (filepath + i, "r") as f:
            l = f.readlines()
            content = []
            p2s = {}
            for j in l:
                k = j.split("\t")
                content.append(k)
                stranname = k[-1][:-1]
                if stranname not in strain:
                    strain.append(stranname)
                strpercent = k[2]
                if stranname not in p2s:
                    p2s.update({stranname: float(strpercent)})
                read = clusters[k[0]]
            p2s.update({"MUM":read})
            percent2strain.append((p2s))
    return (percent2strain, strain)

def sortbyMUMlength(form):
    """sort MUMs by length 
    Args:
        form: pandas dataframe
    Returns:
        The sorted pandas dataframe
    """   
    form['length'] = form.index.str.len()
    form.sort_values('length', ascending=False, inplace=True)
    return(form.drop(["length"], axis=1))

def printTree():
    '''
    Generate and Print Tree
    '''
    print("#############################################################")
    combinedpath = workpath + "combined/"
    os.system("mkdir " + combinedpath)
    os.system("cp " + positivepath + "* " + combinedpath)
    os.system("cp " + negativepath + "* " + combinedpath)
    if "resultparsnp" in os.listdir(workpath):
        print("parsnp result exists, skipping parsnp!")
    else:
        print("running parsnp to get a Newik tree of positive data")
        os.system("python ./parsnp/Parsnp.py -b -c -r ! -d " + positivepath\
            + " -o " + workpath + "resultparsnp/ -p " + threadsnum)
    with open(workpath+"resultparsnp/parsnp.tree", "r") as f:
        l = f.readline()[:-1]
    tree = Tree(l)
    print("#############################################################")
    print("Tree for files, maybe helpful for identifying wrong data")
    print(tree)
    print("rerun without this parameter after change the\
positive/negative dataset!")
    exit()

def dropless90(table):
    """
    when getting tables, drop strains have alignment quality less than 95%
    """
    for i in table.columns:
        if max(list(table[i])) < 95:
            table = table.drop(i ,axis=1)
    return table


def run(argv):
    args = parseArgs(argv)
    # print(args.work_dir)
    # print(args.config_dir)
    global workpath
    global ntpath
    global threadsnum
    global positivepath
    global negativepath
    workpath = args.work_dir
    threadsnum = str(args.processes)
    if workpath[-1] != "/":
        workpath = workpath+"/"
    controlpath = args.config_dir
    try:                                                                        #Create Work Directory
        os.mkdir(workpath)
    except FileExistsError:
        print("File exists: " + workpath)
    ntpath = args.database_path
    '''
    Get Positive/Negative controls
    '''

    print("Work Path:" + workpath)
    print("Config file path: " + controlpath)
    with open(controlpath, "r") as f:
        l = f.readlines()
    [positive,positive_paths,negative,negative_paths,whole,group]=getinfo(l)
    positivepath = workpath + "positive/"
    negativepath =  workpath + "negative/"
    if detectdataexist(workpath):
        print("Found generated positive/negative controls, skipping data\
 gathering process...")
    else:
        print("#############################################################")
        print("Taxon group: " + group)
        print("Positive control group:")
        for i in positive:
            print(i)
        print()
        print("Negative control group:")
        for i in negative:
            print(i)
        download_db(whole, group)                                                       #Download database from NCBI database
        os.system("mkdir " + positivepath)
        os.system("mkdir " + negativepath)
    for i in positive_paths:
        for j in os.listdir(i):
            os.system("cp " + i + "/" + j + " " + positivepath)
    for i in negative_paths:
        for j in os.listdir(i):
            os.system("cp " + i + "/" + j + " " + negativepath)
    '''
    Generate and Print Tree
    '''
    if args.tree_help:
        printTree()
    

    print("#############################################################")
    print("Concatenating Negative Control...")
    randrefindex = np.random.randint(len(os.listdir(positivepath)))
    randreddir = os.listdir(positivepath)[randrefindex]
    gatherdatapath = workpath + "gatherdata/"
    os.system("mkdir " + gatherdatapath)
    with open(gatherdatapath+"concatenatedneg.fasta", "w") as f:
        with open(positivepath+randreddir, "r") as k:
            f.writelines(k.readlines())
        for i in os.listdir(negativepath):
            with open(negativepath+i, "r") as j:
                f.writelines(j.readlines())
    os.system("cp " + positivepath + "*" + " " + gatherdatapath)

    if "MUMs" not in os.listdir(workpath):
        print("#############################################################")
        print("running parsnp")
        print("python ./parsnp/Parsnp.py -c -r ! -d " + gatherdatapath
                    + " -o " + workpath + "MUMs/ -p " + threadsnum)
        os.system("python ./parsnp/Parsnp.py -c -r " + gatherdatapath + randreddir +" -d " + gatherdatapath
                    + " -o " + workpath + "MUMs/ -p " + threadsnum)
    with open(workpath + "MUMs/parsnp.xmfa", "r") as f:
        readlist = f.readlines()
    clusters = generate_read_list(readlist, 40)
#     clusters = {}
#     for i in range(len(nlist)):
#         cluster = nlist[i].split(" ")[2]
#         if cluster not in clusters:                                             #Exclude all LCBs
#             if ('A' in rlist[i]) or ('T' in rlist[i]) or ('C' in rlist[i])\
#                  or ('G' in rlist[i]) or ('-' in rlist[i]):
#                 continue
#             clusters.update({cluster: rlist[i]})
    os.system("mkdir " + workpath + "finalMUMs/")
    for i in clusters:
        with open(workpath + "finalMUMs/" + i + ".fasta", "w") as f:
            f.write(">"+i+"\n")
            f.write(clusters[i])
    clusterl = os.listdir(workpath + "finalMUMs/")
    if args.MUMS_only:
        print('-m: not entering into Blastn step, exiting...')
        exit()
    print("Running Blastn agianst nt database")
    os.system("mkdir " + workpath + "blastresult/")
    for i in pbar(clusterl):
        if ntpath == None:
            os.system("blastn -max_target_seqs 2000 -query -db nt "
            + workpath + "finalMUMs/" + i + " -out "
            + workpath + "blastresult/"+ i
            +".out -outfmt '6 qseqid sseqid pident evalue stitle' -num_threads "
            + threadsnum + " -remote")
        else:
            os.system("blastn -max_target_seqs 2000 -db "+ ntpath +
            " -query " + workpath + "finalMUMs/" + i + " -out "
            + workpath + "blastresult/"+ i
            +".out -outfmt '6 qseqid sseqid pident evalue stitle' -num_threads "
            + threadsnum)
    # os.system("mv *.out " + workpath + "blastresult/")
    os.system('find '+ workpath + 'blastresult/ -name "*" -type\
 f -size 0c | xargs -n 1 rm -f')                                                #Remove all damaged blast results
    blastresultfile = os.listdir(workpath + 'blastresult/')
    [percent2strain, strain] = percenttostrain(workpath + 'blastresult/',\
         blastresultfile, clusters)
    positivelist = []
    others = []
    for i in positive:
        for j in strain:
            if set(j.split(" ")) > set(i.split(" ")):
                positivelist.append(j)
            else:
                others.append(j)
    whole = positivelist+others
    print("generating whole result...")
    res = pd.DataFrame(columns=["MUM"]+whole)
    for i in percent2strain:
        res = res.append(i, ignore_index=True)
    res = res.fillna(0)
    res.set_index(["MUM"], inplace = True)
    print("generating whole complete genome result...")
    rescg = res
    for i in strain:
        keys = i.split(" ")
        if keys[-1] == "genome" and keys[-2] == "complete":
            continue
        rescg = rescg.drop(i, axis=1)
    with open(workpath+"all_strains.csv", "w") as f:
        f.write(res.to_csv())
    with open(workpath+"complete_genomes.csv", "w") as f:
        f.write(rescg.to_csv())

    print("Generating heatmap with all over\
    95% alignment scores strains' complete genomes")
    newcg = pd.read_csv(workpath+"complete_genomes.csv").drop("MUM", axis = 1)
    complete_genomes = dropless90(newcg)
    f, ax = plt.subplots(figsize = (200, 50))
    sns_plot = sns.heatmap(complete_genomes,
            cmap = sns.color_palette("Blues", 500),
            linewidths = 0.1, ax = ax)
    ax.set_title('Blast result for all MUMs')
    ax.set_xlabel('Strains')
    ax.set_ylabel('MUMs')
    sns_plot.figure.savefig(workpath+"competegenomes.png")

run(sys.argv[1:])