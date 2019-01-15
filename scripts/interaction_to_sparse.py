#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Function to convert files from HiC-Pro results to a sparse matrix 
    # add: sanity checking all inputs
    # change temp filename so it's always unique and also do cleanup at the end
    # create directories when needed

#########################################
import scipy
import scipy.sparse
import numpy
import os
import re
import multiprocessing
import subprocess

def HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc):
    """Wrapper function to call individual funcitons"""

    frag_name,frag_prop,frag_amount,valid_chroms,chroms_offsets = Read_resfrag(resfrag,sizes)

    file_valid_pairs, file_self_circle, file_dangling, file_religation = Prepare_files(folder,temporary_loc)

    coo_data = []
    coo_row = []
    coo_col = []
    for current_file in [file_valid_pairs, file_self_circle, file_dangling, file_religation]:
        coo_data, coo_row, coo_col = Update_coo_lists(current_file,coo_data, coo_row, coo_col,chroms_offsets,valid_chroms)

    CSR_mat = scipy.sparse.csr_matrix((coo_data, (coo_row, coo_col)), shape=(len(frag_name), len(frag_name)), dtype = numpy.int32)

    return CSR_mat,frag_name,frag_prop,frag_amount,valid_chroms,chroms_offsets



def Prepare_files(folder,temporary_loc):
    """Find out files in folder remove duplicates and return variables locating files. If present merge files"""


    files = os.listdir(folder)
    if len(files) < 1:
        raise Exception("couldn't find files in specified folder")
    for i in range(len(files)):
        files[i] = os.path.join(folder, files[i])


    regex = re.compile(".*allValidPairs")
    file_valid_pairs = list(filter(regex.match, files))[0]
    

    regex = re.compile(".*SCPairs")
    list_self_circle = list(filter(regex.match, files))
    file_self_circle = os.path.join(temporary_loc, "SCPairs")

    regex = re.compile(".*DEPairs")
    list_dangling = list(filter(regex.match, files))
    file_dangling = os.path.join(temporary_loc, "DEPairs")

    regex = re.compile(".*REPairs")
    list_religation = list(filter(regex.match, files))
    file_religation = os.path.join(temporary_loc, "REPairs")

    if (len(file_self_circle) < 1) or (len(file_dangling) < 1) or (len(file_religation) < 1):
        raise Exception("couldn't find all files in specified folder")

    
    # shcommands = []
    # for files , output in zip((list_self_circle, list_dangling, list_religation),(file_self_circle, file_dangling, file_religation)):
    #     command ="sort -u -k 2,2 -k 3,3 -k 5,5 -k 6,6 " + " ".join(files) + " > " + output
    #     p = subprocess.Popen(command, shell=True)
    #     # runs commands in parallel
    #     shcommands.append(p)
    # for p in shcommands:
    #     #waits for them to finish
    #     p.wait()
    

    return file_valid_pairs, file_self_circle, file_dangling, file_religation



def Read_resfrag(resfrag,sizes):
    """Read the restriction fragment file information and prepare list"""
    # results include
    # the number of fragments in the valid chromosomes for each chromosome
    # the list of fragments with the name of the fragment
    # the list of fragments with tuple with chr, start and end and length of each frag
    with open(sizes,"r") as file_sizes:
        valid_chroms=[]
        for line in file_sizes:
            valid_chroms.append(line.split("\t")[0])

    with open(resfrag,"r") as file_resfrag:
        frag_name=[]
        frag_prop=[]
        frag_amount={}
        for i in valid_chroms:
            frag_amount[i]=0

        for line in file_resfrag:
            splitted_line=line.rstrip().split("\t")
            if splitted_line[0] in valid_chroms:
                frag_name.append(splitted_line[3])
                frag_prop.append((splitted_line[0],int(splitted_line[1]),int(splitted_line[2]),int(splitted_line[2])-int(splitted_line[1])))
                frag_amount[splitted_line[0]] += 1
    
    # quickly generate the offsets for each chromosome. the recorded value is the start of the new chromosome INCLUSIVE. Need to -1 the hic fragment when using the offset
    offsets={}
    for i in range(len(valid_chroms)):
        offsets[valid_chroms[i]] = sum([frag_amount[k] for k in valid_chroms[:i]])
    return frag_name,frag_prop,frag_amount,valid_chroms, offsets



def Update_coo_lists(current_file,data,row,col,offsets,valid_chroms):
    """Takes file and adds reads to the 3 lists containing the data"""
    with open(current_file, "r") as pairs:
        for line in pairs:
            info = line.split()
            end_1 = info[8].split("_")[-1]
            end_2 = info[9].split("_")[-1]
            chr_1 = info[1]
            chr_2 = info[4]
            if (chr_1 in valid_chroms) and (chr_2 in valid_chroms):
                frag_1 = int(end_1) + offsets[chr_1] - 1
                frag_2 = int(end_2) + offsets[chr_2] - 1
                data.append(1)
                data.append(1)
                row.append(frag_1)
                col.append(frag_2)
                row.append(frag_2)
                col.append(frag_1)


    return data, row, col










if __name__=="__main__":
    """test wrapper function, calls the wrapper function for testing. Won't run when imported"""
    import pickle
    # things needed
    
    # hicpro directory containing all the valid ints files
    # restriction fragment definitions
    # sizes of chromosomes just to get the right chromosomes

    folder = os.path.abspath("./testdata/NaiveT_27ac_B1_T1")
    resfrag = os.path.abspath("./testdata/MboI_resfrag_hg38.bed")
    sizes = os.path.abspath("./annotations/hg38.txt")
    temporary_loc = os.path.abspath("./testdata")

    CSR_mat,frag_name,frag_prop,frag_amount,valid_chroms,chroms_offsets = HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc)

    scipy.sparse.save_npz('./testdata/sparse_matrix.npz', CSR_mat)
    with open("./testdata/variables.pi","wb") as picklefile:
        pickle.dump([frag_name,frag_prop,frag_amount,valid_chroms,chroms_offsets],picklefile)