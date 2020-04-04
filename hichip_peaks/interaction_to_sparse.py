#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Function to convert files from HiC-Pro results to a sparse matrix at a restriction site resolution
# restriction site is basically offset by 1. Fragment 0 has restriction site 0 and 1. last fragment will have only 5' but not 3'because that would create an extra index.
# now this will get assigned wrongly the the first and last fragment of the chromosome
# to get the real bp location of the restriction site you can use frag_prop[i][1]
# this way number of fragments is the same as number of sites and it doesn't complicate too much stuff. also logically last site (or first) of a chromosome won't give any meaningful results.


#########################################
import scipy
import scipy.sparse
import numpy
import os
import re
import multiprocessing
import subprocess
import uuid
import logging
import pickle
try:
    #this works only when installed
    from hichip_peaks import helpers
except:
    import helpers 


def HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc,prefix,keeptemp=False,tempcode=str(uuid.uuid4())[0:5]):
    """Wrapper function to call individual funcitons
    CSR_mat: (scipy.sparse.matrix) sparse matrix format of the interaction matrix at a restriction site. All reads are assigned to restriction site based on directionality. the size is the same as the number of restriction FRAGMENTS.
    so it misses the last restriction site
    frag_index: (list) dictionary with all the frag names to frag index. fragment n contains site n and n+1. site n has fragment n-1 and fragment n.
    frag_prop: (list:tuple) the list of fragments with tuple with (chr, start, end , length) of each frag
    frag_amount: (dict) number of fragments in each chromosome as a dictionary
    valid_chroms: (list) list of valid chromosomes
    chroms_offsets: (dict) offsets for each chromosome. the recorded value is the start of the new chromosome INCLUSIVE. n contains referes to site that has last/first site of chromosomes and site 1 of chromosome
    if i call site n that is the site that has both the last of the last chromosome and the start of the new one. it is not really a site because it's just the end and start bp of chromosome. note that whatever you do these sites do not mean anything for further analysis
    """
    # check inputs and make directories
    if not os.path.isdir(temporary_loc):
        os.makedirs(temporary_loc)
    if not os.path.isdir(folder):
        raise Exception("couldn't find folder containing HiC-Pro results")
    if not os.path.isfile(resfrag):
        raise Exception("annotation files couldn't be opened")
    if sizes!=None and not os.path.isfile(sizes) :
        raise Exception("annotation files couldn't be opened")
        
    logging.info("Loading experiment information and read pairs")
    logging.info("#######################################")
    logging.info("Start reading experiment information (restriction fragments and chromosomes)")
    
    # read annotation files and create properties of fragments and matrix
    frag_index,frag_prop,frag_amount,valid_chroms, chroms_offsets = Read_resfrag(resfrag,sizes)

    logging.info("#######################################")
    logging.info("Preparing HiC-Pro output for import")

    # prepare the files from hic-pro and pass filenames back
    file_valid_pairs, file_self_circle, file_dangling, file_religation = Prepare_files(folder,temporary_loc,tempcode,prefix)

    logging.info("#######################################")
    logging.info("Converting HiC-Pro to sparse matrix rappresentation of valid pairs at restriction site resolution")

    # make the sparse matrix. sends a file at a time and adds stuff to the lists that are in coordinate format
    coo_data = []
    coo_row = []
    coo_col = []
    for current_file in [file_valid_pairs, file_self_circle, file_religation]:
        coo_data, coo_row, coo_col = Update_coo_lists_site(current_file,coo_data, coo_row, coo_col,valid_chroms,frag_index)
    coo_data, coo_row, coo_col = Update_coo_lists_site(file_dangling,coo_data, coo_row, coo_col,valid_chroms,frag_index,dangling = True)
    
    # major step that converts the lists into a sparse matrix repressatation
    CSR_mat = scipy.sparse.csr_matrix((coo_data, (coo_row, coo_col)), shape=(len(frag_index)+1, len(frag_index)+1), dtype = numpy.float32)[:len(frag_index),:len(frag_index)]
    
    if keeptemp == False:
        os.remove(file_self_circle)
        os.remove(file_religation)
        os.remove(file_dangling)

    logging.info("#######################################")
    logging.info("Sparse matrix of experiment generated")
    logging.info("Number of read pairs parsed: {}".format(CSR_mat.sum()/2))

    if keeptemp == True:
        logging.info("Saving intermediate files for further use")
        scipy.sparse.save_npz(os.path.join(temporary_loc, prefix + "CSR_matrix.npz"), CSR_mat)
        with open(os.path.join(temporary_loc, prefix + "variables.pickle"),"wb") as picklefile:
            pickle.dump([frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets],picklefile)
    
    return CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets



def Prepare_files(folder,temporary_loc,tempcode,prefix):
    """Find out files in folder remove duplicates and return variables locating files.
    If multiple files are present merge files"""
    
    # check if files are present
    files = os.listdir(folder)
    if len(files) < 1:
        raise Exception("couldn't find files in specified folder")
    for i in range(len(files)):
        files[i] = os.path.join(folder, files[i])

    regex = re.compile(".*allValidPairs")
    try:
        file_valid_pairs = list(filter(regex.match, files))
        if len(file_valid_pairs) > 1:
            raise Exception("Error, There is more than 1 file with .*allValidPairs in the HiC-pro folder, please remove any file that doesn't contain pairs from the folder (eg. stats)")
        file_valid_pairs = file_valid_pairs[0]
    except IndexError:
        raise Exception("couldn't find allValidPairs file")
    # because there can be multiple SC, DE and RE files list all of them.
    regex = re.compile(".*SCPairs")
    list_self_circle = list(filter(regex.match, files))
    file_self_circle = os.path.join(temporary_loc, prefix + tempcode + "SCPairs")

    regex = re.compile(".*DEPairs")
    list_dangling = list(filter(regex.match, files))
    file_dangling = os.path.join(temporary_loc, prefix + tempcode + "DEPairs")

    regex = re.compile(".*REPairs")
    list_religation = list(filter(regex.match, files))
    file_religation = os.path.join(temporary_loc, prefix + tempcode + "REPairs")

    if (len(file_self_circle) < 1) or (len(file_dangling) < 1) or (len(file_religation) < 1):
        raise Exception("couldn't find all files in specified folder")
    
    # runs commands in parallel for sorting and unique the files. also the output is joined to one file
    shcommands = []
    for files , output in zip((list_self_circle, list_dangling, list_religation),(file_self_circle, file_dangling, file_religation)):
        command ="sort -u -k 2,2 -k 3,3 -k 5,5 -k 6,6 " + " ".join(files) + " > " + output
        p = subprocess.Popen(command, shell=True)
        # runs commands in parallel
        shcommands.append(p)
    for p in shcommands:
        #waits for them to finish
        p.wait()
    
    return file_valid_pairs, file_self_circle, file_dangling, file_religation


def Read_resfrag(resfrag,sizes):
    """Read the restriction fragment file information and prepare list with resfrag properties"""
    # assumes resfrag file is sorted as the valid chroms file. valid chroms file only has proper chromosomes, no contigs
    # results include
    # frag_amount the number of fragments in the valid chromosomes for each chromosome
    # frag_index dictionary with all the frag names to frag index
    # frag_prop the list of fragments with tuple with chr, start and end and length of each frag
    # list of valid chromosomes
    # offsets used for the sparse matrix implementation, not used anymore there, bust might still be useful
    
    # note default doesn't have sex chromosomes
    if sizes==None:
        valid_chroms=['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                    'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    else:
        with open(sizes,"r") as file_sizes:
            valid_chroms=[]
            for line in file_sizes:
                valid_chroms.append(line.split("\t")[0])
    
    with helpers.open_by_suffix(resfrag) as file_resfrag:
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
    frag_index=dict()
    for i, item in enumerate(frag_name):
        frag_index[item] = i
    # quickly generate the offsets for each chromosome. the recorded value is the start of the new chromosome INCLUSIVE. 
    offsets={}
    for i in range(len(valid_chroms)):
        offsets[valid_chroms[i]] = sum([frag_amount[k] for k in valid_chroms[:i]])
    return frag_index,frag_prop,frag_amount,valid_chroms, offsets


def Update_coo_lists_site(current_file,data, row, col,valid_chroms,frag_index,dangling = False):
    """Takes file and assigns reads to restriction sites based on direction. 
    Returns the list that is then used to create the sparse matrix"""
    # always appends to lists
    with open(current_file, "r") as pairs:
        for line in pairs:
            info = line.split()
            frag_1 = info[8]
            frag_2 = info[9]
            chr_1 = info[1]
            chr_2 = info[4]
            dir_1 = info[3]
            dir_2 = info[6]
            if (chr_1 in valid_chroms) and (chr_2 in valid_chroms):
                index_frag_1 = frag_index[frag_1]
                index_frag_2 = frag_index[frag_2]
                if dangling:
                    # if dangling end reverse the direction (assign to back instead of front, shouldn't really matter as the two reads compensate themselves)
                    if dir_1 == "-":
                        index_frag_1 += 1
                    if dir_2 == "-":
                        index_frag_2 += 1
                else:
                    # if read is positive assign to next site, otherwise to current site
                    if dir_1 == "+":
                        index_frag_1 += 1
                    if dir_2 == "+":
                        index_frag_2 += 1
                # it's appending twice for the two sides of the matrix. future optimization might include using a upper triangular format to save half the memory
                data.append(1)
                data.append(1)
                row.append(index_frag_1)
                col.append(index_frag_2)
                row.append(index_frag_2)
                col.append(index_frag_1)

    return data, row, col





if __name__=="__main__":
    """test wrapper function, calls the wrapper function for testing. Won't run when imported"""
    import pickle
    # things needed
    
    # hicpro directory containing all the valid ints files
    # restriction fragment definitions
    # sizes of chromosomes just to get the right chromosomes

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[
        logging.StreamHandler()
    ]
    )
    folder = os.path.abspath("./testdata/NaiveT_27ac_B1_T1")
    resfrag = os.path.abspath("./testdata/MboI_resfrag_hg38.bed")
    sizes = None
    temporary_loc = os.path.abspath("./testdata")

    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc,"testdata")

    scipy.sparse.save_npz('./testdata/sparse_matrix.npz', CSR_mat)
    with open("./testdata/variables.pi","wb") as picklefile:
        pickle.dump([frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets],picklefile)