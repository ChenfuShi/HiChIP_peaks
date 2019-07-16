import scipy
import scipy.sparse
import numpy
import pickle
import matplotlib.pyplot

class site_matrix():
    """class for the sparse matrix implementation of the restriction site resolution map for HiC and HiChIP"""
    # note that each chromosome loses the last site and the reads get assigned to the first site of the next chromosome!
    def __init__(self, mat, valid_chroms,frag_prop,frag_index,frag_amount,chroms_offsets):
        """basic initializer using a matrix that has already been created using interaction_to_sparse"""
        self.CSR_mat = mat
        self.valid_chroms = valid_chroms
        self.chroms_offsets = chroms_offsets
        self.frag_prop = frag_prop
        self.frag_index = frag_index
        # chroms_prop start and end of each chromosome
        self.chroms_prop = {}
        for i in chroms_offsets.keys(): 
            self.chroms_prop[i] = (chroms_offsets[i], chroms_offsets[i] + frag_amount[i])
        
        # site prop. list (chr, position, start, end). basically substitute frag_prop. 
        # basically loses half of last fragment of each chromosome. new chromosome starts at 0. but reads include past chromosome
        self.site_prop = [(frag_prop[0][0],0,0,(frag_prop[0][1] + frag_prop[0][2])/2)]
        for i in range(1,len(frag_prop)):
            start = (frag_prop[i-1][1] + frag_prop[i-1][2])/2
            end = (frag_prop[i][1] + frag_prop[i][2])/2
            if frag_prop[i][0] == frag_prop[i-1][0]:
                self.site_prop.append((frag_prop[i][0],frag_prop[i][1],int(start),int(end)))
            else:
                self.site_prop.append((frag_prop[i][0],frag_prop[i][1],0,int(end)))
        # initialize peaks but do not add them
        self.peaks = None
              
        
    def extract_diagonal(self,window):
        """extract the diagonal including the sum of the window in all directions"""
        diagonal = self.CSR_mat.diagonal()#/2 #check if this is making a problem. diagonal gets divided by two, but then the two off diagonals get used twice before and after. 
        # maybe we should consider the diagonal twice? those are religation pairs. don't think it would make a massive difference but yeah check
        if window == 0:
            return numpy.array(diagonal)
        for i in range(1,window+1):
            off_diagonal = CSR_mat.diagonal(k=i).tolist()
            diagonal = [sum(x) for x in zip(diagonal, [0]*i + off_diagonal, off_diagonal + [0]*i)]
        return numpy.array(diagonal)
        
        
    def register_peaks(self,refined_peaks):
        """input refined_peaks and create list of peaks as tuples. Return that list"""
        self.peaks=[]
        i = 0
        while i < len(refined_peaks):
            if refined_peaks[i] == 1:
                start = i
                while refined_peaks[i] == 1:
                    i=i+1
                end = i
                self.peaks.append((start,end))
            i=i+1
        return self.peaks
    

    def find_supporting_interactions(self,peak1,peak2):
        """peak1 and peak2 are two tuples of start-end"""
        return self.CSR_mat[peak1[0]:peak1[1],peak2[0]:peak2[1]].sum()


    def peaks_distance(self,peak1,peak2): 
        """returns distance of the center between two peaks"""   
        pos1 = (self.frag_prop[peak1[1]][1]+self.frag_prop[peak1[0]][1])/2
        pos2 = (self.frag_prop[peak2[1]][1]+self.frag_prop[peak2[0]][1])/2
        return int(abs(pos1-pos2))


    def get_peak_info(self, peak):
        """returns chromosome and start-end of peak"""
        start = (self.frag_prop[peak[0]-1][1] + self.frag_prop[peak[0]-1][2])/2
        end = (self.frag_prop[peak[1]][1] + self.frag_prop[peak[1]][2])/2
        return (self.frag_prop[peak[0]][0],int(start),int(end))


    def get_chromosome(self,chrA,chrB=None):
        """returns a sparse matrix for the chromosome. you can also get a view between two different chromosomes"""
        # this loses the last site of each chromosome
        if chrB == None:
            chrB = chrA
        assert chrA in self.valid_chroms
        assert chrB in self.valid_chroms
        return self.CSR_mat[self.chroms_prop[chrA][0]:self.chroms_prop[chrA][1],self.chroms_prop[chrB][0]:self.chroms_prop[chrB][1]]


    def find_site_index(self,chromosome,bp):
        """find the index of the site that includes the position"""
        index_start, index_end = self.chroms_prop[chromosome]
        for i in range(index_start,index_end,100):
            if self.site_prop[i][3] > bp:
                test_i = i
                break
        while True:
            if self.site_prop[test_i][3] > bp and self.site_prop[test_i][2] <= bp:
                break
            else:
                test_i = test_i - 1
        if test_i < index_start or test_i >= index_end:
            raise Exception("find_site_index error")
        return test_i
    

    def get_region(self,regA,regB = None):
        """returns a sparse matrix for the selected region"""
        # regions are in basepairs. if you have the region in index form there is no point. just slice the matrix
        # (chr, start, end)
        if regB == None:
            regB = regA
        startA = self.find_site_index(regA[0], regA[1])
        endA = self.find_site_index(regA[0], regA[2])
        startB = self.find_site_index(regB[0], regB[1])
        endB = self.find_site_index(regB[0], regB[2])
        return self.CSR_mat[startA:endA,startB,endB]
    

    def viewpoint_extract_indexes(self,reg,distance = 1000):
        # reg is a peak with start and end indexes
        chromosome = self.site_prop[reg[0]][0]
        startChr,endChr = self.chroms_prop[chromosome]
        startView = reg[0]-1000
        endView = reg[1]+1000
        if startView < startChr:
            startView = startChr
        if endView > endChr:
            endView = endChr
        return self.CSR_mat[reg[0]:reg[1],startView:endView].sum(axis=0), [self.site_prop[i][1] for i in range(startView,endView)]
    

    def viewpoint_extract_bp(self,reg,distance = 1000000):
        # reg is a peak with start and end indexes but distance is in bp
        chromosome = self.site_prop[reg[0]][0]
        startChr,endChr = self.chroms_prop[chromosome]
        startBp = self.site_prop[reg[0]][1] - distance
        if startBp < 0:
            startBp = 0
        endBp = self.site_prop[reg[0]][1] + distance
        if endBp > self.site_prop[endChr-1][1]:
            endBp = self.site_prop[endChr-1][1]
        startView = self.find_site_index(chromosome,startBp)
        endView = self.find_site_index(chromosome,endBp)
        return self.CSR_mat[reg[0]:reg[1],startView:endView].sum(axis=0), [self.site_prop[i][1] for i in range(startView,endView)]

