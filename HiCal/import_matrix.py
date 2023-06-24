from scipy.sparse import csr_matrix
import tables
import h5py
import hdf5plugin

import numpy as np
from scipy import signal
import glob

class Heatmap:
    def __init__(self,file_path,symmetric = True):
        self.file_name = get_file_names(file_path)[0]
        ma = make_matrix(self.file_name)
        if symmetric is False:
            self.matrix = ma
        else:
            self.matrix = make_symmetric(ma)
        chr_list = get_intervals(self.file_name)[0]['chr_list']
        unique_chr = reduce_list(chr_list)
        self.chromosomes = [x.decode() for x in unique_chr]

        self.end_idx = get_end_idx(chr_list,unique_chr)
        self.start_idx = get_start_idx(self.end_idx)
        self.idx_dict = make_idx_dict(self.chromosomes,self.start_idx,self.end_idx)

    def chromosome_idx(self,chromosome):
        for key in self.idx_dict.keys():
            if key == chromosome:
                return list(self.idx_dict[key])
            
    def region(self, chromosome):
        idxes = self.chromosome_idx(chromosome)
        return self.matrix[idxes[0]:idxes[1],idxes[0]:idxes[1]]
    
    def inter_region(self, chromosome_x, chromosome_y):
        idxes_x = self.chromosome_idx(chromosome_x)
        idxes_y = self.chromosome_idx(chromosome_y)
        return self.matrix[idxes_x[0]:idxes_x[1],idxes_y[0]:idxes_y[1]]

def get_file_names(file_path = '../data/*sampled_5kb_ICE.h5'):
    "search and store file names in a list"
    file_names = []
    for name in glob.glob(file_path):
        file_names.append(name)
    return file_names


def make_matrix(filename):
    "make the matrix with a given file name"
    with tables.open_file(filename, 'r') as f:
        parts = {}
        try:
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
        except Exception as e:
            log.info('No h5 file. Please check parameters concerning the file type!')
            e
        matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                            shape=parts['shape'])
        matrix_array = matrix.toarray()
    return matrix_array

def get_intervals(filename):
    file = h5py.File(filename)
    key_list = file['intervals'].keys()
    interval_list = {}
    for key in key_list:
        interval_list[key] = file['intervals'][key][()]
    keychr = list(key_list)
    return interval_list, keychr

def reduce_list(full_list):
    unique = []
    for stuff in full_list:
        if stuff not in unique:
            unique.append(stuff)
    return unique

def get_end_idx(chr_list,unique_chr):
    idx = np.zeros(len(unique_chr))
    for j in enumerate(unique_chr):
        for i in enumerate(chr_list):
            if i[1]==j[1]:
                idx[j[0]] = i[0]
    idx = idx.tolist()
    idx = [int(i) for i in idx]
    return idx

def get_start_idx(end_idx):
    idx_start = [0]
    for i in end_idx:
        idx_start.append(int(i+1))
    return idx_start

def make_idx_dict(chromes,start,end):
    idx_dict = {}
    for k in enumerate(chromes):
        idx_dict[k[1]] = [start[k[0]],end[k[0]]]
    return idx_dict

def make_symmetric(ma):
    rotmatrix = np.transpose(ma)
    np.fill_diagonal(rotmatrix,0)
    full_matrix = ma + rotmatrix
    return full_matrix



