import os
import pandas as pd
import numpy as np
from skimage import color
from skimage import io
import matplotlib.pyplot as plt
import json
from collections import namedtuple
import periodictable
import time
import scipy.stats as ss
from numpy import linalg as la
from PIL import Image

#from string import maketrans
from pprint import pprint

threshold = 1

##################################################
# Uses the list() constructor to create a list 4 different types of .tif lists
# based on whether or not they contain certain characters
# -----------------------------------------------------------------------------
# filter(function, iterable) - filters the given iterable with the help of a function 
# os.listdir(path)  		 - Return a list of the entries in the directory given by path.
##################################################

data_dir = "./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/dataset_1_opaques/"
standard_dirs = list(filter(lambda d: d[:14] == "standards_32bt", [d for d in os.listdir('./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/dataset_1_opaques/')]))
mask_dirs = list(filter(lambda d: d[:4] != "obj1" and d[:4] != "obj2" and d[-8:]=="mask.tif",[d for d in os.listdir('./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/dataset_1_opaques/')]))
obj1_dirs = list(filter(lambda d: d[:9] == "obj1_32bt", [d for d in os.listdir('./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/dataset_1_opaques/')]))
obj2_dirs = list(filter(lambda d: d[:9] == "obj2_32bt", [d for d in os.listdir('./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/dataset_1_opaques/')]))

##################################################
# Creates Python Dictionaries based on the list retrieved
# Each dictionary stores an array obtained from the grayscale. The array represents an each pixel point.
# The intensity of each point is assumed to be not absolute.
# -----------------------------------------------------------------------------
# io.imread (fname) - Image reading and writing via imread. The different color bands/channels 
#			  		  are stored in the third dimension, such that a gray-image is MxN
##################################################

standards = {}
for s in standard_dirs:
    tmp = s.split(".")[0].split("_")[-1]
    standards[tmp] = io.imread(data_dir+s)
masks = {}
for m in mask_dirs:
    masks[m.split("std")[0].split("_")[0].split("-")[0]] = io.imread(data_dir+m)
obj1 = {}
for o in obj1_dirs:
    tmp = o.split(".")[0].split("_")[-1]
    obj1[tmp] = io.imread(data_dir+o)
obj2 = {}
for o in obj2_dirs:
    tmp = o.split(".")[0].split("_")[-1]
    obj2[tmp] = io.imread(data_dir+o)
	
##################################################
# Show - takes an array and an int
# -----------------------------------------------------------------------------
# matplotlib.pyplot.imshow(X, cmap=None) - Display an image, i.e. data on a 2D regular raster.
##################################################

def show(array, index):
	plt.imshow(array[index], cmap="gray")

masks.keys() #outputs 'dict_keys(['NiS', 'Ni', 'SCOlv', 'FeS', 'Fe', 'CaTiO3', 'Fe3O4', 'TiO2'])'

standards.keys() #outputs 'dict_keys(['Cr', 'Ni', 'Si', 'Ca', 'S', 'Al', 'P', 'Fe', 'Mg', 'Ti'])'

##################################################
# parsed_weight is a JSON file with parsed Wt% from "standard-element-concentrations.xlsx"
# The standard compounds have a nearly fixed ratio of elements and therfore are compared against the masked
# (unknown) distribution of a known element.
# -----------------------------------------------------------------------------
# periodictable is a library for, obviously, an easy mainpulation of the periodic table
# items() - Method that returns a view object that displays a list of dictionary's (key, value) tuple pairs.
##################################################

f = open('./parsed_weight')
parsed_weight = json.load(f)
weights = {}	
for mask in masks.keys():
    weights[mask] = {}
    if mask == "SCOlv":
        for s in parsed_weight["10"].keys():
            if s in standards.keys():		
                weights[mask][s] = parsed_weight["10"][s]/100
            else:
                weights[mask][standard] = 0
    else:
        for standard in standards.keys():
            # d = dict((str(e), w) for e,w in periodictable.formula(mask).mass_fraction.items())
            # if standard in mask:
                # weights[mask][standard] = d[standard]
            # else:
            weights[mask][standard] = 0
				
# Example of one of the items in the weights dictionary
#{'NiS': {'Cr': 0,
#  'Ni': 0.6466993688738453,
#  'Si': 0,
#  'Ca': 0,
#  'S': 0.35330063112615473,
#  'Al': 0,
#  'P': 0,
#  'Fe': 0,
#  'Mg': 0,
#  'Ti': 0},

##################################################
# mask_maps creates a dict of dict. 
# -----------------------------------------------------------------------------
# numpy.zeros_like - Return an array of zeros with the same shape and type as a given array.
##################################################

def mask(mask, standard):
    t = np.zeros_like(standard)
    for i in range(256):
        for j in range(512):
            if (mask[i][j]):
                t[i][j] = standard[i][j]
    return t

mask_maps = {}
for k in masks.keys():
    mask_maps[k] = {}
    for s in standards.keys():
        mask_maps[k][s] = mask(masks[k], standards[s])

# mask_map example
#{'NiS': {'Cr': array([[0, 0, 0, ..., 0, 0, 0],
#         [0, 0, 0, ..., 0, 0, 0],
#         [0, 0, 0, ..., 0, 0, 0],
#         ...,
#         [0, 0, 0, ..., 0, 0, 0],
#         [0, 0, 0, ..., 0, 0, 0],
#         [0, 0, 0, ..., 0, 0, 0]], dtype=uint32), 'Ni': array([[0, 0, 0, ..., 0, 0, 0]...
         

##################################################
# Function average_std - takes a mask and standard dictionaries. Takes the mean 
# and std of the pixels' intensity values
# -----------------------------------------------------------------------------
# numpy.ndarray.flatten - Return a copy of the array collapsed into one dimension.
##################################################		
	
def average_std(mask, standard):
    t = 0
    d = mask_maps[mask][standard].flatten()
    ret_avg = d[d != 0].mean()
    red_std = d[d != 0].std()
    return ret_avg, red_std
	
##################################################
# Finds s linear model of the intensity of the pixel and the weight of the element
# -----------------------------------------------------------------------------
# append() - method that adds an item to the end of the list.
# dd       - 1d dictionary of all the intensity values of each map 
# numpy.polyfit - Least squares polynomial fit. Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg] of degree deg to points (x, y).
# numpy.linspace - Return evenly spaced numbers over a specified interval.
##################################################		

lin_models = {}
i = 0.00000001
for s in standards.keys():
    x = []
    y = []
    for m in masks.keys():
        dd = mask_maps[m][s].flatten()
        dd = dd[dd != 0]
        for d in dd: 
            if s in weights[m].keys():
                if (weights[m][s] != 0):
                    x.append(weights[m][s])
                else:
                    x.append(i)
                    i = i + 0.01
                    if (i > 0.03):
                        i = 0.01
            y.append(d)
    if (len(x)):
        lin_models[s] = np.poly1d(np.polyfit(y, x, 1))
		
##################################################
# Function read_known_comps() - goes through lines of known compounds, parses them into elements. Each element
# is then again divided to its number of atoms and the element itself. known_comps is a dict of dict that 
# contains each of the compounds' compositions. ie) Fe3O4 : { Fe : 3 , O : 4} 
# ---------------------------------------------------
# endswith() - method that returns True if a string ends with the specified suffix. If not, it returns False.
##################################################	

def read_known_comps():
    known_comps = {}
    file = open("./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/known_compounds.txt")
    for line in file:
        element_dict = {}
        cpd_name = ''
        elements = line.split(',')
        for element in elements:
            dict_entries = element.split('_')
            if (len(dict_entries) == 1):
                if dict_entries[0].endswith('\n'):
                    dict_entries[0] = dict_entries[0][:-1]
                element_dict[dict_entries[0]] = 1
            else:
                if dict_entries[1].endswith('\n'):
                    dict_entries[1] = dict_entries[1][:-1]
                cpd_name += dict_entries[1]
                element_dict[dict_entries[0]] = int(dict_entries[1])
            cpd_name += dict_entries[0]
        known_comps[cpd_name] = element_dict

    return known_comps

##################################################
# read_pt function reads the periodic table
# --------------------------------------------------
####################################################	

def read_pt():
    element_table = {}
    for el in elements:
        element_table[el.symbol] = el.mass
    return element_table

##################################################
# Function weight_calculation - simple stoichiometry to calculate the molecular weight percentage of an element 
# in a compound.
# ---------------------------------------------------
# known compositions : dict of name of compound with a dictionary of elements with num of elements as the values
##################################################	

def weight_calculation(known_compositions, list_of_elements):
    element_table = read_pt()
    calculated_weights = {}
    for compound_key, compound_dict in known_compositions.items():
        cpd_weight = 0
        present_elements = {}
        for element_key, element_num in compound_dict.items():
            cpd_weight += element_num * float(element_table[element_key])
        for element_key, element_num in compound_dict.items():
            if element_key in list_of_elements:
                present_elements[element_key] = element_num * float(element_table[element_key]) / cpd_weight
        calculated_weights[compound_key] = present_elements
    return calculated_weights

##################################################
# Function vectorization_of_known_comps - changes each compound into vectors where each component is the
# %Molecular weight. Since there are 10 elements to be tested, naturally, a 10d vector will be created
# ---------------------------------------------------
# numpy.zeros - Return a new array of given shape and type, filled with zeros.
# enumerate() - method that adds counter to an iterable and returns it (the enumerate object).
###################################################	

def vectorization_of_known_comps(list_of_elements):
    known_comps = read_known_comps()
    known_comps = weight_calculation(known_comps, list_of_elements)
    vec_known_comps = {}
    for compound, composition in known_comps.items():
        vec_compound = np.zeros(len(list_of_elements))
        for ele_index, ele_name in enumerate(list_of_elements):
            if ele_name in composition:
                vec_compound[ele_index] = composition[ele_name]
            else:
                vec_compound[ele_index] = 0
        vec_known_comps[compound] = vec_compound
    return vec_known_comps
	
##################################################
# Function identification - 
# ---------------------------------------------------
# 
# indices are shifted by 1
##################################################	

def identification(percent_comps, list_of_elements):
    known_comps = vectorization_of_known_comps(list_of_elements)

    stddev_of_kc = []  # standard deviation of known compositions

    results_dict = {}

    # This part is purely for visualization
    matrix_size = percent_comps.shape
    discrete_comp = np.chararray((matrix_size[0], matrix_size[1]))
    discrete_comp[:] = ''
	
    for i in range(matrix_size[0]):
        for j in range(matrix_size[1]):
            most_prob = ['Unknown', 10]  # First index is element, Second is percent probability
            for comp_name, composition_vector in known_comps.items():
                norm_of_diff = la.norm(percent_comps[i][j] - composition_vector)
                # if ss.norm(0, stddev_of_kc[most_prob[0]]).pdf(norm_of_diff) < most_prob[1]: # what it should be.
                # also could just make it quadratic would fit better but its like the same thing
                if norm_of_diff < most_prob[1] and norm_of_diff < threshold:
                    most_prob = [comp_name, norm_of_diff]  # change to prob later

            if most_prob[0] in results_dict:
                results_dict[most_prob[0]] += 1
            else:
                results_dict[most_prob[0]] = 1
            # need to add something that takes the most prob ones
            discrete_comp[i][j] = most_prob[0]
            print(most_prob[0])

    for comp_name, comp_count in results_dict.items():
        results_dict[comp_name] = comp_count / (matrix_size[0]*matrix_size[1])

    return results_dict, discrete_comp

##################################################
# # pic_names.txt - contains 32bit tif files of the standard elements
# ---------------------------------------------------
# results - aggregate of the percentage of all elements shown in the image.
# big_boy - The 3rd axis here, idx, is a counter for each of the 10 elements. big_boy elimantes the counter and 
#		    instead creates a doubly nested arrays that sets the relationship between the grayscale and the element.
# shape[] - The shape attribute for numpy arrays returns the dimensions of the array. 
#			If Y has n rows and m columns, then Y.shape is (n,m). So Y.shape[0] is n.
# numpy.dstack(tup) - Stack arrays in sequence depth wise (along third axis). Takes in a sequence of arrays and returns an array.
# zip() - The zip() function take iterables (can be zero or more), makes iterator that aggregates 
#		  elements based on the iterables passed, and returns an iterator of tuples.
# [:]   - Return a shallow copy of the list.
##################################################	
	
list_of_elements = []
gs_matrix = []
f = open("./python/AMNH_Meteorite_Challenge-master/datas/meteorite_mineral_mapper/pic_names.txt")
for line in f:
    try:
        im = Image.open(line[:-1])
        words = line.split('_')
        list_of_elements.append(words[2][:-5])  # This is kind of sloppy
        gs_matrix.append(np.array(im).astype(np.float64))
    except IOError:
        print("oof")
        pass
		
idx = 0
for element, e in zip(gs_matrix, list_of_elements):
    for i in range(element.shape[0]):
        for j in range(element.shape[1]):
            gs_matrix[idx][i][j] = lin_models[e](gs_matrix[idx][i][j])
    idx += 1
big_boy = np.dstack(gs_matrix[:])
results, discrete_comp = identification(big_boy, list_of_elements)
summ = 0
for name, percent in results.items():
    print(name + ": " + str(percent))
	
#list_of_elements
#['Al', 'Ca', 'Cr', 'Fe', 'Mg', 'Ni', 'P', 'S', 'Si', 'Ti']
