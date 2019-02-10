# TODO(123)
# Process image using skimage
from PIL import Image
import numpy as np
import scipy.stats as ss
from numpy import linalg as la
from periodictable import elements
import time

threshold = 1


def add_minerals(known_comps):
    mineral_file = open("standard-element-concentrations.csv")
    element_names = mineral_file.readline()
    element_names = element_names.split(',')

    for mineral in mineral_file:
        element_dict = {}
        element_values = mineral.split(',')
        if element_values[-1].endswith('\n'):
            element_values[-1]=element_values[-1][:-1]
        for i in range(len(element_values)):
            if i == 0:
                pass
            else:
                if element_values[i]!='':
                    element_dict[element_names[i]] = float(element_values[i])/100
                else:
                    element_dict[element_names[i]]=0
        known_comps[element_values[0]] = element_dict
    return known_comps

# input file should follows this format Mg,Si,O_2 for MgSiO2 and a new line for each compound
def read_known_comps():
    known_comps = {}
    file = open("known_compounds.txt")
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

def read_pt():
    element_table = {}
    for el in elements:
        element_table[el.symbol] = el.mass
    return element_table

# known compositions : dict of name of compound with a dictionary of elements with num of elements as the values
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


def vectorization_of_known_comps(list_of_elements):
    known_comps = read_known_comps()
    known_comps = weight_calculation(known_comps, list_of_elements)
    known_comps = add_minerals(known_comps)
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


# indices are shifted by 1
def identification(percent_comps, list_of_elements):
    # lets vectorize the known_comps to the appropriate
    known_comps = vectorization_of_known_comps(list_of_elements)

    stddev_of_kc = []  # standard deviation of known compositions

    results_dict = {}

    # This part is purely for visualization
    matrix_size = percent_comps.shape
    #discrete_comp = np.chararray((matrix_size[0], matrix_size[1]))
    #discrete_comp[:] = ''


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
           #discrete_comp[i][j] = most_prob[0]

    for comp_name, comp_count in results_dict.items():
        results_dict[comp_name] = comp_count / (matrix_size[0]*matrix_size[1])

    return results_dict
