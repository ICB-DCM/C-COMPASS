### SYNTHETIC DATASET

import random
import numpy as np
import pandas as pd
from collections import Counter

filename_fract = 'sim_Fractionation.txt'
filename_marker = 'sim_Markerlist.txt'
filename_tp = 'sim_TotalProteome.txt'

replicates = 3                                                                  # specify number of simulated replicates per condition
conditions = 2                                                                  # specify number of simulated conditions


## SIMULATE FRACTIONATION:
compartments = 5                                                                # specify number of simulated compartments
fractions = 14                                                                  # specify number of fractions per gradient
spread = [3, 4]                                                                 # peak width in fractions (range)
intensities = [5000, 8000]                                                      # peak height as intensities (range)
intensity_variance = 50                                                         # increase for more intensity-variance across a single peaks
peak_accuracy = 200                                                             # increase this value for more accurate peaks

markers = [50, 100]                                                             # number of marker proteins per compartment (range)
unknown_single = [50, 100]                                                      # number of proteins per compartment with one unknown localization (range)
unknown_double = [10, 20]                                                       # number of proteins per compartment with a random second localization (range)
ratios_double = [(75, 25), (50, 50)]                                            # possible mixing ratios for double localizations
unknown_triple = [3, 8]                                                         # number of proteins per compartment with two random additional localizations (range)
ratios_triple = [(50, 25, 25)]                                                  # possible mixing ratios for triple localizations

missing_rep = 0.04                                                              # chance for a protein to be missing in single replicates
missing_cond = 0.02                                                             # chance for a protein to be missing in single conditions
reloc_rate = 0.1                                                                # probability for relocalization


def create_compspecs():
    comp_specs = {}
    for i in range(compartments):
        name = 'Compartment'+str(i+1)
        comp_specs[name] = {}
        
        bins = np.linspace(1, fractions, compartments + 1)
        middle = random.randint(int(bins[i]), int(bins[i + 1] - 1))
        comp_specs[name]['middle'] = middle
        
        sigma = random.uniform(spread[0]/2, spread[1]/2)
        comp_specs[name]['sigma'] = sigma
        
        height = random.uniform(intensities[0], intensities[1])
        comp_specs[name]['height'] = height
        
        comp_specs[name]['number_marker'] = int(random.uniform(markers[0], markers[1]))
        comp_specs[name]['number_single'] = int(random.uniform(unknown_single[0], unknown_single[1]))
        comp_specs[name]['number_double'] = int(random.uniform(unknown_double[0], unknown_double[1]))
        comp_specs[name]['number_triple'] = int(random.uniform(unknown_triple[0], unknown_triple[1]))
    return comp_specs

def reflect_distribution(distribution, lower, upper):
    reflected = np.copy(distribution)
    reflected[distribution < lower] = lower
    reflected[distribution > upper] = upper
    return reflected

def create_profile(middle, sigma, height):
    num_samples = int(np.random.normal(peak_accuracy, peak_accuracy/10, 1))
    random_values = np.random.normal(middle, sigma, num_samples)
    
    reflected_values = reflect_distribution(random_values, 1, fractions)
    discrete_values = np.clip(np.round(reflected_values), 1, fractions).astype(int)
    value_counts = Counter(discrete_values)
    value_counts_x = np.array([value_counts[i] if i in value_counts else 0 for i in range(1, fractions + 1)])
    
    factor = height/max(value_counts_x)
    value_counts_scaled = np.copy(value_counts_x)
    for j in range(fractions):
        factor_rand = float(np.random.normal(factor, height/intensity_variance, 1))
        new_value = value_counts_x[j]*factor_rand
        if new_value >=0:
            value_counts_scaled[j] = new_value
        else:
            while new_value <0:
                factor_rand = float(np.random.normal(factor, height/45, 1))
                new_value = value_counts_x[j]*factor_rand
            value_counts_scaled[j] = new_value
    return value_counts_scaled


comp_specs = {}
for cond in range(conditions):
    comp_specs[cond] = create_compspecs()
comp_list = []
for comp in range(compartments):
    comp_list.append('Compartment'+str(comp+1))

data_columns = ['ProteinName', 'GeneName']
for cond in range(conditions):
    for rep in range(replicates):
        for fract in range(fractions):
            data_columns.append('Con'+str(cond+1)+'_Rep'+str(rep+1)+'_Fr'+str(fract+1).zfill(2))
    for comp in comp_specs[cond]:
        data_columns.append('Amount_'+comp)
data_columns.append('Marker')

# CREATE MARKER PROFILES:
count = 0
all_profiles = []
for comp in comp_list:
    for m in range(comp_specs[cond][comp]['number_marker']):
        count = count +1
        prot_name = 'Prot'+str(count)
        gene_name = 'Gene'+str(count)
        profile_conc = [prot_name, gene_name]
        for cond in range(conditions):
            empty_condition = random.random() < missing_cond
            if empty_condition:
                profile_conc.extend((fractions*replicates + len(comp_list)) * [np.nan])
            else:
                
                specs = comp_specs[cond]
                location = [0 if key != comp else 1 for key in specs]
                middle_1 = specs[comp]['middle']
                sigma_1 = specs[comp]['sigma']
                height_1 = specs[comp]['height']
                for rep in range(replicates):
                    empty_replicate = random.random() < missing_rep
                    if empty_replicate:
                        profile_conc.extend(fractions * [np.nan])
                    else:
                        profile_1 = create_profile(middle_1, sigma_1, height_1)
                        profile_conc.extend(profile_1.astype(float))
                profile_conc.extend(location)
        profile_conc.append(comp)
        all_profiles.append(profile_conc)

# CREATE UNKNOWN SINGLE LOCALIZATIONS:
for comp in comp_list:
    comp_others = [c for c in comp_list if c != comp]
    for s in range(comp_specs[cond][comp]['number_single']):
        count = count +1
        prot_name = 'Prot'+str(count)
        gene_name = 'Gene'+str(count)
        profile_conc = [prot_name, gene_name]
        for cond in range(conditions):
            empty_condition = random.random() < missing_cond
            if empty_condition:
                profile_conc.extend((fractions*replicates + len(comp_list)) * [np.nan])
            else:
            
                specs = comp_specs[cond]
                reloc = random.random() < reloc_rate
                if reloc:
                    c = random.choice(comp_others)
                else:
                    c = comp
                location = [0 if key != c else 1 for key in specs]
                middle_1 = specs[c]['middle']
                sigma_1 = specs[c]['sigma']
                height_1 = specs[c]['height']
                for rep in range(replicates):
                    empty_replicate = random.random() < missing_rep
                    if empty_replicate:
                        profile_conc.extend(fractions * [np.nan])
                    else:
                        profile_1 = create_profile(middle_1, sigma_1, height_1)
                        profile_conc.extend(profile_1.astype(float))
                profile_conc.extend(location)
        profile_conc.append(np.nan)
        all_profiles.append(profile_conc)

# CREATE UNKNOWN DOUBLE LOCALIZATIONS:
for comp in comp_list:
    comp_others = [c for c in comp_list if c != comp]
    for d in range(comp_specs[cond][comp]['number_double']):
        count = count+1
        prot_name = 'Prot'+str(count)
        gene_name = 'Prot'+str(count)
        profile_conc = [prot_name, gene_name]
        for cond in range(conditions):
            empty_condition = random.random() < missing_cond
            if empty_condition:
                profile_conc.extend((fractions*replicates + len(comp_list)) * [np.nan])
            else:
            
                specs = comp_specs[cond]
                reloc = random.random() < reloc_rate
                if reloc:
                    c_1 = random.choice(comp_others)
                    c_others = [co for co in comp_list if co != c_1]
                    c_2 = random.choice(c_others)
                else:
                    c_1 = comp
                    c_2 = random.choice(comp_others)
                    
                middle_1 = specs[c_1]['middle']
                sigma_1 = specs[c_1]['sigma']
                height_1 = specs[c_1]['height']
                middle_2 = specs[c_2]['middle']
                sigma_2 = specs[c_2]['sigma']
                height_2 = specs[c_2]['height']
                
                ratio = random.choice(ratios_double)
                location = [0 if key != c_1 and key != c_2
                            else ratio[0]/100 if key == c_1
                            else ratio[1]/100
                            for key in specs]
                
                for rep in range(replicates):
                    empty_replicate = random.random() < missing_rep
                    if empty_replicate:
                        profile_conc.extend(fractions * [np.nan])
                    else:
                        profile_1 = create_profile(middle_1, sigma_1, height_1).astype(float)
                        profile_2 = create_profile(middle_2, sigma_2, height_2).astype(float)
                        profile_combined = (ratio[0]/100 * profile_1) + (ratio[1]/100 * profile_2)
                        profile_conc.extend(profile_combined.astype(float))
                profile_conc.extend(location)
        profile_conc.append(np.nan)
        all_profiles.append(profile_conc)

# CREATE UNKNOWN TRIPLE LOCALIZATIONS:
for comp in comp_list:
    comp_others = [c for c in comp_list if c != comp]
    for d in range(comp_specs[cond][comp]['number_triple']):
        count = count+1
        prot_name = 'Prot'+str(count)
        gene_name = 'Gene'+str(count)
        profile_conc = [prot_name, gene_name]
        for cond in range(conditions):
            empty_condition = random.random() < missing_cond
            if empty_condition:
                profile_conc.extend((fractions*replicates + len(comp_list)) * [np.nan])
            else:
                
                specs = comp_specs[cond]
                reloc = random.random() < reloc_rate
                if reloc:
                    c_1 = random.choice(comp_others)
                    c_others = [co for co in comp_list if co != c_1]
                    c_2 = random.choice(c_others)
                    c_others = [co for co in comp_list if co != c_1 and co != c_2]
                    c_3 = random.choice(c_others)
                else:
                    c_1 = comp
                    c_2 = random.choice(comp_others)
                    c_others = [co for co in comp_others if co != c_2]
                    c_3 = random.choice(c_others)
                    
                middle_1 = specs[c_1]['middle']
                sigma_1 = specs[c_1]['sigma']
                height_1 = specs[c_1]['height']
                middle_2 = specs[c_2]['middle']
                sigma_2 = specs[c_2]['sigma']
                height_2 = specs[c_2]['height']
                middle_3 = specs[c_3]['middle']
                sigma_3 = specs[c_3]['sigma']
                height_3 = specs[c_3]['height']
                
                ratio = random.choice(ratios_triple)
                location = [0 if key != c_1 and key != c_2 and key != c_3
                            else ratio[0]/100 if key == c_1
                            else ratio[1]/100 if key == c_2
                            else ratio[2]/100
                            for key in specs]
                
                for rep in range(replicates):
                    empty_replicate = random.random() < missing_rep
                    if empty_replicate:
                        profile_conc.extend(fractions * [np.nan])
                    else:
                        profile_1 = create_profile(middle_1, sigma_1, height_1).astype(float)
                        profile_2 = create_profile(middle_2, sigma_2, height_2).astype(float)
                        profile_3 = create_profile(middle_3, sigma_3, height_3).astype(float)
                        profile_combined = (ratio[0]/100 * profile_1) + (ratio[1]/100 * profile_2) + (ratio[2]/100 * profile_3)
                        profile_conc.extend(profile_combined.astype(float))
                profile_conc.extend(location)
        profile_conc.append(np.nan)
        all_profiles.append(profile_conc)

dataset = pd.DataFrame(all_profiles, columns = data_columns)
dataset_shuffled = dataset.sample(frac = 1).reset_index(drop = True)
                
markerset = dataset[['GeneName', 'Marker']].dropna(subset = ['Marker'])


## SIMULATE TOTAL PROTEOME:
tp_intensities = [7000, 10000]
variance = 500
regulated = 20                                                                  # precentage of proteins with changing expression level
changes = [0.1, 10]                                                             # possible fold changes for regulated proteins (range)

proteins = list(dataset['ProteinName'])
      
tp_columns = ['ProteinName']
for cond in range(conditions):
    for rep in range(replicates):
        colname = 'Con'+str(cond+1)+'_Rep'+str(rep+1)
        tp_columns.append(colname)
    tp_columns.append('RelativeRegulation_Con'+str(cond+1))


all_values = []
for prot in proteins:
    height = random.uniform(tp_intensities[0], tp_intensities[1])
    values = [prot]
    changing = random.random() < regulated/100
    
    if changing:
        fc = random.uniform(changes[0], changes[1])
    else:
        fc = np.random.normal(1, 0.2, 1)[0]

    for cond in range(conditions):
        height_cond = height * fc
        variance_cond = variance * fc
        for rep in range(replicates):
            values.append(np.random.normal(height_cond, variance_cond, 1)[0])
        values.append(fc)
    all_values.append(values)

dataset_tp = pd.DataFrame(all_values, columns = tp_columns)

dataset_shuffled.to_csv(filename_fract, sep = '\t', index = False)
markerset.to_csv(filename_marker, sep = '\t', index = False)
dataset_tp.to_csv(filename_tp, sep = '\t', index = False)





















