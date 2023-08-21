import specufex
import os 
import obspy
from specufex import BayesianNonparametricNMF
import h5py
import yaml
import numpy as np
from tqdm import trange
## change this to input arg
yamlPath = "/Users/theresasawi/Documents/11_Manuscripts/Methods_Paper/data/yaml/demo.yaml"




####################################################################################
####################################################################################
###
### Load yaml file settings, creat paths, set parameters
###
####################################################################################
####################################################################################



with open(yamlPath) as stream:
    config = yaml.safe_load(stream)

path_config = config["paths"]
key = path_config["key"]
print("Project key:", key)



projectPath = path_config["projectPath"]
SpecUFEx_H5_name = 'SpecUFEx_' + path_config["h5name"] #f'SpecUFEx_{key}.hdf5'
SpecUFEx_H5_path = projectPath + 'data/H5files/' + SpecUFEx_H5_name


#SpecUFEx parameters
specufex_config = config['specufexParams']
N_patterns_NMF = specufex_config['N_patterns_NMF']
nmf_batchsz = specufex_config['nmf_batchsz']
nmf_nbatch = specufex_config['nmf_nbatch']
N_states_HMM = specufex_config['N_states_HMM']
hmm_batchsz = specufex_config['hmm_batchsz']
hmm_nbatch = specufex_config['hmm_nbatch']


####################################################################################
####################################################################################
###
### Import fingerprints event ids, load and linearize FPs
###
####################################################################################
####################################################################################

ev_IDs = []

with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:
    for evID in fileLoad['fingerprints']:
        ev_IDs.append(evID)





X = linearizeFP(SpecUFEx_H5_path,ev_IDs) #event IDs needed to get fp from H5



col_names = ['fp' + str(a) for a in range(X.shape[1])]
col_names[-1]

fp_df = pd.DataFrame(X,columns=col_names)
fp_df['ev_ID'] = [str(ev) for ev in ev_IDs]
fp_df['event_ID'] = [str(ev) for ev in ev_IDs]

N = len(fp_df)

print(X.shape)


####################################################################################
####################################################################################
###
### Stack spectrograms for SpecUFEx input
###
####################################################################################
####################################################################################

X = []

with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:
    for evID in fileLoad['spectrograms']:
        specMat = fileLoad['spectrograms'].get(evID)[:]
        X.append(specMat)

    X = np.array(X)

# ================
print(np.shape(X))




####################################################################################
####################################################################################
###
### Run SpecUFEx!
###
####################################################################################
####################################################################################

nmf = specufex.BayesianNonparametricNMF(X.shape)

print('fitting NMF')
nmf.fit(X, verbose=0)

print('transforming NMF')
Vs = nmf.transform(X)


hmm = specufex.BayesianHMM(nmf.num_pat, nmf.gain)


print('fitting HMM')
hmm.fit(Vs)


print('transforming HMM')
fingerprints, As, gams = hmm.transform(Vs)


####################################################################################
####################################################################################
###
### save output to H5
###
####################################################################################
####################################################################################

print('writing all output to h5')
with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:


    ##fingerprints are top folder
    if 'fingerprints' in fileLoad.keys():
        del fileLoad["fingerprints"]
    fp_group = fileLoad.create_group('fingerprints')

    if 'SpecUFEX_output' in fileLoad.keys():
        del fileLoad["SpecUFEX_output"]
    out_group = fileLoad.create_group("SpecUFEX_output")

    # write fingerprints: ===============================
    for i, evID in enumerate(fileLoad['spectrograms']):
        fp_group.create_dataset(name= evID, data=fingerprints[i])


    # write the SpecUFEx out: ===========================
    # maybe include these, but they are not yet tested.
    ACM_group = fileLoad.create_group("SpecUFEX_output/ACM")
    STM_group = fileLoad.create_group("SpecUFEX_output/STM")

    for i, evID in enumerate(fileLoad['spectrograms']):
        ACM_group.create_dataset(name=evID,data=Vs[i]) #ACM
        STM_group.create_dataset(name=evID,data=gams[i]) #STM

    gain_group = fileLoad.create_group("SpecUFEX_output/ACM_gain")
    W_group                      = fileLoad.create_group("SpecUFEX_output/W")
    EB_group                     = fileLoad.create_group("SpecUFEX_output/EB")
    ## # # delete probably ! gain_group                   = fileLoad.create_group("SpecUFEX_output/gain")
    #RMM_group                    = fileLoad.create_group("SpecUFEX_output/RMM")

    W_group.create_dataset(name='W',data=nmf.EW)
    EB_group.create_dataset(name=evID,data=hmm.EB)
    gain_group.create_dataset(name='gain',data=nmf.gain) #same for all data
    # RMM_group.create_dataset(name=evID,data=RMM)
