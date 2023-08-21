#1_wf_tosgramH5.py

import pandas as pd
import numpy as np
import os
import glob
import obspy
import specufex
import h5py
import yaml
from tqdm import trange

from f1_spectrogram_functions import wf_to_H5, gen_sgram_QC_noAlias



####################################################################################
####################################################################################
###
### Load yaml file settings, creat paths, set parameters
###
####################################################################################
####################################################################################

# Check if a command line argument is provided
if len(sys.argv) < 2:
    print("Please provide .yaml path")
else:
    # The first command line argument is at index 1 (index 0 is the script name)
    yamlPath = sys.argv[1]
    print("Entered value:", entered_value)



with open(yamlPath) as stream:
    config = yaml.safe_load(stream)

path_config = config["paths"]
key = path_config["key"]
print("Project key:", key)


# build path strings
dataH5_name = f'data_{key}.h5'
projectPath = path_config["projectPath"]
path_waveform = path_config["waveformPath"]
SpecUFEx_H5_name = 'SpecUFEx_' + path_config["h5name"] #f'SpecUFEx_{key}.hdf5'
SpecUFEx_H5_path = projectPath + 'data/H5files/' + SpecUFEx_H5_name


if not os.path.isdir(projectPath + 'data/H5files/'):
    os.mkdir(projectPath + 'data/H5files/')


print("waveform folder path in:", path_waveform)
print(len(glob.glob(path_waveform + "*")), "waveforms in folder")

if not os.path.isdir(projectPath + 'data/H5files/'):
    os.mkdir(projectPath + 'data/H5files/')

dataH5_path = projectPath + 'data/H5files/' + dataH5_name


##spectrogram parameters, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.spectrogram.html


#Spectrogram parameters
sgram_config = config['sgramParams']
fmin = sgram_config['fmin']
fmax = sgram_config['fmax']
winLen_Sec = sgram_config['winLen_Sec']
fracOverlap = sgram_config['fracOverlap']
nfft = sgram_config['nfft']




####################################################################################
####################################################################################
###
### Create catalog with event ids made from filenames, load demo wave
###
####################################################################################
####################################################################################


wf_list = glob.glob(path_waveform + '/*')
wf_list.sort()

filenames = wf_list
ev_ID = [path.split('/')[-1] for path in wf_list]


print('Example evID: ', ev_ID[0])

cat_paths = pd.DataFrame({"ev_ID":ev_ID,
                          "filename":wf_list})


wf_test = obspy.read(cat_paths.filename.iloc[0])

lenData = len(wf_test[0].data)
fs = wf_test[0].stats.sampling_rate

nperseg = int(sgram_config["winLen_Sec"]*fs) #datapoints per window segment
noverlap = int(nperseg*sgram_config["fracOverlap"])  #fraction of window overlapped

#padding must be longer than n per window segment
if nfft < nperseg:
    nfft = nperseg*2
    print("nfft too short; changing to ", nfft)

mode='magnitude'
scaling='spectrum'


# set args for generator
args = {'fs': fs,
        'lenData': lenData,
        'nperseg': nperseg,
        'noverlap': noverlap,
        'nfft': nfft,
        'mode': mode,
        'scaling': scaling,
        'fmin': fmin,
        'fmax': fmax
       }


print(args)



####################################################################################
#############################dd#######################################################
###
### Save waveforms to H5
###
####################################################################################
####################################################################################


evID_keep, wf_example = wf_to_H5(projectPath,dataH5_path,cat_paths,lenData, verbose=1)

###################################
## Save processing information to data H5
###################################
with h5py.File(dataH5_path,'a') as h5file:
    processing_group = h5file.create_group("processing_info")
    processing_group.create_dataset(name= "sampling_rate_Hz", data=fs)#,dtype='S')
    processing_group.create_dataset(name= "lenData", data=lenData)#,dtype='S')







####################################################################################
####################################################################################
###
### Make spectrograms H5, make parameters group in H5
###
####################################################################################
####################################################################################

## Save processing information to spectrogram H5
if os.path.isfile(SpecUFEx_H5_path): ## avoiding this error: ValueError: Unable to create group (name already exists)

    os.remove(SpecUFEx_H5_path)

    with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:

        spec_parameters_group  = fileLoad.create_group(f"spec_parameters")

        spec_parameters_group.clear()

        spec_parameters_group.create_dataset(name= 'fs', data=fs)
        spec_parameters_group.create_dataset(name= 'lenData', data=lenData)
        spec_parameters_group.create_dataset(name= 'nperseg', data=nperseg)
        spec_parameters_group.create_dataset(name= 'noverlap', data=noverlap)
        spec_parameters_group.create_dataset(name= 'nfft', data=nfft)
        spec_parameters_group.create_dataset(name= 'mode', data=mode)
        spec_parameters_group.create_dataset(name= 'scaling', data=scaling)
        spec_parameters_group.create_dataset(name= 'fmin', data=fmin)
        spec_parameters_group.create_dataset(name= 'fmax', data=fmax)






####################################################################################
####################################################################################
###
### Instantiate generator and generate spectrograms
###
####################################################################################
####################################################################################


# put sgrams in h5
gen_sgram = gen_sgram_QC_noAlias(decimation_factor=5,
                                 key = key,
                                evID_list=evID_keep,
                                dataH5_path = dataH5_path,#h5 data file
                                h5File=SpecUFEx_H5_path, #h5 sgram file
                                saveMat=False, #set true to save folder of .mat files
                                sgramOutfile='.', #path to save .mat files
                                **args
                                ) #path to save sgram figures




# def sgramH5

evID_list_QC_sgram = []
spectra_for_avg = []

less10=0
with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:

    n=0
    Nkept=0

    if 'spectrograms' in fileLoad.keys():
        del fileLoad["spectrograms"]

    if 'sgram_normConst' in fileLoad.keys():
        del fileLoad["sgram_normConst"]

    spectrograms_group     = fileLoad.create_group(f"spectrograms")

    sgram_normConst_group  = fileLoad.create_group(f"sgram_normConst")

    lenEv = len(evID_keep)

    while n <= lenEv: ## not sure a better way to execute this? But it works
        try:   #catch generator "stop iteration" error
            evID,sgram,fSTFT,tSTFT, normConstant, Nkept,evID_BADones, i = next(gen_sgram) #next() command updates generator

            n = i+1
            evID = str(evID)


            if not evID in spectrograms_group:


                spectrograms_group.create_dataset(name= evID, data=sgram)
                evID_list_QC_sgram.append(evID)
                spectra_for_avg.append(np.array(sgram))



                if not evID in sgram_normConst_group:

                    sgram_normConst_group.create_dataset(name= evID, data=normConstant)


        except StopIteration: #handle generator error
            break

    print('N events in evID_list_QC_sgram:', len(evID_list_QC_sgram))
    print('N events in evID_BADones:', len(evID_BADones))

    if 'spec_parameters' in fileLoad.keys():
        del fileLoad["spec_parameters"]




####################################################################################
####################################################################################
###
### Save more sgram parameters to H5
###
####################################################################################
####################################################################################


with h5py.File(SpecUFEx_H5_path,'a') as fileLoad:
    try:
        fSTFT_group = fileLoad.create_group(f"fSTFT")
        fSTFT_group.create_dataset(name='fSTFT', data=fSTFT)
    except:
        pass

    try:
        tSTFT_group = fileLoad.create_group(f"tSTFT")
        tSTFT_group.create_dataset(name='tSTFT', data=tSTFT)
    except:
        pass
