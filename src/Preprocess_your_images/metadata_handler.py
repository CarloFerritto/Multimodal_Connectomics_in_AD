import os
import sys
print(sys.path)
import numpy as np
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--dataset_dir', type=str, required=True, help='Dataset Directory')
parser.add_argument('-p', '--patient', type=str, required=True, help='patient id string')
parser.add_argument('-m', '--manufacturer', type=str, required=True, help='Scanner Manufacturer') #hopefully is enough ---> otherwise we have to get the scanner type
parser.add_argument('-dd', '--dwell_dwi', type=float, required=True, help='dwell time dwi')
parser.add_argument('-rotd', '--readout_time_dwi', type=float, required=True, help='REad out time dwi ')
parser.add_argument('-rotf', '--readout_time_fmri', type=float, required=True, help='REad out time fMRi')
parser.add_argument('-ped', '--phase_encoding_dwi', type=str, required=True, help='Phase encoding dwi')
parser.add_argument('-df', '--dwell_fmri', type=float, required=True, help='dwell time fmri')
parser.add_argument('-trf', '--tr_fmri', type=float, required=True, help='Repetition Time fmri')
parser.add_argument('-pef', '--phase_encoding_fmri', type=str, required=True, help='Phase encoding fmri')
args = parser.parse_args()

os.chdir(os.path.join(args.dataset_dir,args.patient))
file_json=os.path.join('func',args.patient+'_task-rest_bold.json')

if ("Siemens" in args.manufacturer) or ("GE" in args.manufacturer):
    f = open(file_json)
    data = json.load(f)
    slicetiming=data['SliceTiming']
elif "Philips" in args.manufacturer:
    # slicetiming for philiphs is not prvided in the .json ---> this slicetiming is a guess based on https://neurostars.org/t/heudiconv-no-extraction-of-slice-timing-data-based-on-philips-dicoms/2201/13 and https://en.wikibooks.org/wiki/SPM/Slice_Timing#Philips_scanners
    slicetiming=[
		0,
		1.5,
		0.0625,
		1.5625,
		0.125,
		1.625,
		0.1875,
		1.6875,
		0.25,
		1.75,
		0.3125,
		1.8125,
		0.375,
		1.875,
		0.4375,
		1.9375,
		0.5,
		2,
		0.5625,
		2.0625,
		0.625,
		2.125,
		0.6875,
		2.1875,
		0.75,
		2.25,
		0.8125,
		2.3125,
		0.875,
		2.375,
		0.9375,
		2.4375,
		1,
		2.5,
		1.0625,
		2.5625,
		1.125,
		2.625,
		1.1875,
		2.6875,
		1.25,
		2.75,
		1.3125,
		2.8125,
		1.375,
		2.875,
		1.4375,
		2.9375]
else:
    raise Exception('Scanner Producer unknown!')


f = open("sliceTimer.txt", "x")
with open('sliceTimer.txt', 'w') as f:
    for line in slicetiming:
        f.write(f"{line}\n")


matrix = np.loadtxt(os.path.join("dwi",args.patient + "_dwi.bvec"), dtype='f', delimiter=' ')
matrix[1,:]=-matrix[1,:]
np.savetxt(os.path.join("dwi",args.patient + "_dwi_real.bvec"),matrix ,fmt='%.8f')

config_dwi=np.zeros((2,4))
#print(args.phase_encoding_dwi)
#print(args.phase_encoding_dwi=='"j"')
#print(args.phase_encoding_dwi=='"j-"')
if args.phase_encoding_dwi=='"j"':
    config_dwi[0,1]=1
    config_dwi[1,1]=-1
elif args.phase_encoding_dwi=='"j-"':
    config_dwi[0,1]=-1
    config_dwi[1,1]=1
elif args.phase_encoding_dwi=='"i"':
    config_dwi[0,0]=1
    config_dwi[1,0]=-1
elif args.phase_encoding_dwi=='"i-"':
    config_dwi[0,0]=-1
    config_dwi[1,0]=1
elif args.phase_encoding_dwi=='"k"':
    config_dwi[0,2]=1
    config_dwi[1,2]=-1
elif args.phase_encoding_dwi=='"k-"':
    config_dwi[0,2]=-1
    config_dwi[1,2]=1
config_dwi[0,3]=args.readout_time_dwi
config_dwi[1,3]=args.readout_time_dwi
np.savetxt("config_dwi.txt", config_dwi, fmt='%.8f', delimiter=' ')

config_fmri=np.zeros((2,4))
if args.phase_encoding_fmri=='"j"':
    config_fmri[0,1]=1
    config_fmri[1,1]=-1
elif args.phase_encoding_fmri=='"j-"':
    config_fmri[0,1]=-1
    config_fmri[1,1]=1
elif args.phase_encoding_fmri=='"i"':
    config_fmri[0,0]=1
    config_fmri[1,0]=-1
elif args.phase_encoding_fmri=='"i-"':
    config_fmri[0,0]=-1
    config_fmri[1,0]=1
elif args.phase_encoding_fmri=='"k"':
    config_fmri[0,2]=1
    config_fmri[1,2]=-1
elif args.phase_encoding_fmri=='"k-"':
    config_fmri[0,2]=-1
    config_fmri[1,2]=1
config_fmri[0,3]=args.readout_time_fmri
config_fmri[1,3]=args.readout_time_fmri
np.savetxt("config_fmri.txt", config_fmri, fmt='%.8f', delimiter=' ')