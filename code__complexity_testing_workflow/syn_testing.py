

import sys
json_fil_name = sys.argv[1]
synaptic_rate_no = int(sys.argv[2])



import os
# json_fil_name = os.path.basename(json_fil_name)


cell_name = (os.path.splitext(os.path.basename(json_fil_name))[0].split('_stage_')[0] )
swc_file_name = cell_name +'_transformed.swc'

from morphology_ult import dendrite_surface_area
swc_file_location = './Required_Files/swc/' + swc_file_name

surface_area = dendrite_surface_area(swc_file_location)

print(int(surface_area))

seed_no = (os.path.splitext(os.path.basename(json_fil_name))[0].split('_stage_')[1].split('_seed_')[1].split('_hof_')[0] )
hof_no = (os.path.splitext(os.path.basename(json_fil_name))[0].split('_stage_')[1].split('_seed_')[1].split('_hof_')[1].split('_')[0] )


simulation_location = './' + cell_name + '_stage_1_seed_' + str(seed_no) + '_hof_' + str(hof_no) + '_rate_' + str(synaptic_rate_no) + '/'


# basic
dt = 0.02
dL = 20
simulation_duration = 16200.0


time_sub_sampling = 10
# space_sub_sampling = 128

sample_rate = 1 / dt
filter_f = 0.01  # 1/ms


sync_per_area = 0.45


plotting_window_0 = int(1000 / dt)
plotting_window = int(2000 / dt)

lcz_window_0 = 10000

synaptic_rate = synaptic_rate_no * 0.001






# Making model 


import cell_functions #Imports csmc_allactive_noaxon
model_processing = 'csmc_allactive_noaxon'


import os
try:
    os.system("rm -r " + simulation_location)  
except:
    pass


from bmtk.builder.networks import NetworkBuilder

# Create network
net = NetworkBuilder('single_neuron')
net.add_nodes(cell_name=cell_name,
              potental='exc',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing=model_processing,
              dynamics_params=json_fil_name,
              morphology=swc_file_name)

net.build()
net.save_nodes(output_dir=simulation_location+'network')



# Create inputs
thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(
    N=int(surface_area * sync_per_area),
    pop_name='tON',
    potential='exc',
    model_type='virtual'
)

thalamus.add_edges(
    source={'pop_name': 'tON'}, target=net.nodes(),
    connection_rule=1,
    syn_weight=0.006,
    delay=2.0,
    weight_function=None,
    target_sections=['basal', 'apical'],
    distance_range=[0.0, 1e10],
    dynamics_params='human_AMPA_ExcToExc.json',
    model_template='exp2syn'
)


thalamus.build()
thalamus.save_nodes(output_dir=simulation_location + 'network')
thalamus.save_edges(output_dir=simulation_location + 'network')





from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

psg = PoissonSpikeGenerator(population='mthalamus')
psg.add(
    node_ids=range( int(surface_area * sync_per_area) ),  # Have 10 nodes to match mthalamus
    firing_rate=synaptic_rate,    # 10 Hz, we can also pass in a nonhomoegenous function/array
    times=(0.0, simulation_duration*0.001)    # Firing starts at 0 s up to 3 s
)
psg.to_sonata(simulation_location + 'inputs/mthalamus_spikes.h5')




from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(
    base_dir= simulation_location,
    config_file='config.json',
    network_dir=simulation_location+'network',
    tstop=simulation_duration,
    dt=dt,
    dL=dL,
    include_examples=False,
    compile_mechanisms=False,
    nsteps_block = 25000,
    spikes_inputs=[('mthalamus', 
                simulation_location + 'inputs/mthalamus_spikes.h5')],

)



# Copy required files
import shutil
import glob

shutil.copyfile(  './Required_Files/json/' + json_fil_name, simulation_location + 'components/biophysical_neuron_models/'+json_fil_name)
shutil.copyfile( './Required_Files/swc/' + swc_file_name, simulation_location + 'components/morphologies/'+swc_file_name)
shutil.copytree('./Required_Files/modfiles', simulation_location + 'components/mechanisms/modfiles')

for file in glob.glob('./Required_Files/synaptic_models/*'):
    shutil.copy(file, simulation_location + 'components/synaptic_models/')





import json

with open(simulation_location +'simulation_config.json') as f:
    sim_conf = json.load(f)


sim_conf.update( 
    {
    "reports": {
    "v_report_all": {
      "variable_name": "v",
      "cells": "all",
      "module": "membrane_report",
      "sections": "all"
            }
        }
    }
)



with open(simulation_location +'simulation_config.json', "w") as outfile:
    outfile.write(json.dumps(sim_conf, indent=4))





from subprocess import call

status=call("nrnivmodl modfiles",cwd = simulation_location + "components/mechanisms",shell=True)
if status:
    print ("NEURON ERROR")
else:
    print ("Compilation done!")



from bmtk.simulator import bionet
import cell_functions #Imports csmc_allactive_noaxon

conf = bionet.Config.from_json(simulation_location +'config.json')
conf.build_env()

net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)



sim.run()


# 
# 
# 
# Post Processing
# 





import h5py
import numpy as np

filename = simulation_location + "output/v_report_all.h5"
f = h5py.File(filename, 'r')
with h5py.File(filename, 'r') as f:
    dd = f['report']['single_neuron']['data'][()]

dd = np.transpose(np.array(dd))

try:
    from reduce_morph import reducr_morphology
    rdc_idx = reducr_morphology(simulation_location)
    dd_rdx = dd[rdc_idx, :]
    from plt_traces import plt_six_traces
    plt_six_traces(dd_rdx, plotting_window_0, plotting_window, sample_rate, simulation_location)
except:
    print("Something else went wrong")



dd = dd[:, lcz_window_0:]

ro,co=dd.shape
TH=np.zeros(ro)


for i in range(ro):
    TH[i]=np.mean(dd[i,:])

dd = dd-np.reshape(TH, (-1, 1))

def lowpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 10):
    import scipy.signal
    sos = scipy.signal.butter(poles, cutoff, 'lowpass', fs=sample_rate, output='sos')
    filtered_data = scipy.signal.sosfiltfilt(sos, data)
    return filtered_data


# def highpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
#     import scipy.signal
#     sos = scipy.signal.butter(poles, cutoff, 'highpass', fs=sample_rate, output='sos')
#     filtered_data = scipy.signal.sosfiltfilt(sos, data)
#     return filtered_data


dd_lp=np.zeros((ro,co))
for i in range(ro):
    dd_lp[i,:] = lowpass(dd[i,:], filter_f, sample_rate)

dd_hp = dd - dd_lp


from plt_traces import plt_three_traces

plt_three_traces(dd, dd_lp, dd_hp, plotting_window_0, plotting_window, sample_rate, filter_f, simulation_location )






def binarization_abs_hilbert(signal_in):
    from scipy.signal import hilbert
    import numpy as np
    ro,co=signal_in.shape
    TH=np.zeros(ro)
    for i in range(ro):
        signal_in[i,:]=abs(hilbert(signal_in[i,:]))
        TH[i]=np.mean(signal_in[i,:])
    for i in range(ro):
        signal_in[i,:]= np.where(signal_in[i,:]>TH[i], 0, 1)
    return signal_in

dd_hp = binarization_abs_hilbert(dd_hp)
dd_lp = binarization_abs_hilbert(dd_lp)



import math
from complexity.lzc.lzc import lz_complexity


ddsp_hp = dd_hp[:, time_sub_sampling::time_sub_sampling]
ddsp_lp = dd_lp[:, time_sub_sampling::time_sub_sampling]


ro,co=ddsp_hp.shape


lzc_local_hp=np.zeros(ro)
lzc_local_lp=np.zeros(ro)

for i in range(ro):
    lzc_local_hp[i] = lz_complexity(ddsp_hp[i,:]) /( len(ddsp_hp[i,:]) / math.log2(len(ddsp_hp[i,:])) )
    lzc_local_lp[i] = lz_complexity(ddsp_lp[i,:]) /( len(ddsp_lp[i,:]) / math.log2(len(ddsp_lp[i,:])) )
    if ((i%10) == 0):
        print('Calculating local LZC: ' + str(i) + ' / ' + str(ro))

print('Soma hf LZC :'  + str(lzc_local_hp[0])+ ' at rate:' + str(synaptic_rate))
print('Soma lf LZC :'  + str(lzc_local_lp[0])+ ' at rate:' + str(synaptic_rate))



lzc_time_hp=np.zeros(co)
lzc_time_lp=np.zeros(co)



for i in range(co):
    temp = np.transpose(ddsp_hp)[i,:]
    lzc_time_hp[i] = lz_complexity(temp) /( len(temp) / math.log2(len(temp)) )
    temp = np.transpose(ddsp_lp)[i,:]
    lzc_time_lp[i] = lz_complexity(temp) /( len(temp) / math.log2(len(temp)) )
    if ((i%10000) == 0):
        print('Calculating time LZC: ' + str(i) + ' / ' + str(co))    

lzc_global_hp = np.mean(lzc_time_hp)
lzc_globla_lp = np.mean(lzc_time_lp)


shape_of_sliced_array = ddsp_hp.shape



# sliced_array = ddsp_hp[0::space_sub_sampling, :]
# shape_of_sliced_array = sliced_array.shape

# ddsp_hp = np.reshape(np.transpose(ddsp_hp[0::space_sub_sampling,:]),-1)
# print('Calculating global hf LZC ')
# lzc_global_hp = lz_complexity(ddsp_hp) /(len(ddsp_hp) / math.log2(len(ddsp_hp)))
# print('Global hf LZC :'  + str(lzc_global_hp)+ ' at rate:' + str(synaptic_rate))


# ddsp_lp = np.reshape(np.transpose(ddsp_lp[0::space_sub_sampling,:]),-1)
# print('Calculating global lf LZC ')
# lzc_globla_lp = lz_complexity(ddsp_lp) /(len(ddsp_lp) / math.log2(len(ddsp_lp)))
# print('Global lf LZC :'  + str(lzc_globla_lp) + ' at rate:' + str(synaptic_rate))




with h5py.File(filename, 'r+') as f:
    del f['report']['single_neuron']['data']
    f['report']['single_neuron']['lzc_local_hp'] = lzc_local_hp
    f['report']['single_neuron']['lzc_local_lp'] = lzc_local_lp
    f['report']['single_neuron']['lzc_hp'] = lzc_global_hp
    f['report']['single_neuron']['lzc_lp'] = lzc_globla_lp
    f['report']['single_neuron']['shape'] = shape_of_sliced_array

