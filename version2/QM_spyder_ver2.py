
# -*- coding: utf-8 -*-

# %% Initialize
from time import sleep
import cmath
import math
#import scipy
#import scipy.optimize,scipy.special,scipy.stats
import time
import lmfit
from resonator_tools import circuit
import matplotlib.pyplot as plt
import numpy as np
import qcodes as qc
from qcodes import (
    Measurement,
    experiments,
    initialise_database,
    initialise_or_create_database_at,
    load_by_run_spec,
    load_or_create_experiment,
)
from qcodes.dataset.plotting import plot_dataset,plot_by_id
from qcodes.logger.logger import start_all_logging
#from resonator_tools import circuit
from qcodes.utils.dataset.doNd import do0d,do1d, do2d
from qcodes.dataset.plotting import plot_dataset
import os
from datetime import date
from qcodes.dataset.data_export import _get_data_from_ds
import qcodes.instrument_drivers.IST_devices.fastDUCK20191120 as IVVI
import qcodes.instrument_drivers.QM_qcodes.version2.configuration as config
import qm
import importlib


today=date.today()
today = today.strftime("%Y%m%d")
#sample=input("Sample name:");
sample='Gatemon_RT_5nm_2'
# define the name of the directory to be created
# path = os.path.abspath(os.getcwd())
path = "D:\Oliver\Resonator_measurements"
path= path +"\\"+today+"_"+sample
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed. Folder already exists?" % path)
else:
    print ("Successfully created the directory %s." % path)

station =qc.Station()
att_res = qc.Parameter('att', label='Attenuation(dB)', unit='dB', set_cmd=None, get_cmd=None)
att_res(40)
station.add_component(att_res)


qc.config["core"]["db_location"] = os.path.join(path, today+'_'+sample+'.db')
initialise_database()
experiment=load_or_create_experiment(experiment_name=today+'_'+sample,sample_name=sample)


dac =IVVI.IST_fastDUCK_20('IVVI', 'COM5')

dac.dac1.step =0.1
dac.dac2.step =0.1
dac.dac3.step =0.1
dac.dac4.step =0.1
dac.dac5.step =0.1
dac.dac6.step =0.1
dac.dac7.step =0.1
dac.dac8.step =0.1
dac.dac9.step =0.1
dac.dac10.step =0.1
dac.dac11.step =0.1
dac.dac12.step =0.1
dac.dac13.step =0.1
dac.dac14.step =0.1
dac.dac15.step =0.1
dac.dac16.step =0.1

station.add_component(dac)

dac.dac1.interdelay =100
dac.dac2.interdelay =100
dac.dac3.interdelay =100
dac.dac4.interdelay =100
dac.dac5.interdelay =100
dac.dac6.interdelay =100
dac.dac7.interdelay =100
dac.dac8.interdelay =100
dac.dac9.interdelay =100
dac.dac10.interdelay =100
dac.dac11.interdelay =100
dac.dac12.interdelay =100
dac.dac13.interdelay =100
dac.dac14.interdelay =100
dac.dac15.interdelay =100
dac.dac16.interdelay =100

gate=dac.dac10

dac.get_dacs()


att_drive = qc.Parameter('att_drive', label='Attenuation(dB)', unit='dB', set_cmd=None, get_cmd=None)
station.add_component(att_drive)
att_drive(40)




# %%%Initialize Octave config -> Where should it be done?

octave_config = qm.octave.QmOctaveConfig()
octave_config.add_device_info('octave1', config.octave_ip , config.octave_port)
octave_config.set_opx_octave_mapping([('con1', 'octave1')])

# %%% Open QM manager

from qcodes.instrument_drivers.QM_qcodes.version2.basic_driver.opx_driver import OPX

qmm = OPX(config=config.config, host=config.opx_ip,
              cluster_name='my_cluster_1', octave=octave_config,
              )


#%% Resonator scan

from qcodes.instrument_drivers.QM_qcodes.version2 import Resonator_scan

res_spec = Resonator_scan.Resonator_scan(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm=qmm.qmm
                                         )
station.add_component(res_spec)

#%%% Code
if_start = 90e6
if_stop = 140e6
npts = 401
lo = 4.6e9

res_spec.n_avg(5e3)
res_spec.if_start(if_start)
res_spec.if_stop(if_stop)
res_spec.n_points(npts)
res_spec.lo(lo)
res_spec.update_freqs()
res_spec.Gain(-10)
Results = res_spec.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Freq_scan_1D"
meas.register_parameter(Results, paramtype='array')
with meas.run() as datasaver:
    res_spec.run_exp()
    datasaver.add_result((Results,
                          Results()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

#%%% Save data

data = _get_data_from_ds(dataset)

freq = data[0][0]['data'].flatten()
R = data[0][1]['data'].flatten()
phi = data[3][1]['data'].flatten()

np.savetxt('background2', [freq,R,phi])

[freq,R,phi] = np.loadtxt('background')
#%%% Close
station.remove_component('ResonatorScan')
res_spec.close()


#%% Qubit Scan


from qcodes.instrument_drivers.QM_qcodes.version2 import Qubit_scan

qubit_spec = Qubit_scan.Qubit_scan(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                   )
station.add_component(qubit_spec)

# %%% Code

importlib.reload(config)

meas_freq = 4.7083e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 5.9e9
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(2e4)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-5)
qubit_spec.wait_time(5000)
Results = qubit_spec.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Qubit_scan_1D"
meas.register_parameter(Results, paramtype='array')
with meas.run() as datasaver:
 
    qubit_spec.run_exp()
    datasaver.add_result((Results,
                          Results()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')


#%%% Close
station.remove_component('QubitScan')
qubit_spec.close()

#%% Qubit scan vs gate

# All instrument uses the same QuantumMachinesManager


if_start_res = 90e6
if_stop_res = 140e6
npts_res = 401
lo_res = 4.6e9
res_spec.n_avg(5e3)
res_spec.if_start(if_start_res)
res_spec.if_stop(if_stop_res)
res_spec.n_points(npts_res)
res_spec.lo(lo_res)
res_spec.update_freqs()
res_spec.Gain(-10)


# Starting parameters
meas_freq = 4.818438e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 5.9e9
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(5e4)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-4)
Results = qubit_spec.get_measurement_parameter()

gates = np.linspace(180,170,11)
f_resonator = qc.Parameter(name='f_resonator',label='f_r',unit='Hz',get_cmd=None,set_cmd=None)
f_qubit = qc.Parameter(name='f_qubit',label='f_q',unit='Hz',get_cmd=None,set_cmd=None)
meas = Measurement(experiment, station)
meas.name = "Qubit_scan_2D"
meas.register_parameter(f_qubit,paramtype='numeric')
meas.register_parameter(f_resonator,paramtype='numeric')
meas.register_parameter(dac.dac10, paramtype='numeric')
meas.register_parameter(Results, setpoints=(dac.dac10,), paramtype='array')
with meas.run() as datasaver:
    for gate in gates:
        dac.dac10(gate)
        time.sleep(120)
        res_spec.update_freqs()
        res_spec.run_exp()
        fr = res_spec.fit_resonator()[0]
        f_resonator(fr)
        print(fr)
        if 4.5e9 < fr < 5.5e9 :
            res_spec.lo(fr-res_spec.if_start()-25e6)
            qubit_spec.lo_resonator(fr-qubit_spec.if_resonator()-2.5e6)
        qubit_spec.update_freqs()
        qubit_spec.run_exp()
        fq = qubit_spec.fit_qubit()
        f_qubit(fq)
        print(fq)
        if 5e9<fq<8e9:
            qubit_spec.lo_qubit(fq+180e6)
        Results = qubit_spec.get_measurement_parameter() #To update setpoints of the results
        datasaver.add_result((dac.dac10,dac.dac10()),
                             (f_resonator,f_resonator()),
                             (f_qubit,f_qubit()),
                             (Results,Results()))
        datasaver.flush_data_to_database()
        
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

#%% Anharmonicity


from qcodes.instrument_drivers.QM_qcodes.version2 import Anharmonicity

anharmonicity = Anharmonicity.Anharmonicity(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                   )

station.add_component(anharmonicity)

# %%% Code
importlib.reload(config)

meas_freq =  4.818438e9
if_start = -350e6
if_stop = -10e6
npts = 341
qubit_freq = 2.302e9
lo_qubit = 2.4e9
if_qubit =  qubit_freq - lo_qubit
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

anharmonicity.n_avg(2e4)
anharmonicity.if_start(if_start)
anharmonicity.if_stop(if_stop)
anharmonicity.n_points(npts)
anharmonicity.lo_qubit(lo_qubit)
anharmonicity.if_qubit(if_qubit)
anharmonicity.lo_resonator(lo_resonator)
anharmonicity.if_resonator(if_resonator)
anharmonicity.update_freqs()
anharmonicity.Gain_resonator(-10)
anharmonicity.Gain_qubit(-8)
Results = anharmonicity.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Anharmonicity_1D"
meas.register_parameter(Results, paramtype='array')
with meas.run() as datasaver:
 
    anharmonicity.run_exp()
    datasaver.add_result((Results,
                          Results()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

#%%% Simulate
anharmonicity.sim_time(10000)
anharmonicity.simulate()
anharmonicity.plot_simulated_wf()
#%%% Close
station.remove_component('Anharmonicity')
anharmonicity.close()

#%% Anharmonicity vs gate

# All instrument uses the same QuantumMachinesManager


importlib.reload(config)

if_start_res = 100e6
if_stop_res = 150e6
npts_res = 401
lo_res = 4.695e9
res_spec.n_avg(5e3)
res_spec.if_start(if_start_res)
res_spec.if_stop(if_stop_res)
res_spec.n_points(npts_res)
res_spec.lo(lo_res)
res_spec.update_freqs()
res_spec.Gain(-10)

# Starting parameters
meas_freq = 4.81819e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 2.4e9
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(5e4)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-11)

anharmonicity.n_avg(5e4)
anharmonicity.if_start(if_start)
anharmonicity.if_stop(if_stop)
anharmonicity.n_points(npts)
anharmonicity.lo_qubit(lo_qubit)
anharmonicity.lo_resonator(lo_resonator)
anharmonicity.if_resonator(if_resonator)
anharmonicity.update_freqs()
anharmonicity.Gain_resonator(-10)
anharmonicity.Gain_qubit(-8)
Results_anharmonicity = anharmonicity.get_measurement_parameter()
meas = Measurement(experiment, station)

gates = np.linspace(550,410,141)

meas.name = "Anharmonicity_2D"
f_qubit = qc.Parameter(name='f_qubit',label='f_q',unit='Hz',get_cmd=None,set_cmd=None)
f_resonator = qc.Parameter(name='f_resonator',label='f_r',unit='Hz',get_cmd=None,set_cmd=None)
meas.register_parameter(dac.dac10, paramtype='numeric')
meas.register_parameter(f_qubit, paramtype='numeric')
meas.register_parameter(f_resonator, paramtype='numeric')
meas.register_parameter(Results_anharmonicity, setpoints=(dac.dac10,), paramtype='array')

with meas.run() as datasaver:
    for gate in gates:
        dac.dac10(gate)
        time.sleep(120)
        res_spec.update_freqs()
        res_spec.run_exp()
        fr = res_spec.fit_resonator()[0]
        f_resonator(fr)
        print(fr)
        if 4.5e9 < fr < 5.5e9 :
            res_spec.lo(fr-res_spec.if_start()-25e6)
            qubit_spec.lo_resonator(fr-qubit_spec.if_resonator()-1e6)
        qubit_spec.update_freqs()
        qubit_spec.run_exp()
        fq = qubit_spec.fit_qubit()
        print(fq)
        f_qubit(fq)
        if 2e9<fq<5e9:
            qubit_spec.lo_qubit(fq+180e6)
            
        anharmonicity.if_qubit(-100e6)
        anharmonicity.lo_qubit(fq+100e6)
        anharmonicity.update_freqs()
        Results_anharmonicity = anharmonicity.get_measurement_parameter() #To update setpoints of the results
        datasaver.add_result((dac.dac10,dac.dac10()),
                             (f_qubit,f_qubit()),
                             (f_resonator,f_resonator()),
                             (Results_anharmonicity,Results_anharmonicity()))
        datasaver.flush_data_to_database()
        
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')







