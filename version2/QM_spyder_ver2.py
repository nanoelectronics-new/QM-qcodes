
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
from qm import generate_qua_script
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

drag_param = qc.Parameter(name='alpha',label='alpha',unit='',get_cmd=None,set_cmd=None)
station.add_component(drag_param)

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
lo = 4.62e9

res_spec.n_avg(5e3)
res_spec.if_start(if_start)
res_spec.if_stop(if_stop)
res_spec.n_points(npts)
res_spec.lo(lo)
res_spec.update_freqs()
res_spec.Gain(-10)
Resonator_Results = res_spec.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Freq_scan_1D"
meas.register_parameter(Resonator_Results, paramtype='array')
with meas.run() as datasaver:
    res_spec.run_exp()
    datasaver.add_result((Resonator_Results,
                          Resonator_Results()))
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


#%% ADC

#%%% Import
from qcodes.instrument_drivers.QM_qcodes.version2 import ADC
importlib.reload(config)

ADC = ADC.ADC(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm=qmm.qmm
                                         )
station.add_component(ADC)

#%%% Code

ADC.readout_operation('readout')

if_freq = 100e6

lo =4.154123e9

ADC.n_avg(1e4)
ADC.IF(if_freq)
ADC.lo(lo)
ADC.Gain(0)
ADC.wait_time(5000)
ADC.readout_pulse_length(68+180)
ADC.update_freqs()
Results = ADC.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "ADC_1D_sliced"
meas.register_parameter(Results, paramtype='array')
with meas.run() as datasaver:
    ADC.run_exp()
    datasaver.add_result((Results,
                          Results()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')


# %%% Simulate
ADC.sim_time(10000)
ADC.simulate()
ADC.plot_simulated_wf()


#%%% Close
station.remove_component('ADC')
ADC.close()


#%% Qubit Scan


from qcodes.instrument_drivers.QM_qcodes.version2 import Qubit_scan
importlib.reload(config)

qubit_spec = Qubit_scan.Qubit_scan(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                   )
station.add_component(qubit_spec)

# %%% Code


meas_freq = 4.735796e9
if_start = -350e6
if_stop = -10e6
npts = 171
lo_qubit = 6.55e9
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(1e5)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-3)
qubit_spec.wait_time(10000)
Qubit_Results = qubit_spec.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Qubit_scan_1D"
meas.register_parameter(Qubit_Results, paramtype='array')
with meas.run() as datasaver:
 
    qubit_spec.run_exp()
    datasaver.add_result((Qubit_Results,
                          Qubit_Results()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

#%%% Plot

from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id=2)
data =_get_data_from_ds(dataset)

time = data[3][0]['data'].flatten()
I = data[0][1]['data'].flatten()
Q = data[1][1]['data'].flatten()

plt.figure()
plt.scatter(I,Q)
#%%% Qubit scan vs Power

importlib.reload(config)

meas_freq = 4.822e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 2.85e9
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
qubit_spec.Gain_qubit(0)
qubit_spec.wait_time(10000)
gains = np.arange(-12,0,0.5)
Results = qubit_spec.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Qubit_scan_2D_power"
meas.register_parameter(qubit_spec.Gain_qubit,paramtype='numeric')
meas.register_parameter(Results, setpoints=(qubit_spec.Gain_qubit,), paramtype='array')
with meas.run() as datasaver:
    for gain in gains:
        qubit_spec.Gain_qubit(gain)
        qubit_spec.run_exp()
        datasaver.add_result((Results,Results()),
                             (qubit_spec.Gain_qubit, qubit_spec.Gain_qubit()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% Simulate
qubit_spec.sim_time(200000)
qubit_spec.simulate()
qubit_spec.plot_simulated_wf()

#%%% Close
station.remove_component('QubitScan')
qubit_spec.close()

#%% Qubit scan vs gate

# All instrument uses the same QuantumMachinesManager


if_start_res = 90e6
if_stop_res = 140e6
npts_res = 401
lo_res = 4.7e9
res_spec.n_avg(5e3)
res_spec.if_start(if_start_res)
res_spec.if_stop(if_stop_res)
res_spec.n_points(npts_res)
res_spec.lo(lo_res)
res_spec.update_freqs()
res_spec.Gain(-10)


# Starting parameters
meas_freq = 4.817438e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 2.4e9
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(1e5)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-8)
qubit_spec.wait_time(10000)
Results = qubit_spec.get_measurement_parameter()

gates = np.linspace(550,410,141)
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
        time.sleep(180)
        res_spec.update_freqs()
        res_spec.run_exp()
        fr = res_spec.fit_resonator()[0]
        f_resonator(fr)
        print(fr)
        if 4.6e9 < fr < 5.2e9 :
            res_spec.lo(fr-res_spec.if_start()-25e6)
            qubit_spec.lo_resonator(fr-qubit_spec.if_resonator()-1e6)
        qubit_spec.update_freqs()
        qubit_spec.run_exp()
        fq = qubit_spec.fit_qubit()
        f_qubit(fq)
        print(fq)
        if 2e9<fq<4.2e9:
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


meas_freq = 4.82187e9
if_start = -350e6
if_stop = -10e6
npts = 341
qubit_freq = 2.68322e9
if_qubit =  -100e6
lo_qubit = qubit_freq-if_qubit
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

anharmonicity.n_avg(1e5)
anharmonicity.if_start(if_start)
anharmonicity.if_stop(if_stop)
anharmonicity.n_points(npts)
anharmonicity.lo_qubit(lo_qubit)
anharmonicity.if_qubit(if_qubit)
anharmonicity.lo_resonator(lo_resonator)
anharmonicity.if_resonator(if_resonator)
anharmonicity.update_freqs()
anharmonicity.Gain_resonator(-10)
anharmonicity.Gain_qubit(-10)
qubit_spec.wait_time(10000)
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

if_start_res = 90e6
if_stop_res = 140e6
npts_res = 401
lo_res = 4.7e9
res_spec.n_avg(5e3)
res_spec.if_start(if_start_res)
res_spec.if_stop(if_stop_res)
res_spec.n_points(npts_res)
res_spec.lo(lo_res)
res_spec.update_freqs()
res_spec.Gain(-10)

# Starting parameters
meas_freq = 4.81919e9
if_start = -350e6
if_stop = -10e6
npts = 341
lo_qubit = 2.4e9
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

qubit_spec.n_avg(1e5)
qubit_spec.if_start(if_start)
qubit_spec.if_stop(if_stop)
qubit_spec.n_points(npts)
qubit_spec.lo_qubit(lo_qubit)
qubit_spec.lo_resonator(lo_resonator)
qubit_spec.if_resonator(if_resonator)
qubit_spec.update_freqs()
qubit_spec.Gain_resonator(-10)
qubit_spec.Gain_qubit(-9)

anharmonicity.n_avg(1e5)
anharmonicity.if_start(if_start)
anharmonicity.if_stop(if_stop)
anharmonicity.n_points(npts)
anharmonicity.lo_qubit(lo_qubit)
anharmonicity.lo_resonator(lo_resonator)
anharmonicity.if_resonator(if_resonator)
anharmonicity.update_freqs()
anharmonicity.Gain_resonator(-10)
anharmonicity.Gain_qubit(-7)
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
        time.sleep(180)
        res_spec.update_freqs()
        res_spec.run_exp()
        fr = res_spec.fit_resonator()[0]
        f_resonator(fr)
        print(fr)
        if 4.5e9 < fr < 5.5e9 :
            res_spec.lo(fr-res_spec.if_start()-25e6)
            qubit_spec.lo_resonator(fr-qubit_spec.if_resonator()-1e6)
            anharmonicity.lo_resonator(fr-qubit_spec.if_resonator()-1e6)
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
        anharmonicity.run_exp()
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


# %% Time Rabi

#%%% Import
from qcodes.instrument_drivers.QM_qcodes.version2 import Rabi
importlib.reload(config)

Rabi = Rabi.Rabi(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                   )

station.add_component(Rabi)

#%%% Code


meas_freq = 4.747967e9
t_start = 5
t_stop = 100
npts = 96
qubit_freq = 6.74163e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator

drag_param(config.alpha_drag)

Rabi.set_config(config.config)

Rabi.n_avg(2e5)
Rabi.t_start(t_start)
Rabi.t_stop(t_stop)
Rabi.n_points(npts)
Rabi.lo_qubit(lo_qubit)
Rabi.if_qubit(if_qubit)
Rabi.lo_resonator(lo_resonator)
Rabi.if_resonator(if_resonator)
Rabi.update_freqs()
Rabi.Gain_resonator(-10)
Rabi.Gain_qubit(13)
Rabi.wait_time(10000)
Rabi.acquisition_mode('full_demodulation')
Rabi_Results = Rabi.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Rabi_1D"
meas.register_parameter(Rabi_Results, paramtype='array')

with meas.run() as datasaver:
    Rabi.run_exp()
    datasaver.add_result((Rabi_Results,
                          Rabi_Results()),)
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% Simulate
Rabi.sim_time(200000)
Rabi.simulate()
Rabi.plot_simulated_wf()


# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id=29)
data =_get_data_from_ds(dataset)

time = data[3][0]['data'].flatten()
mag = data[2][1]['data'].flatten()
phase = data[3][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, freq, phase, c, d):
    return a * np.exp(-b * x)*np.cos(2*np.pi*freq*x+phase) + c + d*x
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[2,0.1,0.05,0,-89,0])
plt.scatter(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{Rabi}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')


# %%% Rabi 2D_IF

meas_freq = 4.822e9
t_start = 5
t_stop = 100
npts = 96
qubit_freq = 2.763e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator

drag_param(config.alpha_drag)

Rabi.set_config(config.config)

Rabi.n_avg(1e5)
Rabi.t_start(t_start)
Rabi.t_stop(t_stop)
Rabi.n_points(npts)
Rabi.lo_qubit(lo_qubit)
Rabi.if_qubit(if_qubit)
Rabi.lo_resonator(lo_resonator)
Rabi.if_resonator(if_resonator)
Rabi.update_freqs()
Rabi.Gain_resonator(-10)
Rabi.Gain_qubit(12)
Rabi.wait_time(10000)
Results = Rabi.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Rabi_2D_IF"
meas.register_parameter(Rabi.if_qubit,paramtype='numeric')
meas.register_parameter(Results, setpoints=(Rabi.if_qubit,), paramtype='array')
freqs = np.arange(-350e6,-50e6,50e6)
with meas.run() as datasaver:
    for freq in freqs:
        Rabi.if_qubit(freq)
        lo =  qubit_freq - freq
        Rabi.lo_qubit(lo)
        Rabi.update_freqs()
        Rabi.Gain_resonator(-10)
        Rabi.Gain_qubit(15)
        Rabi.run_exp()
        datasaver.add_result((Results,Results()),(Rabi.if_qubit,Rabi.if_qubit()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% Rabi sliced

meas_freq = 4.81805e9
t_start = 5
t_stop = 40
npts = 36
qubit_freq = 2.278e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator

drag_param(config.alpha_drag)

Rabi.set_config(config.config)

Rabi.n_avg(2e4)
Rabi.t_start(t_start)
Rabi.t_stop(t_stop)
Rabi.n_points(npts)
Rabi.lo_qubit(lo_qubit)
Rabi.if_qubit(if_qubit)
Rabi.lo_resonator(lo_resonator)
Rabi.if_resonator(if_resonator)
Rabi.readout_pulse_length(68+240)
Rabi.update_freqs()
Rabi.Gain_resonator(-10)
Rabi.Gain_qubit(10)
Rabi.wait_time(10000)

Rabi.acquisition_mode('sliced_demodulation')
Results = Rabi.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "Rabi_1D_sliced_demod"
meas.register_parameter(Results, paramtype='array')

with meas.run() as datasaver:
    Rabi.run_exp()
    datasaver.add_result((Results,Results()),)
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
# %%% Remove

Rabi.close()
station.remove_component('Rabi')


#%% T1

#%%% Import
from qcodes.instrument_drivers.QM_qcodes.version2 import T1

T1 = T1.T1(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                   )

station.add_component(T1)

#%%% Code

importlib.reload(config)


meas_freq = 4.747967e9
t_start = 4 # in clock cycles
t_stop = 60 # in clock cycles
npts = 57
qubit_freq = 6.74163e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator



T1.n_avg(1e5)
T1.t_start(t_start)
T1.t_stop(t_stop)
T1.n_points(npts)
T1.lo_qubit(lo_qubit)
T1.if_qubit(if_qubit)
T1.lo_resonator(lo_resonator)
T1.if_resonator(if_resonator)
T1.update_freqs()
T1.Gain_resonator(-10)
T1.Gain_qubit(15)
T1.wait_time(10000)
Results = T1.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "T1"
meas.register_parameter(Results, paramtype='array')

with meas.run() as datasaver:
    T1.run_exp()
    datasaver.add_result((Results,Results()),)
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% Simulate
T1.sim_time(200000)
T1.simulate()
T1.plot_simulated_wf()

# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id=65)
data =_get_data_from_ds(dataset)

time = data[3][0]['data'].flatten()
mag = data[2][1]['data'].flatten()
phase = data[3][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, c):
    return a * np.exp(-b * x) + c
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[1,0.1,-158])
plt.scatter(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{1}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')
# %%% Remove

T1.close()
station.remove_component('T1')

# %% T2_short

from qcodes.instrument_drivers.QM_qcodes.version2 import T2_short
importlib.reload(T2_short)

T2meas = T2_short.T2_short(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                )

station.add_component(T2meas)


#%%% Code

importlib.reload(config)


meas_freq = 4.747967e9
t_start = 0
t_stop = 100
npts = 101
qubit_freq = 6.74163e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator



T2meas.n_avg(2e5)
T2meas.t_start(t_start)
T2meas.t_stop(t_stop)
T2meas.n_points(npts)
T2meas.lo_qubit(lo_qubit)
T2meas.if_qubit(if_qubit)
T2meas.lo_resonator(lo_resonator)
T2meas.if_resonator(if_resonator)
T2meas.update_freqs()
T2meas.Gain_resonator(-10)
T2meas.Gain_qubit(14)
T2meas.detuning(50e6)
T2meas.alpha_drag(0.15)
T2meas.t_pihalf(6)
T2meas.wait_time(10000)
Results = T2meas.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "T2meas"
meas.register_parameter(Results, paramtype='array')

with meas.run() as datasaver:
    T2meas.run_exp()
    datasaver.add_result((Results,Results()),)
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

 # %%% Simulate
T2meas.sim_time(100000)
T2meas.simulate()
T2meas.plot_simulated_wf()

# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id=66)
data =_get_data_from_ds(dataset)

time = data[3][0]['data'].flatten()
mag = data[2][1]['data'].flatten()
phase = data[3][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, freq, phase, c, d):
    return a * np.exp(-b * x)*np.cos(2*np.pi*freq*x+phase) +c + d*x
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[1,0.01,0.05,1.7,-158,0])
plt.plot(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{2}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')


# %%% Ramsey 2D


meas_freq = 4.829e9
t_start = 0
t_stop = 150
npts = 151
qubit_freq = 3.03474e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator



T2meas.n_avg(1e5)
T2meas.t_start(t_start)
T2meas.t_stop(t_stop)
T2meas.n_points(npts)
T2meas.lo_qubit(lo_qubit)
T2meas.if_qubit(if_qubit)
T2meas.lo_resonator(lo_resonator)
T2meas.if_resonator(if_resonator)
T2meas.update_freqs()
T2meas.Gain_resonator(-10)
T2meas.Gain_qubit(3)
T2meas.detuning(40e6)
T2meas.alpha_drag(0)
T2meas.t_pihalf(9)
T2meas.wait_time(10000)
detunings = np.linspace(-50e6, 50e6,51)
Results = T2meas.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "T2meas"
meas.register_parameter(T2meas.detuning,paramtype='numeric')
meas.register_parameter(Results, setpoints=(T2meas.detuning,), paramtype='array')

with meas.run() as datasaver:
    for value in detunings:
        T2meas.detuning(value)
        T2meas.run_exp()
        datasaver.add_result((T2meas.detuning,T2meas.detuning()), (Results,Results()),)
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% Remove

T2meas.close()
station.remove_component('T2_short')

# %% T2_long

from qcodes.instrument_drivers.QM_qcodes.version2 import T2
importlib.reload(T2)

T2long = T2.T2(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                )

station.add_component(T2long)

#%%% Code


importlib.reload(config)



meas_freq = 4.7469729e9
t_start = 4
t_stop = 400
npts = 199
qubit_freq = 6.73531e9
if_qubit = -250e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 125e6
lo_resonator = meas_freq - if_resonator



T2long.n_avg(2e5)
T2long.t_start(t_start)
T2long.t_stop(t_stop)
T2long.n_points(npts)
T2long.lo_qubit(lo_qubit)
T2long.if_qubit(if_qubit)
T2long.lo_resonator(lo_resonator)
T2long.if_resonator(if_resonator)
T2long.update_freqs()
T2long.Gain_resonator(-10)
T2long.Gain_qubit(15)
T2long.detuning(40e6)
T2long.wait_time(1000)
Results = T2long.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "T2long"
meas.register_parameter(Results, paramtype='array')


# %%% Serialized file
sourceFile = open('debug.py', 'w')
print(generate_qua_script(T2long.get_prog(), T2long.config), file=sourceFile)
sourceFile.close()




# %% T2__Bias


from qcodes.instrument_drivers.QM_qcodes.version2 import T2_short_bias
importlib.reload(T2_short_bias)
importlib.reload(config)

T2bias = T2_short_bias.T2_short_bias(config=config.config, host=config.opx_ip,
                                         cluster_name='my_cluster_1', octave=octave_config, qmm = qmm.qmm
                                )

station.add_component(T2bias)


#%%% Code




meas_freq = 4.829734335e9
t_start = 4
t_stop = 150
npts = 147
qubit_freq = 3.06627e9
if_qubit = -350e6
lo_qubit =  qubit_freq - if_qubit
if_resonator = 250e6
lo_resonator = meas_freq - if_resonator



T2bias.n_avg(1e5)
T2bias.t_start(t_start)
T2bias.t_stop(t_stop)
T2bias.n_points(npts)
T2bias.lo_qubit(lo_qubit)
T2bias.if_qubit(if_qubit)
T2bias.lo_resonator(lo_resonator)
T2bias.if_resonator(if_resonator)
T2bias.update_freqs()
T2bias.Gain_resonator(-10)
T2bias.Gain_qubit(5)
T2bias.bias_amp(0.3)
T2bias.alpha_drag(0)
T2bias.t_pihalf(8)
T2bias.wait_time(10000)
Results = T2bias.get_measurement_parameter()
meas = Measurement(experiment, station)
meas.name = "T2bias"
meas.register_parameter(Results, paramtype='array')

with meas.run() as datasaver:
    T2bias.run_exp()
    datasaver.add_result((Results,Results()),)
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)

[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

 # %%% Simulate
T2bias.sim_time(1000)
T2bias.simulate()
T2bias.plot_simulated_wf()

# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id=17)
data =_get_data_from_ds(dataset)

time = data[3][0]['data'].flatten()
mag = data[2][1]['data'].flatten()
phase = data[3][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, freq, phase, c, d):
    return a * np.exp(-b * x)*np.cos(2*np.pi*freq*x+phase) +c + d*x
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[1,0.01,0.05,1.7,10,143])
plt.plot(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{2}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')
# %%% Remove

T2bias.close()
station.remove_component('T2_short_bias')



# %% Two-probe

from qcodes.parameters import ScaledParameter
from qcodes.instrument_drivers.Keysight.Keysight_34465A import Keysight_34465A   #import the driver

DMM = Keysight_34465A('DMM', 'TCPIP0::10.21.41.141')
station.add_component(DMM)

#%%% parameters
current = ScaledParameter(DMM.volt, gain = 100, name='Current', unit = 'nA') #gain = 1: 1V/nA:1G, 10:100M, 100:10M
bias = ScaledParameter(dac.dac7, division = 1, name = 'Bias', unit = 'mV')
gate = ScaledParameter(dac.dac8, division = 1, name = 'Gate', unit = 'mV')

#%%% Bias sweep
biasstart = -10 
biasstop = 10
npts = 101

biassweep = np.linspace(biasstart, biasstop, npts)
meas = Measurement(experiment, station)
meas.name = "junc2_biassweep_gate-2V"
meas.register_parameter(bias,paramtype='numeric')
meas.register_parameter(current, setpoints=(bias, ), paramtype='numeric')
with meas.run() as datasaver:
    
    for b in biassweep:
        bias(b)
        datasaver.add_result((bias,bias()),
                             (current,current()))
        datasaver.flush_data_to_database()
    bias(0)    
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

#%%% Gate sweep


gatestart = 2000
gatestop = -3000
npts = 2501

gatesweep = np.linspace(gatestart, gatestop, npts)
meas = Measurement(experiment, station)
meas.name = "junc1_gatesweep"
meas.register_parameter(gate,paramtype='numeric')
meas.register_parameter(current, setpoints=(gate, ), paramtype='numeric')
with meas.run() as datasaver:
    bias(10)
    for g in gatesweep:
        gate(g)
        datasaver.add_result((gate,gate()),
                             (current,current()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
    gate(0)
    bias(0)
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s'%(meas.name,dataset.captured_run_id),bbox_inches='tight')


# %% Four-probe

voltage = ScaledParameter(DMM.volt, gain = 1e1, name='Voltage', unit = 'mV') #gain: 1e2 = 10 , 1e1 = 100 , 1=1k
current = ScaledParameter(dac.dac1, gain = 1e-1, name = 'Current', unit = 'nA')# gain: #1e-1=100nA, 1 = 1uA, 1e1=10uA
gate = ScaledParameter(dac.dac8, division = 1, name = 'Gate', unit = 'mV')

#%%% Current sweep


currentstart = -50
currentstop = 50
npts = 101

currentsweep = np.linspace(currentstart, currentstop, npts)
meas = Measurement(experiment, station)
meas.name = "junc1_currentsweep"
meas.register_parameter(current,paramtype='numeric')
meas.register_parameter(voltage, setpoints=(current, ), paramtype='numeric')
with meas.run() as datasaver:
    current(currentstart)
    for c in currentsweep:
        current(c) 
        datasaver.add_result((current,current()),
                             (voltage,voltage()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
    current(0)
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s'%(meas.name,dataset.captured_run_id),bbox_inches='tight')




