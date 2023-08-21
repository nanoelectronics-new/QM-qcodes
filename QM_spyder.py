# -*- coding: utf-8 -*-

# %% Initialize
import os
from time import sleep
import cmath
import math
#import scipy
#import scipy.optimize,scipy.special,scipy.stats
import time
#import lmfit
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
import qcodes.instrument_drivers.rohde_schwarz.SGS100A as sgs
import qcodes.instrument_drivers.IST_devices.fastDUCK20191120 as IVVI
from qcodes_contrib_drivers.drivers.RohdeSchwarz import SMW200A as smw200a

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


# sg_src=sgs.RohdeSchwarz_SGS100A(name='sg_src',address="TCPIP0::rssgs100a111062::inst0::INSTR")
# sg_src.status('OFF')
# sg_src.frequency(4.8e9)
# station.add_component(sg_src)

# sg_src_2=sgs.RohdeSchwarz_SGS100A(name='sg_src_2',address="TCPIP0::rssgs100a111502::inst0::INSTR"),
# sg_src_2.status('OFF')
# sg_src_2.frequency(4.8e9)
# station.add_component(sg_src_2)

# vsg = smw200a.RohdeSchwarz_SMW200A(name='SMW200A5', address='TCPIP::10.21.41.171::hislip0::INSTR' )
# station.add_component(vsg.rfoutput1.frequency)
# station.add_component(vsg.rfoutput1.level)
# station.add_component(vsg.rfoutput1.state)
att_drive = qc.Parameter('att_drive', label='Attenuation(dB)', unit='dB', set_cmd=None, get_cmd=None)
station.add_component(att_drive)
att_drive(40)
# vsg.rfoutput1.frequency(2.2e9)
# vsg.rfoutput1.level(-50)

# vsg.rfoutput1.state('OFF')
# %% Import OPX driver

import qcodes.instrument_drivers.QM_qcodes.opx_driver as OPX
import qcodes.instrument_drivers.QM_qcodes.configuration as config
import importlib
importlib.reload(config)
#opx1 =OPX.OPX(config=config.config,host=config.qop_ip,port=config.qop_port)

#%% Manual output control

from qualang_tools.control_panel import ManualOutputControl
manual_output_control = ManualOutputControl(config.config,host = config.opx_ip,port = config.opx_port)
manual_output_control.turn_on_element('resonator')
manual_output_control.set_amplitude('resonator', 0.1)
manual_output_control.set_frequency('resonator', 50e6)

manual_output_control.close()


#%% Bleufors TC Com

import json 
import requests
DEVICE_IP = '10.21.41.151'
TIMEOUT = 10
# %%% Read temp

def read_temp():
    
    url = 'http://{}:5001/channel/measurement/latest'.format(DEVICE_IP)
    req = requests.get(url, timeout=TIMEOUT)
    data = req.json()
    return data['temperature']


# %%% Set temp
# -*- coding: utf-8 -*-
def set_temp(setpoint):

    url = 'http://{}:5001/heater/update'.format(DEVICE_IP)
    data = {
     'heater_nr': 4,
     'setpoint': setpoint
    }

    requests.post(url, json=data, timeout=TIMEOUT)

    return 0

mxc_temp = qc.Parameter(name='Temperature',label='T',unit='K',get_cmd=read_temp,set_cmd=set_temp)
station.add_component(mxc_temp)

def wait_temp(setpoint):

    while np.abs(mxc_temp()-setpoint) > 5e-3:
        time.sleep(10)


# %% Spectrum analyzer

from qcodes_contrib_drivers.drivers.ThinkRF import R5550

sa = R5550.R5550(name='sa',address ='10.21.41.160')
station.add_component(sa)
# %%% Acquire data


sa.start(5e9-5e3)
sa.stop(5e9+5e3)
sa.rbw(1e2)
sa.vbw(1e0)
sa.averages(1000)
meas = Measurement(experiment, station)
meas.name = "SA"
meas.register_parameter(sa.spectrum,paramtype='array')


with meas.run() as datasaver:
    sa.refresh_faxis()
    datasaver.add_result((sa.spectrum,sa.spectrum()))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%% T sweep

sa.start(5e9-5e3)
sa.stop(5e9+5e3)
sa.rbw(1e2)
sa.vbw(1e0)
sa.averages(2500)
meas = Measurement(experiment, station)
meas.name = "SA_Tsweep"
meas.register_parameter(mxc_temp, paramtype='numeric')
meas.register_parameter(sa.spectrum, setpoints=(mxc_temp,), paramtype='array')
temps= np.linspace(0.8,0.15,14)

with meas.run() as datasaver:
    sa.refresh_faxis()
    for T in temps:
        mxc_temp(T)
        wait_temp(T)
        datasaver.add_result((sa.spectrum,sa.spectrum()),(mxc_temp,mxc_temp()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset=load_by_run_spec(captured_run_id= id)
[axlist, cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s'%(meas.name,dataset.captured_run_id),bbox_inches='tight')

# %%%% Analysis
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id= 13)
data =_get_data_from_ds(dataset)

spectrum = 10**(data[0][2]['data']/10)
temp = np.mean(data[0][0]['data'],axis=1)
spectrum = np.mean(spectrum,axis=1)



#Fit

def func(x, a, b):
    return a*x+b

popt, pcov = curve_fit(func, spectrum*1e10, temp)
plt.scatter(spectrum*1e10,temp)
plt.plot(spectrum*1e10,func(spectrum*1e10, *popt))

# %% Raw ADC

from qcodes.instrument_drivers.QM_qcodes.opx_raw_adc import *
importlib.reload(config)
readout_pulse_length = config.config['pulses']['readout_pulse']['length']

meas_freq = 4.7e9
lo_freq = 4.6e9
resonator_if = meas_freq-lo_freq

lo_qubit =6.15e9
qubit_freq = 6.33e9

if_qubit = qubit_freq - lo_qubit
# sg_src.frequency(lo_freq)
# sg_src.power(18)
# sg_src.status(1)
config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_freq
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_freq
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = resonator_if
config.config['elements']['resonator']['intermediate_frequency'] = resonator_if
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit


opx_raw_adc = OPXRawADC(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
opx_raw_adc.t_meas(readout_pulse_length*1e-9*1) #in 
opx_raw_adc.readout_pulse_length(readout_pulse_length)
opx_raw_adc.smearing(20)
meas = Measurement(experiment, station)
meas.name = "ADC_OPXcalibration"
meas.register_parameter(opx_raw_adc.trace_iq,paramtype='array')

# opx_raw_adc.sim_exp(10000)
with meas.run() as datasaver:
    opx_raw_adc.run_exp()
    get_v = opx_raw_adc.trace_iq.get()
    datasaver.add_result((opx_raw_adc.trace_iq, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
opx_raw_adc.close()

# %%% Plot dataset RAW ADC

# %%%% IQ imbalance

def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1 - g ** 2) * (2 * c ** 2 - 1))
    return [float(N * x) for x in [(1 - g) * c, (1 + g) * s, (1 - g) * s, (1 + g) * c]]
# %%%% Filter

from scipy import signal
from scipy import fft

sos = signal.butter(10, 100e6, 'lowpass', fs=1e9, output='sos')

# %%%% Analysis
dataset =load_by_run_spec(captured_run_id=43)
data =_get_data_from_ds(dataset)
time = data[0][0]['data'].flatten()
I = data[0][1]['data'].flatten()-np.mean(data[0][1]['data'].flatten())
Q = data[1][1]['data'].flatten()-np.mean(data[1][1]['data'].flatten())
#Correcting for the imbalance of the downconversion mixer
# phase = 0 * np.pi/180
# g=-0.00
# c = IQ_imbalance(g,phase)
# c = np.reshape(c, (2,2))
# c2 = np.linalg.inv(np.array([[1,-1],[1,1]]))
# [Inew,Qnew] = np.dot(c2,np.dot(c,[I,Q]))

Inew = signal.sosfilt(sos, I*np.cos(2*np.pi*resonator_if*time*1e-9))
Qnew = signal.sosfilt(sos, Q*np.sin(2*np.pi*resonator_if*time*1e-9))
Pnew = signal.sosfilt(sos, I**2+Q**2)
plt.figure()

plt.plot(time,Pnew)
# plt.plot(time,Q)
# plt.plot(time,np.abs(Inew-1j*Qnew))
# %% Resonator scan 1D
from qcodes.parameters import ElapsedTimeParameter
from qcodes.instrument_drivers.QM_qcodes.opx_resonator_scan import *

def Resonator_scan(LO,IF_start,IF_stop,npts,n_avg,plot=False):
    
    importlib.reload(config)

    lo_freq = LO
    config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_freq
    config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_freq
    opx_freq_scan = OPXSpectrumScan(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
    station.add_component(opx_freq_scan)
    
    readout_pulse_length = config.config['pulses']['readout_pulse']['length']
    fmin= IF_start
    fmax= IF_stop
    npts = npts
    
    opx_freq_scan.f_start(fmin)
    opx_freq_scan.f_stop(fmax)
    opx_freq_scan.n_points(npts)
    opx_freq_scan.t_meas(readout_pulse_length*1e-9*n_avg) #in s
    opx_freq_scan.readout_pulse_length(readout_pulse_length)
    opx_freq_scan.amp(1)
    opx_freq_scan.set_octave()
    # opx_freq_scan.sim_exp(10000)
    
    meas = Measurement(experiment, station)
    meas.name = "Freq_scan_1D"
    meas.register_parameter(opx_freq_scan.trace_mag_phase,paramtype='array')
    with meas.run() as datasaver:
            opx_freq_scan.run_exp()
            get_v = opx_freq_scan.trace_mag_phase.get()
            datasaver.add_result((opx_freq_scan.trace_mag_phase, get_v))
            datasaver.flush_data_to_database()
            id = datasaver.run_id
    dataset =load_by_run_spec(captured_run_id=id)
    
    if plot:        
        [axlist,cbar] = plot_dataset(dataset)
        axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
        axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
    
    station.remove_component('OPXSpectrumScan')
    opx_freq_scan.close_all()
    
    data =_get_data_from_ds(dataset)
    freq = data[0][0]['data'].flatten()
    mag = data[0][1]['data'].flatten()
    phase = data[1][1]['data'].flatten()
    
    return [freq,mag,phase]

# %%% Execute


res_data= Resonator_scan(4.7e9, 95e6, 145e6, 101, 5e3, True)

[freq,mag,phase] = res_data
phase = phase * np.pi/180
mag = 10**(mag/20)
port1 = circuit.notch_port() 
port1.add_data(freq,mag*np.exp(1j*phase))
port1.autofit()



# %% Analysis

dataset=load_by_run_spec(captured_run_id= 17)
data =_get_data_from_ds(dataset)

freq = data[0][0]['data'].flatten()
mag = data[0][1]['data'].flatten()
phase = data[1][1]['data'].flatten()


phase = phase * np.pi/180
mag = 10**(mag/20)
port1 = circuit.notch_port() 
port1.add_data(freq,mag*np.exp(1j*phase))
port1.autofit()
port1.fitresults['fr']
port1.fitresults['Ql']

# %% Resonator scan 2d

from qcodes.instrument_drivers.QM_qcodes.opx_resonator_scan import *
importlib.reload(config)

lo_freq = 4.75e9
# sg_src.frequency(lo_freq)
# sg_src.power(18)
# sg_src.status(1)
config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_freq
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_freq
# print(config.config['mixers']['octave_octave1_1'])
opx_freq_scan = OPXSpectrumScan(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_freq_scan)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
fmin=100e6
fmax= 200e6
npts =201

opx_freq_scan.f_start(fmin)
opx_freq_scan.f_stop(fmax)
opx_freq_scan.n_points(npts)
opx_freq_scan.t_meas(readout_pulse_length*1e-9*5e3) #in s
opx_freq_scan.readout_pulse_length(readout_pulse_length)
# opx_freq_scan.sim_exp(10000)
gains = np.arange(-10,20,1)
meas = Measurement(experiment, station)
meas.name = "Freq_power_scan_2D"
# Gain = qc.Parameter('Gain', label='Gain', unit='dB', set_cmd=None, get_cmd=None)
meas.register_parameter(opx_freq_scan.Gain,paramtype='numeric')
meas.register_parameter(opx_freq_scan.trace_mag_phase,setpoints=(opx_freq_scan.Gain,),paramtype='array')
opx_freq_scan.set_octave()
with meas.run() as datasaver:
    for gain in gains:
        opx_freq_scan.Gain(gain)
        opx_freq_scan.run_exp()
        get_v = opx_freq_scan.trace_mag_phase.get()
        datasaver.add_result((opx_freq_scan.trace_mag_phase, get_v),(opx_freq_scan.Gain,opx_freq_scan.Gain.get()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXSpectrumScan')
opx_freq_scan.close()
  
# %% Qubit scan

from qcodes.instrument_drivers.QM_qcodes.opx_qubit_scan import *

def Qubit_scan(LO,IF_start,IF_stop,df,n_avg,lo_r,if_r,plot = False):
        
    importlib.reload(config)
    
    
    lo_qubit = LO
    resonator_if = if_r
    lo_freq = lo_r
    
    config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_freq
    config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_freq
    config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = resonator_if
    config.config['elements']['resonator']['intermediate_frequency'] = resonator_if
    config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
    config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
    opx_qubit_scan = OPXQubitScan(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
    station.add_component(opx_qubit_scan)
    
    readout_pulse_length = config.config['pulses']['readout_pulse']['length']
    fmin = IF_start
    fmax= IF_stop
    df = df
    freqs = np.arange(fmin,fmax+0.1,df)
    
    npts = len(freqs)
    
    opx_qubit_scan.f_start(fmin)
    opx_qubit_scan.f_stop(fmax)
    opx_qubit_scan.n_points(npts)
    opx_qubit_scan.t_meas(readout_pulse_length*1e-9*n_avg) #in s
    opx_qubit_scan.readout_pulse_length(readout_pulse_length)
    #opx_qubit_scan.amp_resonator(1)
    # opx_qubit_scan.set_qubit_gain(7)
    #opx_qubit_scan.sim_exp(10000)
    
    opx_qubit_scan.set_octave()
    
    meas = Measurement(experiment, station)
    meas.name = "Qubit_scan_1D"
    meas.register_parameter(opx_qubit_scan.trace_mag_phase,paramtype='array')
    with meas.run() as datasaver:
        opx_qubit_scan.run_exp()
        get_v = opx_qubit_scan.trace_mag_phase.get()
        datasaver.add_result((opx_qubit_scan.trace_mag_phase, get_v))
        datasaver.flush_data_to_database()
        id = datasaver.run_id
    dataset =load_by_run_spec(captured_run_id=id)
    
    if plot:
        [axlist,cbar] = plot_dataset(dataset)
        axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
        axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
    station.remove_component('OPXQubitScan')
    opx_qubit_scan.close()
    
    data =_get_data_from_ds(dataset)
    freq = data[0][0]['data'].flatten()
    mag = data[0][1]['data'].flatten()
    phase = data[1][1]['data'].flatten()
    
    return [freq,mag,phase]
    
# %%% Execute


res_data= Resonator_scan(4.7e9, 95e6, 145e6, 101, 5e3, True)
[freq,mag,phase] = res_data
phase = phase * np.pi/180
mag = 10**(mag/20)
port1 = circuit.notch_port() 
port1.add_data(freq,mag*np.exp(1j*phase))
port1.autofit()
fr = port1.fitresults['fr']
lo = 4.7e9
Qubit_scan(2.4e9,-350e6,-10e6,1e6,1e4,lo,fr-lo-2e6,plot = True)

# %% Anharmonicity

importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_anharmonicity import *


meas_freq = 4.822e9
qubit_freq = 2.46e9
if_qubit = -100e6
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
lo_qubit = qubit_freq - if_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit2']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['elements']['qubit2']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][1]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
config.config['mixers']['octave_octave1_2'][1]['intermediate_frequency'] = if_qubit

opx_anharmonicity = OPXAnharmonicity(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_anharmonicity)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
fmin =-350e6
fmax= -10e6
df = 1e6
freqs = np.arange(fmin,fmax+0.1,df)

npts = len(freqs)

opx_anharmonicity.f_start(fmin)
opx_anharmonicity.f_stop(fmax)
opx_anharmonicity.n_points(npts)
opx_anharmonicity.t_meas(readout_pulse_length*1e-9*4e4) #in s
opx_anharmonicity.readout_pulse_length(readout_pulse_length)
opx_anharmonicity.amp_resonator(1)

#opx_anharmonicity.sim_exp(10000)

opx_anharmonicity.set_octave()

meas = Measurement(experiment, station)
meas.name = "Anharmonicity"
meas.register_parameter(opx_anharmonicity.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_anharmonicity.run_exp()
    get_v = opx_anharmonicity.trace_mag_phase.get()
    datasaver.add_result((opx_anharmonicity.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXAnharmonicity')
opx_anharmonicity.close()


# %% Qubit scan 2d
importlib.reload(config)
from qcodes.instrument_drivers.OPX.opx_qubit_scan import *

meas_freq = 4.7044e9
lo_qubit =3.8e9
lo_freq = 4.55e9
resonator_if = meas_freq-lo_freq
# sg_src_2.frequency(lo_qubit)
# sg_src_2.power(-5)
# sg_src_2.IQ_state(1)
# sg_src_2.status(1)
# sg_src.frequency(lo_freq)
# sg_src.power(18)
# sg_src.status(1)

config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_freq
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_freq
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = resonator_if
config.config['elements']['resonator']['intermediate_frequency'] = resonator_if
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
# print(config.config['mixers']['octave_octave1_1'])
opx_qubit_scan = OPXQubitScan(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_qubit_scan)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
fmax = -10e6
fmin= -350e6
df = 1e6
freqs = np.arange(fmin,fmax+0.1,df)

npts = len(freqs)

opx_qubit_scan.f_start(fmin)
opx_qubit_scan.f_stop(fmax)
opx_qubit_scan.n_points(npts)
opx_qubit_scan.t_meas(readout_pulse_length*1e-9*1e4) #in s
opx_qubit_scan.readout_pulse_length(readout_pulse_length)
opx_qubit_scan.set_octave()
# opx_qubit_scan.sim_exp(10000)
gainstart = -10
gainstop = 10
gainstep = 1
powers = np.arange(gainstart,gainstop,gainstep )
meas = Measurement(experiment, station)
meas.name = "Qubit_scan_2D_power_continuouis"
gain = qc.Parameter('Gain', label='Gain(dB)', unit='dB', set_cmd=None, get_cmd=None)
meas.register_parameter(gain,paramtype='numeric')
meas.register_parameter(opx_qubit_scan.trace_mag_phase,setpoints=(gain,),paramtype='array')
with meas.run() as datasaver:
    for power in powers: 
        gain(power)
        opx_qubit_scan.set_qubit_gain(power)
        opx_qubit_scan.run_exp()
        get_v = opx_qubit_scan.trace_mag_phase.get()
        datasaver.add_result((opx_qubit_scan.trace_mag_phase, get_v),(gain,gain()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
opx_qubit_scan.close()
# %% Time Rabi
importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_rabi_time import *

meas_freq = 4.8245e9
qubit_freq = 2.876e9
if_qubit = -250e6
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
lo_qubit = qubit_freq - if_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit

opx_rabi_time = OPXRabiTime(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_rabi_time)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 175
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_rabi_time.t_start(tmin)
opx_rabi_time.t_stop(tmax)
opx_rabi_time.n_points(npts)
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.amp_resonator(1)
opx_rabi_time.amp_qubit(1)
opx_rabi_time.freq_qubit(if_qubit)
opx_rabi_time.wait_time(5000)
#opx_rabi_time.sim_exp(10000)
opx_rabi_time.set_octave()



meas = Measurement(experiment, station)
meas.name = "Rabi_Time_1D"
meas.register_parameter(opx_rabi_time.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_rabi_time.run_exp()
    get_v = opx_rabi_time.trace_mag_phase.get()
    datasaver.add_result((opx_rabi_time.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXRabiTime')
opx_rabi_time.close()

# %%% Analysis

dataset =load_by_run_spec(captured_run_id=15)
data =_get_data_from_ds(dataset)
time = data[0][0]['data'].flatten()
R = data[0][1]['data'].flatten()-np.mean(data[0][1]['data'].flatten())
phase = data[1][1]['data'].flatten()-np.mean(data[1][1]['data'].flatten())

angle =45 *np.pi/180

newdata = np.exp(1j*angle)*(I+1j*Q)

plt.figure()
plt.plot(np.real(newdata))
plt.plot(np.imag(newdata))


# %% Time Rabi - Repetiton rate
importlib.reload(config)
from qcodes.instrument_drivers.OPX.opx_rabi_time import *


qubit_freq = 3.507e9

meas_freq = 5.7613e9
lo_qubit =3.65e9
lo_resonator = 5.65e9

if_qubit = qubit_freq - lo_qubit
if_resonator =meas_freq - lo_resonator

# sg_src_2.frequency(lo_qubit)
# sg_src_2.power(-3)
# sg_src_2.IQ_state(1)
# sg_src_2.status(1)
# sg_src.frequency(lo_freq)
# sg_src.power(18)
# sg_src.status(1)

config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
# print(config.config['mixers']['octave_octave1_1'])
opx_rabi_time = OPXRabiTime(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
# station.add_component(opx_rabi_time)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 1500
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_rabi_time.t_start(tmin)
opx_rabi_time.t_stop(tmax)
opx_rabi_time.n_points(npts)
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.amp_resonator(1)
opx_rabi_time.amp_qubit(1)
opx_rabi_time.freq_qubit(if_qubit)
# opx_rabi_time.sim_exp(10000)

opx_rabi_time.set_octave()

waittimes = [int(1e3),int(1e4),int(1e5),int(1e6)] #wait time in ns
meas = Measurement(experiment, station)
meas.name = "Rabi_Time_2D_reprate"
meas.register_parameter(opx_rabi_time.wait_time,paramtype='numeric')
meas.register_parameter(opx_rabi_time.trace_mag_phase,setpoints=(opx_rabi_time.wait_time,), paramtype='array')
with meas.run() as datasaver:
    for time in waittimes:
        opx_rabi_time.wait_time(time)
        opx_rabi_time.run_exp()
        get_v = opx_rabi_time.trace_mag_phase.get()
        datasaver.add_result((opx_rabi_time.trace_mag_phase, get_v),(opx_rabi_time.wait_time,opx_rabi_time.wait_time()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
opx_rabi_time.close()

  # %% Time Rabi - 2D
importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_rabi_time import *

meas_freq = 4.84956e9
qubit_freq = 3.65593e9
lo_qubit =3.8e9
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
if_qubit = qubit_freq - lo_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
# print(config.config['mixers']['octave_octave1_1'])
opx_rabi_time = OPXRabiTime(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_rabi_time)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 200
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_rabi_time.t_start(tmin)
opx_rabi_time.t_stop(tmax)
opx_rabi_time.n_points(npts)
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.freq_qubit(if_qubit)
opx_rabi_time.wait_time(10000)
opx_rabi_time.set_octave()
amps = np.arange(0.2,1.9,0.05)
#opx_rabi_time.sim_exp(10000)
meas = Measurement(experiment, station)
meas.name = "Rabi_Time_2D"
meas.register_parameter(opx_rabi_time.amp_qubit,paramtype='numeric')
meas.register_parameter(opx_rabi_time.trace_mag_phase,setpoints=(opx_rabi_time.amp_qubit,),paramtype='array')
with meas.run() as datasaver:
    for amp in amps:
        opx_rabi_time.amp_qubit(amp)
        opx_rabi_time.run_exp()
        get_v = opx_rabi_time.trace_mag_phase.get()
        datasaver.add_result((opx_rabi_time.amp_qubit,opx_rabi_time.amp_qubit()),(opx_rabi_time.trace_mag_phase, get_v))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXRabiTime')
opx_rabi_time.close()
# %% T1
importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_T1 import *

meas_freq = 4.8245e9
qubit_freq = 2.876e9
if_qubit = -250e6
if_resonator = 125e6



lo_resonator = meas_freq - if_resonator
lo_qubit = qubit_freq - if_qubit

config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
opx_T1 = OPXT1(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_T1)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 75
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_T1.t_start(tmin)
opx_T1.t_stop(tmax)
opx_T1.n_points(npts)
opx_T1.readout_pulse_length(readout_pulse_length)
opx_T1.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_T1.readout_pulse_length(readout_pulse_length)
opx_T1.set_octave()
#opx_T1.sim_exp(10000)
meas = Measurement(experiment, station)
meas.name = "T1"
meas.register_parameter(opx_T1.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_T1.run_exp()
    get_v = opx_T1.trace_mag_phase.get()
    datasaver.add_result((opx_T1.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXT1')
opx_T1.close()

# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id= 8)
data =_get_data_from_ds(dataset)

time = data[0][0]['data'].flatten()
mag = data[0][1]['data'].flatten()
phase = data[1][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, c):
    return a * np.exp(-b * x) + c
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[4,0.1,29.5])
plt.plot(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{1}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')
# %% T2

importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_T2 import *

meas_freq = 4.8245e9
qubit_freq = 2.876e9
if_qubit = -250e6
if_resonator = 125e6


lo_resonator = meas_freq - if_resonator
lo_qubit = qubit_freq - if_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
opx_T2 = OPXT2(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_T2)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 3
tmax= 200
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_T2.t_start(tmin)
opx_T2.t_stop(tmax)
opx_T2.n_points(npts)
opx_T2.readout_pulse_length(readout_pulse_length)
opx_T2.t_meas(readout_pulse_length*1e-9*2e4) #in s
opx_T2.readout_pulse_length(readout_pulse_length)
opx_T2.freq_qubit(if_qubit)
opx_T2.amp_qubit(0.15)
opx_T2.detuning(20e6)
#opx_T2.sim_exp(10000)
opx_T2.set_octave()

meas = Measurement(experiment, station)
meas.name = "T2"
meas.register_parameter(opx_T2.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_T2.run_exp()
    get_v = opx_T2.trace_mag_phase.get()
    datasaver.add_result((opx_T2.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXT2')
opx_T2.close()

# %%% Fit

from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id= 17)
data =_get_data_from_ds(dataset)

time = data[0][0]['data'].flatten()
mag = data[0][1]['data'].flatten()
phase = data[1][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, c, f, phi):
    return a * np.exp(-b * x)*np.cos(2*np.pi*f*x+phi) + c
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[3,0.1,63,0.02,0])
plt.plot(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{2}^{*}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')
#%% Echo

importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_echo import *

meas_freq = 4.8245e9
qubit_freq = 2.876e9
if_qubit = -250e6
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
lo_qubit = qubit_freq - if_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
opx_echo = OPXEcho(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_echo)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 180
dt = 2
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_echo.t_start(tmin)
opx_echo.t_stop(tmax)
opx_echo.n_points(npts)
opx_echo.readout_pulse_length(readout_pulse_length)
opx_echo.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_echo.readout_pulse_length(readout_pulse_length)
opx_echo.freq_qubit(if_qubit)
opx_echo.amp_qubit(0.15)
#opx_echo.sim_exp(10000)
opx_echo.set_octave()

meas = Measurement(experiment, station)
meas.name = "Echo"
meas.register_parameter(opx_echo.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_echo.run_exp()
    get_v = opx_echo.trace_mag_phase.get()
    datasaver.add_result((opx_echo.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXEcho')
opx_echo.close()

# %%% Fit
from scipy.optimize import curve_fit
dataset=load_by_run_spec(captured_run_id= 16)
data =_get_data_from_ds(dataset)

time = data[0][0]['data'].flatten()
mag = data[0][1]['data'].flatten()
phase = data[1][1]['data'].flatten()

xdata = time
ydata = phase
def func(x, a, b, c):
    return a * (np.exp(-b * x)) + c
plt.figure()
popt, pcov = curve_fit(func, xdata, ydata, p0=[1.5,0.05,31.4])
plt.plot(time,phase)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='$T_{echo}$ ~ %d ns' % (1/popt[1]))
plt.legend()
plt.xlabel('$t_{wait}$ (ns)')
plt.ylabel('$\phi$ (degree)')
#%% Ramsey freq-time

importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_T2 import *



meas_freq = 4.7573e9
qubit_freq = 6.11516e9
lo_qubit =5.9e9
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
if_qubit = qubit_freq - lo_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
opx_T2 = OPXT2(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_T2)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 125
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_T2.t_start(tmin)
opx_T2.t_stop(tmax)
opx_T2.n_points(npts)
opx_T2.readout_pulse_length(readout_pulse_length)
opx_T2.t_meas(readout_pulse_length*1e-9*2e4) #in s
opx_T2.readout_pulse_length(readout_pulse_length)
opx_T2.freq_qubit(if_qubit)
# opx_T2.sim_exp(10000)
opx_T2.set_octave()
detunings = np.arange(-20e6,20e6+0.1e6,1e6)

meas = Measurement(experiment, station)
meas.name = "Ramsey_2D"
meas.register_parameter(opx_T2.detuning,paramtype='numeric')
meas.register_parameter(opx_T2.trace_mag_phase,setpoints=(opx_T2.detuning,),paramtype='array')
with meas.run() as datasaver:
    for detuning in detunings:
        opx_T2.detuning(detuning)
        opx_T2.run_exp()
        get_v = opx_T2.trace_mag_phase.get()
        datasaver.add_result((opx_T2.detuning,opx_T2.detuning()),(opx_T2.trace_mag_phase, get_v))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.add_component('OPXT2')
opx_T2.close()
# %% Rabi -Chevron


importlib.reload(config)
from qcodes.instrument_drivers.QM_qcodes.opx_rabi_time import *

meas_freq = 4.8482e9
qubit_freq = 3.65593e9
lo_qubit =3.8e9
if_resonator = 125e6

lo_resonator = meas_freq - if_resonator
if_qubit = qubit_freq - lo_qubit


config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['octave_octave1_1'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['octave_octave1_1'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['intermediate_frequency'] = if_qubit
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['lo_frequency'] = lo_qubit
config.config['mixers']['octave_octave1_2'][0]['intermediate_frequency'] = if_qubit
# print(config.config['mixers']['octave_octave1_1'])
opx_rabi_time = OPXRabiTime(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)
station.add_component(opx_rabi_time)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
tmin = 4
tmax= 200
dt = 1
taus = np.arange(tmin,tmax+0.1,dt)

npts = len(taus)

opx_rabi_time.t_start(tmin)
opx_rabi_time.t_stop(tmax)
opx_rabi_time.n_points(npts)
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.t_meas(readout_pulse_length*1e-9*5e4) #in s
opx_rabi_time.readout_pulse_length(readout_pulse_length)
opx_rabi_time.freq_qubit(if_qubit)
opx_rabi_time.wait_time(10000)
opx_rabi_time.set_octave()
freqs = np.arange(if_qubit-20e6,if_qubit+20e6+0.1e6,1e6)
#opx_rabi_time.sim_exp(10000)
meas = Measurement(experiment, station)
meas.name = "Rabi_Time_2D"
meas.register_parameter(opx_rabi_time.freq_qubit,paramtype='numeric')
meas.register_parameter(opx_rabi_time.trace_mag_phase,setpoints=(opx_rabi_time.freq_qubit,),paramtype='array')
with meas.run() as datasaver:
    for freq in freqs:
        opx_rabi_time.freq_qubit(freq)
        opx_rabi_time.run_exp()
        get_v = opx_rabi_time.trace_mag_phase.get()
        datasaver.add_result((opx_rabi_time.freq_qubit,opx_rabi_time.freq_qubit()),(opx_rabi_time.trace_mag_phase, get_v))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
station.remove_component('OPXRabiTime')
opx_rabi_time.close()
#%% Qubit - Meas Power


importlib.reload(config)
from qcodes.instrument_drivers.OPX.opx_qubit_measpower import *

if_qubit = 255e6
if_resonator = 104.5e6
lo_qubit = 3e9
lo_resonator = 4.95e9

sg_src.frequency(lo_resonator)
vsg.rfoutput1.frequency(lo_qubit)
vsg.rfoutput1.level(-5)
vsg.rfoutput1.state('ON')

config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['mixer_resonator'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['mixer_resonator'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['mixer_qubit'][0]['lo_frequency'] = lo_qubit

opx_qubit_measpower = OPXMeasPower(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)

readout_pulse_length = config.config['pulses']['readout_pulse']['length']
amin = 0.2
amax= 1.4
da =0.2
amps = np.arange(amin,amax+0.01,da)
fmin = 50e6
fmax = 150e6
df = 1e6
freqs = np.arange(fmin,fmax+0.1,df)
opx_qubit_measpower.df(df)
opx_qubit_measpower.damp(da)
opx_qubit_measpower.amp_start(amin)
opx_qubit_measpower.amp_stop(amax)
opx_qubit_measpower.f_start(fmin)
opx_qubit_measpower.f_stop(fmax)
opx_qubit_measpower.n_points_f(len(freqs))
opx_qubit_measpower.n_points_amp(len(amps))
opx_qubit_measpower.readout_pulse_length(readout_pulse_length)
opx_qubit_measpower.n_avg(1e4) #in s


#opx_qubit_measpower.sim_exp(10000)

meas = Measurement(experiment, station)
meas.name = "Qubit_Meas_power"
meas.register_parameter(opx_qubit_measpower.trace_mag_phase,paramtype='array')
with meas.run() as datasaver:
    opx_qubit_measpower.run_exp()
    get_v = opx_qubit_measpower.trace_mag_phase.get()
    datasaver.add_result((opx_qubit_measpower.trace_mag_phase, get_v))
    datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
opx_qubit_measpower.close()

#%% Rabi-Chevron Amp



importlib.reload(config)
from qcodes.instrument_drivers.OPX.opx_rabi_chevron_amp import *

if_qubit = 255e6
if_resonator = 104.3e6
lo_qubit = 3e9
lo_resonator = 4.95e9

sg_src.frequency(lo_resonator)
vsg.rfoutput1.frequency(lo_qubit)
vsg.rfoutput1.level(-2)
vsg.rfoutput1.state('ON')

config.config['elements']['resonator']['mixInputs']['lo_frequency'] = lo_resonator
config.config['elements']['resonator']['intermediate_frequency'] = if_resonator
config.config['mixers']['mixer_resonator'][0]['lo_frequency'] = lo_resonator
config.config['mixers']['mixer_resonator'][0]['intermediate_frequency'] = if_resonator
config.config['elements']['qubit']['mixInputs']['lo_frequency'] = lo_qubit
config.config['mixers']['mixer_qubit'][0]['lo_frequency'] = lo_qubit

opx_rabi_chevron_amp = OPXRabiChevron_Amp(config=config.config, host=config.opx_ip, port=config.opx_port,host_octave=config.octave_ip,port_octave=config.octave_port)



readout_pulse_length = config.config['pulses']['readout_pulse']['length']
amin = 1
amax= 3
da =0.1
amps = np.arange(amin,amax+0.01,da)
fmin = 50e6
fmax = 150e6
df = 1e6
freqs = np.arange(fmin,fmax+0.1,df)
opx_rabi_chevron_amp.df(df)
opx_rabi_chevron_amp.damp(da)
opx_rabi_chevron_amp.amp_start(amin)
opx_rabi_chevron_amp.amp_stop(amax)
opx_rabi_chevron_amp.f_start(fmin)
opx_rabi_chevron_amp.f_stop(fmax)
opx_rabi_chevron_amp.n_points_f(len(freqs))
opx_rabi_chevron_amp.n_points_amp(len(amps))
opx_rabi_chevron_amp.readout_pulse_length(readout_pulse_length)
opx_rabi_chevron_amp.n_avg(1e4) #in s
opx_rabi_chevron_amp.duration(50)


opx_rabi_chevron_amp.sim_exp(10000)
durations = np.arange(10,1000+0.1,10)
meas = Measurement(experiment, station)
meas.name = "Rabi_Chevron_Amp"
meas.register_parameter(opx_rabi_chevron_amp.duration,paramtype='numeric')
meas.register_parameter(opx_rabi_chevron_amp.trace_mag_phase,setpoints=(opx_rabi_chevron_amp.duration,),paramtype='array')
with meas.run() as datasaver:
    for duration in durations:
        opx_rabi_chevron_amp.duration(duration)
        opx_rabi_chevron_amp.run_exp()
        get_v = opx_rabi_chevron_amp.trace_mag_phase.get()
        datasaver.add_result((opx_rabi_chevron_amp.trace_mag_phase, get_v),(opx_rabi_chevron_amp.duration,opx_rabi_chevron_amp.duration.get()))
        datasaver.flush_data_to_database()
    id = datasaver.run_id
dataset =load_by_run_spec(captured_run_id=id)
[axlist,cbar] = plot_dataset(dataset)
axlist[0].figure.savefig(path+'\\#%s_%s_amp'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
axlist[1].figure.savefig(path+'\\#%s_%s_phase'%(meas.name,dataset.captured_run_id),bbox_inches='tight')
opx_rabi_chevron_amp.close()
