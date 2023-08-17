# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:26:50 2023

@author: Nanoelectronics.ad
"""

from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.QM_qcodes.opx_driver import *
from qm.qua import *
from scipy import signal
from qualang_tools.units import unit
from qualang_tools.bakery import baking
from qualang_tools.bakery.bakery import deterministic_run
import qcodes.instrument_drivers.QM_qcodes.configuration as config
from scipy.signal.windows import gaussian,general_gaussian,cosine
import numpy as np
# noinspection PyAbstractClass

class IQArray(MultiParameter):
    def __init__(self,
        name: str,
        instrument: "OPXEcho",
        start: float,
        stop: float,
        npts: int,):
        # names, labels, and units are the same

        super().__init__('MAg_phase', names=('Mag', 'Phase'), shapes=((npts,), (npts,)),
                         labels=('Mag', 'Phase'),
                         units=('mV', 'degree'),
                         instrument=instrument,
                         # note that EACH item needs a sequence of setpoint arrays
                         # so a 1D item has its setpoints wrapped in a length-1 tuple
                         setpoints = ((instrument.time_axis(),), (instrument.time_axis(),)),
                         setpoint_units = (("ns",), ("ns",)),
                         setpoint_labels = (
                                               ("t",),
                                               ("t",),
                                           ),
                         setpoint_names = (
                                            ("Time",),
                                            ("Time",),
                                                            ),

                         #docstring='param that returns two single values, I and Q')
                         )
        self.set_sweep(start, stop, npts)
#

    def set_sweep(self, start: float, stop: float, npts: int):
        # Needed to update config of the software parameter on sweep change
        # frequency setpoints tuple as needs to be hashable for look up.
        f = tuple(np.linspace(2*int(start), 2*int(stop), num=npts))
        self.setpoints = ((f,), (f,))
        self.shapes = ((npts,), (npts,))

    def get_raw(self):
        R, phase = self.instrument.get_res()
        return R, phase



class OPXEcho(OPX):
    def __init__(self, config: Dict, name: str = "OPXEcho", host=None, port=None,host_octave=None, port_octave=None, **kwargs):
        super().__init__(config, name, host=host, port=port,host_octave=host_octave, port_octave=port_octave, **kwargs)
        self.add_parameter(
            "t_start",
            initial_value=3,
            unit="ns",
            label="t start",
            vals=Numbers(3, 1000),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "t_stop",
            initial_value=3,
            unit="ns",
            label="t stop",
            vals=Numbers(3, 100000),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "n_points",
            initial_value=100,
            unit="",
            vals=Numbers(1, 1e9),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "t_meas",
            unit="s",
            initial_value=0.01,
            vals=Numbers(0, 100),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "time_axis",
            unit="ns",
            label="Time Axis",
            parameter_class=GeneratedSetPoints,
            startparam=self.t_start,
            stopparam=self.t_stop,
            numpointsparam=self.n_points,
            vals=Arrays(shape=(self.n_points,)),
        )
        self.add_parameter(
            "amp_resonator",
            unit="",
            initial_value=1,
            vals=Numbers(-5, 5),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "amp_qubit",
            unit="V",
            initial_value=0.1,
            vals=Numbers(-0.5, 0.5),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "freq_qubit",
            unit="Hz",
            initial_value=100e6,
            vals=Numbers(-350e6, 350e6),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "detuning",
            unit="Hz",
            initial_value=10e6,
            vals=Numbers(-100e6, 100e6),
            get_cmd=None,
            set_cmd=None,
        )


        self.add_parameter(
            "readout_pulse_length",
            unit="ns",
            vals=Numbers(16, 1e7),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "trace_mag_phase",
            start=self.t_start(),
            stop=self.t_stop(),
            npts=self.n_points(),
            parameter_class=IQArray,
        )

    def get_prog(self):
        dt = (self.t_stop() - self.t_start()) / (self.n_points()-1)
        n_avg = round(self.t_meas() * 1e9 / self.readout_pulse_length())

        
        if self.octave:
            lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
            lo_qubit = self.config['mixers']['octave_octave1_2'][0]['lo_frequency']
           
        else:
            lo =self.config['mixers']['mixer_resonator'][0]['lo_frequency']
            lo_qubit = self.config['mixers']['mixer_qubit'][0]['lo_frequency']
            
        self.parameters['trace_mag_phase'].set_sweep(self.t_start(),self.t_stop(),self.n_points())
        
        # Create pulse sequence with baking
        
        t_min = self.t_start()
        t_max = self.t_stop()
        dt = int((t_max-t_min)/(self.n_points()-1))
        t_vec = np.arange(t_min, t_max + dt / 2, dt)
        
        sample_rate = 1e9
        n_samples = len(t_vec)

        baking_list = []  # Stores the baking objects
        t_pihalf = 8
        # Create the different baked sequences, corresponding to the different taus
        mixInput_sample_I = list(self.amp_qubit()*cosine(t_pihalf))
        mixInput_sample_Q = list(np.zeros(t_pihalf))
        mixInput_sample_2I = list(self.amp_qubit()*cosine(2*t_pihalf))
        mixInput_sample_2Q = list(np.zeros(2*t_pihalf))
         
        for i in range(n_samples):
            with baking(config.config, padding_method="left", sampling_rate=sample_rate) as b:
                

              
              # Assign waveforms to quantum element operation
              
                b.add_op("pihalf", "qubit", [mixInput_sample_I, mixInput_sample_Q], digital_marker = 'ON')
                b.add_op("pi", "qubit", [mixInput_sample_2Q, mixInput_sample_2I ], digital_marker = 'ON')
                #init_delay = 2*max_delay+int(3*t_pihalf)# Put initial delay to ensure that all of the pulses will have the same length
                #b.wait(init_delay, 'qubit')  # We first wait the entire duration.
        
                # We add the 2nd pi_half pulse with the phase 'dephasing' (Confusingly, the first pulse will be added later)
                b.wait(2*(t_max-(t_min+i*dt)),'qubit')
                b.play("pihalf", "qubit")
                b.wait(t_min+i*dt,'qubit')
                b.play('pi','qubit')
                b.wait(t_min+i*dt,'qubit')
                b.play('pihalf', 'qubit')
        
              #  b.play_at("pi", "qubit", t=init_delay - i -2*t_pihalf )
              #  b.play_at("pihalf", "qubit", t=init_delay - 2*i -3*t_pihalf )
        
            # Append the baking object in the list to call it from the QUA program
            baking_list.append(b)
            
        self.qm = self.qmm.open_qm(self.config, close_other_machines=True) #need to open once more to update the config file
        with program() as prog:
            n = declare(int)
            j = declare(int)
            t = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            update_frequency('qubit', self.freq_qubit())
            with for_(n, 0, n < n_avg, n + 1):
                with for_(j, 0, j < n_samples, j + 1):
                            
                    deterministic_run(baking_list, j, unsafe=True)
                   # The following wait command is used to align the resonator to happen right after the pulse.
                   # In this specific example, an align() command would have added a slight gap.
                    wait((2*t_max+3*t_pihalf-24)//4,'resonator')
                    reset_phase('resonator')
                    measure("readout", "resonator", None,
                            dual_demod.full("cos", "out1", "sin", "out2", I),
                            dual_demod.full("minus_sin", "out1", "cos", "out2", Q))

                    wait(1250, 'resonator','qubit')
                    save(I, I_st)
                    save(Q, Q_st)
                    


            with stream_processing():
                I_st.buffer(self.n_points()).average().save("I")
                Q_st.buffer(self.n_points()).average().save("Q")

        return prog
    
    def set_octave(self):
        lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
        lo_qubit = self.config['mixers']['octave_octave1_2'][0]['lo_frequency']
        self.qm.octave.set_lo_source('resonator', OctaveLOSource.Internal)
        self.qm.octave.set_lo_frequency('resonator', lo)  # assign the LO inside the octave to element
        self.qm.octave.set_rf_output_gain('resonator', 0)  # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_mode('resonator', RFOutputMode.trig_normal) # set the behaviour of the RF switch
        self.qm.octave.calibrate_element('resonator',[(lo, self.config['elements']['resonator']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
        self.qm =self.qmm.open_qm(self.config)
        self.qm.octave.set_lo_source('qubit', OctaveLOSource.Internal)
        self.qm.octave.set_lo_frequency('qubit', lo_qubit)  # assign the LO inside the octave to element
        self.qm.octave.set_rf_output_gain('qubit', 0)  # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_mode('qubit', RFOutputMode.trig_normal) # set the behaviour of the RF switch
        self.qm.octave.calibrate_element('qubit',[(lo_qubit, self.config['elements']['qubit']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
        self.qm =self.qmm.open_qm(self.config)
        self.qm.octave.set_qua_element_octave_rf_in_port('resonator',"octave1", 1)
        self.qm.octave.set_downconversion('resonator',lo_source=RFInputLOSource.Internal)
        self.qm.octave.set_rf_output_gain('qubit', 10)  # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_gain('resonator', -10)  # can set gain from -10dB to 20dB
    def run_exp(self):
        self.execute_prog(self.get_prog())
        
    def sim_exp(self,duration):
            self.simulate_prog(self.get_prog(),duration)
            
    def get_res(self):
        if (
            self.result_handles is None
        ):
            n = self.n_points()
            return {"I": (0,)*n, "Q": (0,)*n, "R": (0,)*n, "Phi": (0,)*n}
        else:
            self.result_handles.wait_for_all_values()
            u = unit()
            I = u.demod2volts(self.result_handles.get("I").fetch_all(), self.readout_pulse_length())
            Q = u.demod2volts(self.result_handles.get("Q").fetch_all(), self.readout_pulse_length())
            # R = np.sqrt(I ** 2 + Q ** 2)/(self.config['waveforms']['readout_wf']['sample']*self.amp_resonator())
            R = np.sqrt(I ** 2 + Q ** 2)
            phase = np.angle(I + 1j * Q) * 180 / np.pi
             
            return R , phase
