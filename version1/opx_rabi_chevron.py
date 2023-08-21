# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:28:23 2023

@author: Nanoelectronics.ad
"""

from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.OPX.opx_driver import *
from qm.qua import *
from scipy import signal
from qualang_tools.units import unit
# noinspection PyAbstractClass

class IQ2DArray(MultiParameter):
    def __init__(self,
        name: str,
        instrument: "OPXRabiChevron",
        start_f: int,
        stop_f: int,
        npts_f: int,
        start_t: int,
        stop_t : int,
        npts_t: int,):
        # names, labels, and units are the same
        super().__init__('MAg_phase', names=('Mag', 'Phase'), shapes=((npts_t,npts_f), (npts_t,npts_f)),
                         labels=('Mag', 'Phase'),
                         units=('dB', 'degree'),
                         instrument=instrument,
                         # note that EACH item needs a sequence of setpoint arrays
                         # so a 1D item has its setpoints wrapped in a length-1 tuple
                         setpoints = ((instrument.time_axis(),instrument.freq_axis()), (instrument.time_axis(),instrument.freq_axis())),
                         setpoint_units = (("ns",'Hz'), ("ns",'Hz')),
                         setpoint_labels = (
                                               ("t","f",),
                                               ("t","f",),
                                           ),
                         setpoint_names = (
                                            ("Time",'Frequency'),
                                            ("Time",'Frequency'),
                                                            ),

                         #docstring='param that returns two single values, I and Q')
                         )
        self.set_sweep(start_f, stop_f, npts_f, start_t, stop_t, npts_t)
#

    def set_sweep(self, start_f: int, stop_f: int, npts_f: int, start_t: int, stop_t: int, npts_t: int):
        # Needed to update config of the software parameter on sweep change
        # frequency setpoints tuple as needs to be hashable for look up.
        f = tuple(np.linspace(start_f, stop_f, num=npts_f))
        t = tuple(np.linspace(4*start_t, 4*stop_t, num=npts_t))
        self.setpoints = ((t,f), (t,f))
        self.shapes = ((npts_t,npts_f,), (npts_t,npts_f,))


    def get_raw(self):
        R, phase = self.instrument.get_res()
        return R, phase


class OPXRabiChevron(OPX):
    def __init__(self, config: Dict, name: str = "OPXRabiChevron", host=None, port=None,host_octave=None, port_octave=None, **kwargs):
        super().__init__(config, name, host=host, port=port,host_octave=host_octave, port_octave=port_octave, **kwargs)
        self.add_parameter(
            "f_start",
            initial_value=30e6,
            unit="Hz",
            label="f start",
            vals=Numbers(-400e6, 400e6),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "f_stop",
            initial_value=70e6,
            unit="Hz",
            label="f stop",
            vals=Numbers(-400e6, 400e6),
            get_cmd=None,
            set_cmd=None,
        )
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
            "dt",
            initial_value=1,
            unit="",
            vals=Numbers(1, 1e6),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "df",
            initial_value=5e5,
            unit="",
            vals=Numbers(1, 1e9),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "n_avg",
            unit="",
            initial_value=1,
            vals=Numbers(1, 1e5),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "n_points_f",
            initial_value=100,
            unit="",
            vals=Numbers(1, 1e9),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "amp_qubit",
            unit="V",
            initial_value=1,
            vals=Numbers(-2, 2),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "n_points_t",
            initial_value=100,
            unit="",
            vals=Numbers(1, 1e9),
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
            "freq_axis",
            unit="Hz",
            label="Freq Axis",
            parameter_class=GeneratedSetPoints,
            startparam=self.f_start,
            stopparam=self.f_stop,
            numpointsparam=self.n_points_f,
            vals=Arrays(shape=(self.n_points_f,)),
        )
        self.add_parameter(
            "time_axis",
            unit="ns",
            label="Time Axis",
            parameter_class=GeneratedSetPoints,
            startparam=self.t_start,
            stopparam=self.t_stop,
            numpointsparam=self.n_points_t,
            vals=Arrays(shape=(self.n_points_t,)),
        )
        self.add_parameter(
            "trace_mag_phase",
            start_f=self.f_start(),
            stop_f=self.f_stop(),
            npts_f=self.n_points_f(),
            stop_t = self.t_stop(),
            start_t = self.t_start(),
            npts_t = self.n_points_t(),
            parameter_class=IQ2DArray,
        )

    def get_prog(self):
     
        if self.octave:                      
            lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
            lo_qubit = self.config['mixers']['octave_octave1_2'][0]['lo_frequency']

        else:
            lo =self.config['mixers']['mixer_resonator'][0]['lo_frequency']
            lo_qubit = self.config['mixers']['mixer_qubit'][0]['lo_frequency']
    
            
        self.parameters['trace_mag_phase'].set_sweep(lo_qubit+self.f_start(),lo_qubit+self.f_stop(),self.n_points_f(),self.t_start(),self.t_stop(),self.n_points_t())

        with program() as prog:
            n = declare(int)
            f = declare(int)
            t= declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            with for_(n, 0, n < self.n_avg(), n + 1):
             with for_(f, self.f_start(), f < self.f_stop()+self.df()/2.0, f + self.df()):
                 update_frequency('qubit',f)
                 with for_(t, self.t_start(), t<self.t_stop()+self.dt()/2.0, t+self.dt()):                    
                       reset_phase('qubit')
                       play("gauss"*amp(self.amp_qubit()), "qubit",duration=t) #t in clock cycles (4ns)
                       wait(t+17,'resonator')
                       reset_phase('resonator')
                       measure("readout", "resonator", None,
                               dual_demod.full("cos", "out1", "sin", "out2", I),
                               dual_demod.full("minus_sin", "out1", "cos", "out2", Q))
                    
                       wait(20000//4,'resonator','qubit')
                       save(I, I_st)
                       save(Q, Q_st)

            with stream_processing():
                I_st.buffer(self.n_points_t()).buffer(self.n_points_f()).average().save("I")
                Q_st.buffer(self.n_points_t()).buffer(self.n_points_f()).average().save("Q")

        return prog
    
    def set_octave(self):
        lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
        lo_qubit = self.config['mixers']['octave_octave1_2'][0]['lo_frequency']
        self.qm.octave.set_lo_source('resonator', OctaveLOSource.Internal)
        self.qm.octave.set_lo_frequency('resonator', lo)  # assign the LO inside the octave to element
        self.qm.octave.set_rf_output_gain('resonator', 0)  # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_mode('resonator', RFOutputMode.on) # set the behaviour of the RF switch
        self.qm.octave.calibrate_element('resonator',[(lo, self.config['elements']['resonator']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
        self.qm =self.qmm.open_qm(self.config)
        self.qm.octave.set_lo_source('qubit', OctaveLOSource.Internal)
        self.qm.octave.set_lo_frequency('qubit', lo_qubit)  # assign the LO inside the octave to element
        self.qm.octave.set_rf_output_gain('qubit', 0)  # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_mode('qubit', RFOutputMode.on) # set the behaviour of the RF switch
        self.qm.octave.calibrate_element('qubit',[(lo_qubit, self.config['elements']['qubit']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
        self.qm =self.qmm.open_qm(self.config)
        self.qm.octave.set_qua_element_octave_rf_in_port('resonator',"octave1", 1)
        self.qm.octave.set_downconversion('resonator',lo_source=RFInputLOSource.Internal)
        self.qm.octave.set_rf_output_gain('qubit', 10)  # can set gain from -10dB to 20dB
    
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
             R = 20*np.log10(np.sqrt(I ** 2 + Q ** 2)/self.config['waveforms']['readout_wf']['sample'])
             phase = np.angle(I + 1j * Q) * 180 / np.pi
             return R , phase
