from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.QM_qcodes.opx_driver import *
from qm.qua import *
from scipy import signal
from qualang_tools.units import unit

# noinspection PyAbstractClass

class IQArray(MultiParameter):
    def __init__(self,
        name: str,
        instrument: "OPXSpectrumScan",
        start: float,
        stop: float,
        npts: int,):
        # names, labels, and units are the same

        super().__init__('MAg_phase', names=('Mag', 'Phase'), shapes=((npts,), (npts,)),
                         labels=('Mag', 'Phase'),
                         units=('dB', 'degree'),
                         instrument=instrument,
                         # note that EACH item needs a sequence of setpoint arrays
                         # so a 1D item has its setpoints wrapped in a length-1 tuple
                         setpoints = ((instrument.freq_axis(),), (instrument.freq_axis(),)),
                         setpoint_units = (("Hz",), ("Hz",)),
                         setpoint_labels = (
                                               ("f",),
                                               ("f",),
                                           ),
                         setpoint_names = (
                                            ("Frequency",),
                                            ("Frequency",),
                                                            ),

                         #docstring='param that returns two single values, I and Q')
                         )
        self.set_sweep(start, stop, npts)
#

    def set_sweep(self, start: float, stop: float, npts: int):
        # Needed to update config of the software parameter on sweep change
        # frequency setpoints tuple as needs to be hashable for look up.
        f = tuple(np.linspace(int(start), int(stop), num=npts))
        self.setpoints = ((f,), (f,))
        self.shapes = ((npts,), (npts,))

    def get_raw(self):
        R, phase = self.instrument.get_res()
        return R, phase



class OPXSpectrumScan(OPX):
    def __init__(self, config: Dict, name: str = "OPXSpectrumScan", host=None, port=None,host_octave=None, port_octave=None, **kwargs):
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
            "n_points",
            initial_value=100,
            unit="",
            vals=Numbers(1, 1e9),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "Gain",
            initial_value=0,
            unit="dB",
            vals=Numbers(-10, 20),
            get_cmd=None,
            set_cmd=self.set_gain,
        )
        self.add_parameter(
            "t_meas",
            unit="s",
            initial_value=0.01,
            vals=Numbers(0, 1),
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
            numpointsparam=self.n_points,
            vals=Arrays(shape=(self.n_points,)),
        )
        self.add_parameter(
            "amp",
            unit="V",
            initial_value=1,
            vals=Numbers(-5, 5),
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
            start=self.f_start(),
            stop=self.f_stop(),
            npts=self.n_points(),
            parameter_class=IQArray,
        )

    def get_prog(self):
        df = int((self.f_stop() - self.f_start()) / (self.n_points()-1))
        n_avg = round(self.t_meas() * 1e9 / self.readout_pulse_length())
        if self.octave:
            lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
            
        else:
            lo =self.config['mixers']['mixer_resonator'][0]['lo_frequency']
            
        self.parameters['trace_mag_phase'].set_sweep(lo+self.f_start(),lo+self.f_stop(),self.n_points())
        with program() as prog:
            n = declare(int)
            f = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            with for_(n, 0, n < n_avg, n + 1):
                with for_(f, self.f_start(), f <= self.f_stop()+df/2.0, f + df):
                    update_frequency("resonator", f)
                    reset_phase('resonator')
                    measure("readout"*amp(self.amp()), "resonator", None,
                            dual_demod.full("cos", "out1", "sin", "out2", I),
                            dual_demod.full("minus_sin", "out1", "cos", "out2", Q))
                    wait(10000//4,'resonator')
                    save(I, I_st)
                    save(Q, Q_st)

            with stream_processing():
                I_st.buffer(self.n_points()).average().save("I")
                Q_st.buffer(self.n_points()).average().save("Q")

        return prog
    
    def set_octave(self):
        lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
        self.qm.octave.set_lo_source('resonator', OctaveLOSource.Internal)
        self.qm.octave.set_lo_frequency('resonator', lo) # assign the LO inside the octave to element
        self.qm.octave.set_rf_output_gain('resonator', 0) # can set gain from -10dB to 20dB
        self.qm.octave.set_rf_output_mode('resonator', RFOutputMode.trig_normal) # set the behaviour of the RF switch
        self.qm.octave.calibrate_element('resonator',[(lo, self.config['elements']['resonator']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
        self.qm =self.qmm.open_qm(self.config)
        self.qm.octave.set_qua_element_octave_rf_in_port('resonator',"octave1", 1)
        self.qm.octave.set_downconversion('resonator',lo_source=RFInputLOSource.Internal)
        self.qm.octave.set_rf_output_gain('resonator',-10) # can set gain from -10dB to 20dB
        
    def set_gain(self,gain):
        
        self.qm.octave.set_rf_output_gain('resonator', gain)  # can set gain from -10dB to 20dB

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
            df = int((self.f_stop() - self.f_start()) / (self.n_points()-1))
            freqs = np.arange(self.f_start(),self.f_stop()+ df /2, df)
            freq,mag_calib,phase_calib = np.loadtxt("C:\\Users\\Nanoelectronics.ad\\miniconda3\\envs\\qcodes\\Lib\\site-packages\qcodes\\instrument_drivers\\QM_qcodes\\opx_calibration.txt")
            mag_bkg = np.interp(freqs, freq, mag_calib)
            phase_bkg = np.interp(freqs, freq, phase_calib)
            self.result_handles.wait_for_all_values()
            u = unit()
            I = u.demod2volts(self.result_handles.get("I").fetch_all(), self.readout_pulse_length())
            Q = u.demod2volts(self.result_handles.get("Q").fetch_all(), self.readout_pulse_length())
            # R = 20*np.log10(np.sqrt(I ** 2 + Q ** 2)/(self.config['waveforms']['readout_wf']['sample']*self.amp()))
            R = 20*np.log10(np.sqrt(I ** 2 + Q ** 2))-self.Gain()-mag_bkg
            phase = signal.detrend(np.unwrap(np.angle(I + 1j * Q))) * 180 / np.pi - phase_bkg
            return R , phase
