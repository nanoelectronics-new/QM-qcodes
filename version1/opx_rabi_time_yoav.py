from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.OPX.opx_driver import *
from qm.qua import *
from scipy import signal
from qualang_tools.units import unit
# noinspection PyAbstractClass

class IQArray(MultiParameter):
    def __init__(self,
        name: str,
        instrument: "OPXQubitScan",
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
        f = tuple(np.linspace(4*int(start), 4*int(stop), num=npts))
        self.setpoints = ((f,), (f,))
        self.shapes = ((npts,), (npts,))

    def get_raw(self):
        R, phase = self.instrument.get_res()
        return R, phase



class OPXRabiTime(OPX):
    def __init__(self, config: Dict, name: str = "OPXRabiTime", host=None, port=None, **kwargs):
        super().__init__(config, name, host=host, port=port, **kwargs)
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
            vals=Numbers(-2, 2),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "amp_qubit",
            unit="",
            initial_value=1,
            vals=Numbers(-2, 2),
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
        self.parameters['trace_mag_phase'].set_sweep(self.t_start(),self.t_stop(),self.n_points())
        with program() as prog:
            n = declare(int)
            t = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            with for_(n, 0, n < n_avg, n + 1):
                with for_(t, self.t_start(), t <= self.t_stop()+dt/2.0, t + dt):

                    play("gauss"*amp(self.amp_qubit()), "qubit",duration=t) #t in clock cycles (4ns)
                    wait(t, 'resonator')
                    measure("readout"*amp(self.amp_resonator()), "resonator", None,
                            dual_demod.full("cos", "out1", "sin", "out2", I),
                            dual_demod.full("minus_sin", "out1", "cos", "out2", Q))
                    wait(1000//4,'resonator')
                    save(I, I_st)
                    save(Q, Q_st)

            with stream_processing():
                I_st.buffer(self.n_points()).average().save("I")
                Q_st.buffer(self.n_points()).average().save("Q")

        return prog
    
    def run_exp(self):
        self.execute_prog(self.get_prog())

    def get_res(self):
        if (
            self.result_handles is None
        ):
            n = self.n_points()
            return {"I": (0,)*n, "Q": (0,)*n, "R": (0,)*n, "Phi": (0,)*n}
        else:
            u = unit()
            self.result_handles.wait_for_all_values()
            print(self.config['pulses']['readout_pulse']['length'])
            readout_amp = self.config['waveforms']['readout_wf']['sample']
            I = self.result_handles.get("I").fetch_all()
            Q = self.result_handles.get("Q").fetch_all()
            R = 20 * np.log10(np.sqrt((u.demod2volts(I, self.readout_pulse_length()) ** 2) + ( u.demod2volts(Q, self.readout_pulse_length()) ** 2)) / (readout_amp * self.amp_resonator()))
            phase =(np.angle(I + 1j * Q)) * 180 / np.pi
            return R , phase
