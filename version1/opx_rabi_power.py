from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.OPX.opx_driver import *
from qm.qua import *
from scipy import signal

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
                         setpoints = ((instrument.power_axis(),), (instrument.power_axis(),)),
                         setpoint_units = (("V",), ("V",)),
                         setpoint_labels = (
                                               ("Voltage",),
                                               ("Voltage",),
                                           ),
                         setpoint_names = (
                                            ("Voltage",),
                                            ("Voltage",),
                                                            ),

                         #docstring='param that returns two single values, I and Q')
                         )
        self.set_sweep(start, stop, npts)
#

    def set_sweep(self, start: float, stop: float, npts: int):
        # Needed to update config of the software parameter on sweep change
        # frequency setpoints tuple as needs to be hashable for look up.
        f = tuple(np.linspace(start, stop, num=npts))
        self.setpoints = ((f,), (f,))
        self.shapes = ((npts,), (npts,))

    def get_raw(self):
        R, phase = self.instrument.get_res()
        return R, phase



class OPXRabiPower(OPX):
    def __init__(self, config: Dict, name: str = "OPXRabiPower", host=None, port=None, **kwargs):
        super().__init__(config, name, host=host, port=port, **kwargs)
        self.add_parameter(
            "V_start",
            initial_value=0.1,
            unit="V",
            label="V start",
            vals=Numbers(-0.5, 0.5),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "V_stop",
            initial_value=0.1,
            unit="V",
            label="V stop",
            vals=Numbers(-0.5,0.5),
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
            vals=Numbers(0, 1),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "power_axis",
            unit="V",
            label="Power Axis",
            parameter_class=GeneratedSetPoints,
            startparam=self.V_start,
            stopparam=self.V_stop,
            numpointsparam=self.n_points,
            vals=Arrays(shape=(self.n_points,)),
        )
        self.add_parameter(
            "amp_resonator",
            unit="V",
            initial_value=0.1,
            vals=Numbers(-0.5, 0.5),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "drive_length",
            unit="clock cycle(4ns)",
            initial_value=50,
            vals=Numbers(5,1000),
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
            start=self.V_start(),
            stop=self.V_stop(),
            npts=self.n_points(),
            parameter_class=IQArray,
        )

    def get_prog(self):
        dV = (self.V_stop() - self.V_start()) / (self.n_points()-1)
        n_avg = round(self.t_meas() * 1e9 / self.readout_pulse_length())
        self.parameters['trace_mag_phase'].set_sweep(self.V_start(),self.V_stop(),self.n_points())
        with program() as prog:
            n = declare(int)
            V = declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            with for_(n, 0, n < n_avg, n + 1):
                with for_(V, self.V_start(), V < self.V_stop()+dV/2.0, V + dV):
                    play("gauss"*amp(V), "qubit",duration=self.drive_length())
                    align("qubit", "resonator")
                    measure("readout"*amp(self.amp_resonator()), "resonator", None,
                            dual_demod.full("cos", "out1", "sin", "out2", I),
                            dual_demod.full("minus_sin", "out1", "cos", "out2", Q))
                    wait(10000//4,'resonator') #10us repetion rate
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
            self.result_handles.wait_for_all_values()
            I = self.result_handles.get("I").fetch_all()*2*4096/ self.readout_pulse_length()
            Q = self.result_handles.get("Q").fetch_all()*2*4096/self.readout_pulse_length()
            R = 20*np.log(np.sqrt(I ** 2 + Q ** 2)/self.amp_resonator())
            phase = np.angle(I + 1j * Q) * 180 / np.pi
            return R , phase
