import time
from typing import Dict, Optional

import numpy as np
from qcodes import (
    Instrument,
    Parameter, MultiParameter,
)
from qcodes.utils.helpers import abstractmethod
from qcodes.utils.validators import Numbers
from qm import SimulationConfig
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.octave import *
import os
import matplotlib.pyplot as plt

# noinspection PyAbstractClass
class OPX(Instrument):
    """
    Driver for interacting with QM OPX
    """

    def __init__(
            self, config: Dict, name: str = "OPX", host=None, port=None,host_octave = None, port_octave=None, **kwargs) -> None:
        """
        Args:
            name: Name to use internally in QCoDeS
        """
        super().__init__(name, **kwargs)

        self.qm = None
        self.qmm = None
        self.config = None
        self.result_handles = None
        self.job = None
        self.prog_id = None
        self.octave = True
        self.connect(host=host, port=port,host_octave=host_octave,port_octave=port_octave)
        self.set_config(config=config)
        self.add_parameter("results", label="results", get_cmd=self.get_res)

        self.add_parameter(
            "sim_time",
            unit="ns",
            label="sim_time",
            initial_value=100000,
            vals=Numbers(
                4,
            ),
            get_cmd=None,
            set_cmd=None,
        )

    @abstractmethod
    def get_prog(self):
        pass

    def get_res(self):
        return None

    def execute_prog(self, prog):
        self.job = self.qm.execute(prog)
        self.result_handles = self.job.result_handles

    def simulate_prog(self, prog, duration=1000):
        self.job = self.qm.simulate(prog, SimulationConfig(duration))
        samples = self.job.get_simulated_samples()
        plt.figure()
        samples.con1.plot()
        self.result_handles = self.job.result_handles

    def simulate_and_read(self, prog):
        self.simulate_prog(prog, duration=self.sim_time())
        self.result_handles.wait_for_all_values()
        return self.get_res()

    def compile_prog(self, prog):
        self.prog_id = self.qm.compile(prog)

    def execute_compiled_prog(self):
        if self.prog_id is not None:
            pending_job = self.qm.queue.add_compiled(self.prog_id)
            self.job = pending_job.wait_for_execution()
            self.result_handles = self.job.result_handles

    def set_config(self, config):
        self.config = config
        # self.qmm.close_all_quantum_machines()
        self.qm = self.qmm.open_qm(self.config, close_other_machines=True)

    def connect(self, host=None, port=None, host_octave=None, port_octave=None):
        if self.octave:            
            octave_config = QmOctaveConfig()
            # octave_config.set_calibration_db("D:\\Oliver\\Resonator_measurements\\20230316_Farid_qubit_Cu_cavity")
            octave_config.add_device_info('octave1', host_octave, port_octave)
            octave_config.set_opx_octave_mapping([('con1', 'octave1')]) # set default mapping between analog outputs of OPX and the octave       
            self.qmm = QuantumMachinesManager(host=host, port=port,octave= octave_config)
            self.qmm.octave_manager.set_clock("octave1", ClockType.Internal, ClockFrequency.MHZ_10) # internal clock
        else:
            self.qmm = QuantumMachinesManager(host=host, port=port)
            
        self.connect_message()

    def close(self) -> None:
        if self.qm is not None:
            self.qm.close()
        super().close()

    def close(self) -> None:
        self.qmm.close_all_quantum_machines()
        super().close()

    def halt(self) -> None:
        if self.job is not None:
            self.job.halt()

    def connect_message(self, idn_param: str = 'IDN',
                        begin_time: Optional[float] = None) -> None:
        """
        Print a standard message on initial connection to an instrument.

        Args:
            idn_param: Name of parameter that returns ID dict.
                Default ``IDN``.
            begin_time: ``time.time()`` when init started.
                Default is ``self._t0``, set at start of ``Instrument.__init__``.
        """
        idn = {"vendor": "Quantum Machines"}
        idn.update(self.get(idn_param))
        if idn["server"][0] == "1":
            idn["model"] = "OPX"
        elif idn["server"][0] == "2":
            idn["model"] = "OPX+"
        else:
            idn["model"] = ""
        t = time.time() - (begin_time or self._t0)

        con_msg = (
            "Connected to: {vendor} {model} in {t:.2f}s. "
            "QOP Version = {server}, SDK Version = {client}.".format(t=t, **idn)
        )
        print(con_msg)
        self.log.info(f"Connected to instrument: {idn}")

    def get_idn(self) -> Dict[str, Optional[str]]:
        return self.qmm.version()


# noinspection PyAbstractClass
class GeneratedSetPoints(Parameter):
    """
    A parameter that generates a setpoint array from start, stop and num points
    parameters.
    """

    def __init__(self, startparam, stopparam, numpointsparam, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._startparam = startparam
        self._stopparam = stopparam
        self._numpointsparam = numpointsparam

    def get_raw(self):
        return np.linspace(
            self._startparam(), self._stopparam(), self._numpointsparam()
        )


# noinspection PyAbstractClass
class QMDemodParameters(MultiParameter):
    def __init__(self, instr, params, name, names, units, shapes=None, setpoints=None, *args, **kwargs):
       
        if shapes is None:
            shapes = ((),) * len(params)
        if setpoints is None:
            setpoints = ((),) * len(params)
        super().__init__(name=name, names=names, units=units,
            shapes=shapes,
            setpoints=setpoints,
            *args, **kwargs)

        self._instr = instr
        self._params = params

    def get_raw(self):
        vals1 = []
        vals2 = []
        result = self._instr.get_res()
        for param in self._params:
            if param == 'X' or param == 'I':
                vals1.append(result['I'])
            elif param == 'Y' or param == 'Q':
                vals2.append(result['Q'])
            elif param == 'R':
                vals1.append(result['R'])
            elif param == 'Phase' or param == 'Phi':
                vals2.append(result['Phi'])
            else:
                raise NotImplementedError("Only X (I), Y (Q), R or Phase (Phi) are valid inputs")
                
        return np.array([vals1[0], vals2[0]])
