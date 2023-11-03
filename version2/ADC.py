import os

import numpy as np
from qcodes import initialise_or_create_database_at, load_or_create_experiment
from qcodes import Parameter
from qm.qua import *
from qualang_tools.loops import from_array
from qcodes.utils.validators import Numbers
from qcodes.instrument_drivers.QM_qcodes.version2.advanced_driver.advanced_driver import OPXCustomSequence
from qcodes.instrument_drivers.QM_qcodes.version2.basic_driver.opx_driver import ResultParameters
from qm.octave import *
from scipy import signal
from resonator_tools import circuit
from scipy.optimize import curve_fit

def linear(x, a, b):
    return a * x + b

class ADC(OPXCustomSequence):
    """


    Parameters
    ----------
    OPXCustomSequence : TYPE
    if_start
    if_stop
    lo
    n_avg

    Resonator spectrscopy with sweeping the IF frequency of the resoantor element

    Returns
    -------
    self.get_measurement_parameters

    """


    def __init__(
        self,
        config: Dict,
        name: str = "ADC",
        host=None,
        cluster_name =None,
        octave = None,
        close_other_machines=False,
        qmm = None
    ):
        super().__init__(config, name, host=host, cluster_name= cluster_name, octave = octave, qmm=qmm)
        self.counter = 0
        self.measurement_variables = None
        self.close_other_machines = close_other_machines

        self.add_parameter(
            "IF",
            initial_value=50e6,
            unit="Hz",
            label="IF",
            vals=Numbers(-350e6, 350e6),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "lo",
            initial_value=5e9,
            unit="Hz",
            label="LO",
            vals=Numbers(2e9, 18e9),
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
            "wait_time",
            initial_value=5000,
            unit="ns",
            vals=Numbers(100, 1e9),
            get_cmd=None,
            set_cmd=None,
        )

        # Set sweep parameters
        self.opx_scan("0d")
        self.acquisition_mode("raw_adc")

    def pulse_sequence(self):
        """
        OPX pulse sequnce.

        Sweeping the IF frequency of the OPX and measure resonator's resposne

        """
        n = declare(int)
        update_frequency(self.readout_element(), self.IF())
        with for_(n, 0, n < self.n_avg(), n + 1):
            
            reset_phase(self.readout_element())
            self.measurement()
            wait(self.wait_time()//4, self.readout_element())

    def set_gain(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain(self.readout_element(), gain)

    
    def update_freqs(self) -> None:
        """Upload the updates to the OPX and Octave
            based on Instrument Parameters."""

        self.set_sweep_parameters(
            "axis1",
            np.linspace(0, self.readout_pulse_length(), self.readout_pulse_length()),                        
            "ns",
            "time",
        )
        self.set_sweep_parameters(
            "axis2",
            np.linspace(0, self.readout_pulse_length(), self.readout_pulse_length()),                        
            "ns",
            "time2",
        )


        # Update config file
        mixer = self.config['elements'][self.readout_element()
                                        ]['mixInputs']['mixer']
        self.config['elements'][self.readout_element(
        )]['mixInputs']['lo_frequency'] = self.lo()
        self.config['mixers'][mixer][0]['lo_frequency'] = self.lo()
        self.config['elements'][self.readout_element(
        )]['intermediate_frequency'] = self.IF()
        self.config['mixers'][mixer][0]['intermediate_frequency'] = self.IF()

        self.set_config(self.config)

        # Set octave - update LOs and calibrate

        self.set_octave([self.readout_element()]) #reopens qm
        self.set_octave_readout(1)  # input channel
        self.Gain(self.Gain())

    

        

















