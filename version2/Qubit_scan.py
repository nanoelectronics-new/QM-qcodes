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
import lmfit


class Qubit_scan(OPXCustomSequence):
    """


    Parameters
    ----------
    if_resonator:
    if_start:
    if_stop:
    lo_qubit:
    lo_resonator

    Qubit spectrscopy with sweeping the IF frequency of the resoantor element

    Returns
    -------
    self.get_measurement_parameters in the main file fetches the
    results from the OPX

    """

    def __init__(
        self,
        config: Dict,
        name: str = "QubitScan",
        host=None,
        cluster_name=None,
        octave = None,
        close_other_machines=False,
        qmm = None
    ):
        super().__init__(config, name, host=host, cluster_name= cluster_name, octave = octave, qmm = qmm)
        self.counter = 0
        self.measurement_variables = None
        self.close_other_machines = close_other_machines

        self.add_parameter(
            "if_start",
            initial_value=100e6,
            unit="Hz",
            label="f start",
            vals=Numbers(-400e6, 400e6),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "if_stop",
            initial_value=100e6,
            unit="Hz",
            label="f stop",
            vals=Numbers(-400e6, 400e6),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "lo_qubit",
            initial_value=5e9,
            unit="Hz",
            label="LO_qubit",
            vals=Numbers(2e9, 18e9),
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
            "lo_resonator",
            initial_value=5e9,
            unit="Hz",
            label="LO_resonator",
            vals=Numbers(2e9, 18e9),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "if_resonator",
            initial_value=100e6,
            unit="Hz",
            label="if_resonator",
            vals=Numbers(-350e6, 350e6),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "Gain_resonator",
            initial_value=0,
            unit="dB",
            vals=Numbers(-20, 20),
            get_cmd=None,
            set_cmd=self.set_gain_resonator,
        )
        self.add_parameter(
            "Gain_qubit",
            initial_value=0,
            unit="dB",
            vals=Numbers(-20, 20),
            get_cmd=None,
            set_cmd=self.set_gain_qubit,
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
        self.opx_scan("1d")
        self.acquisition_mode("full_demodulation")

    def pulse_sequence(self):
        """
        OPX pulse sequnce.

        Sweeping the IF frequency of the qubit and measure resonator's resposne

        """
        df = int((self.if_stop() - self.if_start()) / (self.n_points()-1))
        n = declare(int)
        f = declare(int)
        with for_(n, 0, n < self.n_avg(), n + 1):
            with for_(f, self.if_start(), f <= self.if_stop()+df/2.0, f + df):
                update_frequency('qubit', f)
                play("saturation", "qubit")
                wait(((2000-32)//4),self.readout_element())
                reset_phase(self.readout_element())
                self.measurement()
                wait(self.wait_time()//4, 'qubit', self.readout_element())

    def set_gain_resonator(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain(self.readout_element(), gain)

    def set_gain_qubit(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain("qubit", gain)


    def update_freqs(self):
        """Upload the updates to the OPX and Octave
            based on Instrument Parameters."""

        self.set_sweep_parameters(
            "axis1",
            np.linspace(self.lo_qubit()+self.if_start(), self.lo_qubit() +
                        self.if_stop(), self.n_points()),
            "Hz",
            "Frequency",
        )

        # Update config file
        mixer_resonator = self.config['elements'][self.readout_element()
                                                  ]['mixInputs']['mixer']
        mixer_qubit = self.config['elements']['qubit'
                                              ]['mixInputs']['mixer']

        self.config['elements'][self.readout_element(
        )]['mixInputs']['lo_frequency'] = self.lo_resonator()
        self.config['elements']['qubit']['mixInputs']['lo_frequency'] = self.lo_qubit()

        self.config['mixers'][mixer_resonator][0]['lo_frequency'] = self.lo_resonator()
        self.config['mixers'][mixer_qubit][0]['lo_frequency'] = self.lo_qubit()

        self.config['elements'][self.readout_element(
        )]['intermediate_frequency'] = self.if_resonator()
        self.config['mixers'][mixer_resonator][0]['intermediate_frequency'] = self.if_resonator()

        self.set_config(self.config)

        # Set octave - update LOs and calibrate

        self.set_octave(['qubit', self.readout_element()])
        self.set_octave_readout(1)  # input channel
        
        self.Gain_qubit(self.Gain_qubit())
        self.Gain_resonator(self.Gain_resonator())
        
    def fit_qubit(self) -> list:
        
        data= self.get_res()
        freq = self.axis1_axis()
        
        phi =data['Phi']
        freq_guess = freq[np.argmax(phi)]
        
        peak1 = lmfit.models.LorentzianModel(prefix='l1_')
      #  peak2 = lmfit.models.LorentzianModel(prefix='l2_')
        background = lmfit.models.LinearModel()
        model = background + peak1 
        model.param_names
        params = model.make_params(slope=1e-5,intercept=np.min(phi),l1_amplitude=1,l1_center=freq_guess,l1_sigma=0.025)
        result = model.fit(phi, params, x=freq/1e9)
        
        return result.best_values['l1_center']
        
      

















