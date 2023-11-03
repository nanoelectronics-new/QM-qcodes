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


class Rabi(OPXCustomSequence):
    """


    Parameters
    ----------
    if_resonator:
    t_start:
    t_stop:
    if_qubit
    lo_qubit:
    lo_resonator

    Rabi measurement: sweeping the time of the drive pulse

    Returns
    -------
    self.get_measurement_parameters in the main file fetches the
    results from the OPX

    """

    def __init__(
        self,
        config: Dict,
        name: str = "Rabi",
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
            "t_start",
            initial_value=10,
            unit="ns",
            label="t_start",
            vals=Numbers(4, 1e4),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "t_stop",
            initial_value=100,
            unit="ns",
            label="t_stop",
            vals=Numbers(4, 1e4),
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
            "if_qubit",
            initial_value=100e6,
            unit="Hz",
            label="IF_qubit",
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
            "amp_qubit",
            initial_value=1,
            unit="",
            vals=Numbers(-2, 2),
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
        dt = int((self.t_stop() - self.t_start()) / (self.n_points()-1))
        n = declare(int)
        t = declare(int)
        with for_(n, 0, n < self.n_avg(), n + 1):
            with for_(t, self.t_start(), t <= self.t_stop()+dt/2.0, t + dt):
                 reset_phase('qubit')
                 play("drag", "qubit",duration=t) #t in clock cycles (4ns)
           
                 # play("gauss"*amp(self.amp_qubit()), "qubit",duration=t) #t in clock cycles (4ns)
                 wait(t+17,'resonator')
                 reset_phase(self.readout_element())
                 self.measurement()
      
                 wait(self.wait_time()//4,'resonator','qubit')

    def set_gain_resonator(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain(self.readout_element(), gain)

    def set_gain_qubit(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain("qubit", gain)


    def update_freqs(self):
        """Upload the updates to the OPX and Octave
            based on Instrument Parameters."""
            
        if self.acquisition_mode() == 'sliced_demodulation':
                      
                self.set_sweep_parameters(
                    "axis2",
                    np.linspace(self.t_start(), self.t_stop(), self.n_points()),
                    "ns",
                    "Drive_Time",
                )
                self.set_sweep_parameters(
                    "axis1",
                    np.arange(0, self.readout_pulse_length(),self.slice_size()*4),
                    "ns",
                    "Demod_Time",
                )
        elif self.acquisition_mode() == 'full_demodulation':
            
                self.set_sweep_parameters(
                    "axis1",
                    np.linspace(self.t_start(), self.t_stop(), self.n_points()),
                    "ns",
                    "Drive_Time",
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
        
        self.config['elements']['qubit']['intermediate_frequency'] = self.if_qubit()
        self.config['mixers'][mixer_qubit][0]['intermediate_frequency'] = self.if_qubit()

        self.set_config(self.config)

        # Set octave - update LOs and calibrate

        self.set_octave(['qubit', self.readout_element()])
        self.set_octave_readout(1)  # input channel
        
        self.Gain_qubit(self.Gain_qubit())
        self.Gain_resonator(self.Gain_resonator())

      

















