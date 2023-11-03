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
from qualang_tools.bakery import baking
from qualang_tools.bakery.bakery import deterministic_run
from scipy import signal
import lmfit
from qm import generate_qua_script
from qualang_tools.config.waveform_tools import *


class T2(OPXCustomSequence):
    """


    Parameters
    ----------
    if_resonator:
    t_start:
    t_stop:
    if_qubit
    lo_qubit:
    lo_resonator

    T2 measurement

    Returns
    -------
    self.get_measurement_parameters in the main file fetches the
    results from the OPX

    """

    def __init__(
        self,
        config: Dict,
        name: str = "T2",
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
            vals=Numbers(0, 1e4),
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
        self.add_parameter(
            "detuning",
            unit="Hz",
            initial_value=10e6,
            vals=Numbers(-100e6, 100e6),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "t_pihalf",
            unit="ns",
            initial_value=5,
            vals=Numbers(4, 100),
            get_cmd=None,
            set_cmd=None,
        )



        # Set sweep parameters
        self.opx_scan("1d")
        self.acquisition_mode("full_demodulation")

    def baking_list(self):

        t_min = self.t_start()
        t_max = self.t_stop()
        dt = int((t_max-t_min)/ (self.n_points()-1))
        t_vec = np.arange(t_min, t_max + dt / 2, dt)

        sample_rate =  1e9
        n_samples = len(t_vec)

        max_delay = t_max  # In samples

        short_baking_list = []  # Stores the baking objects
        long_baking_list = []
        # Create the different baked sequences, corresponding to the different taus
        drag_I,drag_Q  = drag_cosine_pulse_waveforms(0.25, self.t_pihalf(), alpha = 0.25, anharmonicity = -78e6, detuning=0)


        for i in range(16):
            with baking(self.config, padding_method="left") as b:

              # Assign waveforms to quantum element operation

                b.add_op("pihalf", "qubit", [drag_I, drag_Q], digital_marker = 'ON')
               # b.add_op("minus_pihalf", "qubit", [mixInput_sample_Iminus, mixInput_sample_Q], digital_marker = 'ON')
                init_delay = 16# Put initial delay to ensure that all of the pulses will have the same length
                b.wait(init_delay, 'qubit')  # We first wait the entire duration.

                # We add the 2nd pi_half pulse with the phase 'dephasing' (Confusingly, the first pulse will be added later)
                # Play uploads the sample in the original config file (here we use an existing pulse in the config)
                b.frame_rotation_2pi(
                      float(self.detuning() * 1e-9*  i), "qubit"
                   )
                b.play("pihalf", "qubit")

                # We reset frame such that the first pulse will be at zero phase
                # and such that we will not accumulate phase between iterations.
                b.reset_frame("qubit")
                # We add the 1st pi_half pulse. It will be added with the frame at time init_delay - i, which will be 0.
                b.play_at("pihalf", "qubit", t=init_delay - i -self.t_pihalf())

            # Append the baking object in the list to call it from the QUA program
            short_baking_list.append(b)

        self.open_qm(close_other_machines = False)

        for i in range(4):
            with baking(self.config, padding_method="left") as b:
                b.add_op("pihalf", "qubit", [drag_I, drag_Q],
                         digital_marker='ON')
                # We play the first pulse then add a wait of i (which goes from 0 to 3).
                # The padding is on the left, so the extra samples are added there
             
                b.play("pihalf", "qubit")
                b.wait(i, "qubit")
            long_baking_list.append(b)

        self.open_qm(close_other_machines=False)

        return short_baking_list, long_baking_list



    def pulse_sequence(self):
        """
        OPX pulse sequnce.

        Sweeping the IF frequency of the qubit and measure resonator's resposne

        """
        short_baking_list, long_baking_list = self.baking_list()
        dt = int((self.t_stop() - self.t_start()) / (self.n_points()-1))
        n = declare(int)
        t = declare(int)
        t_cycles = declare(int)
        t_left_ns = declare(int)

        with for_(n, 0, n < self.n_avg(), n + 1):
            with for_(t, self.t_start(), t <= self.t_stop() + dt/2.0, t+dt):
                reset_phase('qubit')
                with if_(t < 16):
                    deterministic_run(short_baking_list, t, unsafe=True)
                    wait(4,'resonator')
                with else_():
                    assign(t_cycles, t >> 2)
                    # left shift by 2 is a quick way to multiply by 4
                    assign(t_left_ns, t-(t_cycles << 2))
                    with switch_(t_left_ns, unsafe=True):
                        for j in range(4):
                            with case_(j):
                                long_baking_list[j].run() #4 cycles
                                wait(t_cycles, 'qubit')
                                frame_rotation_2pi(
                                    Cast.mul_fixed_by_int(self.detuning()* 1e-9, t),
                                    "qubit")
                                play("drag_2", 'qubit') #4 cycles
                                reset_frame('qubit')
                                wait(t_cycles+8,'resonator')
                                
                reset_phase(self.readout_element())
                self.measurement()
                wait(self.wait_time()//4, 'resonator', 'qubit')

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
            np.linspace(self.t_start(), self.t_stop(), self.n_points()),
            "ns",
            "Time",
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





















