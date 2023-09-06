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

class Resonator_scan(OPXCustomSequence):
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
        name: str = "ResonatorScan",
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
            "lo",
            initial_value=5e9,
            unit="Hz",
            label="LO",
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
        self.opx_scan("1d")
        self.acquisition_mode("full_demodulation")

    def pulse_sequence(self):
        """
        OPX pulse sequnce.

        Sweeping the IF frequency of the OPX and measure resonator's resposne

        """
        df = int((self.if_stop() - self.if_start()) / (self.n_points()-1))
        n = declare(int)
        f = declare(int)
        with for_(n, 0, n < self.n_avg(), n + 1):
            with for_(f, self.if_start(), f <= self.if_stop()+df/2.0, f + df):
                update_frequency(self.readout_element(), f)
                reset_phase(self.readout_element())
                self.measurement()
                wait(self.wait_time()//4, self.readout_element())

    def set_gain(self, gain):

        # Set gains customized for the elements we use in that specific pusle sequnece
        self.qm.octave.set_rf_output_gain(self.readout_element(), gain)

    def get_res(self):
        """Overwriting parent function to subtract linear
            background.
        """
        df = int((self.if_stop() - self.if_start()) / (self.n_points()-1))
        freqs = np.arange(self.if_start(),self.if_stop()+ df /2, df)
        freq,mag_calib,phase_calib = np.loadtxt("C:\\Users\\Nanoelectronics.ad\\miniconda3\\envs\\qcodes\\Lib\site-packages"
                                                "\\qcodes\\instrument_drivers\\QM_qcodes\\version2\\opx_calibration.txt")
        mag_bkg = 10**(np.interp(freqs, freq, mag_calib)/20)
        mag_bkg = mag_bkg / np.max(mag_bkg)
        phase_bkg = np.interp(freqs, freq, phase_calib)

        output = super().get_res()
        output['R'] = output['R'] / mag_bkg
        output['Phi'] = signal.detrend(
            output['Phi']*np.pi/180)*180/np.pi - phase_bkg
        output['I'] = output['R']*np.cos(output['Phi']*np.pi/180)
        output['Q'] = output['R']*np.sin(output['Phi']*np.pi/180)

        return output
    
    def update_freqs(self) -> None:
        """Upload the updates to the OPX and Octave
            based on Instrument Parameters."""

        self.set_sweep_parameters(
            "axis1",
            np.linspace(self.lo()+self.if_start(), self.lo() +
                        self.if_stop(), self.n_points()),
            "Hz",
            "Frequency",
        )

        # Update config file
        mixer = self.config['elements'][self.readout_element()
                                        ]['mixInputs']['mixer']
        self.config['elements'][self.readout_element(
        )]['mixInputs']['lo_frequency'] = self.lo()
        self.config['mixers'][mixer][0]['lo_frequency'] = self.lo()

        self.set_config(self.config)

        # Set octave - update LOs and calibrate

        self.set_octave([self.readout_element()]) #reopens qm
        self.set_octave_readout(1)  # input channel
        self.Gain(self.Gain())
        
    def fit_resonator(self) -> list:
        
        data = self.get_res()
        
        freqs = self.axis1_axis()
        
        [freq,R_bkg,phi_bkg] = np.loadtxt('background2')
        data['R']  = data['R'] / np.interp(freqs, freq, R_bkg)
        data['Phi']  = data['Phi'] - np.interp(freqs, freq, phi_bkg)
        complex_data = data['R']*np.exp(1j*data['Phi']*np.pi/180)
        
        port1 = circuit.notch_port() 
        port1.add_data(freqs,complex_data)
        port1.autofit()
        port1.plotall()
        
        return [port1.fitresults['fr'], port1.fitresults['Ql']]
    

        

















