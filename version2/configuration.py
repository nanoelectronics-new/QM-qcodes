import numpy as np
from scipy.signal.windows import gaussian,general_gaussian,cosine
from qualang_tools.config.waveform_tools import *


#######################
# AUXILIARY FUNCTIONS #
#######################

# IQ imbalance matrix
def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1 - g ** 2) * (2 * c ** 2 - 1))
    return [float(N * x) for x in [(1 - g) * c, (1 + g) * s, (1 - g) * s, (1 + g) * c]]


#############
# VARIABLES #
#############


qubit_if_freq = 55e6
pump_if_freq = 30e6
pump_lo_freq = 5e9

qubit_lo_freq = 2.7e9
# readout_len = 1e6
calibration_amp = 0.125  # better to be between 0.1V and 0.15V
calibration_pulse_length = 10e3
time_of_flight = 192

opx_ip = '10.21.41.199'
opx_port = 80
octave_ip = '10.21.41.199'
octave_port = 11050


saturation_len = 2000
saturation_amp = 0.2
const_len = 40000
const_amp = 0.125

order = 1
gauss_len = 16
gauss_sigma = gauss_len / 5
gauss_amp = 0.3
alpha_drag = 0.15
#gauss_wf = gauss_amp * general_gaussian(gauss_len, order, gauss_sigma)

anharmonicity = -70e6
drag_I_gen,drag_Q_gen  = drag_cosine_pulse_waveforms(gauss_amp, 5, alpha =alpha_drag, anharmonicity = -80e6, detuning=0)


    
def delayed_wf(amp, length):
    delay = 16 - length
    if delay < 0:
        return amp 
    return np.r_[np.zeros(delay), amp*cosine(length)]


def delayed_wf_2(waveform):
    length =len(waveform)
    delay = 16 - length
    if delay < 0:
        return np.array(waveform)
    return np.r_[np.zeros(delay), waveform]

def delayed_wf_2_rev(waveform):
    length = len(waveform)
    delay = 16 - length
    if delay < 0:
        return waveform
    return np.r_[waveform, np.zeros(delay)]


gauss_wf = delayed_wf(gauss_amp,4)
drag_I = delayed_wf_2(drag_I_gen)
drag_Q = delayed_wf_2(drag_Q_gen)


drag_I_2 = delayed_wf_2_rev(drag_I_gen)
drag_Q_2 = delayed_wf_2_rev(drag_Q_gen)


pi_len = 148
pi_sigma = pi_len / 5
pi_amp = 0.125
# pi_wf = pi_amp * gaussian(pi_len, pi_sigma)



pi_half_len = 8
pi_half_sigma = pi_half_len / 5
pi_half_amp = gauss_amp
pi_half_wf = delayed_wf(pi_half_amp,pi_half_len)


# Resonator



resonator_if_freq = 10e6
#resonator_lo_freq =4.7e9 # Internal synthesizer in the octave
resonator_lo_freq =4.5e9 # Internal synthesizer in the octave

resonator_fake_if_freq = 130e6
resonator_fake_lo_freq =4.7e9 # Internal synthesizer in the octave


short_readout_len = 500
short_readout_amp = 0.4
readout_len = 100
#readout_len = 68
readout_amp = 0.1
long_readout_len = 10000
long_readout_amp = 0.1

# Resonator kick-in

kick_length = 68

kick_amp = 4.5*readout_amp
ringup = 0
clear = np.concatenate([np.full(kick_length, kick_amp),
                        np.full(readout_len, readout_amp)])



# Gate
step_num = 200  # number of dc values for scanning gate
scan_vpp = 0.3  # Max is 0.99
bias_amp = 0.5
bias_len = 2000

##########
# CONFIG #
##########

config = {
    'version': 1,
    'controllers': {
        'con1': {
            'type': 'opx1',
            'analog_outputs': {
                1: {'offset': 0.0,  "shareable": True},  # I qubit
                2: {'offset': 0.0,  "shareable": True},  # Q qubit
                3: {'offset': 0.0,  "shareable": True},  # I resonator
                4: {'offset': 0.0,  "shareable": True},  # Q resonator
                5: {'offset': 0.0,  "shareable": True},  # Qubit Gate
                #6: {'offset': 0.0},  # Qubit Gate
                # 7: {'offset': 0.0},  # Qubit Gate
                # 8: {'offset': 0.0},  # Qubit Gate

            },
            'digital_outputs': {
                1: {"shareable": True},
                2: {"shareable": True},
                3: {"shareable": True},
                4: {"shareable": True},
                5: { "shareable": True},
                6: { "shareable": True}},
            'analog_inputs': {
                1: {'offset': 0.0, 'gain_db': 20,  "shareable": True},  # I from down-conversion
                2: {'offset': 0.0, 'gain_db': 20,  "shareable": True},  # Q from down-conversion
            },
        },
    },
    'elements': {
        'qubit': {
            'mixInputs': {
                'I': ('con1', 3),
                'Q': ('con1', 4),
                'lo_frequency': qubit_lo_freq,
                'mixer': 'octave_octave1_2',
            },
            'intermediate_frequency': qubit_if_freq,
            'operations': {
                'cw': 'const_pulse',
                'saturation': 'saturation_pulse',
                'gauss': 'gaussian_pulse',
                'drag': 'drag_pulse',
                'drag_2': 'drag_pulse_2',
                'pi': 'pi_pulse',
                'pi_half': 'pi_half_pulse',
            },
            "digitalInputs": {
                "switch": {
                      "port": ("con1", 2), #Port number of digital marker that is conencted to the Trig_i input of the octacve
                      "delay": 90,
                      "buffer": 10,
                  },
                },
        },
        'qubit2': { #Defined for three-tone spectroscopy measurements
            'mixInputs': {
                'I': ('con1', 3),
                'Q': ('con1', 4),
                'lo_frequency': qubit_lo_freq,
                'mixer': 'octave_octave1_2',
            },
            'intermediate_frequency': qubit_if_freq,
            'operations': {
                'cw': 'const_pulse',
                'saturation': 'saturation_pulse',
                'gauss': 'gaussian_pulse',
                'pi': 'pi_pulse',
                'pi_half': 'pi_half_pulse',
            },
            "digitalInputs": {
                "switch": {
                      "port": ("con1", 2), #Port number of digital marker that is conencted to the Trig_i input of the octacve
                      "delay": 81,
                      "buffer": 10,
                  },
                },
        },
        
        # 'resonator_fake': {
        #     'mixInputs': {
        #         'I': ('con1', 1),
        #         'Q': ('con1', 2),
        #         'lo_frequency': resonator_fake_lo_freq,
        #         'mixer': 'octave_octave1_1',
        #     },
      #       'intermediate_frequency': resonator_fake_if_freq,
      #       'operations': {
      #           'readout': 'readout_pulse',
      #           'cw': 'const_pulse',
      #           'saturation': 'saturation_pulse',
      #           'gauss': 'gaussian_pulse',
      #           'pi': 'pi_pulse',
      #           'pi_half': 'pi_half_pulse',
      #       },
      #       "digitalInputs": {
      #           "switch": {
      #                 "port": ("con1", 1), #Port number of digital marker that is conencted to the Trig_i input of the octacve
      #                 "delay": 100,
      #                 "buffer": 0,
      #             },
      # },
      #   },
        'resonator': {
            'mixInputs': {
                'I': ('con1', 1),
                'Q': ('con1', 2),
                'lo_frequency': resonator_lo_freq,
                'mixer':'octave_octave1_1',
            },
            'intermediate_frequency': resonator_if_freq,
            'operations': {
                'cw': 'const_pulse',
                'short_readout': 'short_readout_pulse',
                'readout': 'readout_pulse',
                'drag': 'drag_pulse',
                'drag_2': 'drag_pulse_2',
                'long_readout': 'long_readout_pulse',
                'readout_zero': 'readout_pulse_zero',
            },            
             "digitalInputs": {
                 "switch": {
                       "port": ("con1", 1), #Port number of digital marker that is conencted to the Trig_i input of the octacve
                       "delay": 82,
                       "buffer": 12,
                   },
       },
            'outputs': {
                'out1': ('con1', 1),
                'out2': ('con1', 2),
            },
            'time_of_flight': 264, #252
            'smearing': 0,
        },
        # 'gate_sticky': {
        #     "singleInput": {
        #         "port": ("con1", 5),
        #     },
        #     "digitalInputs": {},
        #     'hold_offset': {'duration': 1},
        #     "operations": {
        #         "step": "step_pulse",
        #         'bottom': 'bottom_pulse',
        #     },
        # },
        'gate': {
            "singleInput": {
                "port": ("con1", 5),
            },
            "digitalInputs": {},
            'intermediate_frequency': 0,
            "operations": {
                "bias": "bias_pulse",
            },
        },
    },
    'pulses': {
        'const_pulse': {
            'operation': 'control',
            'length': const_len,
            'waveforms': {
                'I': 'const_wf',
                'Q': 'zero_wf',
            },
            'digital_marker': 'ON',
        },
        'saturation_pulse': {
            'operation': 'control',
            'length': saturation_len,
            'waveforms': {
                'I': 'saturation_drive_wf',
                'Q': 'zero_wf'
            },
            'digital_marker': 'ON',
        },
        'gaussian_pulse': {
            'operation': 'control',
            'length': gauss_len,
            'waveforms': {
                'I': 'gauss_wf',
                'Q': 'zero_wf',
            },
            'digital_marker': 'ON',
        },
        'drag_pulse': {
            'operation': 'control',
            'length': 16,
            'waveforms': {
                'I': 'drag_I_wf',
                'Q': 'drag_Q_wf',
            },
            'digital_marker': 'ON',
        },
        'drag_pulse_2': {
            'operation': 'control',
            'length': 16,
            'waveforms': {
                'I': 'drag_I_wf_2',
                'Q': 'drag_Q_wf_2',
            },
            'digital_marker': 'ON',
        },
        'pi_pulse': {
            'operation': 'control',
            'length': pi_len,
            'waveforms': {
                'I': 'pi_wf',
                'Q': 'zero_wf',
            },
            'digital_marker': 'ON',
        },
        'pi_half_pulse': {
            'operation': 'control',
            'length': 16,
            'waveforms': {
                'I': 'pi_half_wf',
                'Q': 'zero_wf',
            },
            'digital_marker':'ON',
        },        
        'short_readout_pulse': {
            'operation': 'measurement',
            'length': short_readout_len,
            'waveforms': {
                'I': 'short_readout_wf',
                'Q': 'zero_wf',
            },
            'integration_weights': {
                'cos': 'short_cosine_weights',
                'sin': 'short_sine_weights',
                'minus_sin': 'short_minus_sine_weights',
            },
            'digital_marker': 'ON',
        },
        'readout_pulse': {
            'operation': 'measurement',
            'length': readout_len+kick_length,
            'waveforms': {
                'I': 'readout_wf',
                'Q': 'zero_wf',
            },
            'integration_weights': {
                'cos': 'cosine_weights',
                'sin': 'sine_weights',
                'minus_sin': 'minus_sine_weights',
            },
            'digital_marker': 'ON',
            },
        'readout_pulse_zero': {
            'operation': 'measurement',
            'length': readout_len,
            'waveforms': {
                'I': 'zero_wf',
                'Q': 'zero_wf',
            },
            'integration_weights': {
                'cos': 'cosine_weights',
                'sin': 'sine_weights',
                'minus_sin': 'minus_sine_weights',
            },
            'digital_marker': 'ON',
            },
        'long_readout_pulse': {
            'operation': 'measurement',
            'length': long_readout_len,
            'waveforms': {
                'I': 'long_readout_wf',
                'Q': 'zero_wf',
            },
            'integration_weights': {
                'cos': 'long_cosine_weights',
                'sin': 'long_sine_weights',
                'minus_sin': 'long_minus_sine_weights',
            },
            'digital_marker': 'ON',
        },
        'step_pulse': {
            'operation': 'control',
            'length': 16,
            'waveforms': {
                'single': 'step_wf',
            },
        },
        'bottom_pulse': {
            'operation': 'control',
            'length': 16,
            'waveforms': {
                'single': 'bottom_wf',
            },
        },
        'bias_pulse': {
            'operation': 'control',
            'length': bias_len,
            'waveforms': {
                'single': 'bias_wf',
            },
        },
    },

    'waveforms': {
        'const_wf': {'type': 'constant', 'sample': const_amp},
        'saturation_drive_wf': {'type': 'constant', 'sample': saturation_amp},
        'zero_wf': {'type': 'constant', 'sample': 0.0},
        'gauss_wf': {'type': 'arbitrary', 'samples': gauss_wf.tolist()},
        'drag_Q_wf': {'type': 'arbitrary', 'samples': drag_Q.tolist()},
        'drag_I_wf': {'type': 'arbitrary', 'samples': drag_I.tolist()},
        'drag_Q_wf_2': {'type': 'arbitrary', 'samples': drag_Q_2.tolist()},
        'drag_I_wf_2': {'type': 'arbitrary', 'samples': drag_I_2.tolist()},
        'pi_wf': {'type': 'constant', 'sample': pi_amp},
        'pi_half_wf': {'type': 'arbitrary', 'samples': pi_half_wf.tolist()},
        'short_readout_wf': {'type': 'constant', 'sample': short_readout_amp},
        # 'readout_wf': {'type': 'constant', 'sample': readout_amp},
        'readout_wf': {'type': 'arbitrary', 'samples': clear.tolist()},
        'long_readout_wf': {'type': 'constant', 'sample': long_readout_amp},
        'step_wf': {'type': 'constant', 'sample': scan_vpp / step_num},
        'bottom_wf': {'type': 'constant', 'sample': -scan_vpp/2},
        'bias_wf': {'type': 'constant', 'sample': bias_amp},
    },

    'digital_waveforms': {
        'ON': {'samples': [(1, 0)]},
    },

    'integration_weights': {
        'short_cosine_weights': {
            'cosine': [(1.0, short_readout_len)],
            'sine': [(0.0, short_readout_len)],
        },
        'short_sine_weights': {
            'cosine': [(0.0, short_readout_len)],
            'sine': [(1.0, short_readout_len)],
        },
        'short_minus_sine_weights': {
            'cosine': [(0.0, short_readout_len)],
            'sine': [(-1.0, short_readout_len)],
        },
        # 'cosine_weights': {
        #     'cosine': [(0.0,ringup),(1.0, kick_length-ringup),(1.0, readout_len)], 
        #     'sine': [(0.0, readout_len+kick_length)],
        # },
        'cosine_weights': {
            'cosine': [(1.0, readout_len+kick_length)], 
            'sine': [(0.0, readout_len+kick_length)],
        },
        # 'cont_cosine_weights': {
        #     'cosine': [(0.0, readout_len-2500), (1.0, 2500)],
        #     'sine': [(0.0, readout_len)],
        # },
        'sine_weights': {
            'cosine': [(0.0, readout_len+kick_length)],
            'sine': [(1.0, readout_len+ kick_length)],
        },
        # 'cont_sine_weights': {
        #     'cosine': [(0.0, readout_len)],
        #     'sine': [(0.0, readout_len-2500), (1.0, 2500)],
        # },
        'minus_sine_weights': {
            'cosine': [(0.0, readout_len+kick_length)],
            'sine': [(-1.0, readout_len+kick_length)],
        },
        # 'cont_minus_sine_weights': {
        #     'cosine': [(0.0, readout_len)],
        #     'sine': [(0, readout_len-2500), (-1.0, 2500)],
        # },
        'long_cosine_weights': {
            'cosine': [(1.0, long_readout_len)],
            'sine': [(0.0, long_readout_len)],
        },
        'long_sine_weights': {
            'cosine': [(0.0, long_readout_len)],
            'sine': [(1.0, long_readout_len)],
        },
        'long_minus_sine_weights': {
            'cosine': [(0.0, long_readout_len)],
            'sine': [(-1.0, long_readout_len)],
        },
    },
    'mixers': {
         "octave_octave1_1": [
                # {
                #     "intermediate_frequency": resonator_fake_if_freq,
                #     "lo_frequency": resonator_fake_lo_freq,
                #     "correction": (1, 0, 0, 1),
                # },
                {
                    "intermediate_frequency": resonator_if_freq,
                    "lo_frequency": resonator_lo_freq,
                    "correction": (1, 0, 0, 1),
                },
            ],
            "octave_octave1_2": [
                {
                    "intermediate_frequency": qubit_if_freq,
                    "lo_frequency": qubit_lo_freq,
                    "correction": (1, 0, 0, 1),
                },
                {
                    "intermediate_frequency": qubit_if_freq,
                    "lo_frequency": qubit_lo_freq,
                    "correction": (1, 0, 0, 1),
                },
                {
                    "intermediate_frequency": resonator_if_freq,
                    "lo_frequency": resonator_lo_freq,
                    "correction": (1, 0, 0, 1),
                },
            ],
            "octave_octave1_3": [
            {
                "intermediate_frequency": qubit_if_freq,
                "lo_frequency": qubit_lo_freq,
                "correction": (1, 0, 0, 1),
            },
        ],
            "mixer_resonator": [
            {
                "intermediate_frequency": resonator_if_freq,
                "lo_frequency": resonator_lo_freq,
                "correction": IQ_imbalance(0.058, -0.09),
            },
        ],
            "mixer_qubit": [
            {
                "intermediate_frequency": qubit_if_freq,
                "lo_frequency": qubit_lo_freq,
                "correction": (1, 0, 0, 1),
            },
        ],
        },
}
