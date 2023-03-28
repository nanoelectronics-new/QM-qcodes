from qcodes.utils.validators import Arrays
from qcodes.instrument_drivers.OPX.opx_driver import *
from qm.qua import *
from scipy import signal
from qm import generate_qua_script
# noinspection PyAbstractClass

class IQArray(MultiParameter):
    def __init__(self,
        name: str,
        instrument: "OPXRawADC",
        npts: int,):
        # names, labels, and units are the same

        super().__init__('IQ', names=('I', 'Q'), shapes=((npts,), (npts,)),
                         labels=('I', 'Q'),
                         units=('mV', 'mV'),
                         instrument=instrument,
                         # note that EACH item needs a sequence of setpoint arrays
                         # so a 1D item has its setpoints wrapped in a length-1 tuple
                         setpoints = ((instrument.time_axis(),),(instrument.time_axis(),)),
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
        self.set_sweep(npts)
#

    def set_sweep(self, npts: int):
        # Needed to update config of the software parameter on sweep change
        # frequency setpoints tuple as needs to be hashable for look up.
        f = tuple(np.linspace(1,self.instrument.readout_pulse_length()+2*self.instrument.smearing(), num=self.instrument.readout_pulse_length()+
                                                                                                                      2*self.instrument.smearing()))
        self.setpoints = ((f,), (f,))
        self.shapes = ((self.instrument.readout_pulse_length(),), (self.instrument.readout_pulse_length(),))

    def get_raw(self):
        I,Q = self.instrument.get_res()
        return I,Q



class OPXRawADC(OPX):
    def __init__(self, config: Dict, name: str = "OPXRawADC", host=None, port=None,host_octave=None, port_octave=None, **kwargs):
        super().__init__(config, name, host=host, port=port,host_octave=host_octave, port_octave=port_octave, **kwargs)

        self.add_parameter(
            "t_meas",
            unit="s",
            initial_value=0.01,
            vals=Numbers(0, 1),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "t_begin",
            unit="ns",
            initial_value=1,
            vals=Numbers(0, 10),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "readout_pulse_length",
            unit="ns",
            initial_value = 1000,
            vals=Numbers(16, 1e7),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "smearing",
            unit="ns",
            initial_value = 20,
            vals=Numbers(0, 2000),
            get_cmd=None,
            set_cmd=None,
        )
        self.add_parameter(
            "time_axis",
            unit="s",
            label="Time Axis",
            parameter_class=GeneratedSetPoints,
            startparam=self.t_begin,
            stopparam=self.readout_pulse_length,
            numpointsparam=self.readout_pulse_length,
            vals=Arrays(shape=(self.readout_pulse_length,)),
        )
        self.add_parameter(
            "amp",
            unit="V",
            initial_value=0.1,
            vals=Numbers(-0.5, 0.5),
            get_cmd=None,
            set_cmd=None,
        )

        self.add_parameter(
            "trace_iq",
            npts=self.readout_pulse_length(),
            parameter_class=IQArray,
        )

    def get_prog(self):

        n_avg = round(self.t_meas() * 1e9 / self.readout_pulse_length())
        
        if self.octave:
            lo =self.config['mixers']['octave_octave1_1'][0]['lo_frequency']
            self.qm.octave.set_lo_source('resonator', OctaveLOSource.Internal)
            self.qm.octave.set_lo_frequency('resonator', lo) # assign the LO inside the octave to element
            self.qm.octave.set_rf_output_gain('resonator', 0) # can set gain from -10dB to 20dB
            self.qm.octave.set_rf_output_mode('resonator', RFOutputMode.trig_normal) # set the behaviour of the RF switch
            self.qm.octave.calibrate_element('resonator',[(lo, self.config['elements']['resonator']['intermediate_frequency'])]) # can provide many pairs of LO & IFs.
            self.qm =self.qmm.open_qm(self.config)
            self.qm.octave.set_qua_element_octave_rf_in_port('resonator',"octave1", 1)
            self.qm.octave.set_downconversion('resonator',lo_source=RFInputLOSource.Internal)
        else:
            lo =self.config['mixers']['mixer_resonator'][0]['lo_frequency']
            

        self.parameters['trace_iq'].set_sweep(self.readout_pulse_length())
        with program() as prog:
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()
            adc_st = declare_stream(adc_trace=True)
            # with infinite_loop_():
            #     play('readout', 'resonator')

            with for_(n, 0, n < n_avg, n + 1):
                reset_phase("resonator")
                measure("readout", "resonator", adc_st)
                wait(10000 // 4, 'resonator')
       
                # save(I, I_st)
                # save(Q, Q_st)

  
            with stream_processing():
                adc_st.input1().average().save("I")
                adc_st.input2().average().save("Q")
                
        # sourceFile = open('debug.py', 'w')
        # print(generate_qua_script(prog, self.config), file=sourceFile)
        # sourceFile.close()

        return prog
    
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
            self.result_handles.wait_for_all_values()
            print(self.config['pulses']['readout_pulse']['length'])
            I = 1e3*self.result_handles.get("I").fetch_all()/4096
            Q = 1e3*self.result_handles.get("Q").fetch_all()/4096
            return I,Q
