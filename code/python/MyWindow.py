"""
GUI Neuronal Network Simulator

Wonkwon Lee

Modified plot_tools from neurodynex and extended Brian2 to implement neurons and neuronal networks
"""
import sys

from PyQt5.QtWidgets import *
from PyQt5 import uic

import brian2 as b2
import matplotlib.pyplot as plt

import numpy
import random

import LIF_spiking_network as ln
import EIF_spiking_network as en
import AdEx_spiking_network as an


form_class = uic.loadUiType("final_gui.ui")[0]

class MyWindow(QMainWindow, form_class):
    def __init__(self, parent = None):
        super(MyWindow, self).__init__()
        self.setupUi(self)
        self.submit.clicked.connect(self.submit_btn_clicked)
        self.set.clicked.connect(self.set_btn_clicked)
        self.reset.clicked.connect(self.reset_btn_clicked)
 

    def submit_btn_clicked(self):
        # Simple Leaky Integrate-Fire Model
        if self.radio_1.isChecked():           
            rate_monitor, spike_monitor, trace_monitor, monitored_spike_idx = ln.simulate_lif_network( 
                N_Excit = int(self.N_excit.text()),
                N_Inhib = None if self.N_inhib.text() == '' else int(self.N_inhib.text()),
                connection_probability = float(self.cp.text()),
                w0 = float(self.w_excit.text()) * b2.mV,
                g = float(self.w_inhib.text()),
                synaptic_delay = float(self.delay.text())*b2.ms,
                membrane_time_scale = float(self.mtc.text())*b2.ms,
                v_rest = float(self.v_rest.text())*b2.mV,
                v_reset = float(self.v_reset.text())*b2.mV,
                v_spike = float(self.thres.text())*b2.mV,
                abs_refractory_period = float(self.refrac.text())*b2.ms,
                monitored_subset_size = int(self.sim_t_size.text()),
                random_vm_init = self.checkBox.isChecked(),
                external_input=self.checkBox_2.isChecked(),
                sim_time = float(self.sim_t.text())*b2.ms
            )

            # Draw Raster Plot
            self.widget_view1.canvas.ax.clear()
            all_spike_trains = spike_monitor.spike_trains()
            t_max = max(rate_monitor.t / b2.ms)
            t_min = 0.

            count = 0
            for i in monitored_spike_idx:
                ts = all_spike_trains[i]/b2.ms
                idx_spikes = (ts >= t_min) & (ts <= t_max)
                ts_spikes = ts[idx_spikes]      
                self.widget_view1.canvas.ax.scatter(ts_spikes, count * numpy.ones(ts_spikes.shape),
                                  marker=".", c="k", s=15, lw=0)
                count += 1

            # Highlight Raster        
            neuron_num = len(monitored_spike_idx)
            ts = trace_monitor.t/b2.ms
            idx_voltage = (ts >= t_min) & (ts <= t_max)
            fract = numpy.linspace(0, 1, 5)[1:-1]
            highlighted_neu = [int(neuron_num * v) for v in fract]

            for i in range(len(highlighted_neu)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                raster_index = highlighted_neu[i]
                population_index = monitored_spike_idx[raster_index]
                self.widget_view1.canvas.ax.scatter(
                    ts_spikes, raster_index * numpy.ones(ts_spikes.shape),
                    marker=".", c=color, s=100, lw=1)
            traces_i = highlighted_neu

            self.widget_view1.canvas.ax.set_ylim([0, count])
            # self.widget_view1.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view1.canvas.ax.set_ylabel("neuron #")
            # self.widget_view1.canvas.ax.set_title("Raster Plot")
            self.widget_view1.canvas.draw_idle()


            #Draw Voltage Trace
            self.widget_view2.canvas.ax.clear()
            for i in range(len(traces_i)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                population_index = monitored_spike_idx[traces_i[i]]
                self.widget_view2.canvas.ax.plot(
                    ts[idx_voltage], trace_monitor[population_index].v[idx_voltage]/b2.mV,
                    c=color, lw=1.)
            # self.widget_view2.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view2.canvas.ax.set_ylabel("V(t) [mV]")
            # self.widget_view2.canvas.ax.set_title("Voltage Trace")
            self.widget_view2.canvas.draw_idle()
        

        # Exponential Integrate-Fire Model
        elif self.radio_2.isChecked():
            rate_monitor, spike_monitor, trace_monitor, monitored_spike_idx = en.simulate_EIF_network( 
                N_Excit = int(self.N_excit.text()),
                N_Inhib = None if self.N_inhib.text() == '' else int(self.N_inhib.text()),
                connection_probability = float(self.cp.text()),
                w0 = float(self.w_excit.text()) * b2.mV,
                g = float(self.w_inhib.text()),
                synaptic_delay = float(self.delay.text())*b2.ms,
                v_rest = float(self.v_rest.text())*b2.mV,
                v_reset = float(self.v_reset.text())*b2.mV,
                v_spike = float(self.thres.text())*b2.mV,

                v_rheobase = float(self.v_rh.text())*b2.mV,
                delta_T = float(self.delta_t.text())*b2.mV,
                tau = float(self.mtc.text())*b2.ms,                   

                monitored_subset_size = int(self.sim_t_size.text()),
                random_vm_init = self.checkBox.isChecked(),
                external_input=self.checkBox_2.isChecked(),
                sim_time = float(self.sim_t.text())*b2.ms,
            )
            
            # Draw Raster Plot
            self.widget_view1.canvas.ax.clear()
            all_spike_trains = spike_monitor.spike_trains()
            t_max = max(rate_monitor.t / b2.ms)
            t_min = 0.

            count = 0
            for i in monitored_spike_idx:
                ts = all_spike_trains[i]/b2.ms
                idx_spikes = (ts >= t_min) & (ts <= t_max)
                ts_spikes = ts[idx_spikes]      
                self.widget_view1.canvas.ax.scatter(ts_spikes, count * numpy.ones(ts_spikes.shape),
                                  marker=".", c="k", s=15, lw=0)
                count += 1

            # Highlight Raster        
            neuron_num = len(monitored_spike_idx)
            ts = trace_monitor.t/b2.ms
            idx_voltage = (ts >= t_min) & (ts <= t_max)
            fract = numpy.linspace(0, 1, 5)[1:-1]
            highlighted_neu = [int(neuron_num * v) for v in fract]

            for i in range(len(highlighted_neu)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                raster_index = highlighted_neu[i]
                population_index = monitored_spike_idx[raster_index]
                self.widget_view1.canvas.ax.scatter(
                    ts_spikes, raster_index * numpy.ones(ts_spikes.shape),
                    marker=".", c=color, s=100, lw=1)
            traces_i = highlighted_neu

            self.widget_view1.canvas.ax.set_ylim([0, count])
            # self.widget_view1.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view1.canvas.ax.set_ylabel("neuron #")
            # self.widget_view1.canvas.ax.set_title("Raster Plot")
            self.widget_view1.canvas.draw_idle()


            # Draw Voltage Trace
            self.widget_view2.canvas.ax.clear()
            for i in range(len(traces_i)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                population_index = monitored_spike_idx[traces_i[i]]
                self.widget_view2.canvas.ax.plot(
                    ts[idx_voltage], trace_monitor[population_index].v[idx_voltage]/b2.mV,
                    c=color, lw=1.)
            # self.widget_view2.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view2.canvas.ax.set_ylabel("V(t) [mV]")
            # self.widget_view2.canvas.ax.set_title("Voltage Traces")
            self.widget_view2.canvas.draw_idle()


        # Exponential Integrate-Fire Model with Adaptation Current
        elif self.radio_3.isChecked():
            rate_monitor, spike_monitor, trace_monitor, phase_monitor, monitored_spike_idx = an.simulate_AdEx_network(
                N_Excit = int(self.N_excit.text()),
                N_Inhib = None if self.N_inhib.text() == '' else int(self.N_inhib.text()),
                connection_probability = float(self.cp.text()),
                w0 = float(self.w_excit.text()) * b2.mV,
                g = float(self.w_inhib.text()),
                synaptic_delay = float(self.delay.text())*b2.ms,

                v_rest = float(self.v_rest.text())*b2.mV,
                v_reset = float(self.v_reset.text())*b2.mV,
                v_spike = float(self.thres.text())*b2.mV,
                v_rheobase = float(self.v_rh.text())*b2.mV,
                
                C=float(self.cap_c.text())*b2.pF,
                gl=float(self.leak_gl.text())*b2.nS,
                delta_T = float(self.delta_t.text())*b2.mV,
                tau_m = float(self.mtc.text())*b2.ms,
                tau_w = float(self.atc.text())*b2.ms,
                a = float(self.avc.text())*b2.nS,
                b = float(self.sta.text())*b2.nA,

                monitored_subset_size = int(self.sim_t_size.text()),
                random_vm_init = self.checkBox.isChecked(),
                external_input=self.checkBox_2.isChecked(),
                sim_time = float(self.sim_t.text())*b2.ms
            )

            # Draw Raster Plot
            self.widget_view1.canvas.ax.clear()
            all_spike_trains = spike_monitor.spike_trains()
            t_max = max(rate_monitor.t / b2.ms)
            t_min = 0.

            count = 0
            for i in monitored_spike_idx:
                ts = all_spike_trains[i]/b2.ms
                idx_spikes = (ts >= t_min) & (ts <= t_max)
                ts_spikes = ts[idx_spikes]      
                self.widget_view1.canvas.ax.scatter(ts_spikes, count * numpy.ones(ts_spikes.shape),
                                  marker=".", c="k", s=15, lw=0)
                count += 1

            # Highlight Raster
            neuron_num = len(monitored_spike_idx)
            ts = trace_monitor.t/b2.ms
            idx_voltage = (ts >= t_min) & (ts <= t_max)
            fract = numpy.linspace(0, 1, 5)[1:-1]
            highlighted_neu = [int(neuron_num * v) for v in fract]

            for i in range(len(highlighted_neu)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                raster_index = highlighted_neu[i]
                population_index = monitored_spike_idx[raster_index]
                self.widget_view1.canvas.ax.scatter(
                    ts_spikes, raster_index * numpy.ones(ts_spikes.shape),
                    marker=".", c=color, s=100, lw=1)
            traces_i = highlighted_neu

            self.widget_view1.canvas.ax.set_ylim([0, count])
            # self.widget_view1.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view1.canvas.ax.set_ylabel("neuron #")
            # self.widget_view1.canvas.ax.set_title("Raster Plot")
            self.widget_view1.canvas.draw_idle()


            # Draw Voltage Trace
            self.widget_view2.canvas.ax.clear()
            for i in range(len(traces_i)):
                if i == 0:
                    color = "b"
                elif i == 1:
                    color = "r"
                else:
                    color = "k"
                population_index = monitored_spike_idx[traces_i[i]]
                self.widget_view2.canvas.ax.plot(
                    ts[idx_voltage], trace_monitor[population_index].v[idx_voltage]/b2.mV,
                    c=color, lw=1.)
            # self.widget_view2.canvas.ax.set_xlabel("time [ms]")
            # self.widget_view2.canvas.ax.set_ylabel("V(t) [mV]")
            # self.widget_view2.canvas.ax.set_title("Voltage Traces")
            self.widget_view2.canvas.draw_idle()


            # Plot Phase Plane Representation
            self.widget_view3.canvas.ax.clear()
            self.widget_view3.canvas.ax.plot(phase_monitor.v[0] / b2.mV, phase_monitor.w[0] / b2.pA, lw=1)
            # self.widget_view3.canvas.ax.set_xlabel("u [mV]")
            # self.widget_view3.canvas.ax.set_ylabel("w [pAmp]")
            self.widget_view3.canvas.draw_idle()
            
            self.widget_view4.canvas.ax.clear()
            self.widget_view4.canvas.ax.plot(phase_monitor.t / b2.ms, phase_monitor.w[0] / b2.pA, lw=1)
            # self.widget_view4.canvas.ax.set_xlabel("t [ms]")
            # self.widget_view4.canvas.ax.set_ylabel("w [pAmp]")
            self.widget_view4.canvas.draw_idle()
           

        else:
            QMessageBox.about(self, "message", "Choose a neuron model before submit")



    def set_btn_clicked(self):
        if self.radio_1.isChecked():
            self.N_excit.setText('2000')
            self.N_inhib.setText('500')
            self.cp.setText('0.1')
            self.w_excit.setText('0.1') # w0
            self.w_inhib.setText('4')   # g
            self.delay.setText('1.5')
            self.mtc.setText('20')
            self.v_rest.setText('0')
            self.v_reset.setText('10')
            self.thres.setText('20')
            self.refrac.setText('2.0')
            self.sim_t.setText('200')
            self.sim_t_size.setText('100')
            self.checkBox.setCheckState(False)
            self.checkBox_2.setCheckState(True)


        elif self.radio_2.isChecked():
            self.N_excit.setText('2000')
            self.N_inhib.setText('500')
            self.cp.setText('0.1')
            self.w_excit.setText('0.1')
            self.w_inhib.setText('4')
            self.delay.setText('1.5')
            self.mtc.setText('12')
            self.v_rest.setText('-65')
            self.v_reset.setText('-60')
            self.thres.setText('-30')
            self.sim_t.setText('200')
            self.sim_t_size.setText('100')
            self.checkBox.setCheckState(False)
            self.checkBox_2.setCheckState(False)

            # EIF additional parameters
            self.v_rh.setText('-55')
            self.delta_t.setText('2')

            ## Not in EIF
            self.refrac.setText('')


        elif self.radio_3.isChecked():
            def __init_adex__():
                self.N_excit.setText('2000')
                self.N_inhib.setText('500')
                self.cp.setText('0.1')
                self.w_excit.setText('0.1')
                self.w_inhib.setText('4')
                self.delay.setText('1.5')
                self.v_rest.setText('-70.6')
                self.v_reset.setText('-70.6')
                self.thres.setText('-40.4')
                self.sim_t.setText('200')
                self.sim_t_size.setText('100')
                self.checkBox.setCheckState(False)
                self.checkBox_2.setCheckState(False)
                self.mtc.setText('9.36')           
                self.cap_c.setText('281')
                self.leak_gl.setText('30')
                self.v_rh.setText('-50.4')
                self.delta_t.setText('2')
                self.refrac.setText('')
            
            # Regular Spiking, Burst, Transient
            if self.radio_a1.isChecked():
                __init_adex__()
                self.atc.setText('144')         # tauw
                self.avc.setText('4')           # a
                self.sta.setText('0.0805')      # b
                self.v_reset.setText('-70.6')   # Vr

            elif self.radio_a2.isChecked():
                __init_adex__()
                self.atc.setText('20')          # tauw
                self.avc.setText('4')           # a
                self.sta.setText('0.5')         # b
                self.v_reset.setText('-45.4')   # Vr

            elif self.radio_a3.isChecked():
                __init_adex__()
                self.atc.setText('144')         # tauw
                self.avc.setText('3.903')       # a
                self.sta.setText('0')           # b
                self.v_reset.setText('-70.6')   # Vr
                
            else:
                  QMessageBox.about(self, "Error", "Select a firing pattern for AdEx model")
            

        else:
            QMessageBox.about(self, "Error", "Select a network model")



    def reset_btn_clicked(self):
        # Empty all input parameters and toggle off check boxes
        self.N_excit.setText('')
        self.N_inhib.setText('')     
        self.cp.setText('')
        self.w_excit.setText('')
        self.w_inhib.setText('')
        self.delay.setText('')
        self.v_rest.setText('')
        self.v_reset.setText('')
        self.thres.setText('')
        self.sim_t.setText('')
        self.sim_t_size.setText('')
        self.mtc.setText('')            
        self.v_rh.setText('')
        self.delta_t.setText('')
        self.cap_c.setText('')
        self.leak_gl.setText('')
        self.atc.setText('')
        self.avc.setText('')
        self.sta.setText('')
        self.refrac.setText('')
        self.checkBox.setCheckState(False)
        self.checkBox_2.setCheckState(False)

        self.widget_view1.canvas.ax.clear()
        self.widget_view2.canvas.ax.clear()
        self.widget_view3.canvas.ax.clear()
        self.widget_view4.canvas.ax.clear()
        self.widget_view1.canvas.draw_idle()
        self.widget_view2.canvas.draw_idle()
        self.widget_view3.canvas.draw_idle()
        self.widget_view4.canvas.draw_idle()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    myWindow = MyWindow()
    myWindow.show()
    app.exec_()