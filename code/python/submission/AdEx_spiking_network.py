"""
Sparsely connected network of AdEx neurons.

Wonkwon Lee
"""
import brian2 as b2
from brian2 import *
from random import sample
from neurodynex.tools import plot_tools
from numpy import random
import matplotlib.pyplot as plt


# Default Parameters
SYNAPTIC_WEIGHT_W0 = 0.1 * b2.mV
RELATIVE_INHIBITORY_STRENGTH_G = 4.
CONNECTION_PROBABILITY_EPSILON = 0.1
SYNAPTIC_DELAY = 1.5 * b2.ms
# MEMBRANE_RESISTANCE_R = 500 * b2.Mohm
V_REST = -70.6 * b2.mV
RHEOBASE_THRESHOLD_v_rh = -50.4 * b2.mV
SHARPNESS_delta_T = 2.0 * b2.mV
FIRING_THRESHOLD_v_spike = RHEOBASE_THRESHOLD_v_rh + 5 * SHARPNESS_delta_T #-40.4 *b2.mV

CAPACITANCE_C = 281 * b2.pF
LEAK_CONDUCTANCE_gl = 30 * b2.nS
MEMBRANE_TIME_SCALE_tau_m = CAPACITANCE_C/LEAK_CONDUCTANCE_gl #9.36 *b2.ms

# Firing Patterns Parameters
ADAPTATION_TIME_CONSTANT_tau_w = 144.0 * b2.ms 
ADAPTATION_VOLTAGE_COUPLING_a = 4 * b2.nS
SPIKE_TRIGGERED_ADAPTATION_INCREMENT_b = 0.0805 * b2.nA
V_RESET = -70.6 * b2.mV

b2.defaultclock.dt = 0.05 * b2.ms

def simulate_AdEx_network(
        N_Excit=5000,
        N_Inhib=None,
        connection_probability=CONNECTION_PROBABILITY_EPSILON,
        w0=SYNAPTIC_WEIGHT_W0,
        g=RELATIVE_INHIBITORY_STRENGTH_G,
        synaptic_delay=SYNAPTIC_DELAY,

        v_rest=V_REST,
        v_reset=V_RESET,
        v_spike=FIRING_THRESHOLD_v_spike,
        v_rheobase=RHEOBASE_THRESHOLD_v_rh,

        C=CAPACITANCE_C,
        gl=LEAK_CONDUCTANCE_gl,
        delta_T=SHARPNESS_delta_T,
        tau_m=MEMBRANE_TIME_SCALE_tau_m,    
        tau_w=ADAPTATION_TIME_CONSTANT_tau_w,  
        a=ADAPTATION_VOLTAGE_COUPLING_a,
        b=SPIKE_TRIGGERED_ADAPTATION_INCREMENT_b,
        monitored_subset_size=200,
        random_vm_init=False,
        external_input=False,
        sim_time=200.*b2.ms # sim_time need to be more than 10 ms 
        ):

    if N_Inhib is None:
        N_Inhib = int(N_Excit/4)

    J_excit = w0
    J_inhib = -g*w0

    v_spike_str = "v>{:f}*mvolt".format(v_spike / b2.mvolt)
    
    adex_dynamics = """
        dv/dt = (gl*(v_rest - v) + gl*delta_T*exp((v - v_rheobase)/delta_T) + I - w)/C : volt
        dw/dt = (a*(v - v_rest) - w)/tau_w : amp
        I : amp
        """

    network = NeuronGroup(
        N_Excit+N_Inhib, model=adex_dynamics,
        threshold=v_spike_str, reset="v=v_reset;w+=b", method="euler")


    # Intialize each neuron with random membrane voltage between rest voltage and threshold
    if random_vm_init:
        network.v = random.uniform(v_rest/b2.mV, high=v_spike/b2.mV, size=(N_Excit+N_Inhib))*b2.mV
    else:
        network.v = v_rest #network.v = v_rest || network.v[0]

    # network.v = v_rest
    # network.w = 0.0 * b2.pA
    
    # Divide neuron network into excitatory and inhibitory populations    
    exc_neurons = network[:N_Excit]
    inhib_neurons = network[N_Excit:]

    exc_synapses = Synapses(exc_neurons, target=network, on_pre="v += J_excit", delay=synaptic_delay)
    exc_synapses.connect(p=connection_probability)

    inhib_synapses = Synapses(inhib_neurons, target=network, on_pre="v += J_inhib", delay=synaptic_delay)
    inhib_synapses.connect(p=connection_probability)

    if external_input is True:
        external_poisson_input = PoissonInput(target=network, target_var="v", N=1000, rate=13. * b2.Hz, weight=w0)   

    monitored_subset_size = min(monitored_subset_size, (N_Excit+N_Inhib))
    monitored_neurons_subset = sample(range(N_Excit+N_Inhib), monitored_subset_size)
    rate_monitor = PopulationRateMonitor(network)
    spike_monitor = SpikeMonitor(network)
    trace_monitor = StateMonitor(network, "v", record=monitored_neurons_subset)
    phase_monitor = StateMonitor(network, ["v","w"], record=monitored_neurons_subset)


    if external_input is True:
        b2.run(sim_time)
    # No input for initial fixed 10 ms and then after 10ms, network has input step current of 1 namp
    else:   
        b2.run(10*b2.ms)
        network.I = 1*b2.nA
        b2.run(sim_time - 10*b2.ms)
        network.I = 0*b2.nA

    return rate_monitor, spike_monitor, trace_monitor, phase_monitor, monitored_neurons_subset


def plot_adex_state(adex_state_monitor):
    plt.subplot(2, 2, 1)
    plt.plot(adex_state_monitor.t / b2.ms, adex_state_monitor.v[0] / b2.mV, lw=2)
    plt.xlabel("t [ms]")
    plt.ylabel("u [mV]")
    plt.title("Membrane potential")
    plt.subplot(2, 2, 2)
    plt.plot(adex_state_monitor.v[0] / b2.mV, adex_state_monitor.w[0] / b2.pA, lw=2)
    plt.xlabel("u [mV]")
    plt.ylabel("w [pAmp]")
    plt.title("Phase plane representation")
    plt.subplot(2, 2, 3)
    plt.plot(adex_state_monitor.t / b2.ms, adex_state_monitor.w[0] / b2.pA, lw=2)
    plt.xlabel("t [ms]")
    plt.ylabel("w [pAmp]")
    plt.title("Adaptation current")
    plt.show()


def getting_started():
    rate_monitor, spike_monitor, trace_monitor, phase_monitor, monitored_spike_idx = simulate_AdEx_network(N_Excit=1000)
    plot_tools.plot_network_activity(rate_monitor, spike_monitor, trace_monitor, spike_train_idx_list=monitored_spike_idx,
                                     t_min=0.*b2.ms)

    plt.show()
    plot_adex_state(phase_monitor)


if __name__ == "__main__":
    getting_started()