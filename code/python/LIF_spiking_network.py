"""
Implementation of the Brunel 2000 network:
Sparsely connected network of LIF neurons.

Referred to Wulfram Gerstner, Werner M. Kistler, Richard Naud, and Liam Paninski
Neuronal Dynamics (see http://neuronaldynamics.epfl.ch)

Extended existing funtionalities and modified to fit into GUI simulator

Modified by Wonkwon Lee
"""
import brian2 as b2
from brian2 import *
from random import sample
from neurodynex.tools import plot_tools
from numpy import random
import matplotlib.pyplot as plt


# Default Parameters
V_REST = 0. * b2.mV
V_RESET = +10. * b2.mV
FIRING_THRESHOLD = +20. * b2.mV
MEMBRANE_TIME_SCALE = 20. * b2.ms
ABSOLUTE_REFRACTORY_PERIOD = 2.0 * b2.ms
SYNAPTIC_WEIGHT_W0 = 0.1 * b2.mV

RELATIVE_INHIBITORY_STRENGTH_G = 4. 
CONNECTION_PROBABILITY_EPSILON = 0.1
SYNAPTIC_DELAY = 1.5 * b2.ms

b2.defaultclock.dt = 0.05 * b2.ms


def simulate_lif_network(
        N_Excit=5000,
        N_Inhib=None,
        connection_probability=CONNECTION_PROBABILITY_EPSILON,
        w0=SYNAPTIC_WEIGHT_W0,
        g=RELATIVE_INHIBITORY_STRENGTH_G,
        synaptic_delay=SYNAPTIC_DELAY,
        v_rest=V_REST,
        v_reset=V_RESET,
        v_spike=FIRING_THRESHOLD,
        membrane_time_scale=MEMBRANE_TIME_SCALE,   
        abs_refractory_period=ABSOLUTE_REFRACTORY_PERIOD,
        monitored_subset_size=200,
        random_vm_init=False,
        external_input=True,
        sim_time=200.*b2.ms):

    if N_Inhib is None:
        N_Inhib = int(N_Excit/4)

    J_excit = w0
    J_inhib = -g*w0

    lif_dynamics = """
    dv/dt =  (-(v-v_rest) + R * I) / membrane_time_scale : volt (unless refractory)
    I : amp
    """
    R = 10. * b2.Mohm
    
    network = NeuronGroup(
        N_Excit+N_Inhib, model=lif_dynamics,
        threshold="v>v_spike", reset="v=v_reset", refractory=abs_refractory_period,
        method="linear")

    # Intialize each neuron with random membrane voltage between rest voltage and threshold
    if random_vm_init:
        network.v = random.uniform(v_rest/b2.mV, high=v_spike/b2.mV, size=(N_Excit+N_Inhib))*b2.mV
    else:
        network.v = v_rest

    # Divide neuron network into excitatory and inhibitory populations
    excitatory_population = network[:N_Excit]
    inhibitory_population = network[N_Excit:]

    exc_synapses = Synapses(excitatory_population, target=network, on_pre="v += J_excit", delay=synaptic_delay)
    exc_synapses.connect(p=connection_probability)

    inhib_synapses = Synapses(inhibitory_population, target=network, on_pre="v += J_inhib", delay=synaptic_delay)
    inhib_synapses.connect(p=connection_probability)

    # External Neuron Stimuli, imported from Brain2.PoissonInput
    # N : number of external neurons, rate : rate of each inputs, weight : external input strength
    if external_input is True:
        external_poisson_input = PoissonInput(target=network, target_var="v", N=1000, rate=13. * b2.Hz, weight=w0)    

    monitored_subset_size = min(monitored_subset_size, (N_Excit+N_Inhib))
    idx_monitored_neurons = sample(range(N_Excit+N_Inhib), monitored_subset_size)
    rate_monitor = PopulationRateMonitor(network)
    spike_monitor = SpikeMonitor(network, record=idx_monitored_neurons)
    trace_monitor = StateMonitor(network, "v", record=idx_monitored_neurons)

    if external_input is True:
        b2.run(sim_time)
    else:
        b2.run(10*b2.ms)
        network.I = 2.5*b2.nA
        b2.run(sim_time - 10*b2.ms)
        network.I = 0*b2.nA

    return rate_monitor, spike_monitor, trace_monitor, idx_monitored_neurons


def getting_started():
    rate_monitor, spike_monitor, trace_monitor, monitored_spike_idx= simulate_lif_network(N_Excit=2000, external_input=False)
    plot_tools.plot_network_activity(rate_monitor, spike_monitor, trace_monitor, 
        spike_train_idx_list=monitored_spike_idx, t_min=0.*b2.ms)
    plt.show()


if __name__ == "__main__":
    getting_started()
