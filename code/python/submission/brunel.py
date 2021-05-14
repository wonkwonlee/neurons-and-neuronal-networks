"""
Simple Brunel Network with identical LIF neurons 

Wonkwon Lee

Brunel & Hakim Model (1999)

Excitatory and Inhibotry populations are not implemented
"""
from brian2 import *


N = 100
Vr = 10*mV
theta = 20*mV
tau = 20*ms
delta = 2*ms
taurefr = 2*ms
sim_time = .1*second
C = 1000
cp = float(C)/N     #connection probability    
J = .1*mV
muext = 25*mV
sigmaext = 1*mV

eqs = """
dv/dt = (-v+muext + sigmaext * sqrt(tau) * xi)/tau : volt
"""

network = NeuronGroup(N, eqs, threshold='v>theta', reset='v=Vr', refractory=taurefr, method='euler')
network.v = Vr
conn = Synapses(network, network, on_pre='v += -J', delay=delta)
conn.connect(p=cp)
spike_monitor = SpikeMonitor(network)
population_monitor = PopulationRateMonitor(network)

run(sim_time)

subplot(211)
plot(spike_monitor.t/ms, spike_monitor.i, '.')
xlim(0, sim_time/ms)

subplot(212)
plot(population_monitor.t/ms, population_monitor.smooth_rate(window='flat', width=0.5*ms)/Hz)
xlim(0, sim_time/ms)

show()