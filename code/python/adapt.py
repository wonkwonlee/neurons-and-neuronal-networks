"""
Simple AdEx Neuron Network

Wonkwon Lee

Brette R. & Gerstner W. Model (2005)

Synapses and connection probabilty are not implemented
"""
from brian2 import *


# Default Parameters
C = 281 * pF
gL = 30 * nS
taum = C / gL
EL = -70.6 * mV
VT = -50.4 * mV
delta_t = 2 * mV
Vcut = VT + 5 * delta_t

# Pick an electrophysiological behaviour
tauw, a, b, Vr = 144*ms, 4*nS, 0.0805*nA, -70.6*mV # Regular spiking
# tauw,a,b,Vr=20*ms,4*nS,0.5*nA,VT+5*mV # Bursting
# tauw,a,b,Vr=144*ms,2*C/(144*ms),0*nA,-70.6*mV # Fast spiking

eqs = """
dv/dt = (gL*(EL - v) + gL*delta_t*exp((v - VT)/delta_t) + I - w)/C : volt
dw/dt = (a*(v - EL) - w)/tauw : amp
I : amp
"""

neuron = NeuronGroup(1000, model=eqs, threshold='v>Vcut', reset="v=Vr; w+=b", method='euler')
neuron.v = EL
trace = StateMonitor(neuron, 'v', record=0)
spikes = SpikeMonitor(neuron)
raster = PopulationRateMonitor(neuron)

run(20 * ms)
neuron.I = 1*nA
run(100 * ms)
neuron.I = 0*nA
run(20 * ms)

v = trace[0].v[:]
for t in spikes.t:
    i = int(t / defaultclock.dt)
    v[i] = 20*mV

plot(trace.t / ms, v / mV)
xlabel('time (ms)')
ylabel('membrane potential (mV)')
show()