from abc import ABCMeta, abstractmethod
from neuronstate import NeuronState


class EBSynapse(object):
    """
    Abstract base class for a synapse whose plasticity depends on postsynaptic
    events and bursts.
    """

    __metaclass__ = ABCMeta

    def __init__(self, learning_rate, burst_def=0.016, tau_trace=0.020):
        self.weight = 0.
        self.learning_rate = learning_rate
        self.pre = NeuronState(burst_def, tau_trace)
        self.post = NeuronState(burst_def, tau_trace)

    @abstractmethod
    def rule(self, t_int):
        """
        Synaptic plasticity rule.
        :return: Evaluate term on the right-hand side of dw/dt
        """
        pass

    def evolve(self, t_int, timestep):
        """
        Value of weight at next time step
        """
        self.weight += self.rule(t_int)*timestep
        self.pre.evolve_traces(t_int, timestep)
        self.post.evolve_traces(t_int, timestep)

    def reset(self):
        self.weight = 0.
        self.pre.reset()
        self.post.reset()



