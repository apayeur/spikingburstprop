from ebsynapse import EBSynapse
from neuronstate import NeuronStateWithMovingAverage
from neuronstate import NeuronStateWithTwoMovingAverages
import numpy as np

class PreEventSynapse(EBSynapse):
    """
    EB synapse that depends on presynaptic events only (i.e., no homeostatic plasticity)
    """
    def __init__(self, learning_rate, A_plus=1, A_minus=-0.15, burst_def=0.016, tau_trace=0.020):
        super().__init__(learning_rate, burst_def, tau_trace)
        self.A_plus = A_plus
        self.A_minus = A_minus

    def __str__(self):
        return str(self.weight)

    def rule(self, t_int):
        event_trace_pre = self.pre.trace['event']
        event_state_post = self.post.train['event'][t_int]
        burst_state_post = self.post.train['burst'][t_int]

        # note that learning rate is actually rescaled here because
        # event_state_post and burst_state_post should be divided by the size of the timestep
        return self.learning_rate*(self.A_minus*event_state_post + self.A_plus*burst_state_post)*event_trace_pre


class AdaptivePreEventSynapse(EBSynapse):
    """
    EB synapse that depends on presynaptic events only (i.e., no homeostatic plasticity).
    Compared to PreEventSynapse, here A_plus is replaced by avg(E) and A_minus by -avg(B).
    """
    def __init__(self, learning_rate, burst_def=0.016, tau_trace=0.020, tau_ma=10., starting_estimate=(5., 5.*0.15)):
        super().__init__(learning_rate, burst_def, tau_trace)
        self.post_ma = NeuronStateWithMovingAverage(burst_def=burst_def, tau_trace=tau_trace, tau_ma=tau_ma,
                                                    starting_estimate=starting_estimate)
        self.plasticity_on = True

    def __str__(self):
        return str(self.weight)

    def rule(self, t_int):
        event_trace_pre = self.pre.trace['event']
        event_state_post = self.post_ma.train['event'][t_int]
        burst_state_post = self.post_ma.train['burst'][t_int]
        burst_moving_avg = self.post_ma.moving_average['burst']
        event_moving_avg = self.post_ma.moving_average['event']

        # note that learning rate is actually rescaled here because
        # event_state_post and burst_state_post should be divided by the size of the timestep
        return self.learning_rate*(-burst_moving_avg*event_state_post + event_moving_avg*burst_state_post)*event_trace_pre

    def evolve(self, t_int, timestep):
        """
        Value of weight at next time step
        """
        self.pre.evolve_traces(t_int, timestep)
        self.post_ma.evolve_traces(t_int, timestep)
        if self.plasticity_on:
            self.weight += self.rule(t_int)*timestep

    def reset(self):
        self.weight = 0.
        self.pre.reset()
        self.post_ma.reset()


class AdaptivePreEventSynapseTwoTimeConstants(EBSynapse):
    """
    EB synapse that depends on presynaptic events only (i.e., no homeostatic plasticity).
    Compared to PreEventSynapse, here A_plus is replaced by avg(E) and A_minus by -avg(B).
    The time constants for the moving averages are different for E and B.
    """
    def __init__(self, learning_rate, burst_def=0.016, tau_trace=0.020, alpha=10., gamma=10.):
        super().__init__(learning_rate, burst_def, tau_trace)
        self.post_ma = NeuronStateWithTwoMovingAverages(burst_def=burst_def, tau_trace=tau_trace, alpha=alpha, gamma=gamma)
        self.plasticity_on = True

    def __str__(self):
        return str(self.weight)

    def rule(self, t_int):
        event_trace_pre = self.pre.trace['event']
        event_state_post = self.post_ma.train['event'][t_int]
        burst_state_post = self.post_ma.train['burst'][t_int]
        burst_moving_avg = self.post_ma.moving_average['burst']
        event_moving_avg = self.post_ma.moving_average['event']

        # note that learning rate is actually rescaled here because
        # event_state_post and burst_state_post should be divided by the size of the timestep
        return self.learning_rate*(-burst_moving_avg*event_state_post + event_moving_avg*burst_state_post)*event_trace_pre

    def evolve(self, t_int, timestep):
        """
        Value of weight at next time step
        """
        self.pre.evolve_traces(t_int, timestep)
        self.post_ma.evolve_traces(t_int, timestep)
        if self.plasticity_on:
            self.weight += self.rule(t_int)*timestep

    def reset(self):
        self.weight = 0.
        self.pre.reset()
        self.post_ma.reset()