from ebsynapse import EBSynapse


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




