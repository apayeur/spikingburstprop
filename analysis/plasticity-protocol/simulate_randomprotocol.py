from pairingprotocols import RandomProtocol
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt


def simulate(rates, eta, duration, alpha, burst_threshold, tau_pre, starting_estimate=(5., 0.35*5.)):
    """Plasticity protocol consisting of pre and postsynaptic Poisson processes."""
    # Parameters
    nb_reals = 20       # number of realizations
    dt = 0.001          # integration time step
    np.random.seed(3)   # for reproducibility

    # Create synapse
    syn = AdaptivePreEventSynapse(eta,
                                  burst_def=burst_threshold,
                                  tau_trace=tau_pre,
                                  tau_ma=alpha,
                                  starting_estimate=starting_estimate)

    # Main simulation
    weights = np.zeros((nb_reals, len(rates)))
    for n in range(nb_reals):
        if n % 10 == 0:
            print('realization #{}'.format(n+1))
        for idx, r in enumerate(rates):
            syn.reset()
            protocol = RandomProtocol(duration=duration, rate=r)
            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration)

            count = 0
            for x in range(int(duration/dt)):
                count += 1
                syn.evolve(x, dt)
            weights[n, idx] += syn.weight/dt
            # division by dt is because event and burst trains are not divided
            # by dt in the rule (see preeventsynapse.py)

    return np.mean(weights, axis=0), np.std(weights, axis=0)


def analytical(rates, eta, duration, alpha, burst_threshold, tau_pre, starting_estimates=(5, 0.35*5), tau_ref=0.):
    """Analytical result for the weight change."""
    r = rates
    lambda_ = r/(1. - tau_ref*r)
    P_greater = np.exp(-lambda_ * (burst_threshold - tau_ref))
    E = r * P_greater
    B = E * (1 - P_greater)
    E_pre = E
    integrated_Pbar = (B/E)*(aux_fun(duration, starting_estimates[1]/B - 1., starting_estimates[0]/E - 1., alpha)
                             - aux_fun(0, starting_estimates[1]/B - 1., starting_estimates[0]/E - 1., alpha))
    w = eta*tau_pre*E_pre*(B*duration - E*integrated_Pbar)
    return r, w


def aux_fun(x, A, B, alpha):
    return (alpha*(B-A)*np.log(B + np.exp(x/alpha)) + A*x)/B


def plot_cov(synapse, rate, timestep, burst_threshold):
    k, cov_be, cov_eb = synapse.post_ma.compute_eb_corr(timestep)
    tau = np.arange(-0.100, 0.100, 0.0001)
    c_anal = corr_anal(tau, rate, burst_threshold)
    L = int((len(k) + 1) / 2.)
    plt.plot(1000 * timestep * k[L - 1 - int(10. / timestep):L - 1 + int(10. / timestep)],
             cov_be[L - 1 - int(10. / timestep):L - 1 + int(10. / timestep)], label=r'$\langle E(t) B(t+\tau) \rangle$')
    plt.plot(1000 * timestep * k[L - 1 - int(10. / timestep):L - 1 + int(10. / timestep)],
             cov_eb[L - 1 - int(10. / timestep):L - 1 + int(10. / timestep)], label=r'$\langle B(t) E(t+\tau) \rangle$')
    plt.plot(1000*tau, c_anal, label='anal')
    plt.xlim([-100, 100])
    plt.ylabel("Cross-covariance")
    plt.xlabel("Lag [ms]")
    plt.legend()
    plt.tight_layout()
    plt.show()


def corr_anal(tau, rate, burst_threshold):
    P = rate*np.exp(-tau*rate)
    P_greater = np.exp(-rate*burst_threshold)
    E = rate*P_greater
    B = E*(1-P_greater)
    return (tau>0)*(burst_threshold > tau)*E*P + \
           (tau>burst_threshold)*(tau<2*burst_threshold)*E**2*(1-np.exp(rate*(burst_threshold-tau))) + \
           (tau > 2 * burst_threshold) *E*B + \
           (tau<0)*(tau<-burst_threshold)*B*E + (tau<0)*(tau>-burst_threshold)*0.


def burst_fraction(rates, burst_threshold):
    return 1-np.exp(-burst_threshold*rates)


if __name__ == '__main__':
    import sys
    import seaborn as sns
    import os
    sys.path.append('../')

    PLOTPATH = os.path.join('..', '..', 'results', 'learning-rule', 'test')
    if not os.path.exists(PLOTPATH):
        os.makedirs(PLOTPATH)

    plt.style.use('../thesis_mplrc.dms')

    rates = np.linspace(1., 50., 10)
    duration = 100.
    pbar0 = 0.3
    alpha = 15

    for initial_er in [5.]:
        print('simulate initial ER = {}....'.format(initial_er))
        r, w, _, _ = simulate(rates, 0.1, duration, alpha, 0.016, 0.050, starting_estimate=(initial_er, pbar0*initial_er))
        r, w5, _, _ = simulate(rates, 0.1, duration, alpha, 0.016, 0.050, starting_estimate=(initial_er, 0.5*initial_er))

        plt.figure(figsize=(3, 3/1.6))
        plt.plot(r, w, color='black')
        plt.plot(r, w5, color='grey')
        sns.despine()
        plt.plot(r, np.zeros(r.shape), 'k--', lw=1)
        plt.xlabel('Rate [Hz]')
        plt.ylabel('$\Delta W$')
        plt.tight_layout()
        plt.savefig(PLOTPATH + '/DeltaW_AdaptiveLearningRule_InitialER'+str(initial_er)+'.pdf')
        plt.close()


    # testing analytical result
    #rates = np.arange(1, 50)
    plt.figure(figsize=(3, 2))
    r, w3_notauref = analytical(rates, 0.1, 60, 15, 0.016, 0.050, starting_estimates=(5, 0.3 * 5),
                                tau_ref=0.000)
    r, w3 = analytical(rates, 0.1, 60, 15, 0.016, 0.050, starting_estimates=(5, 0.3 * 5), tau_ref=0.002)
    r, w_longref = analytical(rates, 0.1, 60, 15, 0.016, 0.050, starting_estimates=(5, 0.3 * 5), tau_ref=0.004)
    plt.plot(r, w_longref, 'r', label=r'$\tau_\mathrm{ref} = 4$ ms')
    plt.plot(r, w3, 'k', label=r'$\tau_\mathrm{ref} = 2$ ms')
    plt.plot(r, w3_notauref, 'k--', label=r'$\tau_\mathrm{ref} = 0$ ms')
    plt.plot(r, w, '-', color='grey')
    plt.plot(r, np.zeros(r.shape), ':', color='grey', lw=0.5)
    plt.xlabel('Rate [Hz]')
    plt.ylabel('$\Delta W$')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(PLOTPATH + '/Analytical.png')
    plt.close()
