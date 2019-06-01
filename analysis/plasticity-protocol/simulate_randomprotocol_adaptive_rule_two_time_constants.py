from pairingprotocols import RandomProtocol
from preeventsynapse import AdaptivePreEventSynapseTwoTimeConstants
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson processes.
"""

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
    #plt.ylim([np.min(cov_be), 0.0005])
    plt.ylabel("Cross-covariance")
    plt.xlabel("Lag [ms]")
    plt.legend()
    plt.tight_layout()
    plt.show()


def simulate(rates, eta, duration, alpha, burst_threshold, tau_pre):
    np.random.seed(2)

    # parameters for pairing protocol
    tau_moving_average = alpha         # time constant for moving average
    burnin = 3*tau_moving_average      # seconds, relaxation of moving averages
    nb_reals = 1                      # number of realizations

    # learning rate
    learning_rate = eta

    # integration time step
    dt = 0.0005
    sampling_interval = 0.001

    weights = np.zeros(rates.shape)

    # Compute gammas
    gammas = []
    for r in rates:
        gammas.append(get_gamma(r, alpha, burst_threshold))

    # Simulation
    for n in range(nb_reals):
        if n % 1 == 0:
            print('realization #{}'.format(n+1))
        for idx, r in enumerate(rates):
            bp_est = []
            er_est = []
            br_est = []
            weight_trace = []
            protocol = RandomProtocol(duration=duration+burnin, rate=r)
            syn = AdaptivePreEventSynapseTwoTimeConstants(learning_rate, tau_trace=tau_pre, alpha=alpha, gamma=gammas[idx],
                                                          burst_def=burst_threshold)
            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration+burnin)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration+burnin)

            # Relaxation period
            syn.plasticity_on = False
            for t in range(int(burnin/dt)):
                syn.evolve(t, dt)
                if t%int(sampling_interval/dt) == 0:
                    weight_trace.append(syn.weight / dt)
                    er_est.append(syn.post_ma.moving_average['event'])
                    br_est.append(syn.post_ma.moving_average['burst'])
                    if syn.post_ma.moving_average['event'] > 0:
                        bp_est.append(syn.post_ma.moving_average['burst']/syn.post_ma.moving_average['event'])

            # Learning epoch
            syn.plasticity_on = True
            count = 0
            for t in range(int(burnin/dt), int((burnin+duration)/dt)):
                count += 1
                syn.evolve(t, dt)
                if t%int(sampling_interval/dt) == 0:
                    weight_trace.append(syn.weight / dt)
                    er_est.append(syn.post_ma.moving_average['event'])
                    br_est.append(syn.post_ma.moving_average['burst'])
                    if syn.post_ma.moving_average['event'] > 0:
                        bp_est.append(syn.post_ma.moving_average['burst']/syn.post_ma.moving_average['event'])
            weights[idx] += syn.weight/dt
            #  division by dt is because event and burst trains are not divided
            #  by dt in the rule (see preeventsynapse.py)

            P_greater = np.exp(-r * burst_threshold)
            E = r * P_greater
            B = E * (1 - P_greater)

            times = np.arange(len(weight_trace))*sampling_interval
            plt.plot(times, weight_trace)
            plt.show()

            plt.plot(times, er_est, label='ER')
            plt.plot(times, np.ones(len(er_est))*E, 'k')
            plt.plot(times, br_est, label='BR')
            plt.plot(times, np.ones(len(br_est))*B, color='grey')
            plt.legend()
            plt.show()

    weights /= nb_reals

    return rates, weights


def analytical(rates, eta, duration, alpha, burst_threshold, tau_pre):
    weights = []
    for r in rates:
        gamma = get_gamma(r, alpha, burst_threshold)
        P_greater = np.exp(-r*burst_threshold)
        E = r*P_greater
        B = E*(1-P_greater)
        term1 = (1./(1. + alpha*r))\
                *( 1 - np.exp(-burst_threshold*(r + 1./alpha)) )\
                *( E*r - E**2*np.exp(-burst_threshold/alpha) )

        term2 = E**2*np.exp(-burst_threshold/alpha)*(1.-np.exp(-burst_threshold/alpha))

        term3 = -E*B*np.exp(-burst_threshold/gamma)*(1.-np.exp(-burst_threshold/gamma))

        weights.append(duration*eta*tau_pre*E*(term1 + term2 + term3))

    return rates, weights


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


def f_(x, burst_threshold, cst):
    return np.exp(-burst_threshold/x) - np.exp(-2*burst_threshold/x) - cst


def f_alpha(alpha, rates, burst_threshold):
    P_greater = np.exp(-rates * burst_threshold)
    E = rates * P_greater
    B = E * (1 - P_greater)
    r = rates
    term1 = (1. / (1. + alpha * r)) \
            * (1 - np.exp(-burst_threshold * (r + 1. / alpha))) \
            * (r - E * np.exp(-burst_threshold / alpha))

    term2 = (E-B) * np.exp(-burst_threshold / alpha) * (1. - np.exp(-burst_threshold / alpha))

    return term1 + term2


def cst_gamma(rates, alpha, burst_threshold):
    P_greater = np.exp(-rates * burst_threshold)
    E = rates * P_greater
    B = E * (1 - P_greater)
    r =rates
    term1 = (1./ (1. + alpha * r)) \
            * (1 - np.exp(-burst_threshold * (r + 1. / alpha))) \
            * (r - E * np.exp(-burst_threshold / alpha))

    term2 = E * np.exp(-burst_threshold / alpha) * (1. - np.exp(-burst_threshold / alpha))
    return (term1 + term2)/B


def get_gamma(rates, alpha, burst_threshold):
    c = cst_gamma(rates, alpha, burst_threshold)
    sol = fsolve(f_, 1., args=(0.016, c))
    print('gamma = {}'.format(sol[0]))
    return sol[0]


if __name__ == '__main__':
    import sys
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns

    tau = np.arange(0., 0.050, 0.0001)
    rates = np.array([10]) #np.linspace(10, 10, 5)
    duration = 400
    #r, w, w_trace, bp_est = simulate(rates, 0.1, duration, 10., 0.016, 0.020)
    # plt.plot(1000*tau, cov_anal(tau, 10., 0.016))
    # print(burst_fraction(np.linspace(1., 30., 10), 0.016))
    for alpha in np.logspace(0, 2, 5):
        r, w = simulate(rates, 0.1, duration, alpha, 0.016, 0.020)
        r_anal, w_anal = analytical(rates, 0.1, duration, alpha, 0.016, 0.020)

        plt.figure(figsize=(3,3/1.6))
        plt.plot(r, w, color='black', label='simulation')
        plt.plot(r_anal, w_anal, color='grey', lw=1, label='theory')
        sns.despine()
        #plt.plot(r, np.zeros(r.shape), 'k--', lw=1)
        plt.xlabel('Rate [Hz]')
        plt.ylabel('$\Delta W$')
        plt.legend()
        plt.tight_layout()
        plt.savefig('../../results/learning-rule/DeltaW_AdaptiveLearningRuleTwoTimeCsts_Alpha'+str(alpha)+'.pdf')
        plt.close()