from pairingprotocols import RandomProtocol
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson processes.
Here, the plasticity rule is adaptive, i.e., the postsynaptic neuron
estimates its burst probability.
"""


def simulate(rates, eta, duration, alpha, burst_threshold, tau_pre):
    np.random.seed(2)

    # select parameters for pairing protocol
    tau_moving_average = alpha         # time constant for moving average
    burnin = 0*tau_moving_average    # seconds, relaxation of moving averages
    nb_reals = 10

    # create synapse
    learning_rate = eta
    syn = AdaptivePreEventSynapse(learning_rate, tau_trace=tau_pre, tau_ma=tau_moving_average, burst_def=burst_threshold,
                                  starting_estimate=(5., 0.35*5.))

    # integration time step
    dt = 0.001

    weights = np.zeros(rates.shape)
    weight_trace = []
    bp_est = []
    bp = []
    for n in range(nb_reals):
        if n % 10 == 0:
            print('realization #{}'.format(n+1))
        for idx, r in enumerate(rates):
            syn.reset()
            protocol = RandomProtocol(duration=duration+burnin, rate=r)
            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration+burnin)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration+burnin)
            #bp_loc = np.sum(syn.pre.train['burst'])/np.sum(syn.pre.train['event'])
            #bp.append(bp_loc)
            #plot_cov(syn, r, dt, burst_threshold)

            syn.plasticity_on = False
            for t in range(int(burnin/dt)):
                syn.evolve(t, dt)
                weight_trace.append(syn.weight/dt)
                if syn.post_ma.moving_average['event'] > 0:
                    bp_est.append(syn.post_ma.moving_average['burst']/syn.post_ma.moving_average['event'])

            syn.plasticity_on = True
            count = 0
            for x in range(int(burnin/dt), int((burnin+duration)/dt)):
                count += 1
                syn.evolve(x, dt)
                weight_trace.append(syn.weight/dt)
                bp_est.append(syn.post_ma.moving_average['burst']/syn.post_ma.moving_average['event'])
            #print('real bp = {:4.3f}, estimated bp = {4.3f}'.format(bp_loc, bp_est/count))
            #print(r, syn)
            weights[idx] += syn.weight/dt
            #  division by dt is because event and burst trains are not divided
            #  by dt in the rule (see preeventsynapse.py)
    weights /= nb_reals

    return rates, weights, weight_trace, bp_est


def analytical(rates, eta, duration, alpha, burst_threshold, tau_pre):
    # select parameters for pairing protocol
    r = rates
    P_greater = np.exp(-r*burst_threshold)
    E = r*P_greater
    B = E*(1-P_greater)
    term1 = (alpha/(1. + alpha*r))\
            *( 1 - np.exp(-burst_threshold*(r + 1./alpha)) )\
            *( E*r - E**2*np.exp(-burst_threshold/alpha) )

    term2 = alpha*E*(E-B)*np.exp(-burst_threshold/alpha)*(1.-np.exp(-burst_threshold/alpha))

    weights = duration*(eta/alpha)*tau_pre*E*(term1 + term2)

    return r, weights


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
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns

    tau = np.arange(0., 0.050, 0.0001)
    rates = np.linspace(1, 50, 10)
    duration = 60.
    #r, w, w_trace, bp_est = simulate(rates, 0.1, duration, 10., 0.016, 0.020)
    # plt.plot(1000*tau, cov_anal(tau, 10., 0.016))
    # print(burst_fraction(np.linspace(1., 30., 10), 0.016))
    for alpha in np.logspace(0, 2, 5):
        r, w, w_trace, bp_est = simulate(rates, 0.1, duration, alpha, 0.016, 0.020)
        #r_anal, w_anal = analytical(rates, 0.1, duration, alpha, 0.016, 0.020)

        plt.figure(figsize=(3,3/1.6))
        plt.plot(r, w, color='black', label='simulation')
        #plt.plot(r_anal, w_anal, color='grey', lw=1, label='theory')
        sns.despine()
        plt.plot(r, np.zeros(r.shape), 'k--', lw=1)
        plt.xlabel('Rate [Hz]')
        plt.ylabel('$\Delta W$')
        plt.legend()
        plt.tight_layout()
        plt.savefig('../../results/learning-rule/DeltaW_AdaptiveLearningRule_Alpha'+str(alpha)+'.pdf')
        plt.close()