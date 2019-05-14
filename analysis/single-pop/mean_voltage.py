import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../thesis_mplrc.dms')
from scipy.integrate import quad

def meanV(Re, Ri, we, wi, taue, taui):
   """
   :param Re: tau_e * N_e * r_e
   :param we: conductance in unit of leak conductance
   :param alpha_R: R_i/R_e
   :param alpha_w: w_i/w_e
   :return:
   """
   ge = meang(Re, we, taue)
   gi = meang(Ri, wi, taui)

   return (-70. + gi*(-80.))/(1 + ge + gi)

def meanVd(Re, Ri, we, wi, taue, taui, a):
   ge = meang(Re, we, taue)
   gi = meang(Ri, wi, taui)
   return ((1-a)*(-70.) + gi*(-80.))/(1 - a + ge + gi)

def varg(R, w, tau):
   return 0.5*R*tau*w**2

def meang(R, w, tau):
   return R*w*tau

def gT(Re, Ri, we, wi, taue, taui):
   return 1 + meang(Re, we, taue) + meang(Ri, wi, taui)

def varV(Re, Ri, we, wi, taue, taui, taum):
   V_mean = meanV(Re, Ri, we, wi, taue, taui)
   varge = varg(Re, we, taue)
   vargi = varg(Ri, wi, taui)
   g_T = gT(Re, Ri, we, wi, taue, taui)
   taumhat = taum/g_T
   Ve = V_mean**2*varge*taue/(taumhat + taue)/g_T**2
   Vi = (-80. - V_mean)**2*vargi*taui/(taumhat + taui)/g_T**2
   return Ve + Vi


def integrandvar(omega, taux, taum, g_T, tauz, a):
   numerator = (1 + omega**2*tauz**2)**2/np.pi
   denominator = (1 + omega**2*taux**2)*(g_T - a + g_T*tauz**2*omega**2)**2 + \
                   (1 + omega**2*taux**2)*(taum + a*tauz + taum*tauz**2*omega**2)**2*omega**2
   return numerator/denominator


def _varVd(Re, Ri, we, wi, taue, taui, taum, tauz, a):
   g_T = gT(Re, Ri, we, wi, taue, taui)
   V_mean = meanVd(Re, Ri, we, wi, taue, taui, a)
   varge = varg(Re, we, taue)
   vargi = varg(Ri, wi, taui)

   Ve = V_mean**2*varge*taue*2*quad(integrandvar, 0, np.inf, args=(taue, taum, g_T, tauz, a))[0]
   Vi = (-80. - V_mean)**2*vargi*taui*2*quad(integrandvar, 0, np.inf, args=(taui, taum, g_T, tauz, a))[0]
   return Ve + Vi

varVd = np.vectorize(_varVd)


if __name__ == '__main__':
   taue = 5e-3
   taui = 2*taue
   we = np.arange(0, 1, 0.01)
   wi = 2*we
   taum = 16e-3
   taud = 7e-3
   a = -13./170*7
   tauz = 30e-3

   #   Figure 1 - Soma
   plt.figure(1, figsize=(6.5, 7/1.6))
   plt.subplot(221)
   plt.plot(we, meanV(100, 100, we, wi, taue, taui), label=r'$R_e=100$')
   plt.plot(we, meanV(500, 500, we, wi, taue, taui), label=r'$R_e=500$')
   plt.plot(we, meanV(1000, 1000, we, wi, taue, taui), label=r'$R_e=1000$')
   plt.plot(we, meanV(5e3, 5e3, we, wi, taue, taui), label=r'$R_e=5000$')
   plt.text(0.4, 0.3, '$w_i/w_e =$ {:1.0f}'.format(wi[1]/we[1]), horizontalalignment='center', verticalalignment='center',
            transform=plt.gca().transAxes)
   plt.title(r'\textbf{Ai} Effect of $R_e$ on mean voltage', loc='left')
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\bar{V}$ [mV]')
   plt.legend(frameon=False, loc='best')

   plt.subplot(222)
   plt.plot(we, np.sqrt(varV(100, 100, we, wi, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(500, 500, we, wi, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(1000, 1000, we, wi, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(5000, 5000, we, wi, taue, taui, taum)))
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\sigma_{Vm}$ [mV]')
   plt.title(r'\textbf{Aii} Effect of $R_e$ on std($V_m$)', loc='left')

   plt.subplot(223)
   Re = 500
   Ri = Re
   plt.plot(we, meanV(Re, Ri, we, we, taue, taui), label='$w_i=w_e$')
   plt.plot(we, meanV(Re, Ri, we, 2*we, taue, taui), label='$w_i=2w_e')
   plt.plot(we, meanV(Re, Ri, we, 4*we, taue, taui), label='$w_i=4w_e$')
   plt.plot(we, meanV(Re, Ri, we, 8*we, taue, taui), label='$w_i=8w_e$')
   plt.text(0.7, 0.6, '$R_e =$ {}'.format(500), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
   plt.title(r'\textbf{Bi} Effect of $w_i/w_e$ on mean voltage', loc='left')
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\bar{V}$ [mV]')
   plt.legend(frameon=False, loc='best')

   plt.subplot(224)
   plt.plot(we, np.sqrt(varV(Re, Ri, we, we, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(Re, Ri, we, 2*we, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(Re, Ri, we, 4*we, taue, taui, taum)))
   plt.plot(we, np.sqrt(varV(Re, Ri, we, 8*we, taue, taui, taum)))
   plt.xlabel('$w_e$')
   plt.ylabel('$\sigma_{Vm}$ [mV]')
   plt.title(r'\textbf{Bii} Effect of $w_i/w_e$ on std($V_m$)', loc='left')

   plt.tight_layout()
   #plt.show()
   plt.savefig('../../results/single-pop/MeanANDVariance.pdf')
   plt.close()

   #   Figure 2 - Dendrite
   plt.figure(2, figsize=(6.5, 7/1.6))
   plt.subplot(221)
   plt.plot(we, meanVd(100, 100, we, wi, taue, taui, a), label=r'$R_e=100$')
   plt.plot(we, meanVd(500, 500, we, wi, taue, taui, a), label=r'$R_e=500$')
   plt.plot(we, meanVd(1000, 1000, we, wi, taue, taui, a), label=r'$R_e=1000$')
   plt.plot(we, meanVd(5e3, 5e3, we, wi, taue, taui, a), label=r'$R_e=5000$')
   plt.text(0.4, 0.3, '$w_i/w_e =$ {:1.0f}'.format(wi[1]/we[1]), horizontalalignment='center', verticalalignment='center',
            transform=plt.gca().transAxes)
   plt.title(r'\textbf{Ai} Effect of $R_e$ on mean voltage', loc='left')
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\bar{V}_d$ [mV]')
   plt.legend(frameon=False, loc='best')

   plt.subplot(222)
   plt.plot(we, np.sqrt(varVd(100, 100, we, wi, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(500, 500, we, wi, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(1000, 1000, we, wi, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(5000, 5000, we, wi, taue, taui, taud, tauz, a)))
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\sigma_{Vd}$ [mV]')
   plt.title(r'\textbf{Aii} Effect of $R_e$ on std($V_d$)', loc='left')

   plt.subplot(223)
   Re = 500
   Ri = Re
   plt.plot(we, meanVd(Re, Ri, we, we, taue, taui, a), label='$w_i=w_e$')
   plt.plot(we, meanVd(Re, Ri, we, 2*we, taue, taui, a), label='$w_i=2w_e')
   plt.plot(we, meanVd(Re, Ri, we, 4*we, taue, taui, a), label='$w_i=4w_e$')
   plt.plot(we, meanVd(Re, Ri, we, 8*we, taue, taui, a), label='$w_i=8w_e$')
   plt.text(0.7, 0.6, '$R_e =$ {}'.format(500), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
   plt.title(r'\textbf{Bi} Effect of $w_i/w_e$ on mean voltage', loc='left')
   plt.xlabel('$w_e$')
   plt.ylabel(r'$\bar{V}_d$ [mV]')
   plt.legend(frameon=False, loc='best')

   plt.subplot(224)
   plt.plot(we, np.sqrt(varVd(Re, Ri, we, we, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(Re, Ri, we, 2*we, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(Re, Ri, we, 4*we, taue, taui, taud, tauz, a)))
   plt.plot(we, np.sqrt(varVd(Re, Ri, we, 8*we, taue, taui, taud, tauz, a)))
   plt.xlabel('$w_e$')
   plt.ylabel('$\sigma_{Vd}$ [mV]')
   plt.title(r'\textbf{Bii} Effect of $w_i/w_e$ on std($V_d$)', loc='left')

   plt.tight_layout()
   #plt.show()
   plt.savefig('../../results/single-pop/MeanANDVariance_Dend.pdf')
   plt.close()