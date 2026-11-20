# %%
import numpy as np
from combinations import DR_valid_combs
from data_parser import Hadron
import matplotlib.pyplot as plt

# index goes a_{-1} to a_{10} then b_{1} and b_{2}
pn_dg_bestfit_vals = [
    2.30346,
    -9.366346e1,
    2.565390e3,
    -2.5594101e4,
    1.43513109e5,
    -5.0357289e5,
    1.14924802e6,
    -1.72368391e6,
    1.67934876e6,
    -1.01988855e6,
    3.4984035e5,
    -5.1662760e4,
    -5.1885,
    2.9196,
]

pn_pipm_bestfit_vals = [
    6.465e6,
    1.051e1,
    1.979e3,
    5.363e0,
    6.045e5,
    2.549e15,
    1.657e1,
    2.330e7,
    1.119e1,
    2.868e16,
]


def pn_dg_cross(k):
    """
        Computes the cross section for a value k for the:
            pn -> d+gamma

        Parameters
        ----------
        k : `float`
            the difference in momentum between the two interacting particles

        Returns
        -------
        sigma : `float`
            the cross-section, sigma, of the interaction in barns
    """
    val = 0
    if k < 1.28:
        for i in range(12):
            val += pn_dg_bestfit_vals[i]*np.power(k, (i-1))
    else:
        val += np.exp(-pn_dg_bestfit_vals[12]*k-pn_dg_bestfit_vals[13]*np.power(k, 2))
    return val


def hh_pi0_cross(k, types):
    """
        Computes the cross section for a value k for the:
            pp -> dpi^{+}pi^{0}
            pn -> dpi^{0}pi^{0}
            nn -> dpi^{-}pi^{0}

        Parameters
        ----------
        k : `float`
            the difference in momentum between the two interacting particles

        Returns
        -------
        sigma : `float`
            the cross-section, sigma, of the interaction in barns
    """
    if (types == DR_valid_combs[2]) or (types == DR_valid_combs[3]):
        params = [5.099e15, 1.656e1, 2.333e7, 1.133e1, 2.868e16]
    else:
        params = [2.855e6, 1.311e1, 2.961e3, 5.572, 1.461e6]
    a, b, c, d, e = params
    return (a*np.power(k, b))/(np.power((c - np.exp(d*k)), 2) + e)


def hh_dpi_cross(k):
    """
        Computes the cross section for a value k for the:
            pp -> dpi^{+}
            pn -> dpi^{0}
            nn -> dpi^{-}

        Parameters
        ----------
        k : `float`
            the difference in momentum between the two interacting particles

        Returns
        -------
        sigma : `float`
            the cross-section, sigma, of the interaction in barns
    """
    eta = k*2
    if k < 0.7:
        return 0
    params = [170, 1.34, 1.77, 0.38, 0.096]
    a, b, c, d, e = params
    return (a*np.power(eta, b))/(np.power((c - np.exp(d*eta)), 2) + e)


def pn_pipm_cross(k):
    """
        Computes the cross section for a value k for the:
            pn -> dpi^{+}pi^{-}

        Parameters
        ----------
        k : `float`
            the difference in momentum between the two interacting particles

        Returns
        -------
        sigma : `float`
            the cross-section, sigma, of the interaction in barns
    """
    a, b, c, d, e = pn_pipm_bestfit_vals[0:5]
    term1 = (a*np.power(k, b))/(np.power((c - np.exp(d*k)), 2) + e)
    a, b, c, d, e = pn_pipm_bestfit_vals[5:10]
    term2 = (a*np.power(k, b))/(np.power((c - np.exp(d*k)), 2) + e)
    return term1 + term2


# %%
if __name__ == "__main__":
    pp = [Hadron.PROTON, Hadron.PROTON]
    pn = [Hadron.PROTON, Hadron.NEUTRON]
    X = np.linspace(0, 5, 1000)
    Y1 = np.array([pn_dg_cross(val) for val in X])
    Y2 = [hh_pi0_cross(val, pp) for val in X]
    Y3 = [hh_pi0_cross(val, pn) for val in X]
    Y4 = np.array([hh_dpi_cross(val) for val in X])/2
    Y5 = [pn_pipm_cross(val) for val in X]
    plt.plot(X, Y1, label=r"$pn\to d\gamma$")
    plt.plot(X, Y2, label=r"$N N\to d\pi^{0}\pi$")
    plt.plot(X, Y3, label=r"$np\to d\pi^{0}\pi$")
    plt.plot(X, Y4, label=r"$N_{1}N_{2}\to d\pi$")
    plt.plot(X, Y5, label=r"$pn\to d\pi^{+}\pi^{-}$")
    plt.legend()
    # plt.yscale("log")
    plt.loglog()
    plt.ylim([0.01, 1e4])
    plt.xlim([0.1, 10])
    plt.show()

# %%
