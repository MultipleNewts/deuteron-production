# %%
import numpy as np
from scipy.constants import c


class FourVector:
    def __init__(self, p):
        """
            Creates a four-vector that enables Lorentz transformations and finding CoM frames

            Parameters
            ----------
            p : `array`
                the original momentum four-vector of form `[p0, p1, p2, p3] = [E, px*c, py*c, pz*c]`

            Returns
            -------
            ph : `type`
                lorem ipsum dolor sit amet
        """
        self.p = np.array(p, dtype=float)

    def __add__(self, other):
        """
            Defines four-vector addition, allows use of the `+` operator

            Parameters
            ----------
            self : `FourVector`
                A momentum four-vector
            other : `FourVector`
                A momentum four-vector

            Returns
            -------
            sum : `FourVector`
                the sum of the two supplied four-vectors
        """
        return FourVector(self.p + other.p)

    def set_pvectr(self, p):
        """
            Manually sets the four-vector to a new state
        """
        self.p = np.array(p, dtype=float)

    def get_lfactor(self, beta):
        """
            Computes gamma for a given total beta
        """
        return 1/(np.sqrt(1-beta**2))

    def shift_to_v(self, v):
        """
            Shifts the FourVector into the frame defined by the velocity `v`

            Parameters
            ----------
            v : `array`
                the velocity of the frame, of form `v = [vx, vy, vz]`
        """
        if np.count_nonzero(v) == 1:
            self.axis_lorentz(v)
        else:
            self.general_lorentz(v)

    def boost(self, v):
        """
            Boosts the FourVector by some velocity `v`

            Parameters
            ----------
            v : `array`
                the velocity to boost by, of form `v = [vx, vy, vz]`
        """
        self.shift_to_v(-np.array(v))

    def axis_lorentz(self, v):
        """
            Computes a single axis Lorentz transformation of `v`
        """
        idx = np.flatnonzero(v)
        axis = idx[0] + 1
        E = self.p[0]
        p_mag = self.p[axis]
        beta = v[axis-1]/c
        gamma = self.get_lfactor(beta)

        self.p[0] = gamma * (E - beta*p_mag)
        self.p[axis] = gamma * (p_mag - beta*E)

    def general_lorentz(self, v):
        """
            Computes the Lorentz Transformation for a velocity in an arbitrary direction
        """
        bx, by, bz = np.array(v)/c
        beta = np.sqrt(bx**2 + by**2 + bz**2)
        gamma = self.get_lfactor(beta)
        alpha = (gamma**2)/(1+gamma)
        transform = np.array([
            [gamma, -gamma*bx, -gamma*by, -gamma*bz],
            [-gamma*bx, 1+alpha*(bx**2), alpha*bx*by, alpha*bx*bz],
            [-gamma*by, alpha*by*bx, 1+alpha*(by**2), alpha*by*bz],
            [-gamma*bz, alpha*bz*bx, alpha*bz*by, 1+alpha*(bz**2)]
        ])

        self.p = transform @ self.p


# %%
if __name__ == "__main__":
    four_vectors = np.array([
        [938, 0, 0, 0],
        [938, 0, 0, 0],
        [938, 0, 0, 0],
        [1000, 0, 600, 0],
        [1000, 300, 0, 400],
        [1000, 1000, 0, 0]
    ], dtype=float)
    betas = np.array([
        [0, 0, 0],
        [0, 0.8, 0],
        [0.2, 0.8, -0.4],
        [0, 0.5, 0],
        [0, 0.6, 0],
        [0, 0.6, 0]
    ], dtype=float)
    vs = betas*c
    temp = FourVector([0, 0, 0, 0])
    for vect, vel in zip(four_vectors, vs):
        temp.set_pvectr(vect)
        temp.shift_to_v(vel)
        print(f"result: {temp.p}")

    p1 = FourVector([1, 2, -1, 0])
    p2 = FourVector([2, 2, 2, 2])
    CoM = p1 + p2
    print(CoM.p)

# %%
