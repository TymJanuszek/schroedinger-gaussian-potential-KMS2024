import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sc


class WaveFunc:

    def __init__(self, x_steps, dt, dk, k0, minbound, maxbound, v0 = 105, sigma = 0.5):

        self.xsteps = x_steps
        self.dx = (maxbound - minbound) / x_steps
        self.dt = dt
        self.dk = dk
        self.k0 = k0

        self.minbound = minbound
        self.maxbound = maxbound

        self.X = np.linspace(minbound, maxbound, self.xsteps + 1)
        self.psi = 1j * np.zeros(self.xsteps + 1)

        self.V0 = v0
        self.V = v0 * np.exp(-np.power(self.X - 10, 2) / sigma * sigma)

    def count_psi_0(self):
        self.psi = np.sqrt(self.dk) / (np.pi ** 0.25) * np.exp(-0.5 * self.X ** 2 * self.dk ** 2) \
                   * np.exp(1j * self.k0 * self.X)

        self.psi[0] = 0
        self.psi[- 1] = 0

    def count_A(self):
        A = 1j * np.zeros(shape=(self.xsteps - 1, self.xsteps - 1))
        for i in range(self.xsteps - 1):
            A[i, min(i + 1, self.xsteps - 2)] = -1j * self.dt / (2 * self.dx ** 2)
            A[i, max(i - 1, 0)] = -1j * self.dt / (2 * self.dx ** 2)
            A[i, i] = 1 + 1j * self.dt / 2 * (2 / self.dx ** 2 + self.V[i])
        return A

    def count_right_side(self):
        rside = 1j * np.zeros(self.xsteps + 1)
        for i in range(1, self.xsteps - 1):
            rside[i] = self.psi[i] + 1j * self.dt / 2 * (
                        (self.psi[i + 1] - 2 * self.psi[i] + self.psi[i - 1]) / np.power(self.dx, 2) - self.V[i] *
                        self.psi[i])
        return rside

    def simulate(self, duration, speed):
        self.count_psi_0()
        A = self.count_A()

        rng = int(duration / self.dt)
        end = len(self.psi) - 1
        for i in range(rng):
            rside = self.count_right_side()
            self.psi[1:end] = sc.spsolve(A, rside[1:end])

            # print(self.dx*np.sum(np.real(self.psi) ** 2 + np.imag(self.psi) ** 2))

            plt.plot(self.X, np.real(self.psi) ** 2 + np.imag(self.psi) ** 2)
            plt.plot(self.X, self.V / self.V0)
            plt.draw()
            plt.pause(1/speed)
            plt.clf()
