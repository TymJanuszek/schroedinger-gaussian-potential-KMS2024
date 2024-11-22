import scipy.sparse.linalg as sc
import io
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


class WaveFunc:

    def __init__(self, x_steps, dt, dk, k0, minbound, maxbound, v0=105, sigma=0.5):

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

    def count_right_side(self, n):
        rside = 1j * np.zeros(self.xsteps + 1)
        for i in range(1, self.xsteps - 1):
            rside[i] = self.psi[i] + 1j * self.dt / 2 * (
                    (self.psi[i + 1] - 2 * self.psi[i] + self.psi[i - 1]) / np.power(self.dx, 2) - n*self.V[i] *
                    self.psi[i])
        return rside

    def psi2_simulate(self, duration, speed, Vstatus):
        self.count_psi_0()
        A = self.count_A()

        rng = int(duration / self.dt)
        end = len(self.psi) - 1

        fig = plt.figure(figsize=(15, 9))

        for i in range(rng):
            rside = self.count_right_side(Vstatus)
            self.psi[1:end] = sc.spsolve(A, rside[1:end])

            plt.plot(self.X, Vstatus*self.V / self.V0, "g--", label=("Potential V(x)/"+str(self.V0)))
            plt.plot(self.X, np.real(self.psi) ** 2 + np.imag(self.psi) ** 2, "b", label=r"$|\psi(x,t)|^2$")
            plt.legend()
            plt.legend()
            plt.draw()
            plt.pause(1 / speed)
            plt.clf()

    def psi2_simulate_and_save(self, duration, filename, Vstatus):
        self.count_psi_0()
        A = self.count_A()

        rng = int(duration / self.dt)
        end = len(self.psi) - 1

        fig = plt.figure(figsize=(15, 9))
        frames = []
        plt.title(filename)

        for i in range(rng):
            rside = self.count_right_side(Vstatus)
            self.psi[1:end] = sc.spsolve(A, rside[1:end])

            plt.clf()
            plt.plot(self.X, Vstatus*self.V / self.V0, "g--", label=("Potential V(x)/" + str(self.V0)))
            plt.plot(self.X, np.real(self.psi) ** 2 + np.imag(self.psi) ** 2, "b", label=r"$|\psi(x,t)|^2$")
            plt.legend()
            buf = io.BytesIO()
            plt.savefig(buf, format="png")
            buf.seek(0)
            frames.append(Image.open(buf))

        frames[0].save(
            filename,
            save_all=True,
            append_images=frames[1:],
            duration=200,
            loop=0,
        )

    def psi_simulate(self, duration, speed, Vstatus):
        self.count_psi_0()
        A = self.count_A()

        rng = int(duration / self.dt)
        end = len(self.psi) - 1

        fig = plt.figure(figsize=(15, 9))

        for i in range(rng):
            rside = self.count_right_side(Vstatus)
            self.psi[1:end] = sc.spsolve(A, rside[1:end])

            plt.plot(self.X, Vstatus * self.V / self.V0, "g--", label=("Potential V(x)/" + str(self.V0)))
            plt.plot(self.X, np.imag(self.psi), "r", label=r"$Re(\psi(x,t))")
            plt.plot(self.X, np.real(self.psi), "b", label=r"$Im(\psi(x,t))$")
            plt.legend()
            plt.draw()
            plt.pause(1 / speed)
            plt.clf()

    def psi_simulate_and_save(self, duration, filename, Vstatus):
        self.count_psi_0()
        A = self.count_A()

        rng = int(duration / self.dt)
        end = len(self.psi) - 1

        fig = plt.figure(figsize=(15, 9))
        frames = []
        plt.title(filename)

        for i in range(rng):
            rside = self.count_right_side(Vstatus)
            self.psi[1:end] = sc.spsolve(A, rside[1:end])

            plt.clf()
            plt.plot(self.X, Vstatus*self.V / self.V0, "g--", label=("Potential V(x)/" + str(self.V0)))
            plt.plot(self.X, np.imag(self.psi), "r", label=r"$Re(\psi(x,t))")
            plt.plot(self.X, np.real(self.psi), "b", label=r"$Im(\psi(x,t))$")
            plt.legend()
            buf = io.BytesIO()
            plt.savefig(buf, format="png")
            buf.seek(0)
            frames.append(Image.open(buf))

        frames[0].save(
            filename,
            save_all=True,
            append_images=frames[1:],
            duration=200,
            loop=0,
        )