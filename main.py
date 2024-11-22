import schrodinger
import numpy as np

func = schrodinger.WaveFunc(500, 0.002, 1.5, 10, -5, 25)

func.simulate(10, 10000)

