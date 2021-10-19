import kwant
from matplotlib import pyplot

lat = kwant.lattice.square()
syst = kwant.Builder()

syst[lat(1, 0)] = 4
syst[lat(0, 1)] = 4

syst[(lat(1, 0), lat(0, 1))] = 1j

kwant.plot(syst)
