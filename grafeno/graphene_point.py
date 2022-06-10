import kwant
from matplotlib import pyplot
import scipy.sparse.linalg as sla
from math import pi, sqrt, tanh
import numpy as np

graphene = kwant.lattice.honeycomb()
a, b = graphene.sublattices

def make_system(r=10):
    def circle(pos):
        x, y = pos
        return x ** 2 + y ** 2 < r ** 2

    syst = kwant.Builder()

    syst[graphene.shape(circle, (0, 0))] = 4

    hoppings = (((0, 0), a, b), ((0, 1), a, b), ((-1, 1), a, b))
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    # left lead
    sym0 = kwant.TranslationalSymmetry(graphene.vec((-1, 0)))

    def lead0_shape(pos):
        x, y = pos
        return (-0.4 * r < y < 0.4 * r)

    lead0 = kwant.Builder(sym0)
    lead0[graphene.shape(lead0_shape, (0, 0))] = -1
    lead0[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    # The second lead
    sym1 = kwant.TranslationalSymmetry(graphene.vec((1, 0)))

    def lead1_shape(pos):
        x, y = pos
        return (-0.4 * r < y < 0.4 * r)

    lead1 = kwant.Builder(sym1)
    lead1[graphene.shape(lead1_shape, (0, 0))] = -1
    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    return syst, [lead0, lead1]

def plot_conductance(syst, energies):
    # Compute transmission as a function of energy
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(0, 1))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energia [t]")
    pyplot.ylabel("condutância [e^2/h]")
    pyplot.show()


def plot_bandstructure(flead, momenta):
    bands = kwant.physics.Bands(flead)
    energies = [bands(k) for k in momenta]

    pyplot.figure()
    pyplot.plot(momenta, energies)
    pyplot.xlabel("momentum [(a)^-1]")
    pyplot.ylabel("energia [t]")
    pyplot.show()

def main():
    syst, leads = make_system()

    # To highlight the two sublattices of graphene, we plot one with
    # a filled, and the other one with an open circle:
    def family_colors(site):
        return 0 if site.family == a else 1

    # Attach the leads to the system.
    for lead in leads:
        syst.attach_lead(lead)

    # Then, plot the system with leads.
    #kwant.plot(syst, site_color=family_colors, site_lw=0.1,
    #           lead_site_lw=0, colorbar=False)

    # Finalize the system.
    syst = syst.finalized()

    # Compute the band structure of lead 0.
    #pi = np.pi
    #momenta = np.linspace(-pi, pi, 1000)
    #momenta = [0 + 0.02 * pi * i for i in range(101)]
    #plot_bandstructure(syst.leads[1], momenta)

    # Plot conductance.
    energies = np.linspace(0.75, 2.75, 500)
    plot_conductance(syst, energies)

if __name__ == '__main__':
    main()
