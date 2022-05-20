import kwant
import matplotlib.pyplot as plt
import numpy as np

L = 5
W = 5

graphene = kwant.lattice.honeycomb()
syst = kwant.Builder()
a, b = graphene.sublattices

def family_color(site):
    if site.family == a:
        return 'red'
    else:
        return 'green'
#, site_color=family_color
def plot_system(syst):
    kwant.plot(syst, site_lw=0.05)

syst[(a(i, j) for i in range(L) for j in range(W))] = 4
syst[(b(i, j) for i in range(L) for j in range(W))] = 4
syst[graphene.neighbors()] = -1

# lead 1
lat_lead1 = kwant.lattice.square()
sym_lead1 = kwant.TranslationalSymmetry((0, 1))

syst[(lat_lead1(i, 5) for i in range(3, 6))] = 4
syst[((lat_lead1(i+2, 5), b(i, 4)) for i in range(1, 4))] = -1
syst[lat_lead1.neighbors()] = -1

lead1 = kwant.Builder(sym_lead1)
lead1[(lat_lead1(i, 5) for i in range(3, 6))] = 4
lead1[lat_lead1.neighbors()] = -1

# lead 2
lat_lead2 = kwant.lattice.square()
sym_lead2 = kwant.TranslationalSymmetry((0, -1))

syst[(lat_lead2(i, -1) for i in range(1, 4))] = 4
syst[((lat_lead2(i, -1), b(i, 0)) for i in range(1, 4))] = -1
syst[lat_lead2.neighbors()] = -1

lead2 = kwant.Builder(sym_lead2)
lead2[(lat_lead2(i, -1) for i in range(1, 4))] = 4
lead2[lat_lead1.neighbors()] = -1

syst.attach_lead(lead1)
syst.attach_lead(lead2)

sistema = syst.finalized()

def plot_conductance(syst, energies):
    data = []

    for energy in energies:
        smatrix = kwant.smatrix(sistema, energy)
        data.append(smatrix.transmission(1, 0))

    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()

def main():
    plot_system(sistema)
    plot_conductance(sistema, energies=np.linspace(1, 4.25, 200))

if __name__ == "__main__":
    main()

















# graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
#                                  [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites

# def square(pos):
#     return all(-20 < p < 20 for p in pos)

# bulk_graphene = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs))
# bulk_graphene[graphene.shape((lambda pos: True), (0, 0))] = 0
# bulk_graphene[graphene.neighbors(1)] = 1
 
# zigzag_ribbon = kwant.Builder(kwant.TranslationalSymmetry([1, 0]))
# zigzag_ribbon[graphene.shape((lambda pos: abs(pos[1]) < 9), (0, 0))] = 0
# zigzag_ribbon[graphene.neighbors(1)] = 1

# armchair_ribbon = kwant.Builder(kwant.TranslationalSymmetry([0, np.sqrt(3)]))
# armchair_ribbon[graphene.shape((square), (0, 0))] = 0
# armchair_ribbon[graphene.neighbors(1)] = 1
 
# armchair_ribbon.attach_lead(lead)
# armchair_ribbon.attach_lead(lead.reversed())
# 
# kwant.plot(armchair_ribbon)
# syst = armchair_ribbon.finalized()
# 
# kwant.plotter.bands(syst, fig_size=(12, 8));
