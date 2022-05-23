import kwant
import matplotlib.pyplot as plt
import numpy as np

lat = kwant.lattice.honeycomb()
a, b = lat.sublattices

def make_system(L=5, W=5):
    syst = kwant.Builder()
    syst[(a(i, j) for i in range(L) for j in range(W))] = 4
    syst[(b(i, j) for i in range(L) for j in range(W))] = 4
    syst[lat.neighbors()] = -1

    # Create a lead
    lat_lead = kwant.lattice.square()
    sym_lead1 = kwant.TranslationalSymmetry((0, 1))
    
    lead1 = kwant.Builder(sym_lead1)
    lead1[(lat_lead(i, 0) for i in range(2, 7))] = 4
    lead1[lat_lead.neighbors()] = -1
    
    syst[(lat_lead(i, 5) for i in range(2, 7))] = 4
    syst[lat_lead.neighbors()] = -1
    
    # Manually attach sites from graphene to square lattice
    syst[((lat_lead(i+2, 5), b(i, 4)) for i in range(5))] = -1
    
    # lead 2
    lat_lead2 = kwant.lattice.square()
    sym_lead2 = kwant.TranslationalSymmetry((0, -1))
    
    syst[(lat_lead2(i, -1) for i in range(0, 5))] = 4
    syst[((lat_lead2(i, -1), b(i, 0)) for i in range(0, 5))] = -1
    syst[lat_lead2.neighbors()] = -1
    
    lead2 = kwant.Builder(sym_lead2)
    lead2[(lat_lead2(i, -1) for i in range(0, 5))] = 4
    lead2[lat_lead2.neighbors()] = -1
    
    #syst.attach_lead(lead1)
    #syst.attach_lead(lead2)

    return syst, [lead1, lead2]


def family_colors(site):
    return 0 if site.family == a else 1

def plot_system(syst):
    kwant.plot(syst, site_lw=0.05, site_color=family_colors,
            lead_site_lw=0, colorbar=False)

def plot_conductance(syst, energies):
    data = []

    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))

    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energia [t]")
    plt.ylabel("condutância [e^2/h]")
    plt.show()

def plot_bandstructure(flead, momenta):
    bands = kwant.physics.Bands(flead)
    energies = [bands(k) for k in momenta]

    plt.figure()
    plt.plot(momenta, energies)
    plt.xlabel("momentum [(a)^-1]")
    plt.ylabel("energia [t]")
    plt.show()

def main():
    syst, leads = make_system()

    for lead in leads:
        syst.attach_lead(lead)

    plot_system(syst)

    syst = syst.finalized()
    
    plot_conductance(syst, energies=np.linspace(0.75, 2.75, 500))

    pi = np.pi
    #plot_bandstructure(syst.leads[1], momenta=np.linspace(-pi, pi, 1000))

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
