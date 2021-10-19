import kwant
from matplotlib import pyplot

# created lattice general
lat = kwant.lattice.general([(0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)],
                            [(0, 0, 0), (0.25, 0.25, 0.25)])
a, b = lat.sublattices

# make system at cubic form
def make_system(t=1., a=5, b=50, c=5):
    # created cube
    def make_cube(pos):
        x, y, z = pos
        return 0 <= x < a and 0 <= y < b and 0 <= z < c

    syst = kwant.Builder()
    syst[lat.shape(make_cube, (0, 0, 0))] = 0
    syst[lat.neighbors()] = t
    return syst

# call the main function if the script get executed
def main():
    # The standard plotting for 3D.
    syst = make_system()
    kwant.plot(syst)

    syst = make_system(a=1.1, b=1.1, c=1.1)

    # visualize the crystal structure better for a very small system.
    def family_colors(site):
        return 'brown' if site.family == a else 'c'

#    the function has a bug:
#    def family_size(site):
#        return 0.15 if site.family == a else 0.1

    kwant.plot(syst, site_color = family_colors, site_size = 0.1, site_lw = 0.01,
               hop_lw = 0.03)

if __name__ == '__main__':
    main()
