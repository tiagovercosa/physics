import kwant
import matplotlib.pyplot as plt

lat = kwant.lattice.honeycomb()
a, b = lat.sublattices

def make_system(t=1., r=10):

    # circular scattering region
    def circle(pos):
        x, y = pos
        return x**2 + y**2 < r**2

    syst = kwant.Builder()
    syst[lat.shape(circle, (0, 0))] = 0
    syst[lat.neighbors()] = t #To connect neighbors sites
    return syst


def main():
    syst = make_system()

    def family_colors(site):
        return 0 if site.family == a else 1

    kwant.plotter.plot(syst,
                       site_symbol = 'o',
                       site_size = 0.2,
                       site_color = family_colors,
                       site_lw = 0.05,
                       colorbar = False)

if __name__ == '__main__':
    main()
