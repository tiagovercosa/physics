import kwant
from matplotlib import pyplot

lat = kwant.lattice.general([(0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)],
                            [(0, 0, 0), (0, 0, 0.5)])

def make_system(t=1., a=10, b=10, c=10):
    def make_cube(pos):
        x, y, z = pos
        return 0 <= x < a and 0 <= y <= b and 0 <= z < c

    syst = kwant.Builder()
    syst[lat.shape(make_cube, (0, 0, 0))] = None
    syst[lat.neighbors()] = t

    return syst

def main():
    syst = make_system()
    kwant.plot(syst)
    
    syst = make_system(a=1.5, b=1., c=1.5)
    kwant.plot(syst, site_size = 0.1, site_lw = 0.01, hop_lw = 0.02)


if __name__ == '__main__':
    main()
