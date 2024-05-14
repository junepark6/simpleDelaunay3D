import numpy as np
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Triangulation_3 import Delaunay_triangulation_3

def write_cgal_output(coord, outfile, manual=True, verbose=False):
    L = []
    for x,y,z in coord:
        L.append(Point_3(x,y,z))
    T = Delaunay_triangulation_3(L)
    if verbose:
        print('number_of_vertices:',T.number_of_vertices())
        print('number_of_facets:',T.number_of_facets())
        print('number_of_cells:',T.number_of_cells())
        print('number_of_finite_edges:',T.number_of_finite_edges())
        print('number_of_finite_facets:',T.number_of_finite_facets())
        print('number_of_finite_cells:',T.number_of_finite_cells())

    if manual:
        f = open(outfile, 'w')
        natom = coord.shape[0]
        f.write('3\n')
        f.write('%d\n' % natom)
        vertex2index = {}
        XYZ = []
        for i,vertex in enumerate(T.all_vertices()):
            if i == 0: continue
            xyz = vertex.point().__str__()
            XYZ.append(xyz)
        # sort order of atoms

        XYZ.sort(key=lambda x: tuple([float(xx) for xx in x.split()]))
        # write each atom to outfile
        for i,xyz in enumerate(XYZ):
            vertex2index[xyz] = i
            f.write('%s\n' % xyz)

        f.write('%d\n' % T.number_of_finite_cells())
        for cell in T.finite_cells():
            p1 = cell.vertex(0).point().__str__()
            p2 = cell.vertex(1).point().__str__()
            p3 = cell.vertex(2).point().__str__()
            p4 = cell.vertex(3).point().__str__()
            i1 = vertex2index[p1]
            i2 = vertex2index[p2]
            i3 = vertex2index[p3]
            i4 = vertex2index[p4]
            f.write('%d %d %d %d\n' % (i1,i2,i3,i4))
        f.close()
    else:
        T.write_to_file(outfile)

def main():
    natom = 6
    coord = np.random.rand(natom, 3)
    write_cgal_output(coord, "./out_cgal.txt", verbose=True)

if __name__ == "__main__":
   main()