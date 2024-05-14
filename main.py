# This is a sample Python script.
import itertools
import numpy as np
from tetrahedron import Tetrahedron
from valid_with_cgal import write_cgal_output

# visualization
import matplotlib.pyplot as plt
import vedo

class DelaunayTriangulation:
    def __init__(self, file):
        self.tetrahedrons = set()
        self.visit_points = set()
        self.super_tetrahedron = None
        if isinstance(file, str):
            if file.endswith('pdb'):
                self.coord = self.get_coord(file)
            if file.endswith('xyz'):
                self.coord = np.loadtxt(file)
        elif isinstance(file, int):
            self.coord = np.random.rand(file, 3)
        self.coord_wo_numpy = self.coord.tolist()
        self.natom = self.coord.shape[0]

    def get_coord(self, file):
        coord = []
        f = open(file)
        for line in f:
            if line.startswith('ATOM'):
                atom = line[13:16].strip()
                if atom == 'CA':
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord.append([x,y,z])
        f.close()
        return np.array(coord)

    def get_super_tetrahedron(self):
        natom = self.natom
        min_coord = self.coord.min(axis=0)
        max_coord = self.coord.max(axis=0)
        delta = max(max_coord - min_coord)
        p1 = min_coord - delta
        p2 = min_coord + 2 * delta * np.array([1, 0, 0])
        p3 = min_coord + 2 * delta * np.array([0, 1, 0])
        p4 = min_coord + 2 * delta * np.array([0, 0, 1])
        super_tetrahedron = Tetrahedron(p1, p2, p3, p4, index=[natom, natom+1, natom+2, natom+3])
        self.super_tetrahedron = super_tetrahedron
        self.tetrahedrons.add(super_tetrahedron)

    def is_tetrahedron(self, p1, p2, p3, p4):
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        p4 = np.array(p4)
        # Calculate the volume of the tetrahedron
        volume = np.abs(np.dot(p4 - p1, np.cross(p2 - p1, p3 - p1))) / 6.0
        return volume > 0.0

    def add_point(self, point, idx, verbose=True):
        # Select a point one by one and check the conditions of Delaunay triangulation.
        if verbose: print("SELECTED VERTEX:", idx, point)
        bad_tetra = set()
        for tetra in self.tetrahedrons:
            if tetra.check_point_in_sphere(point):
                if verbose: print("DEL TETRA:", tetra.index)
                bad_tetra.add(tetra)

        # Form a new triangulation by connecting a line from
        # the vertex of the erased triangle to
        # the point that did not meet the conditions.
        vertices = set()
        for each in bad_tetra:
            #vertices.update(each.key)
            vertices.update(each.index)
            self.tetrahedrons.remove(each)

        if verbose: print("VERTICES FROM DEL-TETRA:", vertices)
        for i,j,k in itertools.combinations(vertices, 3):
            # point 1
            if i >= self.natom:
                p1 = self.super_tetrahedron.p1234[i - self.natom]
            else:
                p1 = self.coord[i]
            # point 2
            if j >= self.natom:
                p2 = self.super_tetrahedron.p1234[j - self.natom]
            else:
                p2 = self.coord[j]
            # point 3
            if k >= self.natom:
                p3 = self.super_tetrahedron.p1234[k - self.natom]
            else:
                p3 = self.coord[k]

            # check via volume of tetrahedron
            if not self.is_tetrahedron(p1, p2, p3, point): continue

            # build Tetrahedron object
            t = Tetrahedron(p1,p2,p3,point,index=[i,j,k,idx])
            if verbose: print("BUILD TETRA:", [i,j,k,idx])

            pass_all_visitor = True
            for visitor in self.visit_points:
                if visitor in t.p1234: continue
                has_visitor_inside = t.check_point_in_sphere(visitor)
                if has_visitor_inside:
                    pass_all_visitor = False
                    break
            if pass_all_visitor:
                self.tetrahedrons.add(t)

    def triangulation(self):
        for i,coor in enumerate(self.coord):
            self.add_point(coor,i)
            self.visit_points.add(tuple(coor))
            #self.draw_by_vedo(pts=None, faces=None) # if None, use intrinsic data
            #self.draw()
        pass

    def delete_super_tetrahedron(self):
        tetra_contains_super = set()
        superset = set(self.super_tetrahedron.key)
        for tetra in self.tetrahedrons:
            intersec = superset.intersection(set(tetra.key))
            if intersec:
                tetra_contains_super.add(tetra)
        for bad in tetra_contains_super:
            self.tetrahedrons.remove(bad)

    def check_if_all_coords_in_super_tetra(self):
        self.tetrahedrons.clear()
        self.get_super_tetrahedron()
        t = self.tetrahedrons.pop()
        for p1,p2,p3,p4 in itertools.combinations(self.coord, 4):
            if not t.check_point_in_sphere(p1): raise ValueError("has problem in equation! (p1)")
            if not t.check_point_in_sphere(p2): raise ValueError("has problem in equation! (p2)")
            if not t.check_point_in_sphere(p3): raise ValueError("has problem in equation! (p3)")
            if not t.check_point_in_sphere(p4): raise ValueError("has problem in equation! (p4)")

    def check_if_validate_conditions_of_Delaunay(self):
        bad_tetra = set()
        pass_test = True
        for t in self.tetrahedrons:
            for coor in self.coord:
                if tuple(coor) in t.p1234:
                    print("it is angular point!")
                    continue
                if t.check_point_in_sphere(coor):raise ValueError("Do not fulfill the condition!")
                pass_test = False
                #if t.check_point_in_sphere(coor): bad_tetra.add(t)
        #for bad in bad_tetra:
        #    self.tetrahedrons.remove(bad)
        return pass_test

    def draw(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        vertices = None
        for i,tetra in enumerate(self.tetrahedrons):
            point4 = np.array([tetra.p1, tetra.p2, tetra.p3, tetra.p4])
            if i > 0:
                vertices = np.vstack((point4, vertices))
            else:
                vertices = point4

        if len(self.tetrahedrons) > 0:
            ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], linewidth=0, color='skyblue', alpha=0.5)
            ax.scatter(self.coord[:, 0], self.coord[:, 1], self.coord[:, 2], color='red', s=50)
            plt.show()
        else:
            print("can't plot surf because NO tetrahedron!")

    def draw_debug(self, tetra, point):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        vertices = np.array([tetra.p1, tetra.p2, tetra.p3, tetra.p4])
        ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], linewidth=0, color='skyblue', alpha=0.5)
        ax.scatter(point[0], point[1], point[2], color='red', s=50)
        plt.show()

    def draw_by_vedo(self, pts=None, faces=None):
        if not pts:
            pts = np.vstack((self.coord, self.super_tetrahedron.key))
        if not faces:
            faces = [tuple(x.index) for x in self.tetrahedrons]

        print('pts:')
        print(pts)
        print('faces:')
        print(faces)
        # Build the polygonal Mesh object from the vertices and faces
        mesh = vedo.Mesh([pts, faces])
        # Set the backcolor of the mesh to violet
        # and show edges with a linewidth of 2
        mesh.backcolor('violet').linecolor('tomato').alpha(0.3).linewidth(2)
        labels = [vedo.Text3D(i, pos=pts[i], s=.2, c='k') for i in range(len(pts))]
        vedo.show(mesh, labels).close()

    def write_to_file(self, outfile):
        wf = open(outfile, 'w')
        natom = self.natom
        wf.write("3\n") # three-dimension
        wf.write("%d\n" % natom)
        for x,y,z in sorted(self.coord_wo_numpy):
            wf.write('%.5f %.5f %.5f\n' % (x,y,z))
        wf.write("%d\n" % len(self.tetrahedrons)) # should be number_of_finite_cells, but number of tetrahedrons here
        for tetra in self.tetrahedrons:
            i,j,k,l = tetra.index
            wf.write("%d %d %d %d\n" % (i,j,k,l))
        wf.close()

    def read_file(self, outfile):
        f = open(outfile)
        lines = f.readlines()
        f.close()
        XYZ = []
        FACE = []
        ndim = int(lines.pop(0))
        natom = int(lines.pop(0))

        for i in range(natom):
            line = lines.pop(0)
            XYZ.append([float(number) for number in line.strip().split()])
        nface = int(lines.pop(0))

        for x in range(nface):
            line = lines.pop(0)
            if not line.strip(): continue
            i,j,k,l = line.split()
            i = float(i)
            j = float(j)
            k = float(k)
            l = float(l)
            FACE.append((i,j,k,l))
        return XYZ, FACE


if __name__ == '__main__':
    test = True
    show_valid = True
    show_final = True
    valid_file = "./out_cgal.txt"
    inhouse_file = "./out_inhouse.txt"
    #dt = DelaunayTriangulation('3nir.pdb')
    dt = DelaunayTriangulation('point6.xyz')
    #dt = DelaunayTriangulation(12)
    dt.get_super_tetrahedron()
    dt.triangulation()
    dt.delete_super_tetrahedron()
    #dt.draw()

    # write output file
    dt.write_to_file(inhouse_file)
    write_cgal_output(dt.coord, valid_file)

    if show_final:
        print("############################")
        print("## DRAW FINAL 3D DT MESH ###")
        print("############################")
        xyz, face = dt.read_file(inhouse_file)
        dt.draw_by_vedo(pts=xyz, faces=face)

    if show_valid:
        print("############################")
        print("## DRAW VALIDATION FILE  ###")
        print("############################")
        xyz, face = dt.read_file(valid_file)
        dt.draw_by_vedo(pts=xyz, faces=face)



    # test function
    if test:
        dt.check_if_all_coords_in_super_tetra()
        print(dt.check_if_validate_conditions_of_Delaunay())