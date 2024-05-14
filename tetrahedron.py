import numpy as np
from scipy.spatial import distance
import itertools

class Tetrahedron:
    def __init__(self, p1, p2, p3, p4, index=None):
        if len(p1) == 3 and \
                len(p2) == 3 and \
                len(p3) == 3 and \
                len(p4) == 3:
            self.p1 = tuple(p1)
            self.p2 = tuple(p2)
            self.p3 = tuple(p3)
            self.p4 = tuple(p4)
            self.p1234 = [self.p1, self.p2, self.p3, self.p4]
            self.key = sorted([self.p1, self.p2, self.p3, self.p4])
            self.index = index
            self.radius = None
            self.circumcenter = None
            self.set_circumsphere()
            #self.edges = set()
            #for pair in itertools.combinations([self.p1, self.p2, self.p3, self.p4], 2):
            #    self.edges.add(tuple(sorted(pair)))
        else:
            raise ValueError("check 4 given points; size should be 3")

    def __repr__(self):
        return str(self.key)
        #if self.index:
        #    return self.index
        #else:
        #    return self.key

    def __str__(self):
        return("%s|%s|%s|%s" % (self.p1, self.p2, self.p3, self.p4))

    def set_circumsphere(self):
        vertices = np.array(self.p1234)
        squared_dists = distance.pdist(vertices, metric='sqeuclidean')
        squared_dists_mat = distance.squareform(squared_dists)
        with_border = np.insert(np.insert(squared_dists_mat, 0, values=1, axis=1), 0, values=1, axis=0)
        np.fill_diagonal(with_border, 0)
        inv = np.linalg.inv(with_border)
        radius = np.sqrt(inv[0][0] / -2)
        barycentric_coodinates = inv[1:, 0]
        simplex = np.asarray(vertices)
        point = np.asarray(barycentric_coodinates)
        circumcenter = np.dot(simplex.T, point)
        self.radius = radius
        self.circumcenter = circumcenter

    def check_point_in_sphere(self, point):
        if not self.radius:
            print("set radius before running it!")
        point = np.asarray(point)
        if np.sum((self.circumcenter - point) ** 2) <= self.radius ** 2:
            return True
        else:
            return False