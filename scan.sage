import sys

sys.path.append('/home/inaki/polytopes/scan')

import scan


def verify_face(ref_representation, vertices):
    l = LatticePolytope(vertices)
    m = l._sublattice_polytope.normal_form()
    
    ref_representation = [vector(v) for v in ref_representation]
    for v0 in ref_representation:
        rep = [v-v0 for v in ref_representation]
        if LatticePolytope(rep).normal_form() == m:
            return True

    return False

def poly_cb(vertices):
    return Polyhedron(vertices).integral_points()

# We can use our own face point enumeration code if we prefer, just
# pass this as the face_cb parameter to begin_scan below.
def face_cb(vertices):
    l = LatticePolytope(vertices)
    return l.points().columns()

################## Scanning code from here #######################

import cPickle
import os.path

def scan_nvertices(nvertices):
    if isinstance(nvertices, int):
        nvertices = [nvertices]

    # Directory where to save the resulting pickles
    pickledir = "/home/inaki/polytopes/scan/pickles"
    if not os.path.isdir(pickledir):
        os.mkdir(pickledir)
    
    for n in nvertices:
        if n not in range(5,28):
            print "Number of vertices %d should be in the interval [5,27]" % (n,)
            continue
        print "Scanning polytopes with %d vertices" % (n,)
        result = scan.begin_scan(0, -1, "/home/inaki/cvs/palp/class.x",\
                                 "/home/inaki/polytopes/4d/zzdb-partial-%02d" % (n,),\
                                 "/home/inaki/cvs/palp/poly.x")
        fname = os.path.join(pickledir, "v%02d.pickle" % (n,))
        print "Done, saving result to", fname
        fd = open(fname, "w")
        cPickle.dump(result, fd, cPickle.HIGHEST_PROTOCOL)
        fd.close()

def n_entries(a,b,l):
    if (a,b) == (0,0):
        return 0
    if a<b:
        return l[(a,b)]
    return l[(b,a)]

scan_nvertices([5, 27])

# zzz TODO: scan over the last 100k
#result = scan.begin_scan(0, -1, "/home/inaki/cvs/palp/class.x", "/home/inaki/polytopes/4d/zzdb-partial-25", "/home/inaki/cvs/palp/poly.x")

#print str(result)
