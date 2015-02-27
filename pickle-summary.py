#!/usr/bin/python

import cPickle
import os.path

# Directory where the pickles have been saved. Can be an absolute
# path, such as /home/inaki/polytopes/scan/pickles, or a path relative
# to where this script is running.
pickledir = "pickles"

total = 0
Laba = {}
dP0_cube = 0
conifolds = 0
Zn_orbifolds = []
delPezzos = [0]*4
Ypq = []
Histogram = []
very_singular = []
Z2Z2_orbifolds = 0
dP0dP0dP1_models = 0

# Add two lists element-wise. If the lists are of different length,
# pad the shortest list with 0s first.
def add_lists(a,b):
    a = list(a) + [0]*(len(b) - len(a))
    b = list(b) + [0]*(len(a) - len(b))
    return [a[i]+b[i] for i in range(0,len(a))]

for n in range(5,28):
    fname = os.path.join(pickledir, "v%02d.pickle" % (n,))
    fd = open(fname)

    result = cPickle.load(fd)
    
    for key, value in result["Number of L^{aba} cones"].items():
        try:
            Laba[key] += value
        except KeyError:
            Laba[key] = value
    dP0_cube += result["Number of dP0^3 models"]
    conifolds += result["Number of conifolds"]
    Zn_orbifolds = add_lists(Zn_orbifolds,\
                             result["Number of CxC^2/Z_n orbifolds"])
    for n in range(0,4):
        delPezzos[n] += result["Number of del Pezzos"][n]
    partial_Ypq = result["Number of Ypq"]
    if len(partial_Ypq) > len(Ypq):
        for p in range(len(Ypq),len(partial_Ypq)):
            Ypq.append([0]*p)
    for p in range(0,len(partial_Ypq)):
        for q in range(0,p):
            Ypq[p][q]+=partial_Ypq[p][q]
    Histogram = add_lists(Histogram,\
                          result["Interior Point Histogram"])
    for rel_id,n_interior in result["Polytopes with very singular faces"]:
        very_singular.append((rel_id+total,n_interior))
    Z2Z2_orbifolds += result["Number of C^3/(Z2xZ2) orbifolds"]
    dP0dP0dP1_models += result["Number of dP0^2 x dP1 models"]

    n_partial = result["Scan range"][1]+1
    print "Loaded %s: %d polytopes, %d -> %d." %\
          (fname, n_partial, total, total+n_partial-1)

    total += n_partial
    
result = {}
result["Number of L^{aba} cones"] = Laba
result['Number of dP0^3 models'] = dP0_cube
result['Number of conifolds'] = conifolds
result['Number of CxC^2/Z_n orbifolds'] = Zn_orbifolds
result['Number of del Pezzos'] = delPezzos
result['Number of Ypq'] = Ypq
result['Interior Point Histogram'] = Histogram
result['Polytopes with very singular faces'] = very_singular
result['Number of C^3/(Z2xZ2) orbifolds'] = Z2Z2_orbifolds
result['Scan range'] = (0,total-1)
result['Number of dP0^2 x dP1 models'] = dP0dP0dP1_models

fname = os.path.join(pickledir, "summary.pickle")
fd = open(fname, "w")
cPickle.dump(result, fd)
print "Saved summary in", fname
