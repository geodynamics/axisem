#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

layers = [1217.5, 1190., 1160., 1100.]
angles = [[45., 55.], [50., 60.], [55., 65.]]

vp_in = np.ones((len(layers) - 1, len(angles)))

vp_in[0,0] = 0.
vp_in[0,1] = 2.

vp_in[1,0] = 0.
vp_in[1,1] = 5.

vp_in[2,0] = 0.
vp_in[2,1] = 2.


ntheta = 721
dtheta = 180. / (ntheta - 1)

nlayers = len(layers) - 1
nlpl = 2
dr = .01

f = open('model.sph', 'w')

vp = 0.
vs = 0.
rho = 0.

npoints = (nlayers * nlpl + 1)  * ntheta
print >> f, npoints


#r = layers[0] + dr
#for theta in np.linspace(0., 180., ntheta):
#    print >> f, '%7.2f %6.2f %5.2f %5.2f  %5.2f ' % (r, theta, vp, vs, rho)

for l in np.arange(nlayers):
    for r in np.linspace(layers[l] - dr, layers[l+1] + dr, nlpl):
        for theta in np.linspace(0., 180., ntheta):
            if theta < angles[l][0]:
                vp = vp_in[l,0]
            elif theta > angles[l][1]:
                vp = vp_in[l,1]
            else:
                vp = vp_in[l,0] \
                   + (vp_in[l,1] - vp_in[l,0]) / (angles[l][1] - angles[l][0]) \
                        * (theta - angles[l][0])

            print >> f, '%7.2f %6.2f %5.2f %5.2f  %5.2f ' % (r, theta, vp, vs, rho)

vp = 0.
r = layers[-1] - dr
for theta in np.linspace(0., 180., ntheta):
    print >> f, '%7.2f %6.2f %5.2f %5.2f  %5.2f ' % (r, theta, vp, vs, rho)

f.close
