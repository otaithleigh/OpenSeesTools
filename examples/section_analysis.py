import matplotlib.pyplot as plt
import numpy as np

from openseestools import SectionAnalysis
from openseestools.basic import ops

#-------------------------------------------------------------------------------
# Units
#-------------------------------------------------------------------------------
inch = 1.0
kip = 1.0
sec = 1.0

ft = 12*inch
lbf = kip/1000
ksi = kip/inch**2
psi = lbf/inch**2

#-------------------------------------------------------------------------------
# Section definition (W14x159)
#-------------------------------------------------------------------------------
d = 15.0*inch
tw = 0.745*inch
bf = 15.6*inch
tf = 1.19*inch

nfx = 40
nfy = 40

Fy = 50.0*ksi
E = 29000.0*ksi
Eh = E/1000
H = E*Eh/(E - Eh)
G = 11200*ksi
J = 19.7*inch**4


def createSection():
    ops.uniaxialMaterial('Hardening', 1, E, Fy, H, 0.0)
    ops.section('Fiber', 1, '-GJ', G*J)
    # Flanges
    nX = nfx
    nY = int(np.ceil(nfy*tf/d))
    ops.patch('rect', 1, nY, nX, 0.5*d - tf, -0.5*bf, 0.5*d, 0.5*bf)
    ops.patch('rect', 1, nY, nX, -0.5*d, -0.5*bf, -0.5*d + tf, 0.5*bf)

    # Web
    nX = int(np.ceil(nfx*tw/bf))
    nY = int(np.ceil(nfx*(d - 2*tf)/bf))
    ops.patch('rect', 1, nY, nX, -0.5*d + tf, -0.5*tw, 0.5*d - tf, 0.5*tw)


#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------
analysis = SectionAnalysis(createSection, secTag=1)
print('W14x159:')
print('     A = 46.7 in^2')
print('    Ix = 1900 in^4')
print('    Iy =  748 in^4')
print()
analysis.printMaterialInfo()

fig, ax = plt.subplots()
analysis.plotDiscretization(ax)

fig.tight_layout()
fig.savefig('section.svg')
