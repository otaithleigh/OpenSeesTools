import logging
import sys

from openseestools import opensees as ops
from openseestools import captureOutput

logging.basicConfig(format='%(levelname)s:%(name)s.%(funcName)s: %(message)s',
                    level='DEBUG')


@captureOutput
def test_elasticTrussAnalysis():
    # ------------------------------
    # Start of model generation
    # -----------------------------
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)

    # create nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, 144.0, 0.0)
    ops.node(3, 168.0, 0.0)
    ops.node(4, 72.0, 96.0)

    # set boundary condition
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(3, 1, 1)

    # define materials
    ops.uniaxialMaterial('Elastic', 1, 3000.0)

    # define elements
    ops.element('Truss', 1, 1, 4, 10.0, 1)
    ops.element('Truss', 2, 2, 4, 5.0, 1)
    ops.element('Truss', 3, 3, 4, 5.0, 1)

    # create TimeSeries
    ops.timeSeries('Linear', 1)

    # create a plain load pattern
    ops.pattern('Plain', 1, 1)

    # Create the nodal load - command: load nodeID xForce yForce
    ops.load(4, 100.0, -50.0)

    # ------------------------------
    # Start of analysis generation
    # ------------------------------

    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.test('NormDispIncr', 1e-6, 10, 2)
    ops.algorithm('Newton')
    ops.analysis('Static')

    print('AAAAAA')
    ops.analyze(1)

    ux = ops.nodeDisp(4, 1)
    uy = ops.nodeDisp(4, 2)

    assert abs(ux - 0.53009277713228375450) < 1e-12
    assert abs(uy + 0.17789363846931768864) < 1e-12
    return ux, uy


ux, uy = test_elasticTrussAnalysis()
print(f'ux = {ux:+g}')
print(f'uy = {uy:+g}')
print('output:')
print(test_elasticTrussAnalysis.stderr.getvalue())
