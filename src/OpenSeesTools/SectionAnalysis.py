import collections
import json

import numpy as np
import openseespy.opensees as ops

from .basic import OpenSeesAnalysis

__all__ = [
    'SectionAnalysis',
]


class SectionAnalysis(OpenSeesAnalysis):
    def __init__(self, sectionFactory, secTag=1, scratchPath=None,
                 analysisID=0):
        self.sectionFactory = sectionFactory
        self.secTag = int(secTag)
        super().__init__(scratchPath, analysisID)

    def getDiscretization(self):
        # Output files
        file_print = self.scratchFile('print_file', '.json')
        if file_print.exists():
            file_print.unlink()

        # Create model
        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 6)
        self.sectionFactory()
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 0.0, 0.0, 0.0)
        ops.element('zeroLengthSection', 1, 1, 2, self.secTag)
        ops.printModel('-JSON', '-file', str(file_print))
        ops.wipe()

        with open(file_print) as f:
            data = json.load(f)['StructuralAnalysisModel']['properties']

        section = data['sections'][0]
        nFibers = len(section['fibers'])
        fiberMat = np.empty(nFibers, dtype=int)
        fiberLocY = np.empty(nFibers)
        fiberLocZ = np.empty(nFibers)
        fiberArea = np.empty(nFibers)

        for i, fiber in enumerate(section['fibers']):
            fiberMat[i] = fiber['material']
            fiberLocY[i] = fiber['coord'][0]
            fiberLocZ[i] = fiber['coord'][1]
            fiberArea[i] = fiber['area']

        Discretization = collections.namedtuple(
            'Discretization',
            ['fiberMat', 'fiberLocY', 'fiberLocZ', 'fiberArea'])

        if self.deleteFiles:
            file_print.unlink()

        return Discretization(fiberMat, fiberLocY, fiberLocZ, fiberArea)

    def plotDiscretization(self, ax, plotAs2d=False):
        fiberMat, fiberLocY, fiberLocZ, _ = self.getDiscretization()

        C, ic = np.unique(fiberMat, return_inverse=True)
        nMats = len(C)

        if plotAs2d:
            fiberLocZ = np.arange(len(fiberLocY)) + 1

        ax.scatter(fiberLocZ, fiberLocY, 20)

    def printMaterialInfo(self, file=None):
        fprint = lambda *args, **kwargs: print(*args, **kwargs, file=file)
        fiberMat, fiberLocY, fiberLocZ, fiberArea = self.getDiscretization()

        fprint(
            '  Material  |  # Fibers  |    Area    |     Iz     |     Iy     ')
        fprint(
            '------------+------------+------------+------------+------------')
        uniqueFiberMat = np.unique(fiberMat)
        for uMat in uniqueFiberMat:
            ind = np.array(np.nonzero(fiberMat == uMat))
            partSectionArea = np.sum(fiberArea[ind])
            partSectionIy = np.sum(fiberArea[ind]*fiberLocZ[ind]**2)
            partSectionIz = np.sum(fiberArea[ind]*fiberLocY[ind]**2)
            fprint('{:-11d} |{:-11d} |{:-11.6g} |{:-11.6g} |{:-11.6g} '.format(
                uMat, ind.size, partSectionArea, partSectionIz, partSectionIy))
        fprint(
            '------------+------------+------------+------------+------------')
        sectionArea = fiberArea.sum()
        sectionIy = np.sum(fiberArea*fiberLocZ**2)
        sectionIz = np.sum(fiberArea*fiberLocY**2)
        fprint('    Total   |{:-11d} |{:-11.6g} |{:-11.6g} |{:-11.6g} '.format(
            fiberArea.size, sectionArea, sectionIz, sectionIy))
