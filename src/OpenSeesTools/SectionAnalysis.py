import collections
import json

import numpy as np
import openseespy.opensees as ops
from tabulate import tabulate

from .basic import OpenSeesAnalysis

__all__ = [
    'SectionAnalysis',
]

Discretization = collections.namedtuple(
    'Discretization', ['fiberMat', 'fiberLocY', 'fiberLocZ', 'fiberArea'])

_tabulate_header_row_sep_loc = {
    'fancy_grid': None,
    'github': 1,
    'grid': None,
    'html': None,
    'jira': None,
    'latex': 3,
    'latex_booktabs': 3,
    'latex_raw': 3,
    'mediawiki': None,
    'moinmoin': None,
    'orgtbl': 1,
    'pipe': 1,
    'plain': None,
    'presto': 1,
    'psql': 2,
    'rst': None,
    'simple': 1,
    'textile': None,
    'tsv': None,
    'youtrack': None,
}

_tabulate_bottom_row_sep_loc = {
    'fancy_grid': None,
    'github': -1,
    'grid': None,
    'html': None,
    'jira': None,
    'latex': -3,
    'latex_booktabs': -3,
    'latex_raw': -3,
    'mediawiki': None,
    'moinmoin': None,
    'orgtbl': -1,
    'pipe': -1,
    'plain': None,
    'presto': -1,
    'psql': -2,
    'rst': None,
    'simple': -1,
    'textile': None,
    'tsv': None,
    'youtrack': None,
}


class SectionAnalysis(OpenSeesAnalysis):
    """
    Parameters
    ----------
    sectionFactory : function
        Function that creates the section to be analyzed when called with no
        arguments.
    secTag : int, optional
        Tag used by the section created by `sectionFactory`. (default: 1)
    scratchPath : path_like, optional
        Path to the scratch directory. If None, uses the system temporary
        directory. (default: None)
    analysisID : optional
        Unique ID for the analysis. (default: 0)
    """

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

    def printMaterialInfo(self, tablefmt='presto', file=None):
        """Print the material information for the section.

        Parameters
        ----------
        file : optional
            Open file-like descriptor to print to. (default: None)
        """
        fiberMat, fiberLocY, fiberLocZ, fiberArea = self.getDiscretization()

        headers = ['Material', '# Fibers', 'Area', 'Iz', 'Iy']
        rows = []

        uniqueFiberMat = np.unique(fiberMat)
        for uMat in uniqueFiberMat:
            ind = np.array(np.nonzero(fiberMat == uMat))
            partSectionArea = np.sum(fiberArea[ind])
            partSectionIy = np.sum(fiberArea[ind]*fiberLocZ[ind]**2)
            partSectionIz = np.sum(fiberArea[ind]*fiberLocY[ind]**2)
            rows.append(
                [uMat, ind.size, partSectionArea, partSectionIz, partSectionIy])

        sectionArea = fiberArea.sum()
        sectionIy = np.sum(fiberArea*fiberLocZ**2)
        sectionIz = np.sum(fiberArea*fiberLocY**2)
        rows.append(
            ['Total', fiberArea.size, sectionArea, sectionIz, sectionIy])

        table = tabulate(rows, headers, tablefmt=tablefmt, colalign=['right'])
        # Hack in a bottom separator
        header_sep_loc = _tabulate_header_row_sep_loc[tablefmt]
        if header_sep_loc is not None:
            bottom_sep_loc = _tabulate_bottom_row_sep_loc[tablefmt]
            tablelines = table.splitlines()
            table = '\n'.join([
                *tablelines[:bottom_sep_loc], tablelines[header_sep_loc],
                *tablelines[bottom_sep_loc:]
            ])
        print(table, file=file)
