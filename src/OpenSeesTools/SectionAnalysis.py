import typing
import json

import attr
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

from .basic import ops, OpenSeesAnalysis

__all__ = [
    'SectionAnalysis',
]


@attr.s(auto_attribs=True, slots=True)
class SectionDiscretization():
    """Section discretized into constituent fibers.

    Use `SectionAnalysis.getDiscretization()` instead of creating discretized
    sections directly.

    Parameters
    ----------
    fiberMat : np.ndarray(dtype=int)
        Material tags of the section fibers.
    fiberLocY : np.ndarray
        Y-coordinates of the section fibers.
    fiberLocZ : np.ndarray
        Z-coordinates of the section fibers.
    fiberArea : np.ndarray
        Areas of the section fibers.
    """
    fiberMat: np.ndarray
    fiberLocY: np.ndarray
    fiberLocZ: np.ndarray
    fiberArea: np.ndarray

    def getArea(self):
        """Return the total area of the section."""
        return self.fiberArea.sum()

    def getIz(self):
        """Return the total moment of inertia about the Z-axis."""
        return np.sum(self.fiberArea*self.fiberLocY**2)

    def getIy(self):
        """Return the total moment of inertia about the Y-axis."""
        return np.sum(self.fiberArea*self.fiberLocZ**2)


# Index of the header row separator for different tabulate formats.
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

# Index to insert the copied separator into for different tabulate formats.
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
        self._sectionFactory = None
        self.sectionFactory = sectionFactory
        self._secTag = None
        self.secTag = secTag
        self._cachedDiscretization = None
        super().__init__(scratchPath, analysisID)

    @property
    def sectionFactory(self):
        """Function that creates the section to be analyzed when called with no
        arguments."""
        return self._sectionFactory

    @sectionFactory.setter
    def sectionFactory(self, factory):
        # Clear the discretization cache if settings change.
        if self.sectionFactory != factory:
            self._cachedDiscretization = None
        self._sectionFactory = factory

    @property
    def secTag(self):
        """Tag used by the section created by `sectionFactory`."""
        return self._secTag

    @secTag.setter
    def secTag(self, tag):
        # Clear the discretization cache if settings change.
        if self.secTag != int(tag):
            self._cachedDiscretization = None
        self._secTag = int(tag)

    def getDiscretization(self) -> SectionDiscretization:
        """Get the discretization for the section.

        Note that the discretization is cached; this cache is purged when
        changing `sectionFactory` or `secTag`.
        """
        if self._cachedDiscretization is None:
            self.logger.debug('Cache out of date: regenerating discretization')
            self._cachedDiscretization = self._getDiscretization()
        return self._cachedDiscretization

    def _getDiscretization(self) -> SectionDiscretization:
        """Non-caching version of `getDiscretization`."""
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

        return SectionDiscretization(fiberMat, fiberLocY, fiberLocZ, fiberArea)

    def plotDiscretization(self, ax=None, plotAs2d=False):
        """Plot a representation of the discretized section.

        Parameters
        ----------
        ax : optional
            Matplotlib axes to plot on. If None, a new pyplot figure is created.
        plotAs2d : bool, optional
            If True, the x-axis becomes the fiber number instead of the Z-axis
            dimension. (default: False)

        Returns
        -------
        ax
            Matplotlib axes that were plotted.
        """
        disc = self.getDiscretization()

        uniqueMats, _ = np.unique(disc.fiberMat, return_inverse=True)
        nMats = len(uniqueMats)

        if plotAs2d:
            disc.fiberLocZ = np.arange(len(disc.fiberLocY)) + 1
            xlabel = 'Fiber number'
            ylabel = 'Z-coordinate'
        else:
            xlabel = 'Z-coordinate'
            ylabel = 'Y-coordinate'

        if ax is None:
            _, ax = plt.subplots()
        ax.scatter(disc.fiberLocZ, disc.fiberLocY, 20)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def printMaterialInfo(self, file=None, tablefmt='presto', floatfmt='g'):
        """Print the material information for the section.

        Parameters
        ----------
        file : optional
            Open file-like descriptor to print to. (default: None)
        tablefmt : str, optional
            Table format to use. See `tabulate.tabulate_formats` for options.
            (default: 'presto')
        floatfmt : str, optional
            Format to use for floating point numbers. (default: 'g')
        """
        disc = self.getDiscretization()

        headers = ['Material', '# Fibers', 'Area', 'Iz', 'Iy']
        rows = []

        uniqueFiberMat = np.unique(disc.fiberMat)
        for uMat in uniqueFiberMat:
            ind = np.array(np.nonzero(disc.fiberMat == uMat))
            partArea = np.sum(disc.fiberArea[ind])
            partIy = np.sum(disc.fiberArea[ind]*disc.fiberLocZ[ind]**2)
            partIz = np.sum(disc.fiberArea[ind]*disc.fiberLocY[ind]**2)
            rows.append([uMat, ind.size, partArea, partIz, partIy])

        sectionArea = disc.getArea()
        sectionIy = disc.getIy()
        sectionIz = disc.getIz()
        rows.append(
            ['Total', disc.fiberArea.size, sectionArea, sectionIz, sectionIy])

        table = tabulate(rows, headers, tablefmt, floatfmt, colalign=['right'])
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
