"""Wide flange fiber sections.

The fiber sections in this module are constructed using the builder pattern,
with some initial parameters set and then methods called on the returned object
to finish setup::

>>> section = WSection2d.fromName('W14x53', secTag=1, nf=20, axis='strong')
>>> section.setMaterial('Elastic', 1, 29000.0)
<WSection2d(secTag=1, nf=20, axis='strong', d=13.9, tw=0.37, bf=8.06, tf=0.66, k=0.66) [Material: Elastic(tag=1, Es=29000.0)] [Residual stress: None] [Additional stiffness: None]>
"""
from __future__ import annotations

import abc

import attr
import numpy as np

from .basic import ops, patchRect2d, twoFiberSection, fourFiberSectionGJ

#===============================================================================
# Globals
#===============================================================================
__all__ = [
    'WSection2d',
]

WSECTION2D_BENDING_AXES = {'strong', 'weak'}


#===============================================================================
# Exceptions
#===============================================================================
class IncompleteSectionError(Exception):
    """Raised when trying to add an incomplete section to the domain."""
    pass


#===============================================================================
# Validators
#===============================================================================
def is_positive(instance, attribute, value):
    if value <= 0:
        raise ValueError(f'`{attribute.name}` must be positive')


#===============================================================================
# Builder objects
#===============================================================================
@attr.s(auto_attribs=True)
class LehighStressPattern():
    frc: float
    nSectors: int


@attr.s(auto_attribs=True)
class AddedStiffness2d():
    EA: float
    EI: float

    def create(self, matTag):
        ops.uniaxialMaterial('Elastic', matTag, 1.0)
        twoFiberSection(None, matTag, self.EA, self.EI)


@attr.s(auto_attribs=True)
class AddedStiffness3d():
    EA: float
    EIz: float
    EIy: float
    GJ: float

    def create(self, matTag):
        ops.uniaxialMaterial('Elastic', matTag, 1.0)
        fourFiberSectionGJ(None, matTag, self.EA, self.EIy, self.EIz, self.GJ)


#-------------------------------------------------------------------------------
# Materials
#-------------------------------------------------------------------------------
@attr.s(auto_attribs=True)
class MatTag():
    tag: int = attr.ib(converter=int)

    def create(self, fr):
        pass


@attr.s(auto_attribs=True)
class Elastic():
    tag: int = attr.ib(converter=int)
    Es: float = attr.ib(converter=float)

    def create(self, fr):
        ops.uniaxialMaterial('Elastic', tag, self.Es)


@attr.s(auto_attribs=True)
class ElasticPP():
    tag: int = attr.ib(converter=int)
    Es: float = attr.ib(converter=float)
    Fy: float = attr.ib(converter=float)

    def create(self, fr):
        epsyP = self.Fy/self.Es
        epsyN = -self.Fy/self.Es
        eps0 = -fr/self.Es
        ops.uniaxialMaterial('ElasticPP', self.tag, self.Es, epsyP, epsyN, eps0)


@attr.s(auto_attribs=True)
class Steel01():
    tag: int = attr.ib(converter=int)
    Es: float = attr.ib(converter=float)
    Fy: float = attr.ib(converter=float)
    b: float = attr.ib(converter=float)

    def create(self, fr):
        ops.uniaxialMaterial('Steel01', self.tag, self.Fy, self.Es, self.b)


@attr.s(auto_attribs=True)
class Steel02():
    tag: int = attr.ib(converter=int)
    Es: float = attr.ib(converter=float)
    Fy: float = attr.ib(converter=float)
    b: float = attr.ib(converter=float)

    def create(self, fr):
        ops.uniaxialMaterial('Steel02', self.tag, self.Fy, self.Es, self.b,
                             20.0, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0, fr)


WSECTION_MATERIAL_DISPATCH = {
    'MatTag': MatTag,
    'Elastic': Elastic,
    'ElasticPP': ElasticPP,
    'Steel01': Steel01,
    'Steel02': Steel02,
}

MATERIALS_WITH_RESIDUAL_STRESS_SUPPORT = [ElasticPP, Steel02]


#===============================================================================
# Sections
#===============================================================================
@attr.s(repr=False)
class AbstractWSection(abc.ABC):
    def __attrs_post_init__(self):
        self._material = None
        self._residualStress = None
        self._addedStiffness = None
        self._overrideGJ = None

    def __repr__(self):
        init_args = (f'{a.name}={getattr(self, a.name)!r}'
                     for a in self.__attrs_attrs__ if a.init)
        init_args_repr = ', '.join(init_args)

        repr_ = (f'<{self.__class__.__name__}('
                 f'{init_args_repr}'
                 ')'
                 f' [Material: {self._material}]'
                 f' [Residual stress: {self._residualStress}]'
                 f' [Additional stiffness: {self._addedStiffness}]'
                 '>')
        return repr_

    def create(self):
        """Create the section and materials in the current OpenSees domain.

        Only call this after calling all applicable builder functions.

        Returns
        -------
        lastTag : int
            The last material tag used by the section creator.

        Raises
        ------
        IncompleteSectionError
            - If the material for the section has not been set (`setMaterial`)
        """
        if self._material is None:
            raise IncompleteSectionError(f'{self.__class__.__name__}: '
                                         'cannot create section without '
                                         'material information')
        ops.section('Fiber', self.secTag, '-GJ', self.GJ)
        lastTag = self.createFibers()
        if self._addedStiffness is not None:
            lastTag += 1
            self._addedStiffness.create(lastTag)
        return lastTag

    @abc.abstractmethod
    def createFibers(self) -> int:
        pass

    #---------------------------------------------------------------------------
    # Additional section dimensions
    #---------------------------------------------------------------------------
    @property
    def dw(self):
        """Depth of the web."""
        return self.d - 2*self.tf

    @property
    def d1(self):
        """Half the depth of the web."""
        return self.dw/2

    @property
    def d2(self):
        """Half the overall depth."""
        return self.d/2

    @property
    def b1(self):
        """Half the thickness of the web."""
        return self.tw/2

    @property
    def b2(self):
        """Half the width of the flange."""
        return self.bf/2

    @property
    def GJ(self):
        if self._overrideGJ is None:
            G = self._material.Es/(2*(1 + 0.3))
            J = (2*self.bf*self.tf**3 + (self.d - 2*self.tf)*self.tw**2)/3
            self._overrideGJ = G*J
        return float(self._overrideGJ)

    @GJ.setter
    def GJ(self, value):
        self._overrideGJ = value

    #---------------------------------------------------------------------------
    # Builder pattern functions
    #---------------------------------------------------------------------------
    def setMaterial(self, material, *args, **kwargs):
        """Set the material used by the section.

        Parameters
        ----------
        material : {'MatTag', 'Elastic', 'ElasticPP', 'Steel01', 'Steel02'}
            Material to use.
        *args, **kwargs
            Arguments to the material constructor. (See Options section.)

        Options
        -------
        MatTag(tag)
            Use the uniaxial material defined by `tag`.
        Elastic(tag, Es)
            Elastic material.
        ElasticPP(tag, Es, Fy)
            Elastic-perfectly-plastic with optional residual stress pattern.
        Steel01(tag, Es, Fy, b)
            Bilinear steel model without residual stress.
        Steel02(tag, Es, Fy, b)
            GMP steel model with optional residual stress pattern.

        Example
        -------
        >>> section = WSection2d(...).setMaterial('Steel02', 1, 29000, 50, 0.03)
        """
        self._material = WSECTION_MATERIAL_DISPATCH[material](*args, **kwargs)
        return self

    def addLehigh(self, frc, nSectors):
        """Add a Lehigh residual stress pattern to the section.

        Parameters
        ----------
        frc : float
            Maximum residual compressive stress.
        nSectors : int
            Number of sectors used to represent the residual stress pattern. The
            section will create `nSectors` materials to handle the pattern; this
            is done by incrementing the tag given to `setMaterial`. Make sure
            there is enough room!
        """
        self._residualStress = LehighStressPattern(frc, nSectors)
        return self

    @abc.abstractmethod
    def addStiffness(self):
        """Add additional elastic stiffness to the section.

        This method depends on the dimensionality of the problem.
        """
        return self

    def setGJ(self, GJ):
        """Set the torsional stiffness of the section instead of calculating it.

        Parameters
        ----------
        GJ : float
            Torsional stiffness.
        """
        self.GJ = GJ
        return self


@attr.s(auto_attribs=True, repr=False)
class WSection2d(AbstractWSection):
    """Wide flange steel section suitable for two-dimensional analysis.

    Additional options are specified using the builder pattern, e.g.::

    >>> Es = 29000.0
    >>> Fy = 50.0
    >>> section = WSection2d(...).addLehigh(-0.30*Fy, 20).setMaterial('ElasticPP', Es, Fy)

    Call `create` when ready to add the section to the domain.

    Parameters
    ----------
    secTag : int
        Section tag.
    nf : int
        Number of fibers along the primary bending axis.
    axis : {'strong', 'weak'}
        Axis of primary bending. 'strong' creates a strong-axis section, 'weak'
        creates a weak-axis section.
    d : float
        Nominal depth.
    tw : float
        Web thickness.
    bf : float
        Flange width.
    tf : float
        Flange thickness.
    k : float, optional
        Fillet size. If None, uses `tf`. (default: None)

    Builder functions
    -----------------
    setMaterial(kind, *args, **kwargs)
        Define the material for the section.
    addLehigh(frc, nSectors) : optional
        Add a Lehigh residual stress pattern to the section. Only used if the
        material used supports initial stresses.
    addStiffness(EA, EI) : optional
        Add additional elastic stiffness to the section.
    """
    secTag: int = attr.ib(converter=attr.converters.optional(int))
    nf: int = attr.ib(validator=is_positive)
    axis: str = attr.ib(validator=attr.validators.in_(WSECTION2D_BENDING_AXES))
    d: float = attr.ib(validator=is_positive)
    tw: float = attr.ib(validator=is_positive)
    bf: float = attr.ib(validator=is_positive)
    tf: float = attr.ib(validator=is_positive)
    k: float = attr.ib(validator=attr.validators.optional(is_positive),
                       default=None)

    def __attrs_post_init__(self):
        if self.k is None:
            self.k = self.tf
        elif self.k < self.tf:
            raise ValueError('`k` should be larger than `tf`')
        super().__attrs_post_init__()

    #---------------------------------------------------------------------------
    # Other constructors
    #---------------------------------------------------------------------------
    @classmethod
    def fromName(cls, name, secTag, nf, axis, units='US'):
        """Look up shape data to create the section.

        Parameters
        ----------
        name : str
            Name of the shape to look up.
        secTag : int
            Section tag.
        nf : int
            Number of fibers along the primary bending axis.
        axis : {'strong', 'weak'}
            Axis of primary bending. 'strong' creates a strong-axis section,
            'weak' creates a weak-axis section.
        units : {'US', 'SI'}, optional
            Unit system to look up the name and dimensions in. 'US' uses inches
            for length, 'SI' uses millimeters. (default: 'US')
        """
        raise NotImplementedError()

    #---------------------------------------------------------------------------
    # Builder pattern functions
    #---------------------------------------------------------------------------
    def addStiffness(self, EA=0.0, EI=0.0):
        """Add additional elastic stiffness to the section.

        Parameters
        ----------
        EA : float, optional
            Additional elastic axial stiffness. (default: 0.0)
        EI : float, optional
            Additional elastic bending stiffness. (default: 0.0)
        """
        self._addedStiffness = AddedStiffness2d(EA, EI)
        return self

    #---------------------------------------------------------------------------
    # Fiber definition stuff
    #---------------------------------------------------------------------------
    def createFibers(self):
        if self._residualStress is not None and (
                self._material.__class__ in
                MATERIALS_WITH_RESIDUAL_STRESS_SUPPORT):
            lastTag = self.createFibersWithResidualStress()
        else:
            matTag = self._material.tag
            self._material.create(fr=0.0)
            if self.axis == 'strong':
                nff = np.ceil((self.nf/self.d)*self.tf)
                nfw = np.ceil((self.nf/self.d)*(self.d - 2*self.tf))
                patchRect2d(matTag, nff, self.bf, +self.d1, +self.d2)
                patchRect2d(matTag, nfw, self.tw, -self.d1, +self.d1)
                patchRect2d(matTag, nff, self.bf, -self.d2, -self.d1)
            elif self.axis == 'weak':
                nff = np.ceil((self.nf/self.bf)*(self.bf - self.tw)/2)
                nfw = np.ceil((self.nf/self.bf)*self.tw)
                patchRect2d(matTag, nff, 2*self.tf, +self.b1, +self.b2)
                patchRect2d(matTag, nfw, self.d, -self.b1, +self.b1)
                patchRect2d(matTag, nff, 2*self.tf, -self.b2, -self.b1)
            else:
                raise ValueError(f'Unknown bending axis: {axis}')
            lastTag = matTag

        if self.k > self.tf:
            self.createFillets(lastTag)
        return lastTag

    def createFibersWithResidualStress(self):
        frc = self._residualStress.frc
        if self.k > self.tf:
            r = self.k - self.tf
            fillet_area = (1 - 0.25*np.pi)*r**2
            frt = -frc*(self.bf*self.tf)/(self.bf*self.tf + self.tw*self.dw +
                                          4*fillet_area)
        else:
            frt = -frc*(self.bf*self.tf)/(self.bf*self.tf + self.tw*self.dw)
        nSectors = self._residualStress.nSectors

        # Create residual stress materials
        lehighStartTag = self._material.tag
        for i in range(nSectors):
            x = (i + 0.5)/nSectors
            fr = frc + x*(frt - frc)
            self._material.create(fr)
            self._material.tag += 1
        self._material.create(frt)
        self._material.tag += 1

        # Create the actual fibers
        if self.axis == 'strong':
            # Flanges
            bf1 = self.bf/nSectors
            nff = np.ceil(self.nf/self.d)*self.tf
            for i in range(nSectors):
                matTag = lehighStartTag + i
                patchRect2d(matTag, nff, bf1, +self.d1, +self.d2)
                patchRect2d(matTag, nff, bf1, -self.d2, -self.d1)
            # Web
            matTag = lehighStartTag + nSectors
            nfw = np.ceil((self.nf/self.d)*(self.d - 2*self.tf))
            patchRect2d(matTag, nfw, self.tw, -self.d1, +self.d1)
        elif self.axis == 'weak':
            # Flanges
            b21 = self.b2/nSectors
            nff = np.ceil((self.nf/self.bf)*b21)
            for i in range(nSectors):
                bft = self.b2 - i*b21
                matTag = lehighStartTag + i
                patchRect2d(matTag, nff, 2*self.tf, bft - b21, +bft)
                patchRect2d(matTag, nff, 2*self.tf, -bft, b21 - bft)
            # Web
            matTag = lehighStartTag + nSectors
            nfw = np.ceil((self.nf/self.bf)*self.tw)
            patchRect2d(matTag, nfw, self.dw, -self.b1, +self.b1)

        return matTag

    def createFillets(self, matTag):
        r = self.k - self.tf  # Radius of fillet
        area = (1 - 0.25*np.pi)*r**2  # Area of individual fillet
        ybar = 2/(12 - 3*np.pi)*r  # Centroid of fillet from flange face
        if self.axis == 'strong':
            ops.fiber(+self.d1 - ybar, 0.0, 2*area, matTag)
            ops.fiber(-self.d1 + ybar, 0.0, 2*area, matTag)
        elif self.axis == 'weak':
            ops.fiber(+self.b1 + ybar, 0.0, 2*area, matTag)
            ops.fiber(-self.b1 - ybar, 0.0, 2*area, matTag)
        else:
            raise ValueError(f'Unknown bending axis: {axis}')


@attr.s(auto_attribs=True)
class WSection3d(AbstractWSection):
    """Wide flange steel section suitable for two-dimensional analysis.

    Additional options are specified using the builder pattern, e.g.::

    >>> Es = 29000.0
    >>> Fy = 50.0
    >>> section = WSection2d(...).addLehigh(-0.30*Fy, 20).setMaterial('ElasticPP', Es, Fy)

    Parameters
    ----------
    secTag : int
        Section tag.
    nf1 : int
        Number of fibers along the strong axis.
    nf2 : int
        Number of fibers along the weak axis.
    d : float
        Nominal depth.
    tw : float
        Web thickness.
    bf : float
        Flange width.
    tf : float
        Flange thickness.
    k : float, optional
        Fillet size. If None, uses `tf`. (default: None)

    Builder functions
    -----------------
    setMaterial(kind, *args, **kwargs)
        Define the material for the section.
    addLehigh(frc, nSectors) : optional
        Add a Lehigh residual stress pattern to the section. Only used if the
        material used supports initial stresses.
    addStiffness(EA, EIz, EIy, GJ) : optional
        Add additional elastic stiffness to the section.
    """
    secTag: int = attr.ib(converter=attr.converters.optional(int))
    nf1: int = attr.ib(validator=is_positive)
    nf2: str = attr.ib(validator=is_positive)
    d: float = attr.ib(validator=is_positive)
    tw: float = attr.ib(validator=is_positive)
    bf: float = attr.ib(validator=is_positive)
    tf: float = attr.ib(validator=is_positive)
    k: float = attr.ib(validator=attr.validators.optional(is_positive),
                       default=None)

    def __attrs_post_init__(self):
        if self.k is None:
            self.k = self.tf
        elif self.k < self.tf:
            raise ValueError('`k` should be larger than `tf`')
        super().__attrs_post_init__()

    def createFibers(self):
        if self._residualStress is not None and (
                self._material.__class__ in
                MATERIALS_WITH_RESIDUAL_STRESS_SUPPORT):
            lastTag = self.createFibersWithResidualStress()
        else:
            matTag = self._material.tag
            nff1 = np.ceil((self.nf1/self.d)*self.tf)
            nfw1 = np.ceil((self.nf1/self.d)*(self.d - 2*self.tf))
            nff2 = np.ceil((self.nf2/self.bf)*self.bf)
            ops.patch('quad', matTag, )
            lastTag = matTag

        if self.k > self.tf:
            self.createFillets(lastTag)
        return lastTag

    def createFibersWithResidualStress(self):
        frc = self._residualStress.frc
        frt = -frc*(self.bf*self.tf)/(self.bf*self.tf + self.tw*self.dw)
        nSectors = self._residualStress.nSectors

        # Create residual stress materials
        lehighStartTag = self._material.tag
        for i in range(nSectors):
            x = (i + 0.5)/nSectors
            fr = frc + x*(frt - frc)
            self._material.create(fr)
            self._material.tag += 1
        self._material.create(frt)
        self._material.tag += 1

    def createFillets(self, matTag):
        r = self.k - self.tf  # Radius of fillet
        area = (1 - 0.25*np.pi)*r**2  # Area of individual fillet
        ybar = 2/(12 - 3*np.pi)*r  # Centroid of fillet from flange face
        ops.fiber(+self.d1 - (r - ybar), +self.b1 + (r - ybar), area, matTag)
        ops.fiber(-self.d1 + (r - ybar), +self.b1 + (r - ybar), area, matTag)
        ops.fiber(+self.d1 - (r - ybar), -self.b1 - (r - ybar), area, matTag)
        ops.fiber(-self.d1 + (r - ybar), -self.b1 - (r - ybar), area, matTag)
