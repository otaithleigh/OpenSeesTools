import typing

import attr

from .basic import ops, getClassLogger

__all__ = [
    'VariableAnalysis',
]


@attr.s(auto_attribs=True)
class VariableAnalysis():
    """Analysis wrapper with support for multiple tolerances and algorithms.

    Parameters
    ----------
    test : str
        TestMethod to use.
    maxiters : int
        Maximum iterations for the test method.
    tolerances : List[float]
        List of tolerances to try.
    algorithms : List[str]
        List of algorithms to try.
    pflag : int, optional
        Print flag for the test method. (default: 0)
    norm : int
        Norm type for the test method. (default: 2)

    Example
    -------
    >>> analysis = VariableAnalysis('NormDispIncr', 20, tolerances=[1e-6, 1e-5],
                                    algorithms=['Newton', 'BFGS'])
    >>> error = analysis.analyze()
    >>> if error:
    ...     raise ops.OpenSeesError('analyze returned {}'.format(error))
    """
    test: str
    maxiters: int
    tolerances: typing.List[float]
    algorithms: typing.List[str]
    pflag: int = 0
    norm: int = 2

    def __attrs_post_init__(self):
        self.logger = getClassLogger(self.__class__)

    def analyze(self, *args):
        """Run the analysis.

        Currently, the "variable" part of a VariableAnalysis isn't terribly
        smart; it simply iterates through the tolerances and algorithms it has
        been given, resetting the next time `analyze` is called. If the analysis
        step isn't completed successfully after trying all tolerances and
        algorithms, it gives up.

        Parameters
        ----------
        *args
            Arguments to `ops.analyze`.

        Returns
        -------
        error : int
            Return code from `ops.analyze`. 0 indicates success.
        """
        for tol in self.tolerances:
            ops.test(self.test, tol, self.maxiters, self.pflag, self.norm)
            for alg in self.algorithms:
                self.logger.debug(f'Running analysis with tolerance {tol:g}'
                                  f' and algorithm {alg!r}')
                ops.algorithm(alg)
                error = ops.analyze(*args)
                if error == 0:
                    break
            if error == 0:
                break
        return error
