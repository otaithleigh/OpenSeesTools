import matplotlib.pyplot as plt
from tabulate import tabulate

from openseestools import WSection2d, SectionAnalysis

sections = {
    'W14x53': {
        'd': 13.9,
        'bf': 8.06,
        'tw': 0.370,
        'tf': 0.660,
        'k': 1.25,
    },
    'W36x302': {
        'd': 37.3,
        'bf': 16.7,
        'tw': 0.945,
        'tf': 1.68,
        'k': 2.63,
    },
}
expected = {
    'W14x53': {
        'Ag': 15.6,
        'Ix': 541,
        'Iy': 57.7,
    },
    'W36x302': {
        'Ag': 89.0,
        'Ix': 21100,
        'Iy': 1300,
    },
}

section = 'W36x302'

# A992
Es = 29000.0
Fy = 50.0

# Section settings
mattag = 1
sectag = 1
nfibers = 20

frc = -0.30 * Fy
nsectors = 4

floatfmt = '#.3g'

# ===============================================================================
# Create sections
# ===============================================================================
strong_section = (
    WSection2d(sectag, nfibers, 'strong', **sections[section])
    .setMaterial('ElasticPP', mattag, Es, Fy)
    .addLehigh(frc, nsectors)
)
weak_section = (
    WSection2d(sectag, nfibers, 'weak', **sections[section])
    .setMaterial('ElasticPP', mattag, Es, Fy)
    .addLehigh(frc, nsectors)
)


def create_section(section: WSection2d):
    section.create()
    return section.secTag


# ===============================================================================
# Do analysis
# ===============================================================================
A_expected = {'strong': expected[section]['Ag'], 'weak': expected[section]['Ag']}
A_model = {}
Iz_expected = {'strong': expected[section]['Ix'], 'weak': expected[section]['Iy']}
Iz_model = {}
for section in [strong_section, weak_section]:
    analysis = SectionAnalysis(lambda: create_section(section))
    discretization = analysis.getDiscretization()
    print(section.axis.title() + ':')
    analysis.printMaterialInfo(floatfmt=floatfmt)
    analysis.plotDiscretization(plotAs2d=True)
    print()
    print()

    A_model[section.axis] = discretization.getArea()
    Iz_model[section.axis] = discretization.getIz()


# ===============================================================================
# Build table
# ===============================================================================
def error(reference, actual):
    return abs((actual - reference) / reference)


rows = []
headers = [
    'Bending axis',
    'Ag_manual\n(in^2)',
    'Ag_model\n(in^2)',
    'Error\n(%)',
    'Iz_manual\n(in^4)',
    'Iz_model\n(in^4)',
    'Error\n(%)',
]
tablefmt = 'presto'
for axis in ['strong', 'weak']:
    Ae = A_expected[axis]
    Am = A_model[axis]
    Ie = Iz_expected[axis]
    Im = Iz_model[axis]
    eA = 100 * error(Ae, Am)
    eI = 100 * error(Ie, Im)

    rows.append([axis, Ae, Am, eA, Ie, Im, eI])

table = tabulate(rows, headers, tablefmt, floatfmt=floatfmt)
print(table)
plt.show()
