from setuptools import setup

with open('README.rst') as file:
    long_description = file.read()

setup(
    name='openseespy-ext',
    version='0.0.0',

    description='Extensions to OpenSeesPy.',
    long_description=long_description,
    long_description_content_type='text/x-rst',

    package_dir={'': 'src'},
    packages=['openSeesComposite'],

    python_requires='>=3.5',
    install_requires=['openseespy', 'numpy'],

    author='Peter Talley, Mark D. Denavit',
    author_email='ptalley2@vols.utk.edu',
    url='https://github.com/otaithleigh/openseespy-ext',

    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering'
    ]
)
