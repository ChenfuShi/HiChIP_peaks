#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk



#########################################

# to install without copying files
# pip install -e .
# to create distribution packages python setup.py sdist bdist_wheel
# python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# twine upload dist/*

__version__="0.1.2"

# Read contents of long description (same as github readme)
with open("README.md", "r") as fh:
    long_description = fh.read()
# Update version number in __init__.py, This allows people to get the currently installed version from files when they import
with open("./hichip_peaks/__init__.py","w") as init_file:
    init_file.write(f"__version__=\"{__version__}\"")

from setuptools import setup

setup(name='hichip_peaks',
    version=__version__,
    description='A tool to find peaks from hichip data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/ChenfuShi/HiChIP_peaks',
    author='Chenfu Shi',
    author_email='chenfu.shi@postgrad.manchester.ac.uk',
    license='BSD',
    packages=['hichip_peaks'],
    install_requires=[
        'scipy',
        'numpy',
        'statsmodels',
        'matplotlib'
    ],
    entry_points = {
        'console_scripts': 
        ['peak_call=hichip_peaks.main:main',
        'make_bedgraph=hichip_peaks.bedgraph:main',
        'diff_peaks=hichip_peaks.diffpeaks:main',

        ],
    },
    python_requires='>=3.6',
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: BSD License',
        ],
    )


def readme():
    with open('README.md') as f:
        return f.read()