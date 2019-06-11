#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk



#########################################

# to install without copying files
# pip install -e .
# to create distribution packages python setup.py sdist bdist_wheel
# python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# twine upload dist/*

version=0.1

with open("README.md", "r") as fh:
    long_description = fh.read()


from setuptools import setup

setup(name='hichip_peaks',
    version=version,
    description='tool to find peaks from hichip data',
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