#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk



#########################################

version=0.1








from setuptools import setup

setup(name='hichip_tool',
    version=version,
    description='tool to find peaks from hichip data',
    url='https://github.com/ChenfuShi/domain_caller_site',
    author='Chenfu Shi',
    author_email='chenfu.shi@postgrad.manchester.ac.uk',
    license='MIT',
    packages=['scripts','utils'],
    install_requires=[
        'scipy',
        'numpy',
        'os',
    ],
    entry_points = {
        'console_scripts': 
        ['peak_call=scripts.main:main',

        ],
    },
    scripts=['utils/bedgraph.py',
            'utils/IPquality.py',
        ],

    zip_safe=False
    
    )


def readme():
    with open('README.md') as f:
        return f.read()