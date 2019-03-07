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
    packages=['hichip_tool','hichip_utils'],
    install_requires=[
        'scipy',
        'numpy',
    ],
    entry_points = {
        'console_scripts': 
        ['peak_call=hichip_tool.main:main',
        'make_bedgraph=hichip_utils.bedgraph:main',
        'IPquality=hichip_utils.IPquality:main',

        ],
    },
    python_requires='>=3.6',
    zip_safe=False
    
    )


def readme():
    with open('README.md') as f:
        return f.read()