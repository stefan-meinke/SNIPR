from setuptools import setup, find_packages

setup(
    name='ptc_analysis',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'pandas',
        'gffutils'
    ],
    entry_points={
        'console_scripts': [
            'ptc_analysis=ptc_analysis.runner:main'
        ]
    },
    author='Stefan Meinke',
    description='PTC/NMD analysis tool from rMATS splicing events',
    license='MIT',
    python_requires='>=3.6'
)
