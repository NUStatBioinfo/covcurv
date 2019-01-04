from setuptools import setup

setup(
    name='covcurv',
    version='0.1.0',
    packages=['covcurv', 'covcurv.tests'],
    entry_points={
        'console_scripts': ['covcurv=covcurv.__main__:main',
                            'covcurv_test=covcurv.tests.__test__:main',
                            'covcurv_app=covcurv.__main_server__:main'],
    },
    package_data={
        'covcurv': ['resources/*',
                    'resources/shiny/*',
                    'tests/data/*']
    },
    license='MIT',
    url='https://github.com/ffineis/DegNorm',
    author='Frank Fineis',
    author_email='frankfineis2022@u.northwestern.edu',
    description='covcurv: parse RNA-seq samples into gene coverage matrices',
    classifiers=[
        "Programming Language :: Python",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"]
)