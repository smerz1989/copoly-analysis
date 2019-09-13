from setuptools import setup

def readme():
    with open('README.md','r') as read_file:
        return(read_file.read())

setup(name='copoly',
        version='0.1',
        description='Script for launchin and analyzing LAMMPS copolymer simulations',
        long_description=readme(),
        long_description_content_type="text/markdown",
        url='https://github.com/smerz1989/copoly-analysis',
        author='Steven Merz',
        packages=['copoly'],
        setup_requires=[
            'pytest-runner',
            ],
        tests_require=[
            'pytest',
            ],
        install_requires=[
            'matplotlib',
            'seaborn',
            'tqdm',
            'networkx',
            'numpy',
            'paramiko',
            'pandas',
            'matplotlib',
            'python_igraph',
            'mdtraj'
            ],
        include_package_data=True,
        zip_safe=False)
