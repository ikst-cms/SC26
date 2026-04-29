import setuptools
import os

with open("README.md", "r") as fh:
	long_description = fh.read()

# C Module for Monte Carlo Simulation
C_monte_carlo_module = setuptools.Extension("Lacos.CPyMonC.lib.C_MCinterface", ['Lacos/CPyMonC/C_Functions_PyMonC/C_MCinterface.c'],
					extra_compile_args=['-fopenmp', '-O3', '-shared', '-fPIC', '-DOLD_HASH_KEY', '-DOPT_PSITE','-CTIME'],
					extra_link_args=['-lgomp']
					) 

# C Module for Monte Carlo Simulation(MPI)
C_mpi_monte_carlo_module = setuptools.Extension("Lacos.CPyMonC.lib.MPI_C_MCinterface", ['Lacos/CPyMonC/C_Functions_PyMonC/MPI_C_MCinterface.c'],
					extra_compile_args=['-O3', '-shared', '-fPIC', '-DOLD_HASH_KEY','-CTIME']
					)
with open('VERSION') as version_file:
	version = version_file.read().strip()

setuptools.setup(
	name="Lacos",
	version=version,
	author="IKST",
	author_email="info@ikst.res.in",
	description="Lattice Atomic Configuration Software Package",
	long_description=long_description,
	url="www.ikst.res.in",
	include_package_data=True,
	packages=setuptools.find_packages(),		 
        classifiers=(
		"Programming Language :: Python :: 3",
                "Intended Audience :: Developers",
        	"Intended Audience :: Education",
        	"Intended Audience :: Science/Research",
	),
	zip_safe=False,
	entry_points={
		'console_scripts': [
			'CPyMonC = Lacos.CPyMonC.__main__:main',
		],
	},
        install_requires=[
		"numpy==1.20.1",
		"more-itertools==8.4.0",
		"psutil==5.7.0",
		"spglib==2.0.2",		
		"cvxpy==1.2.1",
		"matplotlib==3.2.2",
		"terminaltables==3.1.0",
		"scipy==1.5.0",
		"pymatgen==2023.8.10",
		"phonopy==2.7.1",
		"scikit-learn==0.23.1",
		"click==7.1.2",
		"pandas==1.0.5",
		f"cogue @ file://localhost/{os.getcwd()}/cogue-master",
		"pytz==2020.1",
		"ruamel.yaml==0.17.0",
	],
	setup_requires=['Cython==0.29.24','numpy==1.20.1','setuptools>=18.0'],
	ext_modules=[C_monte_carlo_module, C_mpi_monte_carlo_module]
)