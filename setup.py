from setuptools import setup, find_packages

setup(name='COGfinder',
      version='0.1.0',
      description='Package to find the Cluster of Orthologous Groups (COG) categories for every gene in a fasta files(s)',
      url='https://github.com/jclachance/COGfinder/',
      author='Jean-Christophe Lachance',
      author_email='jelachance@eng.ucsd.edu',
      license='MIT',
      packages=find_packages(),
      package_data={'COGfinder':['data/*.*']},
      include_package_data=True,
      install_requires=['cobra>=0.11.0','numpy>=1.13',
			'BioPython','seaborn','matplotlib'])
