from setuptools import find_packages, setup

setup(
	name = "ASR",
	version = '0.0',
	description = 'Python package containing pipeline for ancestral sequence reconstruction',
	url = 'https://github.com/jgoldford/ASR/',
	author = 'Joshua E. Goldford',
	author_email = 'goldford.joshua@gmail.com',
	packages = find_packages(),
	install_requires = [
		'scipy',
		'numpy',
		'pandas',
		'biopython'],
	include_package_data = True,
)