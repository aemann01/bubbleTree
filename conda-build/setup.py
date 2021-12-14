from setuptools import setup, find_packages

setup(
	name='bubbletree',
	version='1.0',
	packages=find_packages(include=['bubble_tree']),
	install_requires=[
	'pandas>=1.1.3',
	'biopython>=1.79',
	'matplotlib>=3.5.1',
	'seaborn>=0.11.0'])