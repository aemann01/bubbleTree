import setuptools

with open("README.md", "r") as fh:
	long_description=fh.read()

setuptools.setup(
	name="bubbletree",
	version="1.1",
	author="Allison E. Mann",
	author_email="amann3@clemson.edu",
	description="Bubble plot or heatmap ordered by phylogenetic tree",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/aemann01/bubbleTree",
	packages=setuptools.find_packages(),
	scripts=['bubbletree/bubbletree.py'],
	classifiers=(
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License (GPL)",
		"Operating System :: OS Independent",
	),
)
