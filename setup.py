from setuptools import setup, find_packages
#from os.path import basename, splittext
#from glob import glob

def _requirement_packages(filename):
	return open(filename).read().splitlines()

setup(
	name="figp2",
	version="0.1.0",
	license="MIT License",
	description="Symbolic regression with FIGP2",
	author="Raku Shirasawa and Katsushi Takaki",
	url="https://github.com/raku68/FIGP2",
	packages=find_packages("src"), # detect python packages in the codes
	package_dir={"":"src"},
	#py_modules=[splittext(basename(path))[0] for path in glob('src/fgpnls/*.py')]
	zip_safe=False,
	install_requires=_requirement_packages("requirements.txt")
)
