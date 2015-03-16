# -*- encoding: utf-8 -*-

from distutils.core import setup

setup(
	name='MirrorTree',
	version='2.0',
	description='Use different types of correlation along evolution to identify interacting protein pairs.',
	author='Elba Raimúndez, Clàudia Fontserè, Lucas Michel',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Evolutionary Biologists',
		'Programming Language :: Python',
		'Topic :: Bioinformatics :: Phylogenetic trees'
	],
	keywords='Mirror trees',
	packages = ['MirrorTree'],
	package_dir={'MirrorTree': 'MirrorTree'},
	py_modules = ['MirrorTree/modules', 'MirrorTree/functions', 'MirrorTree/visual'],
	scripts = ['MirrorTree/mirrorTree_beta'],
)
