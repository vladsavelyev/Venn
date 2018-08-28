#!/usr/bin/env python
from os.path import join, isfile, abspath, dirname


name = script_name = package_name = 'venn_bed'


from setuptools import setup, find_packages

try:
    from ngs_utils import setup_utils
except:
    version = open('VERSION.txt').read().strip().split('\n')[0]
    setup(version=version)  # For conda-build jinja context to read the version
else:
    version = setup_utils.init(package_name, package_name, __file__)
    from ngs_utils.setup_utils import write_version_py, find_package_files, get_reqs
    from ngs_utils import setup_utils
    setup(
        name=name,
        version=version,
        author='Vlad Saveliev and Alla Mikheenko',
        author_email='vladislav.sav@gmail.com',
        description='Genome capture target coverage evaluation tool',
        long_description=(open('README.md').read()),
        keywords='bioinformatics',
        url='https://github.com/vladsaveliev/TargQC',
        download_url='https://github.com/vladsaveliev/TargQC/releases',
        license='GPLv3',
        packages=find_packages(),
        package_data={
            package_name: [
                '*.html'
            ]
        },
        include_package_data=True,
        zip_safe=False,
        scripts=[
            join('scripts', script_name),
        ],
        install_requires=get_reqs(),
        setup_requires=['numpy'],
        classifiers=[
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Programming Language :: JavaScript',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        test_suite='nose.collector',
        tests_require=['nose'],
    )