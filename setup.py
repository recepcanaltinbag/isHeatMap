from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'isHeatMap',
    version = '1.0.0',
    description = 'heatmap of insertion sequences in a given genome',
    long_description = long_description,
    url = 'https://github.com/recepcanaltinbag/isHeatMap',
    author = 'Recep Can Altınbağ',
    author_email = 'recepcanaltinbag@gmail.com',

    license = 'MIT',

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers'
    ],

    keywords = 'bioinformatics insertion sequences',
    packages = []
    include_package_data = True,
    package_data={
        'data':
        ['test_data/GCA_000982435.1_ASM98243v1_genomic.fna.csv',
        'test_data/GCA_000982435.1_ASM98243v1_genomic.fna.is.fna'],
        'executables':
        ['executables/blastn','executables/makeblastdb']
    },
    entry_points = {
        'console_scripts': [
            'isheatmap = isheatmap.isheatmap:main',
        ],
    }
)



