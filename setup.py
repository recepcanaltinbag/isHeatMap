from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'isHeatMap',
    version = '1.0.9',
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
    packages = find_packages(exclude=['build', '_docs', 'templates']),
    install_requires = ["biopython==1.79","matplotlib","numpy==1.21.4","pandas==1.3.4","scipy==1.10.1","seaborn==0.13.0"],
    include_package_data = True,
    package_data={
        'test_data':
        ['test_data/GCA_000982435.1_ASM98243v1_genomic.fna.csv',
        'test_data/GCA_000982435.1_ASM98243v1_genomic.fna.is.fna'],
        'executables':
        ['executables/blastn','executables/makeblastdb']
    },
    entry_points = {
        'console_scripts': [
            'isheatmap = isheatmap.__main__:main',
        ],
    }
)



