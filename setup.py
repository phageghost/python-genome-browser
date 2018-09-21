import setuptools

VER = '0.3.3'
AUTHOR = 'Dylan Skola'

print('*' * 80)
print('* {:<76} *'.format('python-genome-browser version {} by {}'.format(VER, AUTHOR)))
print('*' * 80)
print()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='pygbrowse',
                 version=VER,
                 description='Tools for making plots of genomic datasets in a genome-browser-like format ',
                 long_description=long_description,
                 url='https://github.com/phageghost/python-genome-browser',
                 author='phageghost',
                 author_email='pygbrowse@phageghost.net',
                 license='MIT',
                 packages=['pygbrowse'],
                 install_requires=['numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn', 'pysam', 'intervaltree'],
                 zip_safe=False,
                 classifiers=(
                     "Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                 ),
                 )
