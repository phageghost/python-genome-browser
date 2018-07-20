import setuptools

VER = '0.0.1'
AUTHOR = 'Dylan Skola'

print('*' * 80)
print('* {:<76} *'.format('pygbrowse version {} by {}'.format(VER, AUTHOR)))
print('*' * 80)
print()

setuptools.setup(name='pygbrowse',
                 version=VER,
                 description='Plotting tools for rendering genome-browser-like views of genomic data',
                 url='http://github.com/phageghost/python-genome-browser',
                 author=AUTHOR,
                 author_email='pygbrowse@phageghost.net',
                 license='MIT',
                 packages=['pygbrowse'],
                 install_requires=['numpy', 'pandas', 'seaborn', 'intervaltree', 'matplotlib'],
                 zip_safe=False)
