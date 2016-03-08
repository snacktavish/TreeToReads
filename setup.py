from setuptools import setup


setup(
    name='TreeToReads',
    version='0.0.3',
    py_modules =['treetoreads',],
    license='No copyright',
    description = 'Tree to Reads - A python script to to read a tree, resolve polytomies, generate mutations and simulate NGS reads.',
    author = 'Emily Jane McTavish  ',
    author_email = 'ejmctavish@gmail.com',
    url = 'https://github.com/snacktavish/TreeToReads', # use the URL to the github repo
    long_description=open('README.md').read(),
    install_requires=[
          'dendropy',
      ],
)
