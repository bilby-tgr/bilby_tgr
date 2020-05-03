#!/usr/bin/env python

from setuptools import setup
import subprocess
import sys
import os

python_version = sys.version_info
if python_version < (3, 5):
    sys.exit("Python < 3.5 is not supported, aborting setup")
print("Confirmed Python version {}.{}.{} >= 3.5.0".format(*python_version[:3]))


def write_version_file(version):
    """ Writes a file with version information to be used at run time

    Parameters
    ----------
    version: str
        A string containing the current version information

    Returns
    -------
    version_file: str
        A path to the version file

    """
    try:
        git_log = subprocess.check_output(
            ['git', 'log', '-1', '--pretty=%h %ai']).decode('utf-8')
        git_diff = (subprocess.check_output(['git', 'diff', '.']) +
                    subprocess.check_output(
                        ['git', 'diff', '--cached', '.'])).decode('utf-8')
        if git_diff == '':
            git_status = '(CLEAN) ' + git_log
        else:
            git_status = '(UNCLEAN) ' + git_log
    except Exception as e:
        print("Unable to obtain git version information, exception: {}"
              .format(e))
        git_status = ''

    version_file = '.version'
    if os.path.isfile(version_file) is False:
        with open('src/' + version_file, 'w+') as f:
            f.write('{}: {}'.format(version, git_status))

    return version_file


def get_long_description():
    """ Finds the README and reads in the description """
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, 'README.rst')) as f:
        long_description = f.read()
    return long_description


# get version info from __init__.py
def readfile(filename):
    with open(filename) as fp:
        filecontents = fp.read()
    return filecontents


VERSION = '0.0.1'
version_file = write_version_file(VERSION)
long_description = get_long_description()

setup(name='bilby_tgr',
      description='Testing GR source models to use with Bilby',
      long_description=long_description,
      long_description_content_type="text/x-rst",
      url='https://git.ligo.org/gregory.ashton/bilby_tgr',
      author='Greg Ashton, Colm Talbot',
      author_email='gregory.ashton@monash.edu',
      license="MIT",
      version=VERSION,
      packages=['bilby_tgr'],
      package_dir={'bilby_tgr': 'src'},
      package_data={'bilby_tgr': [version_file]},
      python_requires='>=3.5',
      install_requires=[
          'bilby',
          'lalsuite',
          'numpy>=1.9',
      ],
      classifiers=[
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent"])
