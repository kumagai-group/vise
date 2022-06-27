import os

from setuptools import setup, find_packages
from distutils.extension import Extension

from vise import __version__

module_dir = os.path.dirname(os.path.abspath(__file__))
reqs_raw = open(os.path.join(module_dir, "requirements.txt")).read()
reqs_list = [r.replace("==", "~=") for r in reqs_raw.split("\n")]
#reqs_list = [r for r in reqs_raw.split("\n")]

# try:
#     from Cython.Distutils import build_ext
# except ImportError:
#     use_cython = False
# else:
#     use_cython = True

cmdclass = {}
ext_modules = []

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='vise',
    version=__version__,
    author="Yu Kumagai",
    author_email="yukumagai@tohoku.ac.jp",
    url='https://github.com/kumagai-group/vise',
    packages=find_packages(),
    license='MIT license',
    description="Package for handling io of vasp package in kumagai group at "
                "IMR, Tohoku university",
    long_description=long_description,
    long_description_content_type="text/markdown",

    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires='>=3.7',
    install_requires=reqs_list,
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'vise = vise.cli.main:main',
            'vise_util = vise.cli.main_util:main',
        ]
    }
)
