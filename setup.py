from setuptools import setup, find_packages
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='vise-yuuukuma',
    version="0.1.dev0",
    author="Yu Kumagai",
    author_email="yuuukuma@gmail.com",
    url='https://github.com/oba-group/vise',
    packages=find_packages(),
    license='MIT license',
    description="Package for construcing the computational materials database " 
                "in Oba group at Tokyo Institute of Technology",
    long_description=long_description,
    long_description_content_type="text/markdown",

    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires='>=3.6',
    # install_requires=['numpy','pymatgen','matplotlib','seekpath','spglib',
    #                   'custodian', 'scipy'],
    cmdclass = cmdclass,
    ext_modules=ext_modules,
    include_package_data=True, install_requires=['PyYAML']
)
