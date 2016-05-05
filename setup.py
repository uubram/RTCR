from setuptools import setup, find_packages
from distutils.extension import Extension
from os import path
from glob import glob
import rtcr

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.rst")) as f:
    long_description= f.read()

setup(
        name = "rtcr",
        version = rtcr.__version__,
        description = "A pipeline for complete and accurate recovery of TCR \
repertoires from high throughput sequencing data",
        long_description = long_description,
        author = "Bram Gerritsen",
        author_email = "b.gerritsen@uu.nl",
        license = "GPLv3",
        entry_points = {
            "console_scripts": ["rtcr = rtcr.__main__:main"]},
        packages = find_packages(),
        include_package_data = True,
        ext_modules = [
            Extension("cseq",
                sources = glob("src/*.c"),
                include_dirs = ["include"],
                language = "c",
                extra_compile_args = ["-std=c99", "-pedantic", "-Wall",
                    "-Wextra", "-O3"],
                )],
        install_requires = ["vtrie>=0.0.1", "editdistance>=0.3.1"])
