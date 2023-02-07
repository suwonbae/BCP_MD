from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
        Pybind11Extension("bcp_md.model.rw",
            sources = ['src/python_ext/bcp_md/model/rw.cpp',
                'src/cpp/bcp_md/model/cpp_rw.cpp'],
            include_dirs = ['src/cpp'],
            ),
        ]

setup(
        name = "bcp_md",
        version = "0.0.0",
        ext_modules = ext_modules,
        package_dir = {'': 'src/python'},
        packages = find_packages(where = 'src/python'),
        cmdclass = {"build_ext": build_ext},
        )
