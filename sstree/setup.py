from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
Cython.Compiler.Options.fast_fail = True
import numpy as np

libSSTreeIncludeDir = "/usr/include/SSTree"

setup(
    ext_modules=[
        Extension(
            name="sstree",
            sources=["sstree.pyx"],
            language="c++",
            include_dirs = ['src', libSSTreeIncludeDir, np.get_include()],
            library_dirs = ["/usr/lib"],
            libraries = ["SSTree"],
        )
    ],
    cmdclass={'build_ext': build_ext}
)
