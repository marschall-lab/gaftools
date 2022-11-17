import os
from setuptools import setup, Extension
from distutils.sysconfig import customize_compiler
import Cython.Build
from Cython.Build import cythonize

MIN_CYTHON_VERSION = "0.29"

def CppExtension(name, sources):
    return Extension(
        name,
        sources=sources,
        language="c++",
        extra_compile_args=["-std=c++11", "-Werror=return-type", "-Werror=narrowing"],
        include_dirs=["src/"],
        undef_macros=["NDEBUG"],
    )

extensions = []

class BuildExt(Cython.Build.build_ext):
    def build_extensions(self):
        # Remove the warning about “-Wstrict-prototypes” not being valid for C++,
        # see http://stackoverflow.com/a/36293331/715090
        customize_compiler(self.compiler)
        if self.compiler.compiler_so[0].endswith("clang"):
            # Clang needs this option in order to find the unordered_set header
            print("detected clang, using option -stdlib=libc++")
            self.compiler.compiler_so.append("-stdlib=libc++")
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        super().build_extensions()

# Avoid compilation if we are being installed within Read The Docs
if os.environ.get("READTHEDOCS") == "True":
    cmdclass = {}
    ext_modules = []
    install_requires = []
else:
    cmdclass = {"build_ext": BuildExt}
    ext_modules = extensions
    install_requires = [
        "pyfaidx>=0.5.5.2",
        "networkx",
    ]

setup(
    use_scm_version={"write_to": "gaftools/_version.py"},
    cmdclass=cmdclass,
    ext_modules=cythonize(ext_modules),
    install_requires=install_requires,
)