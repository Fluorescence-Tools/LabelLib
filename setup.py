import os
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
import re
from subprocess import check_call, check_output, CalledProcessError


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if sys.platform == "win32":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DPYTHON_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if sys.platform == "win32":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j1"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO="{}"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)


def get_version():
    try:
        out = check_output(
            ["git", "show", "-s", "--format=%cd", "--date=short"]
        ).decode()
        out = out.replace("-", ".").strip()
        tag_version = check_output(["git", "describe", "--tags"]).decode()
        if out != tag_version:
            print(
                "WARNING: date-based version ({}) does not match the tag ({})!".format(
                    out, tag_version
                )
            )
        return out
    except CalledProcessError:
        return "1.0.0"


setup(
    name="labellib",
    version=get_version(),
    author="Mykola Dimura",
    author_email="mykola.dimura@gmail.com",
    maintainer="Thomas-Otavio Peulen",
    maintainer_email="thomas.otavio.peulen@gmail.com",
    description="Python bindings for LabelLib",
    long_description="Library for coarse-grained simulations of probes flexibly coupled to biomolecules.",
    url="https://github.com/Fluorescence-Tools/LabelLib",
    license="MPL v2.0",
    ext_modules=[CMakeExtension("labellib", "FlexLabel/python")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=["numpy"],
)
