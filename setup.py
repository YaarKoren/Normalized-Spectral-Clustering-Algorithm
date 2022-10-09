from setuptools import setup, Extension

setup(
    name="spkmeansmodule",
    version="1.0.0",
    description="Python interface for the spkmeansmodule C library",
    author="",
    author_email="",
    ext_modules=[Extension("spkmeansmodule", ["spkmeansmodule.c","spkmeans.c"])]
    )

