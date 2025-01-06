from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext_inplace(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        self.inplace = True  # Set the `inplace` flag

with open("README.md", "r") as fh:
    long_description = fh.read()

cython_extensions = [
    Extension(
        "pair_multiplication.cutils_pair_multiplication",  # Adjust the module path to your package structure
        ["pair_multiplication/cutils_pair_multiplication.pyx"],
        include_dirs=[np.get_include()],
    ),
    #Extension(
    #    "package_name.tableau",  # Add other modules similarly
    #    ["package_name/tableau.pyx"],
    #    include_dirs=[np.get_include()],
    #),
    # Add other Cython files here
]

setup(
    name="pair_multiplication",
    version="0.1",
    author="Bernie Telalovic",
    author_email="bernie.telalovic@gmail.com",
    description="A package for manipulating Young diagrams and related mathematical structures.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BernandaTelalovic/pair_multiplication",  # Replace with your package's URL
    packages=find_packages(),
    ext_modules=cythonize(cython_extensions, language_level=3),  # Ensure extensions are built
    include_dirs=[np.get_include()],  # Include NumPy headers globally
    install_requires=[
        "numpy",
        "scipy",
        "coverage",
    ],
    cmdclass={"build_ext": build_ext_inplace},
    include_package_data=True,  # Include additional files
    setup_requires=["cython", "numpy", "coverage"],  # Add build dependencies
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={
    "pair_multiplication": ["*.so"],
    },
    python_requires='>=3.9',
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'run_tests = run_tests:main',  # Allows running `run_tests` in the command line
        ],
    }
)
