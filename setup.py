from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

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
    install_requires=[
        "numpy",
        "scipy",
        "coverage"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'run_tests = run_tests:main',  # Allows running `run_tests` in the command line
        ],
    }
)
