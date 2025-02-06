from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np
from setuptools.command.build_ext import build_ext as _build_ext
import os
import sys
import subprocess

############################################# MATHEMATICA #####################################################
# Function to check if Mathematica is installed
def find_mathematica():
    """
    Detect if Mathematica is installed and return paths for:
    - WSTP include directory (header files)
    - WSTP library directory (shared libraries)
    
    Returns:
        (bool, str, str): (True, WSTP_INCLUDE, WSTP_LIB) or (False, None, None) if not found.
    """
    print("🔍 Checking for Mathematica installation...", file=sys.stderr)  # Force visible output in pip logs

    # 1️⃣ Check Environment Variables
    env_paths = [
        os.getenv("MATHEMATICA_PATH"),
        os.getenv("WOLFRAM_PATH"),
        os.getenv("WSTP_PATH"),  # Some users might define this manually
    ]
    for env_path in env_paths:
        if env_path and os.path.exists(env_path):
            print(f"✅ Found Mathematica via environment variable: {env_path}", file=sys.stderr)
            return get_wstp_paths(env_path)

    # 2️⃣ Check Common Installation Paths
    possible_paths = {
        "Linux": "/usr/local/Wolfram/Mathematica",
        "MacOS": "/Applications/Mathematica.app",
        "Windows": os.path.join(os.getenv("PROGRAMFILES", ""), "Wolfram Research", "Mathematica"),
    }

    for os_name, path in possible_paths.items():
        if os.path.exists(path):
            print(f"✅ Found Mathematica at {path}", file=sys.stderr)
            return get_wstp_paths(path)

    # 3️⃣ Check if `math` or `wolfram` is in system PATH (Linux/macOS)
    try:
        result = subprocess.run(["which", "math"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stdout.strip():
            print(f"✅ Found `math` at {result.stdout.strip()}", file=sys.stderr)
            return get_wstp_paths(os.path.dirname(result.stdout.strip()))
    except FileNotFoundError:
        pass

    try:
        result = subprocess.run(["which", "wolfram"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stdout.strip():
            print(f"✅ Found `wolfram` at {result.stdout.strip()}", file=sys.stderr)
            return get_wstp_paths(os.path.dirname(result.stdout.strip()))
    except FileNotFoundError:
        pass

    # 4️⃣ On Windows, Check the Registry
    if os.name == "nt":
        try:
            import winreg
            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\Wolfram Research\Mathematica")
            mathematica_path, _ = winreg.QueryValueEx(key, "InstallPath")
            if os.path.exists(mathematica_path):
                print(f"✅ Found Mathematica in Windows registry: {mathematica_path}", file=sys.stderr)
                return get_wstp_paths(mathematica_path)
        except (FileNotFoundError, ImportError):
            pass

    # ❌ If No Installation Found
    print("❌ Mathematica not found. Skipping WSTP setup.", file=sys.stderr)
    return False, None, None

def get_wstp_paths(mathematica_base_path):
    """
    Given Mathematica's base installation path, find WSTP include and library directories.
    """
    if os.name == "nt":  # Windows
        wstp_path = os.path.join(mathematica_base_path, "SystemFiles", "Links", "WSTP", "DeveloperKit", "Windows-x86-64", "CompilerAdditions")
    elif sys.platform == "darwin":  # macOS
        wstp_path = os.path.join(mathematica_base_path, "Contents", "SystemFiles", "Links", "WSTP", "DeveloperKit", "MacOSX-x86-64", "CompilerAdditions")
    else:  # Linux
        wstp_path = os.path.join(mathematica_base_path, "13.1", "SystemFiles", "Links", "WSTP", "DeveloperKit", "Linux-x86-64", "CompilerAdditions")

    if os.path.exists(wstp_path):
        print(f"✅ WSTP paths detected: Include & Lib -> {wstp_path}", file=sys.stderr)
        return True, wstp_path, wstp_path

    print(f"⚠️ Warning: Found Mathematica but could not locate WSTP paths.", file=sys.stderr)
    return False, None, None


# If Mathematica is found, modify init.m
is_mathematica, WSTP_INCLUDE, WSTP_LIB = find_mathematica()

#if (WSTP_INCLUDE is None) or (WSTP_LIB is None):
#    is_mathematica = False

if is_mathematica:

    try:
        print("Mathematica detected! Modifying init.m...", file=sys.stderr)

        # Find Mathematica's init.m path
        mathematica_userbase = os.path.expanduser("~/.Mathematica")  # Linux/mac
        if os.name == "nt":  # Windows
            mathematica_userbase = os.path.join(os.getenv("APPDATA"), "Mathematica")

        init_file = os.path.join(mathematica_userbase, "Kernel", "init.m")

        # Path to the Mathematica package inside your Python package
        mathematica_package_path = os.path.join(os.path.dirname(__file__), "pair_multiplication", "mathematica")

        # Ensure the Kernel directory exists
        os.makedirs(os.path.dirname(init_file), exist_ok=True)

        # Add package path to init.m if not already there
        if os.path.exists(init_file):
            with open(init_file, "r") as f:
                lines = f.readlines()
        else:
            lines = []

        new_line = f'AppendTo[$Path, "{mathematica_package_path}"];\n'

        if new_line not in lines:
            with open(init_file, "a") as f:
                f.write(new_line)
    except Exception as e:
        print('Could not add package to Mathematica path. Skipping Mathematica setup.', file=sys.stderr)
        #is_mathematica = False

else:
    print("Skipping Mathematica setup.", file=sys.stderr)
##############################################################################################################################




class build_ext_inplace(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        self.inplace = True  # Set the `inplace` flag

with open("README.md", "r") as fh:
    long_description = fh.read()
    
    
# Only include WSTP paths if Mathematica is found
wstp_ext = None
#if WSTP_INCLUDE and WSTP_LIB:# and 
if is_mathematica:
    wstp_ext = Extension(
        "pair_multiplication.mathematica_wrapper",
        sources=["pair_multiplication/mathematica/mathematica_wrapper.c"],
        include_dirs=[WSTP_INCLUDE],
        library_dirs=[WSTP_LIB],
        libraries=["WSTP"],  # Link to WSTP
    )


cython_extensions = [
    Extension(
        "pair_multiplication.cutils_pair_multiplication",  # Adjust the module path to your package structure
        ["pair_multiplication/cutils_pair_multiplication.pyx"],
        include_dirs=[np.get_include()],
    ),
]

package_data = ["*.so"]

if is_mathematica:
    cython_extensions.append(wstp_ext)
    package_data.append("mathematica/pair_multiplication.wl")

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
    "pair_multiplication": package_data,
    },
    python_requires='>=3.9',
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'run_tests = run_tests:main',  # Allows running `run_tests` in the command line
        ],
    }
)
