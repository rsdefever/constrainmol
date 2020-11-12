from setuptools import setup, find_packages

requirements = [
    "numpy",
    "parmed",
    "pyomo",
]

__version__ = "0.0.0"

setup(
    name="constrainmol",
    version=__version__,
    packages=find_packages(),
    license="MIT",
    author="Ryan S. DeFever",
    author_email="rdefever@nd.edu",
    url="https://github.com/rsdefever/constrainmol",
    install_requires=requirements,
    python_requires=">=3.6, <4",
    include_package_data=True,
    description="Bond constraints",
)
