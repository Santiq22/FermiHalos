import setuptools

with open('pypi-description.md', 'r') as file:
    long_description = file.read()

setuptools.setup(
    name = "fermihalos",
    version = "0.1.3.dev",
    author = "RAR collaboration",
    author_email = "scollazo@fcaglp.unlp.edu.ar",
    description = "An extended RAR model for dark matter halo astrophysics.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    packages = setuptools.find_packages(include = ["fermihalos"]),
    url = "https://github.com/Santiq22/FermiHalos",
    python_requires = ">=3.10",
    install_requires = ["numpy>=1.23.5", "scipy>=1.15.0"],
    package_data = {"": ["README.md", "pypi-description.md", "LICENSE", "CITATION.bib", "HISTORY.txt", "AUTHORS.txt"]},
    keywords = ["astrophysics", "dark matter", "halo", "RAR"],
    license = "MIT"
)