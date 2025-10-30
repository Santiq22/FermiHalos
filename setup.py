import setuptools

with open('README.md', 'r') as file:
    long_description = file.read()

setuptools.setup(
    name = "fermihalos",
    version = "0.1.0",
    author = "RAR collaboration",
    author_email = "scollazo@fcaglp.unlp.edu.ar",
    description = "An extended RAR model for dark matter halo astrophysics.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    packages = setuptools.find_packages(include = ["fermihalos"]),
    url = "https://github.com/Santiq22/FermiHalos",
    python_requires = ">=3",
    install_requires = ["numpy", "scipy"],
    package_data = {"": ["README.md", "LICENSE", "CITATION.bib"]},
    keywords = ["astrophysics", "dark matter", "halo", "RAR"],
    license = "MIT"
)