from setuptools import find_packages, setup

setup(
    name="pyrangeyes",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,  # Ensure package data is included
    package_data={
        "pyrangeyes": ["data/*"],  # Specify the path to include data folder contents
    },
)
