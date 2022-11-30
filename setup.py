from setuptools import setup, find_packages

with open("./assets/requirements.txt") as requirements_file:
    requirements = requirements_file.readlines()

setup(
    name='assets',
    version='0.0.1',
    packages=find_packages(include=['assets']),
    install_requires=requirements,
)
