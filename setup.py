from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = [
        line.strip() for line in f if line.strip() and not line.startswith("#")
    ]

setup(
    name="psrmatch",
    version="1.0.0",
    packages=find_packages(exclude=["legacy", "legacy.*"]),
    install_requires=requirements,
)
