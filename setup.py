from setuptools import setup, find_packages

with open("README.md", "r") as fh:
  long_description = fh.read()

setup(
    name='py-moc',
    version='0.0.1rc5',
    description="A simple model suite for the MOC",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://pymoc.github.io',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
