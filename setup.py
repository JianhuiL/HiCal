import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HiCal",
    version="0.0.1",
    author="Jianhui Li",
    author_email="jianhui.li@columbia.edu",
    description="HiC heatmap analysis tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JianhuiL/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)