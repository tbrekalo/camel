import sys
from skbuild import setup

setup(
    name="camelpy",
    version="0.2.0",
    author="Tvrtko Brekalo",
    author_email="brekalo.tvrtko@gmail.com",
    description="Sequencing coverage static library",
    url="https://github.com/tbrekalo/camel",
    license="MIT",
    packages=['camelpy'],
    package_dir={'': ''},
    cmake_install_dir="camelpy",
    include_package_data=True,
    python_requires=">=3.8")
