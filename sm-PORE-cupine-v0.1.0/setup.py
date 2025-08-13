from setuptools import setup, find_packages
import glob
import os

with open("requirements.txt") as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from smPOREcupine import __version__, _program

setup(name="smPOREcupine",
      version=__version__,
      packages=["smPOREcupine"],
      description="Pipeline to perform signal alignment of nanopore signal traces to a signal reference",
      entry_points="""
      [console_scripts]
      {program} = smPOREcupine.execute_workflows:main
      """.format(program = _program),
      install_requires=required,
      include_package_data=True,
      package_data={
          "smPOREcupine": ['Snakefile','config/cluster_pluto.json', 'config/config_env_pluto.yaml']
      },
      keywords=[])
