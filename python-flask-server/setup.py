# coding: utf-8

import sys
from setuptools import setup, find_packages

NAME = "swagger_server"
VERSION = "1.0.0"

# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

REQUIRES = ["connexion"]

setup(
    name=NAME,
    version=VERSION,
    description="Cancer Dependency Map API",
    author_email="translator@broadinstitute.org",
    url="",
    keywords=["Swagger", "Cancer Dependency Map API"],
    install_requires=REQUIRES,
    packages=find_packages(),
    package_data={'': ['swagger/swagger.yaml']},
    include_package_data=True,
    entry_points={
        'console_scripts': ['swagger_server=swagger_server.__main__:main']},
    long_description="""\
    This site provides an API access to [Cancer Dependency Map (DepMap) data](https://depmap.org/portal/download/). The goal of the Cancer Dependency Map is to create a comprehensive preclinical reference map connecting  tumor features with tumor dependencies to accelerate the development of precision treatments.  By integrating data beyond those collected at the Broad, DepMap hopes to develop a complete  understanding of the vulnerabilities of cancer, identify targets for therapeutic development,  and design strategies to optimize patient responses to those therapies. By using this site, you agree to DepMap&#39;s [Terms and Conditions](https://depmap.org/portal/terms/).
    """
)

