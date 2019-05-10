import connexion
import six

from swagger_server.models.cell_line import CellLine  # noqa: E501
from swagger_server.models.gene import Gene  # noqa: E501
from swagger_server.models.protein import Protein  # noqa: E501
from swagger_server import util


def cell_lines_get():  # noqa: E501
    """Retrieve list of cell lines

     # noqa: E501


    :rtype: List[CellLine]
    """
    return 'do some magic!'


def genes_get():  # noqa: E501
    """Retrieve list of genes

     # noqa: E501


    :rtype: List[Gene]
    """
    return 'do some magic!'


def proteins_get():  # noqa: E501
    """Retrieve list of protein antibodies

     # noqa: E501


    :rtype: List[Protein]
    """
    return 'do some magic!'
