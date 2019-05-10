import connexion
import six

from swagger_server.models.cell_line import CellLine  # noqa: E501
from swagger_server.models.gene import Gene  # noqa: E501
from swagger_server.models.protein import Protein  # noqa: E501
from swagger_server import util

from swagger_server.db.queries import select_all_genes
from swagger_server.db.queries import select_all_cell_lines
from swagger_server.db.queries import select_all_proteins

def cell_lines_get():  # noqa: E501
    """Retrieve list of cell lines

     # noqa: E501


    :rtype: List[CellLine]
    """
    return select_all_cell_lines()



def genes_get():  # noqa: E501
    """Retrieve list of genes

     # noqa: E501


    :rtype: List[Gene]
    """
    return select_all_genes()


def proteins_get():  # noqa: E501
    """Retrieve list of protein antibodies

     # noqa: E501


    :rtype: List[Protein]
    """
    return select_all_proteins()
