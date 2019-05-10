import connexion
import six

from swagger_server.models.copy_number_value import CopyNumberValue  # noqa: E501
from swagger_server.models.dependency_value import DependencyValue  # noqa: E501
from swagger_server.models.expression_value import ExpressionValue  # noqa: E501
from swagger_server.models.mutation import Mutation  # noqa: E501
from swagger_server.models.protein_array_value import ProteinArrayValue  # noqa: E501
from swagger_server.models.rnai_value import RnaiValue  # noqa: E501
from swagger_server import util


def copy_number_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve copy-number values by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[CopyNumberValue]
    """
    return 'do some magic!'


def copy_number_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve copy-number values by gene

     # noqa: E501

    :param entrez_gene_id: 
    :type entrez_gene_id: int

    :rtype: List[CopyNumberValue]
    """
    return 'do some magic!'


def gene_dependency_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve gene dependency by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[DependencyValue]
    """
    return 'do some magic!'


def gene_dependency_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve gene dependency by gene

     # noqa: E501

    :param entrez_gene_id: 
    :type entrez_gene_id: int

    :rtype: List[DependencyValue]
    """
    return 'do some magic!'


def gene_expression_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve gene expression by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[ExpressionValue]
    """
    return 'do some magic!'


def gene_expression_by_gene_ensembl_gene_get(ensembl_gene):  # noqa: E501
    """Retrieve gene expression by gene

     # noqa: E501

    :param ensembl_gene: 
    :type ensembl_gene: str

    :rtype: List[ExpressionValue]
    """
    return 'do some magic!'


def mutations_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve mutations by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[Mutation]
    """
    return 'do some magic!'


def mutations_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve mutations by gene

     # noqa: E501

    :param entrez_gene_id: 
    :type entrez_gene_id: int

    :rtype: List[Mutation]
    """
    return 'do some magic!'


def protein_array_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve protein array values by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[ProteinArrayValue]
    """
    return 'do some magic!'


def protein_array_by_protein_antibody_name_get(antibody_name):  # noqa: E501
    """Retrieve protein array values by gene

     # noqa: E501

    :param antibody_name: 
    :type antibody_name: str

    :rtype: List[ProteinArrayValue]
    """
    return 'do some magic!'


def rnai_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve rnai values by cell line

     # noqa: E501

    :param depmap_id: 
    :type depmap_id: str

    :rtype: List[RnaiValue]
    """
    return 'do some magic!'


def rnai_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve rnai values by gene

     # noqa: E501

    :param entrez_gene_id: 
    :type entrez_gene_id: int

    :rtype: List[RnaiValue]
    """
    return 'do some magic!'
