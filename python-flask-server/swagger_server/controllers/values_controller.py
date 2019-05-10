import connexion
import six

from swagger_server.models.copy_number_value import CopyNumberValue  # noqa: E501
from swagger_server.models.dependency_value import DependencyValue  # noqa: E501
from swagger_server.models.expression_value import ExpressionValue  # noqa: E501
from swagger_server.models.mutation import Mutation  # noqa: E501
from swagger_server.models.protein_array_value import ProteinArrayValue  # noqa: E501
from swagger_server.models.rnai_value import RnaiValue  # noqa: E501
from swagger_server import util

from swagger_server.db.queries import select_copy_number_by_cell_line
from swagger_server.db.queries import select_copy_number_by_gene
from swagger_server.db.queries import select_gene_dependency_by_cell_line
from swagger_server.db.queries import select_gene_dependency_by_gene
from swagger_server.db.queries import select_gene_expression_by_cell_line
from swagger_server.db.queries import select_gene_expression_by_gene
from swagger_server.db.queries import select_mutations_brief_by_cell_line
from swagger_server.db.queries import select_mutations_brief_by_gene
from swagger_server.db.queries import select_mutations_by_cell_line
from swagger_server.db.queries import select_mutations_by_gene
from swagger_server.db.queries import select_protein_array_by_antibody_name
from swagger_server.db.queries import select_protein_array_by_cell_line
from swagger_server.db.queries import select_rnai_by_cell_line
from swagger_server.db.queries import select_rnai_by_gene


def copy_number_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve gene expression by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[CopyNumberValue]
    """
    return select_copy_number_by_cell_line(depmap_id)


def copy_number_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve gene expression by gene

     # noqa: E501

    :param entrez_gene_id:
    :type entrez_gene_id: int

    :rtype: List[CopyNumberValue]
    """
    return select_copy_number_by_gene(entrez_gene_id)


def gene_dependency_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve gene expression by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[DependencyValue]
    """
    return select_gene_dependency_by_cell_line(depmap_id)


def gene_dependency_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve gene expression by gene

     # noqa: E501

    :param entrez_gene_id:
    :type entrez_gene_id: int

    :rtype: List[DependencyValue]
    """
    return select_gene_dependency_by_gene(entrez_gene_id)


def gene_expression_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve gene expression by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[ExpressionValue]
    """
    return select_gene_expression_by_cell_line(depmap_id)


def gene_expression_by_gene_ensembl_gene_get(ensembl_gene):  # noqa: E501
    """Retrieve gene expression by gene

     # noqa: E501

    :param ensembl_gene:
    :type ensembl_gene: str

    :rtype: List[ExpressionValue]
    """
    return select_gene_expression_by_gene(ensembl_gene)


def mutations_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve mutations by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[Mutation]
    """
    return select_mutations_by_cell_line(depmap_id)


def mutations_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve mutations by gene

     # noqa: E501

    :param entrez_gene_id:
    :type entrez_gene_id: int

    :rtype: List[Mutation]
    """
    return select_mutations_by_gene(entrez_gene_id)


def protein_array_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve protein array values by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[ProteinArrayValue]
    """
    return select_protein_array_by_cell_line(depmap_id)


def protein_array_by_protein_antibody_name_get(antibody_name):  # noqa: E501
    """Retrieve protein array values by gene

     # noqa: E501

    :param antibody_name:
    :type antibody_name: str

    :rtype: List[ProteinArrayValue]
    """
    return select_protein_array_by_antibody_name(antibody_name)


def rnai_by_cell_line_depmap_id_get(depmap_id):  # noqa: E501
    """Retrieve rnai values by cell line

     # noqa: E501

    :param depmap_id:
    :type depmap_id: str

    :rtype: List[RnaiValue]
    """
    return select_rnai_by_cell_line(depmap_id)


def rnai_by_gene_entrez_gene_id_get(entrez_gene_id):  # noqa: E501
    """Retrieve rnai values by gene

     # noqa: E501

    :param entrez_gene_id:
    :type entrez_gene_id: int

    :rtype: List[RnaiValue]
    """
    return select_rnai_by_gene(entrez_gene_id)
