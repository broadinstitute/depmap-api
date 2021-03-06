# coding: utf-8

from __future__ import absolute_import
from datetime import date, datetime  # noqa: F401

from typing import List, Dict  # noqa: F401

from swagger_server.models.base_model_ import Model
from swagger_server import util


class Gene(Model):
    """NOTE: This class is auto generated by the swagger code generator program.

    Do not edit the class manually.
    """

    def __init__(self, hgnc_id: str=None, gene_symbol: str=None, ensembl_gene: str=None, entrez_gene_id: int=None, omim: List[str]=None):  # noqa: E501
        """Gene - a model defined in Swagger

        :param hgnc_id: The hgnc_id of this Gene.  # noqa: E501
        :type hgnc_id: str
        :param gene_symbol: The gene_symbol of this Gene.  # noqa: E501
        :type gene_symbol: str
        :param ensembl_gene: The ensembl_gene of this Gene.  # noqa: E501
        :type ensembl_gene: str
        :param entrez_gene_id: The entrez_gene_id of this Gene.  # noqa: E501
        :type entrez_gene_id: int
        :param omim: The omim of this Gene.  # noqa: E501
        :type omim: List[str]
        """
        self.swagger_types = {
            'hgnc_id': str,
            'gene_symbol': str,
            'ensembl_gene': str,
            'entrez_gene_id': int,
            'omim': List[str]
        }

        self.attribute_map = {
            'hgnc_id': 'hgnc_id',
            'gene_symbol': 'gene_symbol',
            'ensembl_gene': 'ensembl_gene',
            'entrez_gene_id': 'entrez_gene_id',
            'omim': 'omim'
        }

        self._hgnc_id = hgnc_id
        self._gene_symbol = gene_symbol
        self._ensembl_gene = ensembl_gene
        self._entrez_gene_id = entrez_gene_id
        self._omim = omim

    @classmethod
    def from_dict(cls, dikt) -> 'Gene':
        """Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The gene of this Gene.  # noqa: E501
        :rtype: Gene
        """
        return util.deserialize_model(dikt, cls)

    @property
    def hgnc_id(self) -> str:
        """Gets the hgnc_id of this Gene.


        :return: The hgnc_id of this Gene.
        :rtype: str
        """
        return self._hgnc_id

    @hgnc_id.setter
    def hgnc_id(self, hgnc_id: str):
        """Sets the hgnc_id of this Gene.


        :param hgnc_id: The hgnc_id of this Gene.
        :type hgnc_id: str
        """

        self._hgnc_id = hgnc_id

    @property
    def gene_symbol(self) -> str:
        """Gets the gene_symbol of this Gene.


        :return: The gene_symbol of this Gene.
        :rtype: str
        """
        return self._gene_symbol

    @gene_symbol.setter
    def gene_symbol(self, gene_symbol: str):
        """Sets the gene_symbol of this Gene.


        :param gene_symbol: The gene_symbol of this Gene.
        :type gene_symbol: str
        """

        self._gene_symbol = gene_symbol

    @property
    def ensembl_gene(self) -> str:
        """Gets the ensembl_gene of this Gene.


        :return: The ensembl_gene of this Gene.
        :rtype: str
        """
        return self._ensembl_gene

    @ensembl_gene.setter
    def ensembl_gene(self, ensembl_gene: str):
        """Sets the ensembl_gene of this Gene.


        :param ensembl_gene: The ensembl_gene of this Gene.
        :type ensembl_gene: str
        """

        self._ensembl_gene = ensembl_gene

    @property
    def entrez_gene_id(self) -> int:
        """Gets the entrez_gene_id of this Gene.


        :return: The entrez_gene_id of this Gene.
        :rtype: int
        """
        return self._entrez_gene_id

    @entrez_gene_id.setter
    def entrez_gene_id(self, entrez_gene_id: int):
        """Sets the entrez_gene_id of this Gene.


        :param entrez_gene_id: The entrez_gene_id of this Gene.
        :type entrez_gene_id: int
        """

        self._entrez_gene_id = entrez_gene_id

    @property
    def omim(self) -> List[str]:
        """Gets the omim of this Gene.


        :return: The omim of this Gene.
        :rtype: List[str]
        """
        return self._omim

    @omim.setter
    def omim(self, omim: List[str]):
        """Sets the omim of this Gene.


        :param omim: The omim of this Gene.
        :type omim: List[str]
        """

        self._omim = omim
