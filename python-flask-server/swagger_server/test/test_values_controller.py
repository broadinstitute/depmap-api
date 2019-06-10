# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.models.copy_number_value import CopyNumberValue  # noqa: E501
from swagger_server.models.dependency_value import DependencyValue  # noqa: E501
from swagger_server.models.expression_value import ExpressionValue  # noqa: E501
from swagger_server.models.mutation import Mutation  # noqa: E501
from swagger_server.models.protein_array_value import ProteinArrayValue  # noqa: E501
from swagger_server.models.rnai_value import RnaiValue  # noqa: E501
from swagger_server.test import BaseTestCase


class TestValuesController(BaseTestCase):
    """ValuesController integration test stubs"""

    def test_copy_number_by_cell_line_depmap_id_get(self):
        """Test case for copy_number_by_cell_line_depmap_id_get

        Retrieve copy-number values by cell line
        """
        response = self.client.open(
            '/depmap/copy_number/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_copy_number_by_gene_entrez_gene_id_get(self):
        """Test case for copy_number_by_gene_entrez_gene_id_get

        Retrieve copy-number values by gene
        """
        response = self.client.open(
            '/depmap/copy_number/by_gene/{entrez_gene_id}'.format(entrez_gene_id=56),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_gene_dependency_by_cell_line_depmap_id_get(self):
        """Test case for gene_dependency_by_cell_line_depmap_id_get

        Retrieve gene dependency by cell line
        """
        response = self.client.open(
            '/depmap/gene_dependency/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_gene_dependency_by_gene_entrez_gene_id_get(self):
        """Test case for gene_dependency_by_gene_entrez_gene_id_get

        Retrieve gene dependency by gene
        """
        response = self.client.open(
            '/depmap/gene_dependency/by_gene/{entrez_gene_id}'.format(entrez_gene_id=56),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_gene_expression_by_cell_line_depmap_id_get(self):
        """Test case for gene_expression_by_cell_line_depmap_id_get

        Retrieve gene expression by cell line
        """
        response = self.client.open(
            '/depmap/gene_expression/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_gene_expression_by_gene_ensembl_gene_get(self):
        """Test case for gene_expression_by_gene_ensembl_gene_get

        Retrieve gene expression by gene
        """
        response = self.client.open(
            '/depmap/gene_expression/by_gene/{ensembl_gene}'.format(ensembl_gene='ensembl_gene_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_mutations_by_cell_line_depmap_id_get(self):
        """Test case for mutations_by_cell_line_depmap_id_get

        Retrieve mutations by cell line
        """
        response = self.client.open(
            '/depmap/mutations/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_mutations_by_gene_entrez_gene_id_get(self):
        """Test case for mutations_by_gene_entrez_gene_id_get

        Retrieve mutations by gene
        """
        response = self.client.open(
            '/depmap/mutations/by_gene/{entrez_gene_id}'.format(entrez_gene_id=56),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_protein_array_by_cell_line_depmap_id_get(self):
        """Test case for protein_array_by_cell_line_depmap_id_get

        Retrieve protein array values by cell line
        """
        response = self.client.open(
            '/depmap/protein_array/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_protein_array_by_protein_antibody_name_get(self):
        """Test case for protein_array_by_protein_antibody_name_get

        Retrieve protein array values by gene
        """
        response = self.client.open(
            '/depmap/protein_array/by_protein/{antibody_name}'.format(antibody_name='antibody_name_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_rnai_by_cell_line_depmap_id_get(self):
        """Test case for rnai_by_cell_line_depmap_id_get

        Retrieve rnai values by cell line
        """
        response = self.client.open(
            '/depmap/rnai/by_cell_line/{depmap_id}'.format(depmap_id='depmap_id_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_rnai_by_gene_entrez_gene_id_get(self):
        """Test case for rnai_by_gene_entrez_gene_id_get

        Retrieve rnai values by gene
        """
        response = self.client.open(
            '/depmap/rnai/by_gene/{entrez_gene_id}'.format(entrez_gene_id=56),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
