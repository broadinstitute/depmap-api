# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.models.cell_line import CellLine  # noqa: E501
from swagger_server.models.gene import Gene  # noqa: E501
from swagger_server.models.protein import Protein  # noqa: E501
from swagger_server.test import BaseTestCase


class TestTypesController(BaseTestCase):
    """TypesController integration test stubs"""

    def test_cell_lines_get(self):
        """Test case for cell_lines_get

        Retrieve list of cell lines
        """
        response = self.client.open(
            '/depmap/cell_lines',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_genes_get(self):
        """Test case for genes_get

        Retrieve list of genes
        """
        response = self.client.open(
            '/depmap/genes',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_proteins_get(self):
        """Test case for proteins_get

        Retrieve list of protein antibodies
        """
        response = self.client.open(
            '/depmap/proteins',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
