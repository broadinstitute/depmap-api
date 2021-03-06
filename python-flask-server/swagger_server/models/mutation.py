# coding: utf-8

from __future__ import absolute_import
from datetime import date, datetime  # noqa: F401

from typing import List, Dict  # noqa: F401

from swagger_server.models.base_model_ import Model
from swagger_server import util


class Mutation(Model):
    """NOTE: This class is auto generated by the swagger code generator program.

    Do not edit the class manually.
    """

    def __init__(self, entrez_gene_id: int=None, gene_symbol: str=None, depmap_id: str=None, ncbi_build: int=None, chromosome: str=None, start_position: int=None, end_position: int=None, strand: str=None, variant_classification: str=None, variant_type: str=None, reference_allele: str=None, tumor_seq_allele: str=None, dbsnp_rs: List[str]=None, dbsnp_val_status: str=None, genome_change: str=None, annotation_transcript: str=None, cdna_change: str=None, codon_change: str=None, protein_change: str=None, is_deleterious: str=None, tcga_hotspot_count: str=None, cosmic_hotspot_count: str=None, exac_af: str=None, cga_wes_ac: str=None, sanger_wes_ac: str=None, sanger_recalibwes_ac: str=None, rnaseq_ac: str=None, hc_ac: str=None, rd_ac: str=None, wgs_ac: str=None, variant_annotation: str=None):  # noqa: E501
        """Mutation - a model defined in Swagger

        :param entrez_gene_id: The entrez_gene_id of this Mutation.  # noqa: E501
        :type entrez_gene_id: int
        :param gene_symbol: The gene_symbol of this Mutation.  # noqa: E501
        :type gene_symbol: str
        :param depmap_id: The depmap_id of this Mutation.  # noqa: E501
        :type depmap_id: str
        :param ncbi_build: The ncbi_build of this Mutation.  # noqa: E501
        :type ncbi_build: int
        :param chromosome: The chromosome of this Mutation.  # noqa: E501
        :type chromosome: str
        :param start_position: The start_position of this Mutation.  # noqa: E501
        :type start_position: int
        :param end_position: The end_position of this Mutation.  # noqa: E501
        :type end_position: int
        :param strand: The strand of this Mutation.  # noqa: E501
        :type strand: str
        :param variant_classification: The variant_classification of this Mutation.  # noqa: E501
        :type variant_classification: str
        :param variant_type: The variant_type of this Mutation.  # noqa: E501
        :type variant_type: str
        :param reference_allele: The reference_allele of this Mutation.  # noqa: E501
        :type reference_allele: str
        :param tumor_seq_allele: The tumor_seq_allele of this Mutation.  # noqa: E501
        :type tumor_seq_allele: str
        :param dbsnp_rs: The dbsnp_rs of this Mutation.  # noqa: E501
        :type dbsnp_rs: List[str]
        :param dbsnp_val_status: The dbsnp_val_status of this Mutation.  # noqa: E501
        :type dbsnp_val_status: str
        :param genome_change: The genome_change of this Mutation.  # noqa: E501
        :type genome_change: str
        :param annotation_transcript: The annotation_transcript of this Mutation.  # noqa: E501
        :type annotation_transcript: str
        :param cdna_change: The cdna_change of this Mutation.  # noqa: E501
        :type cdna_change: str
        :param codon_change: The codon_change of this Mutation.  # noqa: E501
        :type codon_change: str
        :param protein_change: The protein_change of this Mutation.  # noqa: E501
        :type protein_change: str
        :param is_deleterious: The is_deleterious of this Mutation.  # noqa: E501
        :type is_deleterious: str
        :param tcga_hotspot_count: The tcga_hotspot_count of this Mutation.  # noqa: E501
        :type tcga_hotspot_count: str
        :param cosmic_hotspot_count: The cosmic_hotspot_count of this Mutation.  # noqa: E501
        :type cosmic_hotspot_count: str
        :param exac_af: The exac_af of this Mutation.  # noqa: E501
        :type exac_af: str
        :param cga_wes_ac: The cga_wes_ac of this Mutation.  # noqa: E501
        :type cga_wes_ac: str
        :param sanger_wes_ac: The sanger_wes_ac of this Mutation.  # noqa: E501
        :type sanger_wes_ac: str
        :param sanger_recalibwes_ac: The sanger_recalibwes_ac of this Mutation.  # noqa: E501
        :type sanger_recalibwes_ac: str
        :param rnaseq_ac: The rnaseq_ac of this Mutation.  # noqa: E501
        :type rnaseq_ac: str
        :param hc_ac: The hc_ac of this Mutation.  # noqa: E501
        :type hc_ac: str
        :param rd_ac: The rd_ac of this Mutation.  # noqa: E501
        :type rd_ac: str
        :param wgs_ac: The wgs_ac of this Mutation.  # noqa: E501
        :type wgs_ac: str
        :param variant_annotation: The variant_annotation of this Mutation.  # noqa: E501
        :type variant_annotation: str
        """
        self.swagger_types = {
            'entrez_gene_id': int,
            'gene_symbol': str,
            'depmap_id': str,
            'ncbi_build': int,
            'chromosome': str,
            'start_position': int,
            'end_position': int,
            'strand': str,
            'variant_classification': str,
            'variant_type': str,
            'reference_allele': str,
            'tumor_seq_allele': str,
            'dbsnp_rs': List[str],
            'dbsnp_val_status': str,
            'genome_change': str,
            'annotation_transcript': str,
            'cdna_change': str,
            'codon_change': str,
            'protein_change': str,
            'is_deleterious': str,
            'tcga_hotspot_count': str,
            'cosmic_hotspot_count': str,
            'exac_af': str,
            'cga_wes_ac': str,
            'sanger_wes_ac': str,
            'sanger_recalibwes_ac': str,
            'rnaseq_ac': str,
            'hc_ac': str,
            'rd_ac': str,
            'wgs_ac': str,
            'variant_annotation': str
        }

        self.attribute_map = {
            'entrez_gene_id': 'entrez_gene_id',
            'gene_symbol': 'gene_symbol',
            'depmap_id': 'depmap_id',
            'ncbi_build': 'ncbi_build',
            'chromosome': 'chromosome',
            'start_position': 'start_position',
            'end_position': 'end_position',
            'strand': 'strand',
            'variant_classification': 'variant_classification',
            'variant_type': 'variant_type',
            'reference_allele': 'reference_allele',
            'tumor_seq_allele': 'tumor_seq_allele',
            'dbsnp_rs': 'dbsnp_rs',
            'dbsnp_val_status': 'dbsnp_val_status',
            'genome_change': 'genome_change',
            'annotation_transcript': 'annotation_transcript',
            'cdna_change': 'cdna_change',
            'codon_change': 'codon_change',
            'protein_change': 'protein_change',
            'is_deleterious': 'is_deleterious',
            'tcga_hotspot_count': 'tcga_hotspot_count',
            'cosmic_hotspot_count': 'cosmic_hotspot_count',
            'exac_af': 'exac_af',
            'cga_wes_ac': 'cga_wes_ac',
            'sanger_wes_ac': 'sanger_wes_ac',
            'sanger_recalibwes_ac': 'sanger_recalibwes_ac',
            'rnaseq_ac': 'rnaseq_ac',
            'hc_ac': 'hc_ac',
            'rd_ac': 'rd_ac',
            'wgs_ac': 'wgs_ac',
            'variant_annotation': 'variant_annotation'
        }

        self._entrez_gene_id = entrez_gene_id
        self._gene_symbol = gene_symbol
        self._depmap_id = depmap_id
        self._ncbi_build = ncbi_build
        self._chromosome = chromosome
        self._start_position = start_position
        self._end_position = end_position
        self._strand = strand
        self._variant_classification = variant_classification
        self._variant_type = variant_type
        self._reference_allele = reference_allele
        self._tumor_seq_allele = tumor_seq_allele
        self._dbsnp_rs = dbsnp_rs
        self._dbsnp_val_status = dbsnp_val_status
        self._genome_change = genome_change
        self._annotation_transcript = annotation_transcript
        self._cdna_change = cdna_change
        self._codon_change = codon_change
        self._protein_change = protein_change
        self._is_deleterious = is_deleterious
        self._tcga_hotspot_count = tcga_hotspot_count
        self._cosmic_hotspot_count = cosmic_hotspot_count
        self._exac_af = exac_af
        self._cga_wes_ac = cga_wes_ac
        self._sanger_wes_ac = sanger_wes_ac
        self._sanger_recalibwes_ac = sanger_recalibwes_ac
        self._rnaseq_ac = rnaseq_ac
        self._hc_ac = hc_ac
        self._rd_ac = rd_ac
        self._wgs_ac = wgs_ac
        self._variant_annotation = variant_annotation

    @classmethod
    def from_dict(cls, dikt) -> 'Mutation':
        """Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The mutation of this Mutation.  # noqa: E501
        :rtype: Mutation
        """
        return util.deserialize_model(dikt, cls)

    @property
    def entrez_gene_id(self) -> int:
        """Gets the entrez_gene_id of this Mutation.


        :return: The entrez_gene_id of this Mutation.
        :rtype: int
        """
        return self._entrez_gene_id

    @entrez_gene_id.setter
    def entrez_gene_id(self, entrez_gene_id: int):
        """Sets the entrez_gene_id of this Mutation.


        :param entrez_gene_id: The entrez_gene_id of this Mutation.
        :type entrez_gene_id: int
        """

        self._entrez_gene_id = entrez_gene_id

    @property
    def gene_symbol(self) -> str:
        """Gets the gene_symbol of this Mutation.


        :return: The gene_symbol of this Mutation.
        :rtype: str
        """
        return self._gene_symbol

    @gene_symbol.setter
    def gene_symbol(self, gene_symbol: str):
        """Sets the gene_symbol of this Mutation.


        :param gene_symbol: The gene_symbol of this Mutation.
        :type gene_symbol: str
        """

        self._gene_symbol = gene_symbol

    @property
    def depmap_id(self) -> str:
        """Gets the depmap_id of this Mutation.


        :return: The depmap_id of this Mutation.
        :rtype: str
        """
        return self._depmap_id

    @depmap_id.setter
    def depmap_id(self, depmap_id: str):
        """Sets the depmap_id of this Mutation.


        :param depmap_id: The depmap_id of this Mutation.
        :type depmap_id: str
        """

        self._depmap_id = depmap_id

    @property
    def ncbi_build(self) -> int:
        """Gets the ncbi_build of this Mutation.


        :return: The ncbi_build of this Mutation.
        :rtype: int
        """
        return self._ncbi_build

    @ncbi_build.setter
    def ncbi_build(self, ncbi_build: int):
        """Sets the ncbi_build of this Mutation.


        :param ncbi_build: The ncbi_build of this Mutation.
        :type ncbi_build: int
        """

        self._ncbi_build = ncbi_build

    @property
    def chromosome(self) -> str:
        """Gets the chromosome of this Mutation.


        :return: The chromosome of this Mutation.
        :rtype: str
        """
        return self._chromosome

    @chromosome.setter
    def chromosome(self, chromosome: str):
        """Sets the chromosome of this Mutation.


        :param chromosome: The chromosome of this Mutation.
        :type chromosome: str
        """

        self._chromosome = chromosome

    @property
    def start_position(self) -> int:
        """Gets the start_position of this Mutation.


        :return: The start_position of this Mutation.
        :rtype: int
        """
        return self._start_position

    @start_position.setter
    def start_position(self, start_position: int):
        """Sets the start_position of this Mutation.


        :param start_position: The start_position of this Mutation.
        :type start_position: int
        """

        self._start_position = start_position

    @property
    def end_position(self) -> int:
        """Gets the end_position of this Mutation.


        :return: The end_position of this Mutation.
        :rtype: int
        """
        return self._end_position

    @end_position.setter
    def end_position(self, end_position: int):
        """Sets the end_position of this Mutation.


        :param end_position: The end_position of this Mutation.
        :type end_position: int
        """

        self._end_position = end_position

    @property
    def strand(self) -> str:
        """Gets the strand of this Mutation.


        :return: The strand of this Mutation.
        :rtype: str
        """
        return self._strand

    @strand.setter
    def strand(self, strand: str):
        """Sets the strand of this Mutation.


        :param strand: The strand of this Mutation.
        :type strand: str
        """

        self._strand = strand

    @property
    def variant_classification(self) -> str:
        """Gets the variant_classification of this Mutation.


        :return: The variant_classification of this Mutation.
        :rtype: str
        """
        return self._variant_classification

    @variant_classification.setter
    def variant_classification(self, variant_classification: str):
        """Sets the variant_classification of this Mutation.


        :param variant_classification: The variant_classification of this Mutation.
        :type variant_classification: str
        """

        self._variant_classification = variant_classification

    @property
    def variant_type(self) -> str:
        """Gets the variant_type of this Mutation.


        :return: The variant_type of this Mutation.
        :rtype: str
        """
        return self._variant_type

    @variant_type.setter
    def variant_type(self, variant_type: str):
        """Sets the variant_type of this Mutation.


        :param variant_type: The variant_type of this Mutation.
        :type variant_type: str
        """

        self._variant_type = variant_type

    @property
    def reference_allele(self) -> str:
        """Gets the reference_allele of this Mutation.


        :return: The reference_allele of this Mutation.
        :rtype: str
        """
        return self._reference_allele

    @reference_allele.setter
    def reference_allele(self, reference_allele: str):
        """Sets the reference_allele of this Mutation.


        :param reference_allele: The reference_allele of this Mutation.
        :type reference_allele: str
        """

        self._reference_allele = reference_allele

    @property
    def tumor_seq_allele(self) -> str:
        """Gets the tumor_seq_allele of this Mutation.


        :return: The tumor_seq_allele of this Mutation.
        :rtype: str
        """
        return self._tumor_seq_allele

    @tumor_seq_allele.setter
    def tumor_seq_allele(self, tumor_seq_allele: str):
        """Sets the tumor_seq_allele of this Mutation.


        :param tumor_seq_allele: The tumor_seq_allele of this Mutation.
        :type tumor_seq_allele: str
        """

        self._tumor_seq_allele = tumor_seq_allele

    @property
    def dbsnp_rs(self) -> List[str]:
        """Gets the dbsnp_rs of this Mutation.


        :return: The dbsnp_rs of this Mutation.
        :rtype: List[str]
        """
        return self._dbsnp_rs

    @dbsnp_rs.setter
    def dbsnp_rs(self, dbsnp_rs: List[str]):
        """Sets the dbsnp_rs of this Mutation.


        :param dbsnp_rs: The dbsnp_rs of this Mutation.
        :type dbsnp_rs: List[str]
        """

        self._dbsnp_rs = dbsnp_rs

    @property
    def dbsnp_val_status(self) -> str:
        """Gets the dbsnp_val_status of this Mutation.


        :return: The dbsnp_val_status of this Mutation.
        :rtype: str
        """
        return self._dbsnp_val_status

    @dbsnp_val_status.setter
    def dbsnp_val_status(self, dbsnp_val_status: str):
        """Sets the dbsnp_val_status of this Mutation.


        :param dbsnp_val_status: The dbsnp_val_status of this Mutation.
        :type dbsnp_val_status: str
        """

        self._dbsnp_val_status = dbsnp_val_status

    @property
    def genome_change(self) -> str:
        """Gets the genome_change of this Mutation.


        :return: The genome_change of this Mutation.
        :rtype: str
        """
        return self._genome_change

    @genome_change.setter
    def genome_change(self, genome_change: str):
        """Sets the genome_change of this Mutation.


        :param genome_change: The genome_change of this Mutation.
        :type genome_change: str
        """

        self._genome_change = genome_change

    @property
    def annotation_transcript(self) -> str:
        """Gets the annotation_transcript of this Mutation.


        :return: The annotation_transcript of this Mutation.
        :rtype: str
        """
        return self._annotation_transcript

    @annotation_transcript.setter
    def annotation_transcript(self, annotation_transcript: str):
        """Sets the annotation_transcript of this Mutation.


        :param annotation_transcript: The annotation_transcript of this Mutation.
        :type annotation_transcript: str
        """

        self._annotation_transcript = annotation_transcript

    @property
    def cdna_change(self) -> str:
        """Gets the cdna_change of this Mutation.


        :return: The cdna_change of this Mutation.
        :rtype: str
        """
        return self._cdna_change

    @cdna_change.setter
    def cdna_change(self, cdna_change: str):
        """Sets the cdna_change of this Mutation.


        :param cdna_change: The cdna_change of this Mutation.
        :type cdna_change: str
        """

        self._cdna_change = cdna_change

    @property
    def codon_change(self) -> str:
        """Gets the codon_change of this Mutation.


        :return: The codon_change of this Mutation.
        :rtype: str
        """
        return self._codon_change

    @codon_change.setter
    def codon_change(self, codon_change: str):
        """Sets the codon_change of this Mutation.


        :param codon_change: The codon_change of this Mutation.
        :type codon_change: str
        """

        self._codon_change = codon_change

    @property
    def protein_change(self) -> str:
        """Gets the protein_change of this Mutation.


        :return: The protein_change of this Mutation.
        :rtype: str
        """
        return self._protein_change

    @protein_change.setter
    def protein_change(self, protein_change: str):
        """Sets the protein_change of this Mutation.


        :param protein_change: The protein_change of this Mutation.
        :type protein_change: str
        """

        self._protein_change = protein_change

    @property
    def is_deleterious(self) -> str:
        """Gets the is_deleterious of this Mutation.


        :return: The is_deleterious of this Mutation.
        :rtype: str
        """
        return self._is_deleterious

    @is_deleterious.setter
    def is_deleterious(self, is_deleterious: str):
        """Sets the is_deleterious of this Mutation.


        :param is_deleterious: The is_deleterious of this Mutation.
        :type is_deleterious: str
        """

        self._is_deleterious = is_deleterious

    @property
    def tcga_hotspot_count(self) -> str:
        """Gets the tcga_hotspot_count of this Mutation.


        :return: The tcga_hotspot_count of this Mutation.
        :rtype: str
        """
        return self._tcga_hotspot_count

    @tcga_hotspot_count.setter
    def tcga_hotspot_count(self, tcga_hotspot_count: str):
        """Sets the tcga_hotspot_count of this Mutation.


        :param tcga_hotspot_count: The tcga_hotspot_count of this Mutation.
        :type tcga_hotspot_count: str
        """

        self._tcga_hotspot_count = tcga_hotspot_count

    @property
    def cosmic_hotspot_count(self) -> str:
        """Gets the cosmic_hotspot_count of this Mutation.


        :return: The cosmic_hotspot_count of this Mutation.
        :rtype: str
        """
        return self._cosmic_hotspot_count

    @cosmic_hotspot_count.setter
    def cosmic_hotspot_count(self, cosmic_hotspot_count: str):
        """Sets the cosmic_hotspot_count of this Mutation.


        :param cosmic_hotspot_count: The cosmic_hotspot_count of this Mutation.
        :type cosmic_hotspot_count: str
        """

        self._cosmic_hotspot_count = cosmic_hotspot_count

    @property
    def exac_af(self) -> str:
        """Gets the exac_af of this Mutation.


        :return: The exac_af of this Mutation.
        :rtype: str
        """
        return self._exac_af

    @exac_af.setter
    def exac_af(self, exac_af: str):
        """Sets the exac_af of this Mutation.


        :param exac_af: The exac_af of this Mutation.
        :type exac_af: str
        """

        self._exac_af = exac_af

    @property
    def cga_wes_ac(self) -> str:
        """Gets the cga_wes_ac of this Mutation.


        :return: The cga_wes_ac of this Mutation.
        :rtype: str
        """
        return self._cga_wes_ac

    @cga_wes_ac.setter
    def cga_wes_ac(self, cga_wes_ac: str):
        """Sets the cga_wes_ac of this Mutation.


        :param cga_wes_ac: The cga_wes_ac of this Mutation.
        :type cga_wes_ac: str
        """

        self._cga_wes_ac = cga_wes_ac

    @property
    def sanger_wes_ac(self) -> str:
        """Gets the sanger_wes_ac of this Mutation.


        :return: The sanger_wes_ac of this Mutation.
        :rtype: str
        """
        return self._sanger_wes_ac

    @sanger_wes_ac.setter
    def sanger_wes_ac(self, sanger_wes_ac: str):
        """Sets the sanger_wes_ac of this Mutation.


        :param sanger_wes_ac: The sanger_wes_ac of this Mutation.
        :type sanger_wes_ac: str
        """

        self._sanger_wes_ac = sanger_wes_ac

    @property
    def sanger_recalibwes_ac(self) -> str:
        """Gets the sanger_recalibwes_ac of this Mutation.


        :return: The sanger_recalibwes_ac of this Mutation.
        :rtype: str
        """
        return self._sanger_recalibwes_ac

    @sanger_recalibwes_ac.setter
    def sanger_recalibwes_ac(self, sanger_recalibwes_ac: str):
        """Sets the sanger_recalibwes_ac of this Mutation.


        :param sanger_recalibwes_ac: The sanger_recalibwes_ac of this Mutation.
        :type sanger_recalibwes_ac: str
        """

        self._sanger_recalibwes_ac = sanger_recalibwes_ac

    @property
    def rnaseq_ac(self) -> str:
        """Gets the rnaseq_ac of this Mutation.


        :return: The rnaseq_ac of this Mutation.
        :rtype: str
        """
        return self._rnaseq_ac

    @rnaseq_ac.setter
    def rnaseq_ac(self, rnaseq_ac: str):
        """Sets the rnaseq_ac of this Mutation.


        :param rnaseq_ac: The rnaseq_ac of this Mutation.
        :type rnaseq_ac: str
        """

        self._rnaseq_ac = rnaseq_ac

    @property
    def hc_ac(self) -> str:
        """Gets the hc_ac of this Mutation.


        :return: The hc_ac of this Mutation.
        :rtype: str
        """
        return self._hc_ac

    @hc_ac.setter
    def hc_ac(self, hc_ac: str):
        """Sets the hc_ac of this Mutation.


        :param hc_ac: The hc_ac of this Mutation.
        :type hc_ac: str
        """

        self._hc_ac = hc_ac

    @property
    def rd_ac(self) -> str:
        """Gets the rd_ac of this Mutation.


        :return: The rd_ac of this Mutation.
        :rtype: str
        """
        return self._rd_ac

    @rd_ac.setter
    def rd_ac(self, rd_ac: str):
        """Sets the rd_ac of this Mutation.


        :param rd_ac: The rd_ac of this Mutation.
        :type rd_ac: str
        """

        self._rd_ac = rd_ac

    @property
    def wgs_ac(self) -> str:
        """Gets the wgs_ac of this Mutation.


        :return: The wgs_ac of this Mutation.
        :rtype: str
        """
        return self._wgs_ac

    @wgs_ac.setter
    def wgs_ac(self, wgs_ac: str):
        """Sets the wgs_ac of this Mutation.


        :param wgs_ac: The wgs_ac of this Mutation.
        :type wgs_ac: str
        """

        self._wgs_ac = wgs_ac

    @property
    def variant_annotation(self) -> str:
        """Gets the variant_annotation of this Mutation.


        :return: The variant_annotation of this Mutation.
        :rtype: str
        """
        return self._variant_annotation

    @variant_annotation.setter
    def variant_annotation(self, variant_annotation: str):
        """Sets the variant_annotation of this Mutation.


        :param variant_annotation: The variant_annotation of this Mutation.
        :type variant_annotation: str
        """

        self._variant_annotation = variant_annotation
