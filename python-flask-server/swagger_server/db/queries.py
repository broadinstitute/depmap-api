
import sys
import sqlite3

from swagger_server.models.gene import Gene
from swagger_server.models.cell_line import CellLine
from swagger_server.models.protein import Protein

from swagger_server.models.copy_number_value import CopyNumberValue
from swagger_server.models.dependency_value import DependencyValue
from swagger_server.models.expression_value import ExpressionValue
from swagger_server.models.rnai_value import RnaiValue
from swagger_server.models.mutation import Mutation
from swagger_server.models.protein_array_value import ProteinArrayValue
from swagger_server.models.protein_array_value_protein import ProteinArrayValueProtein



connection = sqlite3.connect("DepMap.db", check_same_thread=False)

def select_all_genes():
    """
        Select all genes from SQLite database
    """

    query = """
        SELECT GENE_SYMBOL, ENTREZ_GENE_ID, ENSEMBL_GENE, MIM_NUMBER FROM GENE
    """
    cur = connection.cursor()
    cur.execute(query)
    rows = cur.fetchall()

    genes = []
    for row in rows:
        gene = Gene(gene_symbol=row[0], entrez_gene_id=row[1], ensembl_gene=row[2], omim=row[3])
        genes.append(gene)
    cur.close()
    return genes


def select_all_cell_lines():
    """
        Select all cell lines from SQLite database
    """

    query = """
        SELECT
            DEPMAP_ID,
            CCLE_NAME,
            ALIASES,
            COSMIC_ID,
            SANGER_ID,
            PRIMARY_DISEASE,
            SUBTYPE_DISEASE,
            GENDER,
            SOURCE
        FROM CELL_LINE
    """
    cur = connection.cursor()
    cur.execute(query)
    rows = cur.fetchall()

    cell_lines = []
    for row in rows:
        cell_line = CellLine(
            depmap_id = row[0],
            ccle_name = row[1],
            aliases = (row[2].split(';') if row[2] != "" else []),
            cosmic_id = row[3],
            sanger_id = row[4],
            primary_disease = row[5],
            subtype_disease = row[6],
            gender = row[7],
            source = row[8]
        )
        cell_lines.append(cell_line)
    cur.close()
    return cell_lines


def select_all_proteins():
    """
        Select all proteins from SQLite database
    """

    query = """
        SELECT
            ANTIBODY_NAME, TARGET_GENES, VALIDATION_STATUS, COMPANY, CATALOG_NUMBER
        FROM PROTEIN
    """
    cur = connection.cursor()
    cur.execute(query)
    rows = cur.fetchall()

    proteins = []
    for row in rows:
        target_genes = row[1].split(' ')
        protein = Protein(
            antibody_name=row[0],
            target_genes=target_genes,
            validation_status = row[2],
            company = row[3],
            catalog_number = row[4])
        proteins.append(protein)
    cur.close()
    return proteins


def select_copy_number_by_cell_line(depmap_id):
    """
        Select copy number values by cell line from SQLite database
    """

    query = """
        SELECT ENTREZ_GENE_ID, DEPMAP_ID, VALUE
        FROM COPY_NUMBER
        INNER JOIN GENE ON COPY_NUMBER.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON COPY_NUMBER.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE DEPMAP_ID = ?
    """
    return select_copy_number(query, depmap_id)


def select_copy_number_by_gene(entrez_gene_id):
    """
        Select copy number values by gene from SQLite database
    """

    query = """
        SELECT ENTREZ_GENE_ID, DEPMAP_ID, VALUE
        FROM COPY_NUMBER
        INNER JOIN GENE ON COPY_NUMBER.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON COPY_NUMBER.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE ENTREZ_GENE_ID = ?
    """
    return select_copy_number(query, entrez_gene_id)


def select_copy_number(query, query_id):
    """
        Execute copy-number query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    values = []
    for row in rows:
        value = CopyNumberValue(entrez_gene_id=row[0], depmap_id=row[1], value=row[2])
        values.append(value)
    cur.close()
    return values


def select_gene_dependency_by_cell_line(depmap_id):
    """
        Select gene-dependency values by cell line from SQLite database
    """

    query = """
        SELECT ENTREZ_GENE_ID, DEPMAP_ID, VALUE
        FROM GENE_KNOCKDOWN
        INNER JOIN GENE ON GENE_KNOCKDOWN.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON GENE_KNOCKDOWN.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE DEPMAP_ID = ?
    """
    return select_gene_dependency(query, depmap_id)


def select_gene_dependency_by_gene(entrez_gene_id):
    """
        Select gene-dependency values by gene from SQLite database
    """

    query = """
        SELECT ENTREZ_GENE_ID, DEPMAP_ID, VALUE
        FROM GENE_KNOCKDOWN
        INNER JOIN GENE ON GENE_KNOCKDOWN.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON GENE_KNOCKDOWN.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE ENTREZ_GENE_ID = ?
    """
    return select_gene_dependency(query, entrez_gene_id)


def select_gene_dependency(query, query_id):
    """
        Execute gene-dependency query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    values = []
    for row in rows:
        value = DependencyValue(entrez_gene_id=row[0], depmap_id=row[1], value=row[2])
        values.append(value)
    cur.close()
    return values


def select_gene_expression_by_cell_line(depmap_id):
    """
        Select gene-expression values by cell line from SQLite database
    """

    query = """
        SELECT ENSEMBL_GENE, DEPMAP_ID, VALUE
        FROM GENE_EXPRESSION
        INNER JOIN GENE ON GENE_EXPRESSION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON GENE_EXPRESSION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE DEPMAP_ID = ?
    """
    return select_gene_expression(query, depmap_id)


def select_gene_expression_by_gene(ensembl_gene):
    """
        Select gene-expression values by gene from SQLite database
    """

    query = """
        SELECT ENSEMBL_GENE, DEPMAP_ID, VALUE
        FROM GENE_EXPRESSION
        INNER JOIN GENE ON GENE_EXPRESSION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON GENE_EXPRESSION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE ENSEMBL_GENE = ?
    """
    return select_gene_expression(query, ensembl_gene)


def select_gene_expression(query, query_id):
    """
        Execute gene-expression query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    values = []
    for row in rows:
        value = ExpressionValue(ensembl_gene=row[0], depmap_id=row[1], value=row[2])
        values.append(value)
    cur.close()
    return values


def select_rnai_by_cell_line(depmap_id):
    """
        Select rnai values by cell line from SQLite database
    """

    query = """
        SELECT GENES, DEPMAP_ID, VALUE
        FROM RNAI
        INNER JOIN CELL_LINE ON RNAI.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        INNER JOIN (
          SELECT GENE_COMBO_ID, GROUP_CONCAT(ENTREZ_GENE_ID,',') AS GENES
          FROM GENE_COMBO
          INNER JOIN GENE ON GENE.GENE_ID = GENE_COMBO.GENE_ID
          GROUP BY GENE_COMBO_ID
        ) AS X ON X.GENE_COMBO_ID = RNAI.GENE_COMBO_ID
        WHERE DEPMAP_ID = ?
    """
    return select_rnai(query, depmap_id)


def select_rnai_by_gene(entrez_gene_id):
    """
        Select rnai values by gene from SQLite database
    """

    query = """
        SELECT GENES, DEPMAP_ID, VALUE
        FROM (
          SELECT DISTINCT GENE_COMBO_ID
          FROM GENE_COMBO
          INNER JOIN GENE ON GENE.GENE_ID = GENE_COMBO.GENE_ID
          WHERE ENTREZ_GENE_ID = ?
        ) AS GC
        INNER JOIN RNAI ON RNAI.GENE_COMBO_ID = GC.GENE_COMBO_ID
        INNER JOIN CELL_LINE ON RNAI.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        INNER JOIN (
          SELECT GENE_COMBO_ID, GROUP_CONCAT(ENTREZ_GENE_ID,',') AS GENES
          FROM GENE_COMBO
          INNER JOIN GENE ON GENE.GENE_ID = GENE_COMBO.GENE_ID
          GROUP BY GENE_COMBO_ID
        ) AS X ON X.GENE_COMBO_ID = RNAI.GENE_COMBO_ID
    """
    return select_rnai(query, entrez_gene_id)


def select_rnai(query, query_id):
    """
        Execute rnai query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    values = []
    for row in rows:
        entrez_gene_ids = row[0].split(',')
        value = RnaiValue(entrez_gene_ids=entrez_gene_ids, depmap_id=row[1], value=row[2])
        values.append(value)
    cur.close()
    return values


def select_protein_array_by_cell_line(depmap_id):
    """
        Select protein-array values by cell line from SQLite database
    """

    query = """
        SELECT ANTIBODY_NAME, TARGET_GENES, DEPMAP_ID, VALUE
        FROM PROTEIN_ARRAY
        INNER JOIN PROTEIN ON PROTEIN_ARRAY.PROTEIN_ID = PROTEIN.PROTEIN_ID
        INNER JOIN CELL_LINE ON PROTEIN_ARRAY.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE DEPMAP_ID = ?
    """
    return select_protein_array(query, depmap_id)


def select_protein_array_by_antibody_name(antibody_name):
    """
        Select protein-array values by antibody from SQLite database
    """

    query = """
        SELECT ANTIBODY_NAME, TARGET_GENES, DEPMAP_ID, VALUE
        FROM PROTEIN_ARRAY
        INNER JOIN PROTEIN ON PROTEIN_ARRAY.PROTEIN_ID = PROTEIN.PROTEIN_ID
        INNER JOIN CELL_LINE ON PROTEIN_ARRAY.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE ANTIBODY_NAME = ?
    """
    return select_protein_array(query, antibody_name)


def select_protein_array(query, query_id):
    """
        Execute protein-array query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    values = []
    for row in rows:
        protein = ProteinArrayValueProtein(antibody_name=row[0], target_genes=row[1].split(' '))
        value = ProteinArrayValue(protein=protein, depmap_id=row[2], value=row[3])
        values.append(value)
    cur.close()
    return values


def select_mutations_brief_by_cell_line(depmap_id):
    """
        Select mutations by cell line from SQLite database
    """

    query = """
        SELECT
          ENTREZ_GENE_ID, GENE_SYMBOL, DEPMAP_ID,
          NCBI_BUILD, CHROMOSOME, START_POSITION, END_POSITION,
          STRAND, GENOME_CHANGE, ANNOTATION_TRANSCRIPT, CDNA_CHANGE,
          CODON_CHANGE, PROTEIN_CHANGE
        FROM MUTATION
        INNER JOIN GENE ON MUTATION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON MUTATION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE DEPMAP_ID = ?
    """
    return select_mutations_brief(query, depmap_id)


def select_mutations_brief_by_gene(entrez_gene_id):
    """
        Select mutations by gene from SQLite database
    """

    query = """
        SELECT
          ENTREZ_GENE_ID, GENE_SYMBOL, DEPMAP_ID,
          NCBI_BUILD, CHROMOSOME, START_POSITION, END_POSITION,
          STRAND, GENOME_CHANGE, ANNOTATION_TRANSCRIPT, CDNA_CHANGE,
          CODON_CHANGE, PROTEIN_CHANGE
        FROM MUTATION
        INNER JOIN GENE ON MUTATION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON MUTATION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        WHERE ENTREZ_GENE_ID = ?
    """
    return select_mutations_brief(query, entrez_gene_id)


def select_mutations_brief(query, query_id):
    """
        Execute mutations query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    mutations = []
    for row in rows:
        mutation = Mutation(
            entrez_gene_id = row[0],
            gene_symbol = row[1],
            depmap_id = row[2],
            ncbi_build = row[3],
            chromosome = row[4],
            start_position = row[5],
            end_position = row[6],
            strand = row[7],
            genome_change = row[8],
            annotation_transcript = row[9],
            cdna_change = row[10],
            codon_change = row[11],
            protein_change = row[12]
        )
        mutations.append(mutation)
    cur.close()
    return mutations


def select_mutations_by_cell_line(depmap_id):
    """
        Select mutations by cell line from SQLite database
    """

    query = """
        SELECT
          ENTREZ_GENE_ID, GENE_SYMBOL, DEPMAP_ID,
          NCBI_BUILD, CHROMOSOME, START_POSITION, END_POSITION,
          STRAND, GENOME_CHANGE, ANNOTATION_TRANSCRIPT, CDNA_CHANGE,
          CODON_CHANGE, PROTEIN_CHANGE,
          VARIANT_CLASSIFICATION, VARIANT_TYPE, VARIANT_ANNOTATION,
          REFERENCE_ALLELE, TUMOR_SEQ_ALLELE,
          MUTATION_ID
        FROM MUTATION
        INNER JOIN GENE ON MUTATION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON MUTATION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        INNER JOIN VARIANT_TYPE ON VARIANT_TYPE.VARIANT_TYPE_ID = MUTATION.VARIANT_TYPE_ID
        INNER JOIN ALLELE ON ALLELE.ALLELE_ID = MUTATION.ALLELE_ID
        WHERE DEPMAP_ID = ?
    """
    mutation_info = select_mutation_info_by_cell_line(depmap_id)
    return select_mutations(query, depmap_id, mutation_info)


def select_mutations_by_gene(entrez_gene_id):
    """
        Select mutations by gene from SQLite database
    """

    query = """
        SELECT
          ENTREZ_GENE_ID, GENE_SYMBOL, DEPMAP_ID,
          NCBI_BUILD, CHROMOSOME, START_POSITION, END_POSITION,
          STRAND, GENOME_CHANGE, ANNOTATION_TRANSCRIPT, CDNA_CHANGE,
          CODON_CHANGE, PROTEIN_CHANGE,
          VARIANT_CLASSIFICATION, VARIANT_TYPE, VARIANT_ANNOTATION,
          REFERENCE_ALLELE, TUMOR_SEQ_ALLELE,
          MUTATION_ID
        FROM MUTATION
        INNER JOIN GENE ON MUTATION.GENE_ID = GENE.GENE_ID
        INNER JOIN CELL_LINE ON MUTATION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        INNER JOIN VARIANT_TYPE ON VARIANT_TYPE.VARIANT_TYPE_ID = MUTATION.VARIANT_TYPE_ID
        INNER JOIN ALLELE ON ALLELE.ALLELE_ID = MUTATION.ALLELE_ID
        WHERE ENTREZ_GENE_ID = ?
    """
    mutation_info = select_mutation_info_by_gene(entrez_gene_id)
    return select_mutations(query, entrez_gene_id, mutation_info)


def select_mutations(query, query_id, mutation_info):
    """
        Execute mutation query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    mutations = []
    for row in rows:
        mut_id = row[18]
        mut_info = {}
        if mut_id in mutation_info:
            mut_info = mutation_info[mut_id]
        dbsnp_rs = None
        if 'dbSNP_RS' in mut_info:
            dbsnp_rs = mut_info['dbSNP_RS'].split('|')

        mutation = Mutation(
            entrez_gene_id = row[0],
            gene_symbol = row[1],
            depmap_id = row[2],
            ncbi_build = row[3],
            chromosome = row[4],
            start_position = row[5],
            end_position = row[6],
            strand = row[7],
            genome_change = row[8],
            annotation_transcript = row[9],
            cdna_change = row[10],
            codon_change = row[11],
            protein_change = row[12],
            variant_classification = row[13],
            variant_type = row[14],
            variant_annotation = row[15],
            reference_allele = row[16],
            tumor_seq_allele = row[17],
            dbsnp_rs = dbsnp_rs,
            dbsnp_val_status = mut_info.get('dbSNP_Val_Status'),
            is_deleterious = mut_info.get('isDeleterious'),
            tcga_hotspot_count = mut_info.get('TCGAhsCnt'),
            cosmic_hotspot_count = mut_info.get('COSMIChsCnt'),
            exac_af = mut_info.get('ExAC_AF'),
            va_wes_ac = mut_info.get('VA_WES_AC'),
            cga_wes_ac = mut_info.get('CGA_WES_AC'),
            sanger_wes_ac = mut_info.get('SangerWES_AC'),
            sanger_recalibwes_ac = mut_info.get('SangerRecalibWES_AC'),
            rnaseq_ac = mut_info.get('RNAseq_AC'),
            hc_ac = mut_info.get('HC_AC'),
            rd_ac = mut_info.get('RD_AC'),
            wgs_ac = mut_info.get('WGS_AC')
        )
        mutations.append(mutation)
    cur.close()
    return mutations


def select_mutation_info_by_gene(entrez_gene_id):
    """
        Select mutation info by gene from SQLite database
    """

    query = """
        SELECT
          MUTATION.MUTATION_ID, PROPERTY_NAME, VALUE
        FROM MUTATION
        INNER JOIN GENE ON MUTATION.GENE_ID = GENE.GENE_ID
        INNER JOIN MUTATION_INFO ON MUTATION.MUTATION_ID = MUTATION_INFO.MUTATION_ID
        INNER JOIN MUTATION_PROPERTY ON MUTATION_INFO.PROPERTY_ID = MUTATION_PROPERTY.PROPERTY_ID
        WHERE ENTREZ_GENE_ID = ?
    """
    return select_mutation_info(query, entrez_gene_id)


def select_mutation_info_by_cell_line(depmap_id):
    """
        Select mutation info by cell line from SQLite database
    """

    query = """
        SELECT
          MUTATION.MUTATION_ID, PROPERTY_NAME, VALUE
        FROM MUTATION
        INNER JOIN CELL_LINE ON MUTATION.CELL_LINE_ID = CELL_LINE.CELL_LINE_ID
        INNER JOIN MUTATION_INFO ON MUTATION.MUTATION_ID = MUTATION_INFO.MUTATION_ID
        INNER JOIN MUTATION_PROPERTY ON MUTATION_INFO.PROPERTY_ID = MUTATION_PROPERTY.PROPERTY_ID
        WHERE DEPMAP_ID = ?
    """
    return select_mutation_info(query, depmap_id)


def select_mutation_info(query, query_id):
    """
        Execute mutation-info query and return results
    """
    cur = connection.cursor()
    cur.execute(query, (query_id,))
    rows = cur.fetchall()

    mutation_info = {}
    for row in rows:
        mut_id = row[0]
        property_name = row[1]
        value = row[2]
        if mut_id not in mutation_info:
            mutation_info[mut_id] = {}
        mutation_info[mut_id][property_name] = value
    cur.close()
    return mutation_info

def main():
    pass

if __name__ == '__main__':
    main()
