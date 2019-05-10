
import csv
import sys
import sqlite3
import requests
from collections import namedtuple

def db_connection(database):
    """
        Connect to an SQL database
    """
    return sqlite3.connect(database)


def create_metadata_table(connection, table_name, headers):
    """
        Create metadata table in SQLite database, use headers from a file as column names
    """
    cur = connection.cursor()
    statement = "CREATE TABLE "+table_name+" ( "+table_name+"_ID INTEGER PRIMARY KEY NOT NULL, "
    statement += ", ".join([header+" TEXT" for header in headers])+");"
    cur.execute(statement)
    cur.close()


def load_metadata(connection, filename, table_name, id_index=[0]):
    """
        Load metadata table from a file to the SQLite database
    """
    id_map = {}

    with open(sys.argv[1]+'/'+filename,'r') as cell_line_file:
        input = csv.reader(cell_line_file)
        headers = [header.replace(' ','_').upper() for header in next(input)]

        create_metadata_table(connection, table_name, headers)
        id = 1
        cur = connection.cursor()
        for row in input:
            statement = "INSERT INTO "+table_name+" ("+table_name+"_ID, "+", ".join(headers)+")\n"
            statement += "VALUES ("+str(id)+", '"+"', '".join(row)+"');"
            cur.execute(statement)
            for index in id_index:
                id_map[row[index]]=id
            id += 1
        cur.close()
        connection.commit()
    return id_map


def load_cell_lines(connection):
    """
        Load cell lines from a file to the SQLite database
    """
    cell_line_ids = load_metadata(connection, "DepMap-2019q1-celllines_v2.csv", "CELL_LINE", [0,1])
    cell_line_ids['KE97_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'] = cell_line_ids['ACH-000167']
    cell_line_ids['NCIH684_LIVER'] = cell_line_ids['ACH-000089']
    cell_line_ids['S117_SOFT_TISSUE'] = cell_line_ids['ACH-000037']
    cell_line_ids['AZ521_STOMACH'] = cell_line_ids['ACH-001015']
    return cell_line_ids


def create_gene_table(connection):
    """
        Create a gene table in the SQLite database
    """
    statement = """
        CREATE TABLE GENE (
          GENE_ID        INTEGER PRIMARY KEY NOT NULL,
          GENE_SYMBOL    TEXT,
          ENTREZ_GENE_ID INTEGER,
          ENSEMBL_GENE   TEXT,
          MIM_NUMBER     TEXT
        );
    """
    cur = connection.cursor()
    cur.execute(statement)
    cur.close()


def load_omim_genes():
    """
        Load gene mapping file from OMIM
    """
    omim_genes = {}
    with open(sys.argv[1]+"/mim2gene.txt",'r') as file:
        input = csv.reader(file,delimiter='\t')
        for row in input:
            if len(row) == 5 and row[1] == "gene":
                omim_gene = {
                    "MIM": "MIM:"+row[0],
                    "ENTREZ_GENE_ID": row[2],
                    "SYMBOL": row[3],
                    "ENSEMBL_GENE": row[4]
                }
                omim_genes[omim_gene['MIM']] = omim_gene
                omim_genes[omim_gene['ENSEMBL_GENE']] = omim_gene
                omim_genes[omim_gene['ENTREZ_GENE_ID']] = omim_gene
    return omim_genes


def load_ensembl_genes(cur, omim_genes, genes):
    """
        Load genes and their ensembl ids from gene-expression file to the SQLite database
    """
    gene_id = 1
    with open(sys.argv[1]+"/CCLE_depMap_19Q1_TPM.csv",'r') as file:
        input = csv.reader(file)
        headers = next(input)
        for header in headers:
            gene_info = header.split(' ')
            if len(gene_info) == 2:
                symbol = gene_info[0]
                ensembl = gene_info[1].strip('( )')
                genes.headers[header] = gene_id
                genes.ensembl[ensembl] = gene_id
                if symbol not in genes.symbols:
                    genes.symbols[symbol] = gene_id
                statement = """
                    INSERT INTO GENE (GENE_ID, GENE_SYMBOL, ENSEMBL_GENE)
                    VALUES ({gene_id}, '{symbol}', '{ensembl}');
                """.format(gene_id=gene_id, symbol=symbol, ensembl=ensembl)
                cur.execute(statement)
                if ensembl in omim_genes:
                    omim_gene = omim_genes[ensembl]
                    statement = """
                        UPDATE GENE SET MIM_NUMBER = '{mim}'
                        WHERE GENE_ID = {gene_id};
                    """.format(mim=omim_gene['MIM'], gene_id=gene_id)
                    cur.execute(statement)
                gene_id += 1
    return gene_id


MY_GENE_INFO_URL = "http://mygene.info/v3/gene/{}?fields=ensembl"


def find_entrez_gene(symbol, entrez, omim_genes, genes):
    """
        Match entrez gene to an existing entrez gene or an ensembl gene
    """
    if entrez in genes.entrez:
        return genes.entrez[entrez]
    if entrez in omim_genes:
        omim_gene = omim_genes[entrez]
        if omim_gene['ENSEMBL_GENE'] in genes.ensembl:
            return genes.ensembl.get(omim_gene['ENSEMBL_GENE'])
    if symbol in genes.symbols:
        return genes.symbols[symbol]
    my_gene_info = requests.get(MY_GENE_INFO_URL.format(entrez)).json()
    if 'ensembl' in my_gene_info and 'gene' in my_gene_info['ensembl']:
        ensembl = my_gene_info['ensembl']
        return genes.ensembl.get(ensembl['gene'])
    return None


def add_entrez_gene(cur, header, symbol, entrez, db_gene_id, omim_genes, genes, gene_id):
    """
        Add entez gene to the SQLite database
    """
    if db_gene_id != None:
        if header != None:
            genes.headers[header] = db_gene_id
        if entrez not in genes.entrez and entrez != '0':
            genes.entrez[entrez] = db_gene_id
            statement = """
                UPDATE GENE SET ENTREZ_GENE_ID = {entrez}
                WHERE GENE_ID = {gene_id};
            """.format(entrez=entrez, gene_id=db_gene_id)
            cur.execute(statement)
    else:
        if header != None:
            genes.headers[header] = gene_id
        if symbol not in genes.symbols:
            genes.symbols[symbol] = gene_id
        if entrez == '0':
            statement = """
                INSERT INTO GENE (GENE_ID, GENE_SYMBOL)
                VALUES ({gene_id}, '{symbol}');
            """.format(gene_id=gene_id, symbol=symbol)
        else:
            genes.entrez[entrez] = gene_id
            statement = """
                INSERT INTO GENE (GENE_ID, GENE_SYMBOL, ENTREZ_GENE_ID)
                VALUES ({gene_id}, '{symbol}', {entrez});
            """.format(gene_id=gene_id, symbol=symbol, entrez=entrez)
        cur.execute(statement)
        if entrez in omim_genes:
            omim_gene = omim_genes[entrez]
            statement = """
                UPDATE GENE SET MIM_NUMBER = '{mim}'
                WHERE GENE_ID = {gene_id};
            """.format(mim=omim_gene['MIM'], gene_id=gene_id)
            cur.execute(statement)
        gene_id += 1
    return gene_id


def load_entrez_genes(cur, omim_genes, genes, gene_id, filename):
    """
        Load entrez genes from a data file; genes as columns
    """
    with open(sys.argv[1]+'/'+filename,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        for header in headers:
            gene_info = header.split(' ')
            if header not in genes.headers and len(gene_info) == 2:
                symbol = gene_info[0]
                entrez = gene_info[1].strip('( )')
                db_gene_id = find_entrez_gene(symbol, entrez, omim_genes, genes)
                gene_id = add_entrez_gene(cur, header, symbol, entrez, db_gene_id, omim_genes, genes, gene_id)
    return gene_id


def load_rnai_genes(cur, omim_genes, genes, gene_id, filename):
    """
        Load entrez genes from a data file; gene combinations as rows
    """
    gene_combos = []
    with open(sys.argv[1]+'/'+filename,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        for row in input:
            header = row[0]
            gene_info = header.split(' ')
            if len(gene_info) == 2:
                symbols = gene_info[0].split('&')
                entrez_list = gene_info[1].strip('( )').split('&')
                if len(symbols) == len(entrez_list):
                    combo = []
                    for symbol, entrez in zip(symbols, entrez_list):
                        db_gene_id = find_entrez_gene(symbol, entrez, omim_genes, genes)
                        gene_id = add_entrez_gene(cur, None, symbol, entrez, db_gene_id, omim_genes, genes, gene_id)
                        combo.append(genes.entrez[entrez])
                    gene_combos.append({"header": header, "genes":combo})
                else:
                    print("WARNING: "+header+" uneven sizes")
    print(str(gene_id-1)+" genes loaded")
    return gene_combos


def create_gene_combo_table(connection):
    """
        Create a gene combinations table in the SQLite database
    """
    statement = """
        CREATE TABLE GENE_COMBO (
          GENE_COMBO_ID  INTEGER,
          GENE_ID        INTEGER
        );
    """
    cur = connection.cursor()
    cur.execute(statement)
    cur.close()


def load_gene_combos(connection, gene_combos, genes):
    """
        Load gene combinations to the SQLite database
    """
    create_gene_combo_table(connection)
    cur = connection.cursor()
    gene_combo_id = 1
    for combo in gene_combos:
        genes.combos[combo['header']] = gene_combo_id
        for db_gene_id in combo['genes']:
            statement = """
                INSERT INTO GENE_COMBO (GENE_COMBO_ID, GENE_ID)
                VALUES ({gene_combo_id}, {gene_id});
            """.format(gene_combo_id=gene_combo_id, gene_id=db_gene_id)
            cur.execute(statement)
        gene_combo_id += 1
    cur.close()
    connection.commit()


def load_genes(connection, omim_genes):
    """
        Load genes from data files to the SQLite database
    """
    create_gene_table(connection)
    Genes = namedtuple('Genes',['headers','ensembl','symbols','entrez','combos'])
    genes = Genes({},{},{},{},{})
    cur = connection.cursor()
    gene_id = load_ensembl_genes(cur, omim_genes, genes)
    print(str(gene_id-1)+" genes loaded")
    gene_id = load_entrez_genes(cur, omim_genes, genes, gene_id, "gene_effect_corrected.csv")
    print(str(gene_id-1)+" genes loaded")
    gene_id = load_entrez_genes(cur, omim_genes, genes, gene_id, "public_19Q1_gene_cn.csv")
    print(str(gene_id-1)+" genes loaded")
    gene_combos = load_rnai_genes(cur, omim_genes, genes, gene_id, "D2_combined_gene_dep_scores.csv")
    cur.close()
    connection.commit()
    load_gene_combos(connection, gene_combos, genes)
    return genes


def create_value_table(connection, table):
    """
        Create a value table in the SQLite database
    """
    statement = """
        CREATE TABLE {0} (
          {0}_ID        INTEGER PRIMARY KEY,
          {1}_ID  INTEGER NOT NULL,
          {2}_ID       INTEGER NOT NULL,
          VALUE         REAL,
          FOREIGN KEY({1}_ID) REFERENCES {1}({1}_ID),
          FOREIGN KEY({2}_ID) REFERENCES {2}({2}_ID)
          );
    """ .format(table.name, table.dim1, table.dim2)
    cur = connection.cursor()
    cur.execute(statement)
    cur.close()


def insert_value(cur, table, id1, id2, value):
    statement = """
        INSERT INTO {table.name} ({table.dim1}_ID, {table.dim2}_ID, VALUE)
        VALUES ({id1}, {id2}, {value})
    """.format(table=table, id1=id1, id2=id2, value=value)
    cur.execute(statement)


def load_values(file, connection, table, cell_line_map, gene_id_map):
    """
        Load values from a file to the SQLite database
    """
    create_value_table(connection, table)
    with open(sys.argv[1]+'/'+file,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        gene_ids = [ gene_id_map.get(header) for header in headers]

        for row in input:
            if row[0] in cell_line_map:
                cell_line_id = cell_line_map[row[0]]
                cur = connection.cursor()
                for i in range(1,len(row)):
                    gene_id = gene_ids[i]
                    value = row[i]
                    if value != 'NA':
                        insert_value(cur, table, cell_line_id, gene_id, value)
                cur.close()
                connection.commit()
            else:
                print(row[0]+" NOT FOUND")


def load_rnai_values(file, connection, table, cell_line_map, gene_combo_map):
    """
        Load rnai values from a file to the SQLite database
    """
    create_value_table(connection, table)
    with open(sys.argv[1]+'/'+file,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        for header in headers:
            if header not in cell_line_map and header != '':
                print(header+" NOT FOUND")
        cell_lines = [cell_line_map.get(header) for header in headers]

        for row in input:
            if row[0] in gene_combo_map:
                gene_combo_id = gene_combo_map[row[0]]
                cur = connection.cursor()
                for i in range(1,len(row)):
                    cell_line_id = cell_lines[i]
                    value = row[i]
                    if cell_line_id != None and value != 'NA':
                        insert_value(cur, table, cell_line_id, gene_combo_id, value)
                cur.close()
                connection.commit()


def create_mutation_tables(connection):
    """
        Create tables for mutation information
    """
    statements = [
        """
        CREATE TABLE VARIANT_TYPE (
            VARIANT_TYPE_ID         INTEGER PRIMARY KEY NOT NULL,
            VARIANT_TYPE            TEXT NOT NULL,
            VARIANT_CLASSIFICATION  TEXT,
            VARIANT_ANNOTATION      TEXT
        );
        """,

        """
        CREATE TABLE ALLELE (
            ALLELE_ID         INTEGER PRIMARY KEY NOT NULL,
            REFERENCE_ALLELE  TEXT,
            TUMOR_SEQ_ALLELE  TEXT
        );
        """,

        """
        CREATE TABLE MUTATION (
          MUTATION_ID            INTEGER PRIMARY KEY,
          CELL_LINE_ID           INTEGER NOT NULL,
          GENE_ID                INTEGER NOT NULL,
          VARIANT_TYPE_ID        INTEGER NOT NULL,
          ALLELE_ID              INTEGER NOT NULL,
          NCBI_BUILD             INTEGER,
          CHROMOSOME             TEXT,
          START_POSITION         INTEGER,
          END_POSITION           INTEGER,
          STRAND                 TEXT,
          GENOME_CHANGE          TEXT,
          ANNOTATION_TRANSCRIPT  TEXT,
          CDNA_CHANGE            TEXT,
          CODON_CHANGE           TEXT,
          PROTEIN_CHANGE         TEXT,

          FOREIGN KEY(GENE_ID) REFERENCES GENE(GENE_ID),
          FOREIGN KEY(CELL_LINE_ID) REFERENCES CELL_LINE(CELL_LINE_ID),
          FOREIGN KEY(VARIANT_TYPE_ID) REFERENCES VARIANT_TYPE(VARIANT_TYPE_ID),
          FOREIGN KEY(ALLELE_ID) REFERENCES ALLELE(ALLELE_ID)
        );
        """,

        """
        CREATE TABLE MUTATION_PROPERTY (
            PROPERTY_ID   INTEGER NOT NULL,
            PROPERTY_NAME  TEXT NOT NULL
        );
        """,

        """
        CREATE TABLE MUTATION_INFO (
            MUTATION_INFO_ID  INTEGER PRIMARY KEY,
            MUTATION_ID       INTEGER NOT NULL,
            PROPERTY_ID       INTEGER NOT NULL,
            VALUE             TEXT NOT NULL,

            FOREIGN KEY(MUTATION_ID) REFERENCES MUTATION(MUTATION_ID)
            FOREIGN KEY(PROPERTY_ID) REFERENCES PROPERTY(PROPERTY_ID)
        );
        """
    ]
    for statement in statements:
        cur = connection.cursor()
        cur.execute(statement)
        cur.close()


def rm_na(string):
    """
        Replace NA with an empty string
    """
    if string == 'NA':
        return ''
    return string.replace("'","''")


def load_variants(connection, variant_set, mutation_metadata):
    """
        Load variant info in the SQLite database
    """
    cur = connection.cursor()
    for i,(var_class, var_type, var_anno) in enumerate(variant_set, 1):
        mutation_metadata.variants[(var_class, var_type, var_anno)] = i
        statement = """
            INSERT INTO VARIANT_TYPE (VARIANT_TYPE_ID, VARIANT_TYPE, VARIANT_CLASSIFICATION, VARIANT_ANNOTATION)
            VALUES ({id}, '{var_type}', '{var_class}', '{var_anno}');
        """.format(id=i, var_type=var_type, var_class=var_class, var_anno=var_anno)
        cur.execute(statement)

    cur.close()
    connection.commit()
    print(str(len(mutation_metadata.variants))+" variants loaded")


def load_alleles(connection, allele_set, mutation_metadata):
    """
        Load allele info in the SQLite database
    """
    cur = connection.cursor()
    for i,(ref_allele,tumor_allele) in enumerate(allele_set, 1):
        mutation_metadata.alleles[(ref_allele,tumor_allele)] = i
        statement = """
            INSERT INTO ALLELE (ALLELE_ID, REFERENCE_ALLELE, TUMOR_SEQ_ALLELE)
            VALUES ({id}, '{ref_allele}', '{tumor_allele}');
        """.format(id=i, ref_allele=ref_allele, tumor_allele=tumor_allele)
        cur.execute(statement)

    cur.close()
    connection.commit()
    print(str(len(mutation_metadata.alleles))+" alleles loaded")


def load_mutation_metadata(connection, filename, omim_genes, genes):
    """
        Load mutation metadata from a file
    """
    MutationMetadata = namedtuple('MutationMetadata',['genes','variants','alleles'])
    mutation_metadata = MutationMetadata({},{},{})
    variant_set = set()
    allele_set = set()
    gene_id = max(genes.entrez.values())+1
    with open(sys.argv[1]+'/'+filename,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        symbol_index = headers.index('Hugo_Symbol')
        entrez_index = headers.index('Entrez_Gene_Id')
        var_class_index = headers.index('Variant_Classification')
        var_type_index = headers.index('Variant_Type')
        ref_allele_index = headers.index('Reference_Allele')
        tumor_allele_index = headers.index('Tumor_Seq_Allele1')
        var_anno_index = headers.index('Variant_annotation')

        cur = connection.cursor()
        for row in input:
            symbol = row[symbol_index]
            entrez = row[entrez_index]
            key = symbol + ' ('+entrez+')'
            if key not in mutation_metadata.genes:
                db_gene_id = find_entrez_gene(symbol, entrez, omim_genes, genes)
                if db_gene_id == None and entrez == '0':
                    URL = "http://mygene.info/v3/query?q={}&species=human"
                    my_gene_info = requests.get(URL.format(symbol)).json()
                    if my_gene_info['total'] == 1:
                        hit = my_gene_info['hits'][0]
                        symbol = hit['symbol']
                        if 'entrezgene' in hit:
                            entrez = str(hit['entrezgene'])
                        db_gene_id = find_entrez_gene(symbol, entrez, omim_genes, genes)
                gene_id = add_entrez_gene(cur, None, symbol, entrez, db_gene_id, omim_genes, genes, gene_id)
                if db_gene_id == None:
                    db_gene_id = gene_id-1
                mutation_metadata.genes[key] = db_gene_id

            var_class = rm_na(row[var_class_index])
            var_type = rm_na(row[var_type_index])
            var_anno = rm_na(row[var_anno_index])
            variant_set.add((var_class, var_type, var_anno))

            ref_allele = rm_na(row[ref_allele_index])
            tumor_allele = rm_na(row[tumor_allele_index])
            allele_set.add((ref_allele,tumor_allele))

        cur.close()
        connection.commit()
        print(str(gene_id-1)+" genes loaded")
        load_variants(connection, variant_set, mutation_metadata)
        load_alleles(connection, allele_set, mutation_metadata)

    return mutation_metadata


def load_property_names(connection, property_names, count_properties):
    """
        Load property names from a mutation file in the SQLite database
    """
    cur = connection.cursor()
    for property_id, property_name in enumerate(property_names,1):
        statement = """
            INSERT INTO MUTATION_PROPERTY (PROPERTY_ID, PROPERTY_NAME)
            VALUES ({id}, '{name}')
        """.format(id=property_id, name=property_name)
        cur.execute(statement)
    for property_id, (_,property_name) in enumerate(count_properties,property_id+1):
        statement = """
            INSERT INTO MUTATION_PROPERTY (PROPERTY_ID, PROPERTY_NAME)
            VALUES ({id}, '{name}')
        """.format(id=property_id, name=property_name)
        cur.execute(statement)
    cur.close()
    connection.commit()


def load_mutations(connection, filename, mutation_metadata, cell_line_map):
    """
        Load mutations from a file in the SQLite database
    """
    with open(sys.argv[1]+'/'+filename,'r') as file:
        input = csv.reader(file)
        headers = next(input)
        symbol_index = headers.index('Hugo_Symbol')
        entrez_index = headers.index('Entrez_Gene_Id')
        build_index = headers.index('NCBI_Build')
        chromosome_index = headers.index('Chromosome')
        start_index = headers.index('Start_position')
        end_index = headers.index('End_position')
        strand_index = headers.index('Strand')
        var_class_index = headers.index('Variant_Classification')
        var_type_index = headers.index('Variant_Type')
        ref_allele_index = headers.index('Reference_Allele')
        tumor_allele_index = headers.index('Tumor_Seq_Allele1')
        genome_change_index = headers.index('Genome_Change')
        transcript_index = headers.index('Annotation_Transcript')
        cell_line_index = headers.index('Tumor_Sample_Barcode')
        cdna_change_index = headers.index('cDNA_Change')
        codon_change_index = headers.index('Codon_Change')
        protein_change_index = headers.index('Protein_Change')
        var_anno_index = headers.index('Variant_annotation')

        property_names = ['dbSNP_RS','dbSNP_Val_Status','isDeleterious','ExAC_AF',
            'VA_WES_AC','CGA_WES_AC','SangerWES_AC','SangerRecalibWES_AC','RNAseq_AC',
            'HC_AC','RD_AC','WGS_AC']
        property_index = [headers.index(name) for name in property_names]
        count_properties = [('isTCGAhotspot','TCGAhsCnt'),('isCOSMIChotspot','COSMIChsCnt')]
        count_properties_index = [(headers.index(presence), headers.index(name))
            for presence, name in count_properties
        ]
        load_property_names(connection, property_names, count_properties)

        cur = connection.cursor()
        mut_id = 1
        for row in input:
            cell_line = row[cell_line_index]
            if cell_line in cell_line_map:
                cell_line_id = cell_line_map[cell_line]

                symbol = row[symbol_index]
                entrez = row[entrez_index]
                key = symbol + ' ('+entrez+')'
                gene_id = mutation_metadata.genes[key]

                var_class = rm_na(row[var_class_index])
                var_type = rm_na(row[var_type_index])
                var_anno = rm_na(row[var_anno_index])
                var_id = mutation_metadata.variants[(var_class, var_type, var_anno)]

                ref_allele = rm_na(row[ref_allele_index])
                tumor_allele = rm_na(row[tumor_allele_index])
                allele_id = mutation_metadata.alleles[(ref_allele,tumor_allele)]

                build = row[build_index]
                chromosome = row[chromosome_index]
                start_pos = row[start_index]
                end_pos = row[end_index]
                strand = row[strand_index]

                genome_change = rm_na(row[genome_change_index])
                transcript = rm_na(row[transcript_index])
                cdna_change = rm_na(row[cdna_change_index])
                codon_change = rm_na(row[codon_change_index])
                protein_change = rm_na(row[protein_change_index])

                statement = """
                    INSERT INTO MUTATION (MUTATION_ID, CELL_LINE_ID, GENE_ID,
                        VARIANT_TYPE_ID, ALLELE_ID, NCBI_BUILD, CHROMOSOME,
                        START_POSITION, END_POSITION, STRAND, GENOME_CHANGE,
                        ANNOTATION_TRANSCRIPT, CDNA_CHANGE, CODON_CHANGE, PROTEIN_CHANGE)
                    VALUES ({mut_id}, {cell_line_id}, {gene_id},
                        {var_id}, {allele_id}, {build}, '{chromosome}',
                        {start_pos}, {end_pos}, '{strand}', '{genome_change}',
                        '{transcript}', '{cdna_change}', '{codon_change}', '{protein_change}')
                """.format(mut_id=mut_id, cell_line_id=cell_line_id, gene_id=gene_id,
                        var_id=var_id, allele_id=allele_id, build=build, chromosome=chromosome,
                        start_pos=start_pos, end_pos=end_pos, strand=strand, genome_change=genome_change,
                        transcript=transcript, cdna_change=cdna_change, codon_change=codon_change,
                        protein_change=protein_change)
                cur.execute(statement)

                for property_id, index in enumerate(property_index,1):
                    property_value = rm_na(row[index])
                    if property_value != '':
                        statement = """
                            INSERT INTO MUTATION_INFO (MUTATION_ID, PROPERTY_ID, VALUE)
                            VALUES ({mut_id}, {prop_id}, '{value}')
                        """.format(mut_id=mut_id, prop_id=property_id, value=property_value)
                        cur.execute(statement)

                for property_id, (property_presence_index,property_name_index) in enumerate(count_properties_index,property_id+1):
                    property_present = rm_na(row[property_presence_index])
                    property_value = rm_na(row[property_name_index])
                    if property_present != 'FALSE' and property_value != '0' and property_value != '':
                        statement = """
                            INSERT INTO MUTATION_INFO (MUTATION_ID, PROPERTY_ID, VALUE)
                            VALUES ({mut_id}, {prop_id}, '{value}')
                        """.format(mut_id=mut_id, prop_id=property_id, value=property_value)
                        cur.execute(statement)

                mut_id += 1
            else:
                print(row[0]+" NOT FOUND")
        cur.close()
        connection.commit()


def build_indexes(connection):
    """
        Build indexes for gene and cell-line tables
    """
    statements = [
        "CREATE INDEX IDX_CL_ID ON CELL_LINE(CELL_LINE_ID);",
        "CREATE INDEX IDX_CL_NAME ON CELL_LINE(CCLE_NAME);",
        "CREATE INDEX IDX_CL_DMID ON CELL_LINE(DEPMAP_ID);",
        "CREATE INDEX IDX_GENE_ID ON GENE(GENE_ID);",
        "CREATE INDEX IDX_GENE_SYMBOL ON GENE(GENE_SYMBOL);",
        "CREATE INDEX IDX_GENE_ENTREZ ON GENE(ENTREZ_GENE_ID);",
        "CREATE INDEX IDX_GENE_ENSEMBL ON GENE(ENSEMBL_GENE);",
        "CREATE INDEX IDX_GENE_MIM ON GENE(MIM_NUMBER);",
        "CREATE INDEX IDX_COMBO_GCID ON GENE_COMBO(GENE_COMBO_ID);",
        "CREATE INDEX IDX_COMBO_GID ON GENE_COMBO(GENE_ID);",
        "CREATE INDEX IDX_PROD_ID ON PROTEIN(PROTEIN_ID);",
        "CREATE INDEX IDX_PROT_NAME ON PROTEIN(ANTIBODY_NAME);"
    ]
    for statement in statements:
        cur = connection.cursor()
        print(statement)
        cur.execute(statement)
        cur.close()


def build_value_indexes(connection, table):
    """
        Build indexes for value tables
    """
    statements = [
        "CREATE INDEX IDX_{table.abbrev}_{table.dim1} ON {table.name}({table.dim1}_ID);".format(table=table),
        "CREATE INDEX IDX_{table.abbrev}_{table.dim2} ON {table.name}({table.dim2}_ID);".format(table=table)
    ]
    for statement in statements:
        cur = connection.cursor()
        print(statement)
        cur.execute(statement)
        cur.close()


def build_mutation_indexes(connection):
    """
        Build indexes for mutation tables
    """
    statements = [
        "CREATE INDEX IDX_VAR_ID ON VARIANT_TYPE(VARIANT_TYPE_ID);",
        "CREATE INDEX IDX_ALLELE_ID ON ALLELE(ALLELE_ID);",
        "CREATE INDEX IDX_MUT_ID ON MUTATION(MUTATION_ID);",
        "CREATE INDEX IDX_MUT_CL ON MUTATION(CELL_LINE_ID);",
        "CREATE INDEX IDX_MUT_GENE ON MUTATION(GENE_ID);",
        "CREATE INDEX IDX_PROP_ID ON MUTATION_PROPERTY(PROPERTY_ID);",
        "CREATE INDEX IDX_MUT_INFO_ID ON MUTATION_INFO(MUTATION_ID);"
    ]
    for statement in statements:
        cur = connection.cursor()
        print(statement)
        cur.execute(statement)
        cur.close()


def main():
    """
        Load depMap data from files and load them into an SQLite database
    """
    print("loading "+sys.argv[1])
    connection = db_connection(sys.argv[2])
    print("loading cell lines")
    cell_line_map = load_cell_lines(connection)
    print("loading genes")
    omim_genes = load_omim_genes()
    genes = load_genes(connection, omim_genes)
    create_mutation_tables(connection)

    mutations_file = "depmap_19Q1_mutation_calls.csv"
    mutation_metadata=load_mutation_metadata(connection, mutations_file, omim_genes, genes)
    print(str(len(mutation_metadata.genes))+" mut genes")
    print(str(len(mutation_metadata.variants))+" mut variants")
    print(str(len(mutation_metadata.alleles))+" mut alleles")

    protein_map = load_metadata(connection, "CCLE_RPPA_Ab_info_20180123.csv", "PROTEIN")
    print(str(len(protein_map))+" proteins loaded")

    gene_id_map = genes.headers
    ValueTable = namedtuple('ValueTable',['name','abbrev','dim1','dim2'])

    print("loading gene expression")
    ge_table = ValueTable("GENE_EXPRESSION","GEX","CELL_LINE","GENE")
    gene_expr_file = "CCLE_depMap_19Q1_TPM.csv"
    load_values(gene_expr_file, connection, ge_table, cell_line_map, gene_id_map)

    print("loading gene knockdown")
    gk_table = ValueTable("GENE_KNOCKDOWN","GKD","CELL_LINE","GENE")
    gene_knockdown_file = "gene_effect_corrected.csv"
    load_values(gene_knockdown_file, connection, gk_table, cell_line_map, gene_id_map)

    print("loading copy number")
    cn_table = ValueTable("COPY_NUMBER","CN","CELL_LINE","GENE")
    gene_knockdown_file = "public_19Q1_gene_cn.csv"
    load_values(gene_knockdown_file, connection, cn_table, cell_line_map, gene_id_map)

    print("loading RNAi dependencies")
    rnai_table = ValueTable("RNAI","RNAI","CELL_LINE","GENE_COMBO")
    gene_knockdown_file = "D2_combined_gene_dep_scores.csv"
    load_rnai_values(gene_knockdown_file, connection, rnai_table, cell_line_map, genes.combos)

    print("loading mutations")
    load_mutations(connection, mutations_file, mutation_metadata, cell_line_map)

    print("loading protein array data")
    prot_table = ValueTable("PROTEIN_ARRAY","RPPA","CELL_LINE","PROTEIN")
    file = "CCLE_RPPA_20180123.csv"
    load_values(file, connection, prot_table, cell_line_map, protein_map)


    print("building indexes")
    build_indexes(connection)
    build_value_indexes(connection, ge_table)
    build_value_indexes(connection, gk_table)
    build_value_indexes(connection, cn_table)
    build_value_indexes(connection, rnai_table)
    build_value_indexes(connection, prot_table)
    build_mutation_indexes(connection)
    connection.close()
    print("done")


if __name__ == '__main__':
    main()
