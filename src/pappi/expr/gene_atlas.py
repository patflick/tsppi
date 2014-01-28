from .expr import GeneExpression
from .. import id_mapping
from .. import sql
import csv


# TODO remove this one a proper handleing of id_mapping tables is implemented
from .. import data_config


class GeneAtlas(GeneExpression):
    def __init__(self, filename, sql_connection):
        """
        """
        # firstly initialize the superclass
        GeneExpression.__init__(self, filename, sql_connection, "gene_atlas")

        # set the raw file properties
        self.file_field_seperator = ","

    def linearize_table(self):
        """
        The GeneAtlas table is in a [Gene Annotation]x[Tissue] format, this has
        to be linearized before continuing.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('linear')
        linear_cols = ['Gene_ID']
        sql.linearize_table(src_table, linear_cols, "Tissue",
                            "ExpressionValue", dst_table, self.sql_conn)

    def id_mapping(self):
        """
        Maps the gene identifier to the unifying ID defined in `id_mapping`.
        """
        # TODO do properly once we have a generalized id mapping function
        #      with generalized mapping tables

        # import the annotation tables and create the mapping table
        id_mapping.import_u133a_annot_file(data_config.U133A_ANNOT_FILE,
                                           self.sql_conn)
        id_mapping.import_gnf1h_annot_file(data_config.GNF1H_ANNOT_FILE,
                                           self.sql_conn)
        id_mapping.create_mapping_table('chip_annot', 'hgnc', self.sql_conn,
                                        True,
                                        [id_mapping.GNF1H_MAPPING_TABLE_NAME,
                                         id_mapping.U133A_MAPPING_TABLE_NAME])

        # first map to hgnc
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["Gene_ID"], "chip_annot",
                                  dst_table, 'hgnc',
                                  self.sql_conn, True)
        # then to the unifying ID
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["Gene_ID"], 'hgnc',
                                  dst_table, id_mapping.UNIFYING_ID,
                                  self.sql_conn, True)

    # GeneAtlas doesn't supply reliablilites, thus no filter operation will
    # be performed

    def normalize_table(self):
        """
        Puts the table into the format: [Gene], [Tissue], [ExpressionValue]
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('normalized')
        sqlquery = ('SELECT Cast(Gene_ID AS TEXT) as Gene, Tissue AS Type, '
                    'Tissue, CAST(ExpressionValue AS REAL) AS ExpressionValue '
                    'FROM ' + src_table + ' '
                    'WHERE ExpressionValue != ""')
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def classify(self):
        """
        Classifies the ExpressionValue into a binary value Expressed.
        """
        # set the classification condition
        self.classify_cond = ">= 100.0"
        # call the super-class function
        GeneExpression.classify(self)
