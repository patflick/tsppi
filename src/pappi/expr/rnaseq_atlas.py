from .expr import GeneExpression
from .. import id_mapping
from .. import sql
import csv


class RnaSeqAtlas(GeneExpression):

    def __init__(self, filename, sql_connection):
        """
        """
        # firstly initialize the superclass
        GeneExpression.__init__(self, filename, sql_connection, "rnaseq_atlas")

        # set the raw file properties
        self.file_field_seperator = "\t"

    def linearize_table(self):
        """
        The RNA Seq Atlas table is in a [Gene]x[Tissue] format, this has to
        be linearized before continuing.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('linear')
        linear_cols = ["entrez_gene_id", "ensembl_gene_id", "hgnc_symbol",
                       "transcript", "transcript_length"]
        sql.linearize_table(src_table, linear_cols, "Tissue",
                            "ExpressionValue", dst_table, self.sql_conn)

    def id_mapping(self):
        """
        Maps the gene identifier to the unifying ID defined in `id_mapping`.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["hgnc_symbol"], "hgnc",
                                  dst_table, id_mapping.UNIFYING_ID,
                                  self.sql_conn, True)

    # RNA Seq Atlas doesn't supply reliablilites, thus no filter operation will
    # be performed

    def normalize_table(self):
        """
        Puts the table into the format: [Gene], [Tissue], [ExpressionValue]
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('normalized')
        sqlquery = ('SELECT hgnc_symbol as Gene, Tissue AS Type, Tissue, '
                    'ExpressionValue '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def classify(self):
        """
        Classifies the ExpressionValue into a binary value Expressed.
        """
        # set the classification condition
        # RPKM >= 1.0
        self.classify_cond = ">= 1.0"
        # call the super-class function
        GeneExpression.classify(self)
