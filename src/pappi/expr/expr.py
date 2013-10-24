from .. import sql
from ..table_manager import TableManager
import csv


class GeneExpression(TableManager):

    file_field_seperator = '\t'
    file_has_header = True
    file_quoting = csv.QUOTE_NONE

    """
    An interface/super-class to all Gene/Protein-Expression data sets.
    """
    def __init__(self, filename, sql_connection, dataset_name):
        """
        Constructs a new GeneExpression class. Should not be used
        directly but only from overwriting subclasses.
        """
        self.filename = filename
        self.sql_conn = sql_connection
        self.name = dataset_name

    def import_raw_file(self):
        """
        Imports the expression file in raw format into the SQL database
        If the specific datasets need some kind of raw pre-processing,
        then this function MUST be overwritten.
        """
        table_name = self.name + "_raw"
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header, csv_quoting=self.file_quoting,
                       sql_conn=self.sql_conn)

    def filter(self):
        """
        Filters the expression data for genes with certain
        reliablilities/confidences.
        """
        pass

    def normalize_table(self):
        """
        Normalizes the table structure to the format
        [Gene],[Type] (,[Specific Type],...), [ExpressionValue]
        and especially gets rid of the [Gene] x [Type] format.
        """
        pass

    def id_mapping(self):
        """
        Maps the gene/protein identifier to a unifying identifier.
        """
        pass

    def classify(self):
        """
        Thresholds or classifies genes into either expressed ( = 1)
        or non-expressed (= 0).
        """
        raise NotImplementedError("The classify method has to be "
                                  "implemented by all subclasses")

