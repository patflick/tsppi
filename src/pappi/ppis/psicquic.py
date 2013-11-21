from . import ppi
from .. import sql
import re


class Psicquic(ppi.PPI):
    """
    Imports and filters the string-db PPI network.
    """

    def __init__(self, filename, sql_connection, name):
        """
        Creates a new PPI network, with the given filename and sql
        connection object.

        Parameter:
            - filename:         Filename of the ppi network file.
            - sql_connection:   SQL connection object for the database,
                                into which the PPI network is to be
                                imported.
            - name:             Name of the particular PSICQUIC service
                                (i.e. MINT, DIP, etc)
        """
        # copy references into the new instance
        ppi.PPI.__init__(self, filename, sql_connection, "psicquic_" + name,
                         "Gene1", "Gene2", "uniprot")

        # init custom fields
        self.has_header = False
        self.file_field_seperator = '\t'

    def import_raw_file(self):
        """
        Import the PSI-MITAB v 2.5 file, but only the first two columns
        ( the interacting proteins) and the reliablility column.
        """
        # function to extract a uniprot id
        def psicqic_get_uniprot_id(s):
            p = re.compile(r'uniprotkb:(\w+)')
            m = p.match(s)
            if m:
                return m.group(1)
            else:
                return None

        # custom row iterator wrapper for the PSI-MITAB 2.5 format
        def psicquic_row_iter(base_iter):
            for row in base_iter:
                # simply "parses" uniprot and confidence from PSI MITAB 2.5:
                # https://code.google.com/p/psicquic/wiki/MITAB25Format
                # TODO maybe look in alt names if none found here
                gene1 = psicqic_get_uniprot_id(row[0])
                gene2 = psicqic_get_uniprot_id(row[1])
                confidence = row[14]
                yield [gene1, gene2, confidence]
            return
        # import the psicquic ppi
        table_name = self.next_tmp_table("raw")
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header,
                       column_names=["Gene1", "Gene2", "Confidence"],
                       column_types=["varchar(16)"]*3,
                       row_iterator_wrapper=psicquic_row_iter,
                       csv_quoting=self.file_quoting, sql_conn=self.sql_conn)

    def filter(self):
        # don't filter yet
        # TODO (maybe) implement other SW design to handle differently filtered
        #              PPIs
        pass
