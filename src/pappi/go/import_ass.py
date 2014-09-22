from .. import sql

GO_ASSOC_TABLE = "go_gene_assoc"


def import_go_association(filename, sql_conn, table=GO_ASSOC_TABLE):
    """
    Imports a GO gene association file.

    @param filename     The name of the file to import.
    @param sql_conn     The SQL connection to be used for import.
    @param table        The name of the table to be created from the CSV
    """
    sql.import_csv(filename, table, '\t', False, import_columns=[2, 4],
                   column_names=["Gene", "GOTerm"], sql_conn=sql_conn,
                   skip_rows=9)
