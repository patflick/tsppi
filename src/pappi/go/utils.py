
import sqlite3


def name2id(name):
    """
    Returns the integer ID representation of a GO Term given as string.

    @param name The name of the GO-Term given as string
    @returns    The integer ID of the GO-Term
    """
    return int(name.split(":")[1])


# loads GO associations from SQL
def load_go_associations_sql(sql_conn, only_genes=None):
    """
    Loads GO associations from SQL using the given SQL connection.

    @param sql_conn     The SQL connection object.
    @param only_genes   Restrict the loaded genes to this set (default: None).
    @returns    A dict() object mapping gene names to sets of GO terms
    """
    cur = sql_conn.cursor()
    cur.execute("SELECT * FROM go_gene_assoc")
    assoc = {}
    for row in cur.fetchall():
        gene = row[0]
        go_term = row[1]
        go_term = name2id(go_term)
        if only_genes and gene not in only_genes:
            continue
        if type(go_term) is str:
            go_term = self.name2id(go_term)
        if gene in assoc:
            assoc[gene].add(go_term)
        else:
            assoc[gene] = set([go_term])
    return assoc
