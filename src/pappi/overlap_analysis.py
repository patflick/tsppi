from . import sql


PPIS_TO_ANALYZE = ['bossi', 'ccsb', 'havu', 'string', 'psicquic_all']
EXPRS_TO_ANALYZE = ['emtab', 'gene_atlas', 'hpa', 'rnaseq_atlas']


def create_ppi_all_ids_table(ppi, sql_conn):
    """
    Creates a table named `ppi`_ids that holds all the distinct IDs used
    in the ppi network.
    """
    sqlquery = ('SELECT Gene1 as Gene FROM ' + ppi + ' UNION '
                'SELECT Gene2 as Gene FROM ' + ppi)
    sql.new_table_from_query(ppi + '_ids', sqlquery, sql_conn)


def create_expr_all_ids_table(expr, sql_conn):
    """
    Creates a table named `expr`_ids that holds all the distinct IDs used
    in the expression data set.
    """
    sqlquery = ('SELECT DISTINCT Gene FROM ' + expr)
    sql.new_table_from_query(expr + '_ids', sqlquery, sql_conn)


def get_binary_fields(fields, src_field):
    """
    Returns a SQL style list of all PPIs with the src_ppi
    given set to 1 and all other set to 0.
    """
    return ', '.join(('1' if p == src_field else '0') + ' AS ' + p
                     for p in fields)


def calc_ppi_edge_overlap(sql_conn, ppis=PPIS_TO_ANALYZE):
    """
    Calculates the overlap of the edges of all combinations of PPIs
    in order to generate the data needed for a Venn-diagram.
    """
    union = ' UNION '.join('SELECT Gene1 || "-" || Gene2 as Interaction, '
                           + get_binary_fields(ppis, p) + ' '
                           'FROM ' + p for p in ppis)
    sql.new_table_from_query('ppi_edge_overlap_1', union, sql_conn)
    # create second table by grouping by the edge and summing the bit patterns
    mask_sum = ', '.join('SUM(' + p + ') AS ' + p for p in ppis)
    comb_int = ('SELECT Interaction, ' + mask_sum + ' FROM '
                'ppi_edge_overlap_1 GROUP BY Interaction')
    sql.new_table_from_query('ppi_edge_overlap_2', comb_int, sql_conn)
    # create final table by grouping by the bit pattern and
    # counting edges
    ppis_list = ', '.join(ppis)
    sqlquery = ('SELECT ' + ppis_list + ', COUNT(*) AS count '
                'FROM ppi_edge_overlap_2 GROUP BY ' + ppis_list)
    sql.new_table_from_query('ppi_edge_overlap', sqlquery, sql_conn)


def calc_ppi_id_overlap(sql_conn, ppis=PPIS_TO_ANALYZE):
    """
    Calculates the overlap of the gene identifiers of all combinations of PPIs
    in order to calculate the data needed for a Venn-diagram.
    """
    union = ' UNION '.join('SELECT DISTINCT Gene1 as Gene, '
                           + get_binary_fields(ppis, p) + ' '
                           'FROM ' + p + ' '
                           ' UNION '
                           'SELECT DISTINCT Gene2 as Gene, '
                           + get_binary_fields(ppis, p) + ' '
                           'FROM ' + p for p in ppis)
    sql.new_table_from_query('ppi_id_overlap_1', union, sql_conn)
    # create second table by grouping by the edge and summing the bit patterns
    mask_sum = ', '.join('SUM(' + p + ') AS ' + p for p in ppis)
    comb_int = ('SELECT Gene, ' + mask_sum + ' FROM '
                'ppi_id_overlap_1 GROUP BY Gene')
    sql.new_table_from_query('ppi_id_overlap_2', comb_int, sql_conn)
    # create final table by grouping by the bit pattern and
    # counting edges
    ppis_list = ', '.join(ppis)
    sqlquery = ('SELECT ' + ppis_list + ', COUNT(*) AS count '
                'FROM ppi_id_overlap_2 GROUP BY ' + ppis_list)
    sql.new_table_from_query('ppi_id_overlap', sqlquery, sql_conn)


def calc_expr_overlap(sql_conn, exprs=EXPRS_TO_ANALYZE):
    """
    Calculates the gene overlap of all combinations of expression data sets
    in order to generate the data needed for a Venn-diagram.
    """
    union = ' UNION '.join('SELECT Gene, '
                           + get_binary_fields(exprs, p) + ' '
                           'FROM ' + p for p in exprs)
    sql.new_table_from_query('expr_overlap_1', union, sql_conn)
    # create second table by grouping by the gene and summing the bit patterns
    mask_sum = ', '.join('SUM(' + p + ') AS ' + p for p in exprs)
    comb_int = ('SELECT Gene, ' + mask_sum + ' FROM '
                'expr_overlap_1 GROUP BY Gene')
    sql.new_table_from_query('expr_overlap_2', comb_int, sql_conn)
    # create final table by grouping by the bit pattern and
    # counting genes
    exprs_list = ', '.join(exprs)
    sqlquery = ('SELECT ' + exprs_list + ', COUNT(*) AS count '
                'FROM expr_overlap_2 GROUP BY ' + exprs_list)
    sql.new_table_from_query('expr_overlap', sqlquery, sql_conn)


def calc_pairwise_expr_ppi_id_overlap(sql_conn, exprs=EXPRS_TO_ANALYZE,
                                      ppis=PPIS_TO_ANALYZE, verbose=True):
    """
    Calculates the pairwise ID overlap of the expression data sets given
    by `exprs` and the PPI networks given by `ppis`.
    """
    # create table for overlap results
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS overlap_pairwise_expr_ppi')
    cur.execute('CREATE TABLE overlap_pairwise_expr_ppi (`expr` varchar(16), '
                '`ppi` varchar(16), `expr_size` int, `ppi_size` int, '
                '`overlap_size` int)')

    # loop through all combinations
    for p in ppis:
        # create list of ids (if doesn't exist yet)
        if not sql.table_exists(p + '_ids', sql_conn):
            if verbose:
                print("Creating ID table for " + p)
            create_ppi_all_ids_table(p, sql_conn)
        # get count of ppi ids
        if verbose:
            print("Getting ID count from " + p)
        cur.execute('SELECT COUNT(*) FROM ' + p + '_ids')
        ppi_size = cur.fetchone()[0]
        # loop through all expression data sets
        for e in exprs:
            # get expr size
            # TODO: this could be done outside of for p in ppis loop
            #       (more efficient)
            if verbose:
                print("Getting ID count from " + e)
            cur.execute('SELECT COUNT(DISTINCT Gene) FROM ' + e)
            expr_size = cur.fetchone()[0]

            # get overlap size
            if verbose:
                print("Getting overlap size between " + e + " and " + p)
            cur.execute('SELECT COUNT(*) FROM ('
                        'SELECT Gene FROM ' + p + '_ids INTERSECT '
                        'SELECT DISTINCT Gene FROM ' + e + ')')
            overlap_size = cur.fetchone()[0]

            # save results to database
            cur.execute('INSERT INTO overlap_pairwise_expr_ppi '
                        '(expr, ppi, expr_size, ppi_size, overlap_size) '
                        'VALUES (?,?,?,?,?)',
                        [e, p, expr_size, ppi_size, overlap_size])

    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()


def calc_pairwise_expr_ppi_edge_overlap(sql_conn, exprs=EXPRS_TO_ANALYZE,
                                        ppis=PPIS_TO_ANALYZE, verbose=True):
    """
    Calculates the pairwise overlap of the expression data sets given
    by `exprs` and the edges of the PPI networks given by `ppis`.

    Contrary to the function `calc_pairwise_expr_ppi_id_overlap`, this returns
    the count of edges of the PPIS (instead of count of distinct IDs)
    and the number of edges where _both_ IDs are in the expression data set
    as the overlap_size.
    """
    # create table for overlap results
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS overlap_pairwise_expr_ppi_edges')
    cur.execute('CREATE TABLE overlap_pairwise_expr_ppi_edges '
                '(`expr` varchar(16), `ppi` varchar(16), `expr_size` int, '
                '`ppi_size` int, `overlap_size` int)')

    # loop through all combinations
    for e in exprs:
        # create distinct id table
        if not sql.table_exists(e + '_ids', sql_conn):
            create_expr_all_ids_table(e, sql_conn)

        if verbose:
            print("Getting ID count from " + e)
        # get size of the expr data set
        cur.execute('SELECT COUNT(*) FROM ' + e + '_ids')
        expr_size = cur.fetchone()[0]

        for p in ppis:
            if verbose:
                print("Getting ID count from " + p)
            # get size of the ppi
            cur.execute('SELECT COUNT(*) FROM ' + p)
            ppi_size = cur.fetchone()[0]

            # inner join both ppi columns with the expr genes
            if verbose:
                print("Getting overlap size for " + p + " and " + e)
            cur.execute('SELECT COUNT(*) FROM ('
                        'SELECT a.Gene1, a.Gene2 '
                        'FROM ' + p + ' AS a '
                        'INNER JOIN ' + e + '_ids AS b ON a.Gene1 = b.Gene '
                        'INNER JOIN ' + e + '_ids AS c ON a.Gene2 = c.Gene)')
            overlap_size = cur.fetchone()[0]

            # save results to database
            cur.execute('INSERT INTO overlap_pairwise_expr_ppi_edges '
                        '(expr, ppi, expr_size, ppi_size, overlap_size) '
                        'VALUES (?,?,?,?,?)',
                        [e, p, expr_size, ppi_size, overlap_size])

    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()
