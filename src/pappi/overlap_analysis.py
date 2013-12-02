import itertools

from . import sql

PPIS_TO_ANALYZE = ['bossi', 'ccsb', 'havu', 'string', 'psicquic_all']
#PPIS_TO_ANALYZE = ['psicquic_dip', 'psicquic_i2d_imex', 'psicquic_innatedb_imex', 'psicquic_intact', 'psicquic_matrixdb', 'psicquic_mbinfo', 'psicquic_mint', 'psicquic_molcon', 'psicquic_mpidb', 'psicquic_uniprot']
EXPRS_TO_ANALYZE = ['emtab', 'gene_atlas', 'hpa', 'hpa_all', 'rnaseq_atlas']


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


def get_and_save_id_overlap(table1, table2, id_field, result_table, sql_conn,
                            verbose=False):
    # check parameters
    for t in [table1, table2]:
        if not sql.table_exists(t + '_ids', sql_conn):
            raise ValueError("table1_ids and table2_ids must be existing"
                             "SQL tables")

    # get sql cursor
    cur = sql_conn.cursor()

    # get sizes of the datasets (number of distinct IDs)
    id_counts = []
    for t in [table1, table2]:
        if verbose:
            print("Getting ID count from " + t)
        cur.execute('SELECT COUNT(DISTINCT ' + id_field + ') FROM '
                    + t + '_ids')
        c = cur.fetchone()[0]
        id_counts.append(c)

    # get overlap size
    if verbose:
        print("Getting size of ID overlap between " + table1 + " and "
              + table2)
    cur.execute('SELECT COUNT(*) FROM ('
                'SELECT DISTINCT ' + id_field + ' FROM ' + table1 + '_ids '
                'INTERSECT '
                'SELECT DISTINCT ' + id_field + ' FROM ' + table2 + '_ids)')
    overlap_size = cur.fetchone()[0]

    # insert results into the given table
    if verbose:
        print("Saving results from overlap analysis of "
              + table1 + " and " + table2)
    cur.execute('INSERT INTO ' + result_table + ' VALUES (?,?,?,?,?)',
                [table1, table2, id_counts[0], id_counts[1], overlap_size])

    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()


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

    # create id tables
    for p in ppis:
        if not sql.table_exists(p + '_ids', sql_conn):
            if verbose:
                print("Creating ID table for " + p)
            create_ppi_all_ids_table(p, sql_conn)
    for e in exprs:
        if not sql.table_exists(e + '_ids', sql_conn):
            if verbose:
                print("Creating ID table for " + e)
            create_expr_all_ids_table(e, sql_conn)

    for p in ppis:
        for e in exprs:
            get_and_save_id_overlap(e, p, 'Gene',
                                    'overlap_pairwise_expr_ppi',
                                    sql_conn, verbose=verbose)

#    # loop through all combina    # close cursor and commit changes to DB
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
            overlap_size = get_subset_ppi_size(p, e + '_ids', sql_conn)

            # save results to database
            cur.execute('INSERT INTO overlap_pairwise_expr_ppi_edges '
                        '(expr, ppi, expr_size, ppi_size, overlap_size) '
                        'VALUES (?,?,?,?,?)',
                        [e, p, expr_size, ppi_size, overlap_size])

    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()


def calc_pairwise_ppi_id_overlap(sql_conn, ppis=PPIS_TO_ANALYZE,
                                 verbose=False):
    """
    Calculates the pairwise overlap of the gene/protein IDs of all PPI
    networks given by `ppis`.

    @param sql_conn:    The SQL connection to be used.
    @param ppis:        A list of PPI names that correspond to tables in the
                        given database.
    """
    # first create the result table
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS overlap_pairwise_ppi_ids')
    cur.execute('CREATE TABLE overlap_pairwise_ppi_ids'
                '(`ppi1` varchar(16), `ppi2` varchar(16), `ppi1_size` int, '
                '`ppi2_size` int, `overlap_size` int)')
    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()

    # create all_ids tables if not yet existing
    for p in ppis:
        if not sql.table_exists(p + '_ids', sql_conn):
            if verbose:
                print("Creating ID table for " + p)
            create_ppi_all_ids_table(p, sql_conn)

    # get and save the id overlaps of all PPI combinations
    for x, y in itertools.combinations(ppis, 2):
        get_and_save_id_overlap(x, y, 'Gene', 'overlap_pairwise_ppi_ids',
                                sql_conn, verbose=verbose)


def get_subset_ppi_size(ppi, id_subset, sql_conn):
    """
    Returns the size of the subset of the given PPI where only
    edges are considered where both protein/gene IDs are within
    the given set of IDs `id_subset`.

    @param ppi:         The SQL table defining the PPI.
    @param id_subset:   The SQL table holding a set of IDs from which the
                        subset of the given PPI is constructed.
    @param sql_conn:    The SQL connection object.
    @returns:           The number of edges in the subset of the PPI.
    """
    # get SQL cursor
    cur = sql_conn.cursor()

    # execute SQL statement calculating the size of the subset PPI
    cur.execute('SELECT COUNT(*) FROM ('
                'SELECT a.Gene1, a.Gene2 '
                'FROM ' + ppi + ' AS a '
                'INNER JOIN ' + id_subset + ' AS b ON a.Gene1 = b.Gene '
                'INNER JOIN ' + id_subset + ' AS c ON a.Gene2 = c.Gene)')
    result = cur.fetchone()[0]
    cur.close()

    return result


def calc_pairwise_ppi_edge_overlap(sql_conn, ppis=PPIS_TO_ANALYZE,
                                   verbose=False):
    """
    Calculates the pairwise overlap of the edges of all PPI
    networks given by `ppis` (i.e. the shared edges between
    those PPI networks).

    @param sql_conn:    The SQL connection to be used.
    @param ppis:        A list of PPI names that correspond to tables in the
                        given database.
    """
    # first create the result table
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS overlap_pairwise_ppi_edges')
    cur.execute('CREATE TABLE overlap_pairwise_ppi_edges'
                '(`ppi1` varchar(16), `ppi2` varchar(16), '
                '`ppi1_size` int, `ppi2_size` int, '
                '`ppi1_shared_ids_size` int, `ppi2_shared_ids_size`, '
                '`overlap_size` int)')

    # create all_ids tables if not yet existing
    for p in ppis:
        if not sql.table_exists(p + '_ids', sql_conn):
            if verbose:
                print("Creating ID table for " + p)
            create_ppi_all_ids_table(p, sql_conn)

    for x, y in itertools.combinations(ppis, 2):
        # create table of shared ids
        sqlquery = ('SELECT Gene FROM ' + x + '_ids '
                    'INTERSECT '
                    'SELECT Gene FROM ' + y + '_ids')
        sql.new_table_from_query('overlap_ids_tmp', sqlquery, sql_conn)

        # get sizes of ppis and subset ppis
        ppi_sizes = []
        subset_ppi_sizes = []
        for p in [x, y]:
            # whole PPi size
            cur.execute('SELECT COUNT(*) FROM ' + p)
            ppi_size = cur.fetchone()[0]
            ppi_sizes.append(ppi_size)

            # subset (of shared IDs) sizes
            ppis_sub_size = get_subset_ppi_size(p, 'overlap_ids_tmp', sql_conn)
            subset_ppi_sizes.append(ppis_sub_size)

        cur.execute('DROP TABLE overlap_ids_tmp')

        # get overlap size
        cur.execute('SELECT COUNT(*) FROM ('
                    'SELECT Gene1, Gene2 FROM ' + x + ' '
                    'INTERSECT '
                    'SELECT Gene1, Gene2 FROM ' + y + ')')
        overlap_size = cur.fetchone()[0]

        # save results to db
        cur.execute('INSERT INTO overlap_pairwise_ppi_edges '
                    'VALUES (?,?,?,?,?,?,?)',
                    [x, y, ppi_sizes[0], ppi_sizes[1], subset_ppi_sizes[0],
                     subset_ppi_sizes[1], overlap_size])

    # close cursor and commit changes to DB
    cur.close()
    sql_conn.commit()
