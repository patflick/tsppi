from . import sql


PPIS_TO_ANALYZE = ['bossi', 'ccsb', 'havu']
EXPRS_TO_ANALYZE = ['emtab', 'gene_atlas', 'hpa', 'rnaseq_atlas']


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
