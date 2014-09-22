
# load sql config and get connection
source("sql_config.R")
con <- get_sql_conn()

table <- "overlap_pairwise_ppi_ids"
ppi_id_data <- dbGetQuery(con, paste("SELECT * FROM ", table))

table <- "overlap_pairwise_ppi_edges"
ppi_edge_data <- dbGetQuery(con, paste("SELECT * FROM ", table))

for (row in 1:nrow(ppi_edge_data))
{
    ppi1 <- ppi_edge_data[row, 'ppi1']
    ppi2 <- ppi_edge_data[row, 'ppi2']
    id_overlap <- ppi_id_data[row, 'overlap_size']
    edge_overlap <- ppi_edge_data[row, 'overlap_size']
    ppi1_subset_size <- ppi_edge_data[row, 'ppi1_shared_ids_size']
    ppi2_subset_size <- ppi_edge_data[row, 'ppi2_shared_ids_size']

    n_max <- id_overlap
    m_max <- n_max * (n_max - 1) / 2

    # using the hypergeometric distrubution, get if actual overlap is
    # bigger or smaller than randomly expected overlap and the according p value
    expected_overlap <- ppi1_subset_size * 1.0 * ppi2_subset_size / m_max

    # put the variables into the notation of:
    # http://en.wikipedia.org/wiki/Hypergeometric_distribution
    N <- m_max
    # the order of these two don't matter (it's a symmetrical problem)
    K <- ppi1_subset_size
    n <- ppi2_subset_size
#    print('expected')
#    print(expected_overlap)
#    print('real')
#    print(edge_overlap)
#    print('phyper of mean')
#    print(phyper(expected_overlap, K, N - K, n))
#    print(phyper(expected_overlap, n, N - n, K))
    print(expected_overlap < edge_overlap)
    print(1 - phyper(edge_overlap, K, N - K, n))

}
