/**
 * Benchmarks for subgraph algorithms
 *  - 3 methods for clustering coefficients
 *  - 2 methods for betweenness
 *  - 2 methods for PLP clustering
 */

// STL includes
#include <iostream>
#include <vector>

// own includes
#include <SQLiteIO.h>
#include <PpiGraph.h>
#include <cputimer.h>
#include <ClusterTools.h>
#include <subgraph_algos.h>
#include "ppi_utils.h"

// include auxiliary networKit
#include <auxiliary/Log.h>


//#define LOG(msg) std::cerr << msg << std::endl;
#define LOG(msg)

/*********************************************************************
 *                       Benchmark Betweenness                       *
 *********************************************************************/


void benchmark_ts_betweenness(tsppi::TsPpiGraph& tsppi)
{
    CPUTimer timer;

    double naive_time, fast_time;

    LOG("Start Benchmark: Naive");
    timer.start();
    std::vector<std::vector<double> > naive_ts_bw = tsppi::algo::subgraph_betweenness(tsppi.subgraphs);
    timer.stop();
    naive_time = timer.getTime();
    LOG("Time for Naive: " << timer.getTime() << " s");

    LOG("Start Benchmark: Fast");
    timer.start();
    std::vector<std::vector<double> > fast_ts_bw = tsppi::algo::subgraph_betweenness_fast(tsppi.subgraphs);
    timer.stop();
    fast_time = timer.getTime();
    LOG("Time for Fast: " << timer.getTime() << " s");

    std::cout << naive_time << ";" << fast_time;
}

void benchmark_all_bw(tsppi::SQLiteIO& sqlio)
{
    foreach_ppi_expr([&](std::string p, std::string e){
        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(p, e);
        std::cout << p << ";" << e << ";";
        benchmark_ts_betweenness(ts_graph);
        std::cout << std::endl;
    });
}


/*********************************************************************
 *                 Benchmark clustering coefficients                 *
 *********************************************************************/

void benchmark_ts_clustercoeff(tsppi::TsPpiGraph& tsppi)
{
    CPUTimer timer;

    double naive_time, fast_time, faster_time;

    LOG("Start Benchmark: Naive");
    timer.start();
    std::vector<std::vector<double> > naive_ts_bw = tsppi::algo::subgraph_cc(tsppi.subgraphs);
    timer.stop();
    naive_time = timer.getTime();
    LOG("Time for Naive: " << timer.getTime() << " s");

    LOG("Start Benchmark: Fast");
    timer.start();
    std::vector<std::vector<double> > fast_ts_bw = tsppi::algo::subgraph_cc_neighbor_comb(tsppi.subgraphs);
    timer.stop();
    fast_time = timer.getTime();
    LOG("Time for Fast: " << timer.getTime() << " s");

    LOG("Start Benchmark: Fast");
    timer.start();
    std::vector<std::vector<double> > faster_ts_bw = tsppi::algo::subgraph_cc_neighbor_comb_vec(tsppi.subgraphs);
    timer.stop();
    faster_time = timer.getTime();
    LOG("Time for Fast: " << timer.getTime() << " s");

    std::cout << naive_time << ";" << fast_time << ";" << faster_time;
}

void benchmark_all_cc(tsppi::SQLiteIO& sqlio)
{
    std::cout << "ppi;expr;N;M;T;naive;neighcomb;tsvector" << std::endl;
    foreach_ppi_expr([&](std::string p, std::string e){
        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(p, e);
        std::cout << p << ";" << e << ";" << ts_graph.subgraphs.numberOfNodes() << ";" << ts_graph.subgraphs.numberOfEdges() << ";" <<  ts_graph.numberOfTissues() << ";";
        benchmark_ts_clustercoeff(ts_graph);
        std::cout << std::endl;
    });
}


/*********************************************************************
 *                     Benchmark PLP clustering                      *
 *********************************************************************/

void benchmark_ts_PLP(tsppi::TsPpiGraph& tsppi, bool printHist=false)
{
    CPUTimer timer;

    double naive_time, fast_time;

    LOG("Start Benchmark: Naive");
    timer.start();
    std::vector<NetworKit::Partition > partitions = tsppi::algo::subgraph_PLP(tsppi.subgraphs);
    timer.stop();
    naive_time = timer.getTime();
    LOG("Time for Naive: " << timer.getTime() << " s");

    LOG("Start Benchmark: Fast");
    timer.start();
    std::vector<NetworKit::Partition > partitions_2 = tsppi::algo::subgraph_PLP_vec(tsppi.subgraphs);
    timer.stop();
    fast_time = timer.getTime();
    LOG("Time for Fast: " << timer.getTime() << " s");


    // print histogram
    if (printHist)
    {
        std::cout << "Histogram of cluster sizes for Naive:" << std::endl;
        clusterSizeHist(partitions);
        std::cout << "Histogram of cluster sizes for fast:" << std::endl;
        clusterSizeHist(partitions_2);
    }

    // print timings
    std::cout << naive_time << ";" << fast_time;
}

void benchmark_all_PLP(tsppi::SQLiteIO& sqlio)
{
    std::cout << "ppi;expr;N;M;T;naive;ts" << std::endl;
    foreach_ppi_expr([&](std::string p, std::string e){
        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(p, e);
        std::cout << p << ";" << e << ";" << ts_graph.subgraphs.numberOfNodes() << ";" << ts_graph.subgraphs.numberOfEdges() << ";" <<  ts_graph.numberOfTissues() << ";";
        benchmark_ts_PLP(ts_graph);
        std::cout << std::endl;
    });
}


/*********************************************************************
 *                        main() entry point                         *
 *********************************************************************/

void print_usage()
{
    std::cerr << "Usage: ./benchmark [OPTIONS] <sqlite_db_filename>" << std::endl;
    std::cerr << "  Options:" << std::endl;
    std::cerr << "      -a  The algorithm to test, can be: `cc`, `bw`, or `plp`" << std::endl;
    std::cerr << "      -p  (optional) The PPI to test this on. `-e` has to be set as well" << std::endl;
    std::cerr << "      -e  (optional) The Expression data to test this on. `-p` has to be set as well" << std::endl;
    std::cerr << "           If neigher -p or -e is set, the given algorithm is tested" << std::endl;
    std::cerr << "           on all PPIs and expression datasets." << std::endl;
}


int main(int argc, char *argv[])
{
    Aux::Log::setLogLevel("ERROR");

    --argc;
    ++argv;
    if (argc < 1){
        print_usage();
        exit(1);
    }

    std::string algo;
    std::string ppi = "";
    std::string expr = "";
    std::string sqlite_db_file = "";

    // parse optional arguments
    while(argc > 0){
        if (argv[0][0] == '-')
        {
            switch(argv[0][1])
            {
                case 'a':
                    // set algorithm to benchmark
                    algo = argv[1];
                    ++argv;
                    --argc;
                    break;
                case 'p':
                    // set PPI to load
                    ppi = argv[1];
                    ++argv;
                    --argc;
                    break;
                case 'e':
                    // set expr to load
                    expr = argv[1];
                    ++argv;
                    --argc;
                    break;
                default:
                    print_usage();
                    exit(1);
            }
            --argc;
            ++argv;
        }
        else
        {
            // this is the database file
            sqlite_db_file = argv[0];
            ++argv;
            --argc;
        }
    }

    if (sqlite_db_file == "" || algo == "")
    {
        print_usage();
        exit(1);
    }

    // open Sqlite database file
    tsppi::SQLiteIO sqlio(sqlite_db_file);

    if (ppi != "")
    {
        assert(expr != "");

        tsppi::TsPpiGraph ts_graph = sqlio.load_tsppi_graph(ppi, expr);
        if (algo == "bw")
        {
            benchmark_ts_betweenness(ts_graph);
        } else if (algo == "cc")
        {
            benchmark_ts_clustercoeff(ts_graph);
        } else if (algo == "plp")
        {
            benchmark_ts_PLP(ts_graph);
        } else {
            print_usage();
            exit(1);
        }
    }
    else
    {
        // benchmark all ppi and expr
        if (algo == "bw")
        {
            benchmark_all_bw(sqlio);
        } else if (algo == "cc")
        {
            benchmark_all_cc(sqlio);
        } else if (algo == "plp")
        {
            benchmark_all_PLP(sqlio);
        } else {
            print_usage();
            exit(1);
        }
    }

    return 0;
}
