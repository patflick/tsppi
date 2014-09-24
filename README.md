Analysis of tissue-specific protein-protein interaction networks
================================================================

The code in this repository implements an analysis pipeline for tissue-specific
protein-protein interaction networks.

In this text we use the following abbreviations:

- **PP**: protein-protein interaction network
- **tsPPI**: tissue-specific PPI

## Code organization

Different parts of the pipeline are accomplished by different programming
languages. Python is used mainly for data import and processing. Most data comes
in text format (tab separated or comma separated) from various sources. This
data is imported and processed into a unified format and saved into a **SQLite**
table.

Graph analysis is done with the **NetworKit** parallel graph analysis library in
order to speed up graph analysis tasks. To this end we implement our own
extensions and extended algorithms for analyzing multiple Subgraphs
simultaneously. This is useful, since all **tsPPIs** are represented as a set of
subgraphs or a common **PPI** *parent* network. We implement this `Subgraphs`
data structure and the extended algorithms using `C++` and make them accessible
inside our python based pipeline by extending **NetworKit's** cython API.

Final data analysis and visualization is done by a set of **R** scripts, making
heavy use of the **ggplot2** visualizations library.

The code is organized into different parts:

- `src` contains the python pipeline and other python scripts
- `src/pappi` contains all custom python modules used for the python pipeline
- `analysis` contains **R** scripts for analyzing data and generating figures
  and plots
- `figs` contains generated figures which are saved by the **R** scripts in
  `analysis`
- `data` contains the raw data and scripts to automatically download it
- `ppi_networkit` contains the Subgraph algorithms and the custom cython module,
  both of which extend the `NetworKit` graph library and it's cython python
  module.
- TODO: add submodule for fastSemSim

## Installation

### Dependencies

- `python3` (sqlite3)
- `g++` (version >= 4.7)
- TODO: fastSemSim
- `scons` (for NetworKit compilation)
- `cython`
- `sqlite3-dev`
- `mercurial`


### Compiling

To compile the NetworKit interacting cython module (including the extended
algorithms for Subgraph datastructures), first install the before mentioned
dependencies.

Then you can simply run the `build.sh` script in the main directory of the
repository:

```sh
./build.sh
```


### Data

Most of the datasets can can be automatically downloaded. Only a few need
to be manually procured. See the [`data/README.md`](data/README.md) for details.


## Running

### Analysis Pipeline

Make sure all needed data is available in the `data` folder before proceeding.

Running the python pipeline consists of running 3 steps:
(approximate run times are on a 4 core 3.4Ghz Intel CPU with 8 GiB RAM)

1. Data import and preprocessing (~10-20 minutes)
2. Running the graph analysis on all networks (~4-8 hours)
3. Running graph clustering/community-detection on all graphs (~6 hours)

To run this pipeline, execute the following:

```sh
cd src
python3 init_data.py
python3 graph_properties.py
python3 networkit_clustering.py
```

The final analysis and data visualization is implemented in R, all scripts are
available in the `analysis` folder. Run these scripts to get the individual
figures and graphs.

### Algorithm benchmarks and tests

#### BPScore

The different *BPScore* algorithms can be benchmarked using the
`bpscore_benchmark.py` script. To run this execute:

```sh
cd src
python3 bpscore_benchmark.py
```

#### Subgraph algoritms (ppi_networkit)

In order to execute the benchmarks and tests for the modified Subgraph
algorithms implemented in the `ppi_networkit` cython module,
run in the `ppi_networkit` folder:

For the tests:

```sh
./build/tests
```

And for the benchmarks:

```sh
./build/benchmark
```


