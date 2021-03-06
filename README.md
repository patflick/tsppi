Analysis of tissue-specific protein-protein interaction networks
================================================================

The code in this repository implements an analysis pipeline for tissue-specific
protein-protein interaction networks.

In this text we use the following abbreviations:

- **PPI**: protein-protein interaction network
- **tsPPI**: tissue-specific PPI

This is the Master thesis project of *Patrick Flick*. The Master thesis is
available
[here](https://github.com/r4d2/tsppi/blob/master/docs/master_thesis_tsppi.pdf?raw=true).
(See also [`docs`](docs)). The thesis is licensed under the **CC-BY 4.0**
license.

## Code organization

Different parts of the pipeline are accomplished by different programming
languages. Python is used mainly for data import and processing. Most data comes
in text format (tab separated or comma separated) from various sources. This
data is imported and processed into a unified format and saved into a **SQLite**
table.

Graph analysis is performed with the
[**NetworKit**](https://networkit.iti.kit.edu/) parallel graph analysis library
in order to speed up graph analysis tasks. To this end we implement our own
extensions and extended algorithms for analyzing multiple Subgraphs
simultaneously. This is useful, since all **tsPPIs** are represented as a set of
subgraphs or a common **PPI** *parent* network. We implement this `Subgraphs`
data structure and the extended algorithms using `C++` and make them accessible
inside our python based pipeline by extending **NetworKit's** cython API.

Final data analysis and visualization is done by a set of **R** scripts, making
heavy use of the **ggplot2** visualizations library.

The code is organized into different parts:

- [`src`](src) contains the python pipeline and other python scripts
- [`src/pappi`](src/pappi) contains all custom python modules used for the
  python pipeline
- [`analysis`](analysis) contains **R** scripts for analyzing data and
  generating figures and plots
- [`figs`](figs) contains generated figures which are saved by the **R** scripts
  in [`analysis`](analysis)
- [`data`](data) contains the raw data and scripts to automatically download it
- [`ppi_networkit`](ppi_networkit) contains the Subgraph algorithms and the
  custom cython module, both of which extend the
  [**NetworKit**](https://networkit.iti.kit.edu/) graph library and it's cython
  python module.
- [`docs`](docs) contains the master thesis document. This documents the
  algorithms and analysis implemented, and shows and interprets the results
  achieved.
- TODO: add submodule for fastSemSim

## Installation

### Dependencies

- `python3` (sqlite3)
- `g++` (version >= 4.7)
- `scons` (for NetworKit compilation)
- `cython`
- `sqlite3-dev`
- `mercurial`
- [`fastSemSim`](http://sourceforge.net/projects/fastsemsim) (only needed for
  running the BPScore benchmarks with
  [`bpscore_benchmark.py`](src/bpscore_benchmark.py) )


### Compiling

To compile the NetworKit interacting cython module (including the extended
algorithms for Subgraph datastructures), first install the before mentioned
dependencies.

Then you can simply run the [`build.sh`](build.sh) script in the main directory
of the repository:

```sh
./build.sh
```


### Data

Most of the datasets can can be automatically downloaded. Only a few need
to be manually procured. See the [`data/README.md`](data/README.md) for details.


## Running

### Analysis Pipeline

Make sure all needed data is available in the [`data`](data) folder before
proceeding.

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
available in the [`analysis`](analysis) folder. Run these scripts to get the
individual figures and graphs.

### Algorithm benchmarks and tests

#### BPScore

The different *BPScore* algorithms can be benchmarked using the
[`bpscore_benchmark.py`](src/bpscore_benchmark.py) script. To run this execute:

```sh
cd src
python3 bpscore_benchmark.py
```

#### Subgraph algoritms (ppi_networkit)

In order to execute the benchmarks and tests for the modified Subgraph
algorithms implemented in the `ppi_networkit` cython module,
run in the [`ppi_networkit`](ppi_networkit) folder:

For the tests:

```sh
./build/tests
```

And for the benchmarks:

```sh
./build/benchmark
```

## Copyright Notice / Licensing

The code (`src`, `ppi_networkit`, `analysis`) is licensed under the **MIT**
license (see [LICENSE](LICENSE)).

The Master Thesis and all figures (folder `docs` and `figs`) are licensed under
the **CC BY 4.0** (Creative Commons Attribution 4.0 International) license
(see [docs/LICENSE](docs/LICENSE)).
