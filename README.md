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
  and oplots
- `figs` contains generated figures which are saved by the **R** scripts in
  `analysis`
- `data` contains the raw data and scripts to automatically download it
- TODO: add other repos (ppi_networkit, networkit and fastSemSim) and their
  dependencies

## Installation

### Dependencies

- python3 (sqlite3)
- NetworKit
- C++11 compiler
- fastSemSim
- 

### Data

(see data folder readme)

### Compiling

TODO (building python modules



## Running

executing python pipeline
executing benchmarks and tests of ppi_networkit
plot all results

