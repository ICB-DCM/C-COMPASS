# C-COMPASS

[![PyPI](https://badge.fury.io/py/ccompass.svg)](https://badge.fury.io/py/ccompass)
[![Documentation](https://readthedocs.org/projects/c-compass/badge/?version=latest)](https://c-compass.readthedocs.io)
[![DOI](https://zenodo.org/badge/916143374.svg)](https://doi.org/10.5281/zenodo.14712134)


**C-COMPASS** (Cellular COMPartmentclASSifier) is an open-source software tool designed to predict the spatial distribution of proteins across cellular compartments. It uses a neural network-based regression model to analyze multilocalization patterns and integrate protein abundance data while considering different biological conditions. C-COMPASS is designed to be accessible to users without extensive computational expertise, featuring an intuitive graphical user interface.

The data analyzed by C-COMPASS typically derives from proteomics fractionation samples that result in compartment-specific protein profiles. Our tool can be used to analyze datasets derived from various experimental techniques.

## Key Features

- **Protein Localization Prediction**: Use a neural network to predict the spatial distribution of proteins within cellular compartments.
- **Dynamic Compartment Composition Analysis**: Model changes in compartment composition based on protein abundance data under various conditions.
- **Comparison of Biological Conditions**: Compare different biological conditions to identify and quantify relocalization of proteins and re-organization of cellular compartments.
- **Multi-Omics Support**: Combine your proteomics experiment with different omics measurements such as lipidomics to bring your project to the spacial multi-omics level.
- **User-Friendly Interface**: No coding skills required; the tool features a simple GUI for conducting analysis.

## Documentation

Further documentation is available at https://c-compass.readthedocs.io/en/latest/.

## Installation

### Single-file executables

Single-file executables that don't require a Python installation are available
on the release page for Linux, Windows, and MacOS. Download the appropriate
file for your operating system and run it.

On Windows, make sure to install the Microsoft C and C++ (MSVC) runtime
libraries before ([further information](ttps://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170),
[direct download](https://aka.ms/vs/17/release/vc_redist.x64.exe)).

### Via pip

```bash
# install
pip install --pre ccompass

# launch the GUI
ccompass
# or alternatively: `python -m ccompass`
```

Note that C-COMPASS currently requires Python>=3.10, and due to its
`tensorflow` dependency Python<=3.12.

On Ubuntu linux, installing the `python3-tk` package is required:

```bash
sudo apt-get install python3-tk
```

### Troubleshooting

If you encounter any issues during installation, please refer to the
[troubleshooting guide](https://c-compass.readthedocs.io/en/latest/troubleshooting.html).

## Usage

See also https://c-compass.readthedocs.io/en/latest/usage.html.

### Graphical User Interface (GUI)

* The GUI will guide you through the process of loading and analyzing your
  proteomics dataset, including fractionation samples and Total Proteome
  samples.
* Follow the on-screen instructions to perform the analysis and configure
  settings only if required
* Standard parameters should fit for the majority of experiments.
  You **don't need to change the default settings!**

### Computation time (on a standard desktop computer)

* Preprocessing of Gradient and TotalProteome Data takes only up to a few
  minutes.
* Neural Network training for a dataset with three conditions and four
  replicates takes around 1-2h.
* Calculation of static predictions (per condition) takes a few minutes.
* Calculation of conditional comparisons (global comparison) takes up to
  30min. (for the above-mentioned dataset)
* Calculation of class-centric statistics and comparison takes up to 10 min.
  (for the above-mentioned dataset)

## Contributing

Contributions to C-COMPASS are welcome!

For further information, please refer to
[https://c-compass.readthedocs.io/en/latest/contributing.html](https://c-compass.readthedocs.io/en/latest/contributing.html).

## License

C-COMPASS is licensed under the BSD 3-Clause License.

## Citation

If you use C-COMPASS in your research, please cite the following publication:

```bibtex
@Article{HaasTra2024,
  author           = {Haas, Daniel Thomas and Trautmann, Eva-Maria and Mao, Xia and Gerl, Mathias J. and Klose, Christian and Cheng, Xiping and Hasenauer, Jan and Krahmer, Natalie},
  journal          = {bioRxiv},
  title            = {{C-COMPASS}: a neural network tool for multi-omic classification of cell compartments},
  year             = {2024},
  doi              = {10.1101/2024.08.05.606647},
  elocation-id     = {2024.08.05.606647},
  eprint           = {https://www.biorxiv.org/content/early/2024/08/08/2024.08.05.606647.full.pdf},
  publisher        = {Cold Spring Harbor Laboratory},
  url              = {https://www.biorxiv.org/content/early/2024/08/08/2024.08.05.606647},
}
```

## Contact

For any questions, contact `daniel.haas@helmholtz-munich.de` or post an
issue at https://github.com/ICB-DCM/C-COMPASS/issues/.
