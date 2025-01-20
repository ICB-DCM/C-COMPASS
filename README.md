# C-COMPASS

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
pip install ccompass
# launch the GUI
ccompass    # or alternatively: `python -m ccompass`
```

Note that C-COMPASS currently requires Python>=3.10, and due to its
`tensorflow` dependency Python<=3.12.

On Ubuntu linux, installing the `python3-tk` package is required:

```bash
sudo apt-get install python3-tk
```

## Usage

See also https://c-compass.readthedocs.io/en/latest/usage.html.

### Graphical User Interface (GUI)

* The GUI will guide you through the process of loading and analyzing your proteomics dataset, including fractionation samples and Total Proteome samples.
*Follow the on-screen instructions to perform the analysis and configure settings only if required
* Standard parameters should fit for the majority of experiments.
  You **don't need to change the default settings!**

### Command-Line Usage (Optional)
You can also run the software via the command line:
> python CCMPS.py

### Computation time (on a standard desktop computer)
- Preprocessing of Gradient and TotalProteome Data takes only up to a few minutes.
- Neural Network training for a dataset with three conditions and four replicates needs around 1-2h.
- Calculation of static predictions (per condition) takes a few minutes.
- Calculation of conditional comparisons (global comparison) takes up to 30 min. (for the above-mentioned dataset)
- Calculation of class-centric statistics and comparison takes up to 10 min. (for the above-mentioned dataset)

### Status:
- The appearance of the GUI will be improved in the near future. Progress bars will be included, as well as some help sections.
- Computation time will be optimized in the near future.
- Principal analysis steps and calculations will be kept as they are in version 1.0 unless changes are suggested by the reviewers.

### Contributing

Contributions to C-COMPASS are welcome! To contribute:

1. **Fork the repository** on GitHub.
2. **Create a new branch** for your changes.
3. **Commit your changes**.
4. **Submit a pull request**.

### License

C-COMPASS is licensed under the BSD 3-Clause License.

### Trouble-Shooting

* **SmartScreen Warning**: If Windows blocks the application via SmartScreen, this is due to the software being unsigned. Please consult your IT department to bypass this restriction if necessary.
* **Long Path Issues on Windows**: If your system encounters long path errors, you can activate them in your registry under 'HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem' by setting the value for **LongPathsEnabled* from 0 to 1.

### Contact

For any questions, contact `daniel.haas@helmholtz-munich.de` or post an
issue at https://github.com/ICB-DCM/C-COMPASS/issues/.
