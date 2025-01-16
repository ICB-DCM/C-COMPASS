===========
Usage Guide
===========

0. Data Preparation
===================

* To analyze your spatial proteomics datasets, you need the proteomics report file(s) derived from your spectral search software, such as MaxQuant, Spectronaut, DIANN, or others. Your data must be reported as a **pivot report table**, meaning that your table includes one column per sample, as well as additional columns for further information. The necessary columns are:

  * One column per sample (fraction).
  * One column containing a **unique identifier** (e.g., protein groups, protein ID, etc.).
  * One column containing key names that match the key names in your marker list (usually gene names). Ensure these keys are compatible, including case sensitivity.

* Furthermore, you need a file containing your marker proteins. C-COMPASS provides prepared marker lists from previous publications, or you can use a custom export from a database relevant to your project. This file must include at least two columns:

  * A column containing key names matching those in your dataset (usually gene names, see above).
  * A column containing **class annotations** (for spatial proteomics experiments, this should represent the compartments where the marker proteins are located).

* An additional dataset containing the total proteomes of the fractionation samples (proteomes derived from whole cell/tissue lysate) can be provided for **class-centric analysis** of compartments. This file should contain:

  * One column per total proteome sample.
  * One column containing the **same unique identifier** as used in the fractionation samples (see above).


Additional Notes
----------------

* All input files must be **tab-delimited** (.tsv or .txt).
* If using an export file from **Perseus**, ensure that the file does not contain a second-layer header.
* Input datasets (for both fractionation and total proteome) can be stored in the same file or split across different files. If they are split, ensure that the **identifiers** are consistent.


1. Graphical User Interface (GUI)
=================================

C-COMPASS allows you to save and load your sessions via the main toolbar.

A session can be saved as a NumPy (``.npy``) file, which includes all datasets,
marker lists, settings, analyses, trainings, and statistics. These will be
fully restored upon loading.

2. Before training
==================

#. **Data Import**

   #. There are two tabs for data import: `Fractionation` and `TotalProteome`.

   #. Fractionation data can be analyzed independently, but TotalProteome is
      required for final class-centric statistics.

   #. Use the `Add file...` button to import datasets.
      Multiple datasets can be imported and will appear in the dropdown menu.
      To remove a dataset, select it from the dropdown and click `Remove.`

   #. The table will display all column names found in the selected dataset.

#. **Sample Annotation**

   #. For Fractionation data: Assign the condition, replicate number, and
      fraction numbers by selecting the relevant column names and clicking the
      appropriate button.

   #. For TotalProteome data: Follow the same steps as Fractionation data,
      using consistent condition names.

   #. Set the identifier column (e.g., `ProteinGroups`) for both Fractionation and
      TotalProteome datasets using the "Set Identifier" button.
      Ensure compatibility between these columns.

   #. For other columns, either remove them or mark them as `Keep.`
      Data marked as `Keep` will not be used in the analysis but will be
      available for export.

   #. **IMPORTANT**: Ensure that the column matching the marker list's naming
      (usually the gene name column) is kept.

#. **Pre-Processing**

   #. Once columns are annotated, click `Process Fract.` or `Process TP`
      to import the data.

   #. Fractionation and TotalProteome data can be processed independently.

#. **Marker List Import**

   #. In the `Marker Selection` frame, load marker lists via the `Add...`
      button. Multiple marker lists can be imported, and individual lists can
      be removed using the `Remove` button.

   #. Imported marker lists will be displayed in the box.

   #. For each marker list, specify the key column (e.g., gene names)
      and the class column (e.g., compartment).

   #. In the `Fract. Key` section, select the column from the fractionation dataset that contains the compatible key naming. If the identifier and key column are the same, select `[IDENTIFIER].`

#. **Marker Check & Matching**

   #. Click `Manage...` to view all class annotations from the marker lists.
      Unselect any classes you do not want in the analysis or rename them.

   #. Classes with different nomenclatures (e.g., ``ER`` vs. ``Endoplasmic Reticulum``) can be merged by giving them the same name.

   #. Median profiles of marker proteins and Pearson correlation matrices
      can be displayed via the corresponding buttons.
      Export options for plots and tables are available.

   #. Confirm your marker selection by clicking `Match!`.

3. Training
===========

#. Start the training process by clicking `Train C-COMPASS`.

#. Various network architectures will be trained and evaluated for optimal results. This process may take over an hour, depending on dataset size.

#. Progress will be shown in the background console window.

#. **Hint**: Save your session after training to avoid repeating the process.

#. **Note**: Future versions will optimize training time while maintaining calculation accuracy.

4. After training
=================

#. **Statistics**

   #. After training, create `Static Statistics` via `Predict Proteome`
      to generate quantitative classifications for each condition.

   #. Predictions can be exported or imported for comparison across sessions,
      ensuring compatible identifiers.

   #. Use the `Report` button to export results.

   #. Create simple plots and export them, along with the corresponding data tables.

#. **Conditional Comparison - Global Changes**

   #. `Calculate Global Changes` compares localization across conditions,
      providing relocalization results.

   #. Results can be displayed and exported similarly to the statistics.

#. **Conditional Comparison - Class-centric Changes**

   #. **CPA (Class-centric Protein Amount)**: The amount of protein within a compartment, normalized by total proteome data. This is a relative value that requires comparison across conditions.

   #. **CFC (Class-centric Fold-Change)**: The fold change of proteins across conditions within a compartment, based on CPA values. Only proteins with valid fractionation and total proteome data for both conditions will have CFC values.

5. Spatial Lipidomics
======================

#. C-COMPASS has been used for spatial lipidomics analysis, though no dedicated feature currently exists for multi-omics analysis.

#. You can concatenate proteomics and lipidomics datasets into one file before importing into C-COMPASS. Lipids will be treated like proteins, and spatial information can be derived similarly.

#. Future versions of C-COMPASS will include features specifically designed for lipidomics.

6. Parameters
=============

#. All parameters are set to default values used in our publication. It is not recommended to change them unless you are familiar with the procedure and its impact on results.

#. **Parameters - Fractionation**

   #. Parameters for analysis and visualization can be adjusted independently.

   #. **Min. valid fractions**: Profiles with fewer valid values across fractions can be filtered out.

   #. **Found in at least X Replicates**: Proteins found in fewer replicates than specified will be removed.

   #. **Pre-scaling**: Options include MinMax scaling or Area scaling.

   #. **Exclude Proteins from Worst Correlated Replicate**: Removes the replicate with the lowest Pearson correlation.

   #. **Post-scaling**: Same options as Pre-scaling, useful for median profiles.

   #. **Remove Baseline Profiles**: Removes profiles with only 0 values after processing.

#. **Parameters - TotalProteome**

   #. **Found in at least X**: Similar to Fractionation data, this filters proteins found in fewer replicates.

   #. **Imputation**: Missing values can be replaced by 0 or other values.

#. **Parameters - Marker Selection**

   #. Discrepancies across marker lists can be handled by excluding markers or taking the majority annotation.

#. **Parameters - Spatial Prediction**

   #. **WARNING**: Changes here are not recommended!

   #. Various upsampling, noise, and SVM filtering methods are available for marker prediction.

#. **Other parameters** for network training and optimization can be configured, including dense layer activation, output activation, loss function, optimizers, and number of epochs.
