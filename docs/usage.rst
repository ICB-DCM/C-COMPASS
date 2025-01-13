III. Usage Guide
================

1. Graphical User Interface (GUI)  
-------------------------------------
   
   a. C-COMPASS allows you to save and load your sessions via the main toolbar.
   
   b. A session can be saved as a NumPy (.npy) file, which includes all datasets, marker lists, settings, analyses, trainings, and statistics. These will be fully restored upon loading.

2. Pre-Training  
-------------------

   a. **Data Import**  

      i. There are two tabs for data import: Fractionation and TotalProteome.

      ii. Fractionation data can be analyzed independently, but TotalProteome is required for final class-centric statistics.

      iii. Use the "Add file..." button to import datasets. Multiple datasets can be imported and will appear in the dropdown menu. To remove a dataset, select it from the dropdown and click "Remove."

      iv. The table will display all column names found in the selected dataset.

   b. **Sample Annotation**

      i. For Fractionation data: Assign the condition, replicate number, and fraction numbers by selecting the relevant column names and clicking the appropriate button.

      ii. For TotalProteome data: Follow the same steps as Fractionation data, using consistent condition names.

      iii. Set the identifier column (e.g., ProteinGroups) for both Fractionation and TotalProteome datasets using the "Set Identifier" button. Ensure compatibility between these columns.

      iv. For other columns, either remove them or mark them as "Keep." Data marked as "Keep" will not be used in the analysis but will be available for export.

      v. **IMPORTANT**: Ensure that the column matching the marker list's naming (usually the gene name column) is kept.

   c. **Pre-Processing**

      i. Once columns are annotated, click "Process Fract." or "Process TP" to import the data.

      ii. Fractionation and TotalProteome data can be processed independently.

   d. **Marker List Import**

      i. In the "Marker Selection" frame, load marker lists via the "Add..." button. Multiple marker lists can be imported, and individual lists can be removed using the "Remove" button.

      ii. Imported marker lists will be displayed in the box.

      iii. For each marker list, specify the key column (e.g., gene names) and the class column (e.g., compartment).

      iv. In the "Fract. Key" section, select the column from the fractionation dataset that contains the compatible key naming. If the identifier and key column are the same, select "[IDENTIFIER]."

   e. **Marker Check & Matching**

      i. Click "Manage..." to view all class annotations from the marker lists. Unselect any classes you do not want in the analysis or rename them.

      ii. Classes with different nomenclatures (e.g., "ER" vs. "Endoplasmic Reticulum") can be merged by giving them the same name.

      iii. Median profiles of marker proteins and Pearson correlation matrices can be displayed via the corresponding buttons. Export options for plots and tables are available.

      iv. Confirm your marker selection by clicking "Match!."

3. Training  
---------------

   a. Start the training process by clicking "Train C-COMPASS."

   b. Various network architectures will be trained and evaluated for optimal results. This process may take over an hour, depending on dataset size.

   c. Progress will be shown in the background console window.

   d. **Hint**: Save your session after training to avoid repeating the process.

   e. **Note**: Future versions will optimize training time while maintaining calculation accuracy.

4. Post-Training  
--------------------

   a. **Statistics**

      i. After training, create "Static Statistics" via "Predict Proteome" to generate quantitative classifications for each condition.

      ii. Predictions can be exported or imported for comparison across sessions, ensuring compatible identifiers.

      iii. Use the "Report" button to export results.

      iv. Create simple plots and export them, along with the corresponding data tables.

   b. **Conditional Comparison - Global Changes**

      i. "Calculate Global Changes" compares localization across conditions, providing relocalization results.

      ii. Results can be displayed and exported similarly to the statistics.

   c. **Conditional Comparison - Class-centric Changes**

      i. **CPA (Class-centric Protein Amount)**: The amount of protein within a compartment, normalized by total proteome data. This is a relative value that requires comparison across conditions.

      ii. **CFC (Class-centric Fold-Change)**: The fold change of proteins across conditions within a compartment, based on CPA values. Only proteins with valid fractionation and total proteome data for both conditions will have CFC values.

5. Spatial Lipidomics  
-------------------------

   a. C-COMPASS has been used for spatial lipidomics analysis, though no dedicated feature currently exists for multi-omics analysis.

   b. You can concatenate proteomics and lipidomics datasets into one file before importing into C-COMPASS. Lipids will be treated like proteins, and spatial information can be derived similarly.

   c. Future versions of C-COMPASS will include features specifically designed for lipidomics.

6. Parameters
-----------------

   a. All parameters are set to default values used in our publication. It is not recommended to change them unless you are familiar with the procedure and its impact on results.

   b. **Parameters - Fractionation**

      i. Parameters for analysis and visualization can be adjusted independently.

      ii. **Min. valid fractions**: Profiles with fewer valid values across fractions can be filtered out.

      iii. **Found in at least X Replicates**: Proteins found in fewer replicates than specified will be removed.

      iv. **Pre-scaling**: Options include MinMax scaling or Area scaling.

      v. **Exclude Proteins from Worst Correlated Replicate**: Removes the replicate with the lowest Pearson correlation.

      vi. **Post-scaling**: Same options as Pre-scaling, useful for median profiles.

      vii. **Remove Baseline Profiles**: Removes profiles with only 0 values after processing.

   c. **Parameters - TotalProteome**

      i. **Found in at least X**: Similar to Fractionation data, this filters proteins found in fewer replicates.

      ii. **Imputation**: Missing values can be replaced by 0 or other values.

   d. **Parameters - Marker Selection**

      i. Discrepancies across marker lists can be handled by excluding markers or taking the majority annotation.

   e. **Parameters - Spatial Prediction**

      i. **WARNING**: Changes here are not recommended!

      ii. Various upsampling, noise, and SVM filtering methods are available for marker prediction.

   f. **Other parameters** for network training and optimization can be configured, including dense layer activation, output activation, loss function, optimizers, and number of epochs.
