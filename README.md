### **Primary Outcome**
The product of this project is "Supplementary Dataset S1.csv" which is located in the root directory. This is the dataset used to benchmark the performance and behavior of the different base pairing annotation tools.
The pipeline for the project is contained within main.py, which finds the consensus annotations to construct the supplementary dataset.

### **Description of Data**
Within the root directory, there is a subdirectory titled "Data". This directory is home to all the data used in the project which includes PDB, PDBx, tool annotation, and results analysis files. All of the raw and processed data used throughout the workflow is already present within this directory (thus there is no need to download and re-execute the main.py script). The Data directory is currently consistent with version 3.95 of RNA 3D Motif Atlas, and holds 2.28 GB of data.

### **Instructions for Running**
**Reason for downloading and running.** Data consistent with version 3.95 of RNA 3D Motif Atlas is already present and processed within the "Data" subdirectory. If a more recent version of RNA 3D Motif Atlas is desired, then the entire project should be downloaded, and the main.py script run. The benchmark annotation produced from this will not be titled "Supplementary Dataset S1.csv" but be called "benchmark.csv" and its path is '\Data\ResultsAnalysis\BenchmarkDatasetInfo\benchmark.csv'.

**Usage Warnings.** Running the script to generate a new benchmark file will take between 1 and 2 days. This time constraint cannot be fully reduced by using a more powerful computer. This is due to the pipeline having a bottlenecks in Steps 4. 

In Step 4, the PDBs are sent to the ClaRNA web service to be annotated. This step alone is likely to take at least a day and consists of waiting for the ClaRNA server to respond. 

Step 5 prevents the pipeline from running automatically through a single sitting. In this step, the same PDBs are sent to DSSR for annotation, but this step must be done manually due to DSSR requiring a paid license. Thus a user will have to pause or stop the pipeline from running, get the DSSR tool annotations, place them in the directory '\Data\Raw\AnnotationTools\DSSR_Annotations\', and then rerun main.py while specifying to start at Step 5 (instead of Step 1) when prompted by the console. Currently the pipeline expects the DSSR annotations to be present and thus may break without them.

**Steps for running.** 
NOTE: If the pipeline was interrupted, run main.py again and enter an integer representing the pipeline step that had been running when the error occurred.
1) Run main.py.
2) Enter the integer "1".
2.5) Step 4 might prompt the user to retry. It is up to your discretion to retry or continue. Continuing merely means the size of the available dataset is smaller, but the pipeline will continue to function.
3) Step 5 prompts the console, don't continue the pipeline until DSSR annotations have been manually added in '\Data\Raw\AnnotationTools\DSSR_Annotations\'. Continue by entering 'y' or 'Y'.
3.5) If the pipeline was interrupted during Step 5, rerun main.py and enter the integer 5 when prompted and then the character 'y' or 'Y' when prompted.
4) Look for '\Data\ResultsAnalysis\BenchmarkDatasetInfo\benchmark.csv' which is the updated version of "Supplementary Dataset S1.csv".

### **Developer Notes**
- Pipeline needs to be adjusted to allow working with versions of R3DMA that aren't 3.95 (introduce a config.txt file).
- Pipeline needs to be streamlined such that the user no longer needs to be conscious of which pipeline step is being executed.
