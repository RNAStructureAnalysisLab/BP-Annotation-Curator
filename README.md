The product of this project is "Supplementary Dataset S1.csv".
The pipeline for the project is contained within main.py, which finds the consensus annotations to construct the supplementary dataset.

To use, merely run main.py and a prompt will appear in the console. If running for the first time, use the prompt to start at Step 1. The script will then execute the specified step (step 1) and all other remaining steps. When the script finishes, the benchmark file will be stored as \Data\ResultsAnalysis\BenchmarkDatasetInfo\benchmark.csv (this is the same as 'Supplementary Dataset S1.csv').

Running the entire pipeline is expected to take more than one day. This is primarily due to Step 4 and Step 5. Step 4 is dependent on the ClaRNA web service which can take more than a day annotating the PDB files. Step 5 requires receiving annotations from DSSR.

NOTE: Step 5, the DSSR annotations, must be done manually. It is not integrated within the pipeline. Thus if DSSR annotations are desired, cancel the script once Step 5 starts, add the DSSR annotations in the \Data\Raw\AnnotationTools\DSSR_Annotations directory. Then run main.py again, but start at Step 5 instead of Step 1 and tell the console to proceed when notified.
