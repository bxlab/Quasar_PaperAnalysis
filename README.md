# Quasar_PaperAnalysis
All scripts and information needed to perform all analyses and generate all figures from the HiFive QuASAR paper.

In order to rerun these analyses, simply download or clone this repository, then run the analysis.sh script. There are two lines that should be edited to match user preferences. Line 2 should be set to the number of CPUs available for MPI to use (1 will not use mpirun). Line 13 should be uncommented if the user wishes to run the complete analysis. Otherwise the included results files will be used and figures will be generated from those. Please note that the complete analysis may take up to 48 hours on 100 CPUs and will take significant (~250Gb) storage space. All of the Hi-C data may be found and downloaded from https://bx.bio.jhu.edu/data/quasar.
