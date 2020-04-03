# Breast cancer classifier based on silver-seq data

Please find here a tutorial together with some example data for a breast cancer recurrence classification practice (Zhou, Z., Wu, Q., Yan, Z., Zheng, H., Chen, C.J., Liu, Y., Qi, Z., Calandrelli, R., Chen, Z., Chien, S. and Su, H.I., 2019. Extracellular RNA in a single droplet of human serum reflects physiologic and disease states. Proceedings of the National Academy of Sciences, 116(38), pp.19200-19208.) 

The tutorial is the R program below and you'll be guided to construct a classifier for breast cancer recurrence status classification.

You also need some data input to run the R program. Please find the data here: https://drive.google.com/open?id=1i-8qJifMFZ-71xtoaGCYlDuRY1-YWywg. You'll be able to find the following files,
* tpm_96_nodup.txt, TPM values for the 96 breast cancer samples. Data for the samples are organized in columns (C1 to C96)
* readcounts_96_nodup.txt, read counts for the 96 cancer samples. Data for the samples are organized in columns (C1 to C96)
* patient_info.csv, meta info for all samples. The column recurStatus contains recurrence status of the corresponding donor, "R" for recurrence and "N" for non-recurrence.
* preselectedList, the 750 pre-selected breast cancer biomarker genes in Ensembl gene ID
