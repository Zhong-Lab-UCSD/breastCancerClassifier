# Breast cancer classifier based on silver-seq data

This is a data platform and tutorial on breast cancer recurrence classification in paper of Zhou et al.(Zhou, Z., Wu, Q., Yan, Z., Zheng, H., Chen, C.J., Liu, Y., Qi, Z., Calandrelli, R., Chen, Z., Chien, S. and Su, H.I., 2019. Extracellular RNA in a single droplet of human serum reflects physiologic and disease states. Proceedings of the National Academy of Sciences, 116(38), pp.19200-19208.) The tutorial builds an example classifier of breast cancer recurrence status in R.

The required data is at https://drive.google.com/open?id=1i-8qJifMFZ-71xtoaGCYlDuRY1-YWywg, which including the following items,
* tpm_96_nodup.txt, TPM values for 96 breast cancer samples, in which the columns are in the order of C1 to C96
* readcounts_96_nodup.txt, read counts for 96 cancer samples, in which the columns are in the order of C1 to C96
* patient_info.csv, contains meta info for all samples, in which column recurStatus represents the recurrence status of each sample, with R denotes recurrence and N non-recurrence.
* preselectedList, 750 breast cancer biomarker genes shown in Ensembl gene ID

