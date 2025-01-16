# Detection of mutational signatures in panel sequencing data
This project aims to address whether ERCC2 mutation-associated mutational signatures can be reliably identified from panel sequencing data. Using whole exome sequencing (WES) data from bladder cancer samples, mutational signatures will be extracted, and a simulated panel sequencing data will be generated. A computational model inspired by the SigMA framework will be developed to identify mutational signatures in WES and simulated data, followed by testing on real-world panel sequencing datasets. By enabling the identification of ERCC2 mutations and their associated mutational signatures in routine panel sequencing, this study has the potential to enhance clinical diagnostics and increase our understanding of ERCC2-driven mutagenesis in cancer. 

# Background
Cancer genomes carry a complex history of somatic mutations that are the result of various mutational processes that have operated during the cell’s lifetime. Depending on their origin, these mutational processes can be endogenous, such as errors in DNA replication, or exogenous like UV radiation and chemotherapy. Regardless, both types of mutational processes leave distinct patterns known as mutational signatures, that reflect the underlying biological mechanisms and exposures associated with cancer development and progression [1, 2]. An important part of cancer etiology has been identifying and characterizing mutational signatures, such as the comprehensive cataloging efforts of mutational signatures in the COSMIC database [3]. Most studies use whole genome sequencing (WGS) or whole exome sequencing (WES) for signature analysis, but panel sequencing is another alternative for detecting mutational signatures in a focused genomic context, especially considering its cost–effectiveness and prevalence in clinical diagnostics [4, 5]. 

The nucleotide excision repair (NER) pathway plays an important role in maintaining genomic integrity by repairing bulky DNA lesions induced by environmental or chemical agents, such as UV radiation and platinum-based chemotherapies. NER can be initiated in two different ways, depending on the lesion type: transcription coupled repair (TC-NER), activated by RNA polymerase stalling at lesions in transcribed regions, and global genome repair (GC-NER), which is able to operate and recognize lesions across the entire genome. Following lesion recognition, TC-NER and GG-NER converge on a common NER pathway that excises and replaces the damaged DNA strand in an error-free manner [6, 7]. Mutations affecting components of the NER pathway have great implications for cancer risk, progression, and treatment response [8]. 

Excision repair cross-complementation group 2 (ERCC2), a key NER gene that encodes the xeroderma pigmentosum group D (XPD) protein, a 5’-3’ helicase that is a component of the transcription factor II H (TFIIH) complex, is essential for unwinding DNA at damaged sites and facilitating repair. Point mutations in ERCC2 vary in location but are known to impair XPD-s helicase activity, stability, and interaction with other TFIIH components [9]. Somatic missense mutations in ERCC2 have been identified as putative drivers in approximately 10% of muscle-invasive bladder cancer cases, emphasizing its clinical significance. These mutations, predominantly localized to the helicase domain of XPD,  impair the functionality of the NER pathway to repair DNA lesions, leading to increased sensitivity to cisplatin-based chemotherapy, a key treatment for muscle-invasive bladder cancer [10, 11]. Additionally, ERCC2 mutations are associated with specific mutational signatures, such as Signature 2 and Signature 5, along with insertion-deletion (ID) signatures like ID8 [12, 13]. The ID8 signature, characterized by ≥ 5 bp deletions with short (≤ 2 bp) or no micro homologies and no repeats at deletion boundaries, has been determined to have the largest weight in a logistic regression-based ERCC2-mutant classifier trained on WES data [12]. The ability to detect the presence of these signatures offers insights into the biological consequences of ERCC2 deficiency. 

The increasing interest in applying mutational signatures in clinical diagnostics has led to the exploration of panel sequencing data as an alternative to WGS and WES. Panel sequencing targets only selected genomic regions, making it a more practical approach for routine diagnostics. Computational tools, such as SigMA (Signature Multivariate Analysis) which accurately determines the presence of  mutational signatures associated with homologous recombination (HR) deficiency from targeted gene panels, have demonstrated the potential for accurately identifying mutational signatures from panel sequencing data, providing a way to integrate mutation analysis into clinical workflows [14]. However, the ability of effectively detecting ERCC2 mutation-associated mutational signatures, such as ID8 and SBS5 from panel sequencing data remains unexplored. 

# Methods
## Data 
Whole-exome sequencing (WES) data from The Cancer Genome Atlas (TCGA) bladder cancer (BLCA) dataset which contains 414 tumor samples was used for the majority of the project. This dataset includes somatic mutation data, allowing us to explore mutational signatures associated with bladder cancer. Additionally, panel sequencing data (MSK-IMPACT505) from AACR Project GENIE, v16.1, was used to assess the performance of the developed computational model [15]. 

## Extraction of Mutational Signatures
To extract and identify the mutational signatures present in the dataset, the MutationalPatterns R package was used [16]. Aside from extracting the mutational signatures, the objective at this stage was to compare the samples with mutations in the ERCC2 helicase domain to those without ERCC2 mutations, focusing on the differences in mutational patterns. The R package facilitated the analysis of somatic mutations across different mutation types, including single base substitutions (SBS), double base substitutions (DBS), and insertions/deletions (indels). The workflow included the following steps:

### Construction of Mutational Matrices
Mutational count matrices which summarize the frequency of each mutation type across samples, were generated for SBS, DBS and indels. Furthermore the SBS count matrix was filtered to retain only samples with more than 50 mutations, ensuring sufficient data quality for accurate signature extraction. To account for biases introduced by genomics context, the SBS matric was normalized based on the trinucleotide frequencies within the exome
	 
### Visualization of Mutational Profiles
Mutational profiles for SBS, DBS and indels were plotted to identify the types of mutations that are most prevalent within the dataset. These offered a detailed view of the mutation spectrum and contextual biases, enabling comparison between ERCC2-mutant and wild-type samples. 

### Signature Fitting 
To extract the mutational signatures, signature fitting was performed using the COSMIC mutational signatures, focusing on those that have been previously associated with bladder cancer, ensuring biological relevance and reducing the influence of irrelevant signatures. To visualize the signature composition within the samples, the relative signature contribution was plotted, both for ERCC2-mutated and wild-type samples. The fitting process was validated by comparing the original mutational profiles with the reconstructed profiles generated for SBS, DBS, and indels. These plots evaluate the accuracy of the signature fitting, illustrating how well the extracted signatures explain the observed mutation data.

## Developing a Computational Model for ERCC2-Associated Mutational Signature Identification
As part of this project, I developed a computational model similar to SigMA to identify mutational signatures associated with ERCC2 mutations. The primary goal of the model was to distinguish ERCC2-mutant from ERCC2-wildype samples using the signatures SBS5 and ID8, and including features such as cosine similarity, likelihood measures and signature exposures. 

### Data Processing and Signature Fitting 
The initial step involved processing the count matrices by removing samples that were associated with microsatellite instability (MSI; TCGA-XF-A8HG and TCGA-ZF-AA4W) and POLE mutations (TCGA-DK-A6AW), as identified in previous studies [17, 18]. This filtering step ensures that confounding mutational processes will not influence downstream analyses. Additionally, the SBS count matrix was filtered so that it only contains samples where the total number of mutations is 50 or higher.

Signature fitting using non-negative least squares (NNLS) approach was performed next. The first fitting round included all COSMIC mutational signatures associated with bladder cancer. Additionally, a second fitting round was performed, restricting the analysis only to signatures with non-zero contribution in most samples from the first round. This refinement improved the specificity of the signature identification, and the results from this step were used for downstream analysis. For both rounds of fitting, the mutational count matrix used was not normalized in terms of trinucleotide frequency.

### Hierarchical Clustering 
Hierarchical clustering was performed on all WES samples with the goal to group them according to their mutational signature compositions, as determined in the previous step. SBS and INDEL samples were clustered separately using average linkage. The optimal number of clusters was determined by evaluating within-cluster sum of squares (WCSS) and manual exploration of cluster numbers to ensure biological relevance. A centroid, representing the average mutational spectrum of all samples in that cluster, was calculated for each SBS and INDEL cluster. These centroids were further normalized, making them interpretable as probability distributions in likelihood calculations. Lastly, the ERCC2-mutant samples were identified and evaluated whether they cluster together or exhibit distinct clustering patterns compared to wild-type samples. 

### Simulation of Artificial Panel Sequencing Data
Artificial WES samples were simulated to generate 250 ERCC2-mutant and 250 ERCC2-wildtype samples. The total mutation counts for each artificial sample were randomly chosen from the mutation count distribution observed in ERCC2-mutant and wildtype samples in the WES data. Each simulated sample included both SBS and INDEL mutations. The simulated samples were then filtered so they only include the regions covered by the MSK-IMPACT505 panel, which targets 505 genes. Furthermore, the simulated SBS count matrix was filtered, so that only samples that have more than 5 mutations in total were kept. The distribution of the total mutation count per sample of the artificial panel samples was compared to that of the actual MSK-IMPACT505 panel sequencing data for both ERCC2-mutant and wildtype groups to validate the simulation process.

### Likelihood and Posterior Probability Calculations
Using the normalized cluster centroids as probability distributions, the likelihood for a mutational spectrum St  = (n1, n2, …, n96) to belong to the ith cluster. The average probability distribution of the ith cluster is represented as Di = (pi1, pi2, …, pi96). In the case of INDELS the number of mutation types is 83 instead of 96. The probability of St being generated by Di was defined as the product of probabilities of 96 (or 83) mutation types, where mutation probabilities are denoted as pk, and nk denotes the counts of each mutation type. Separate calculations were performed for SBS and INDEL samples, using the following equation:

![Equation](https://latex.codecogs.com/png.latex?P%28S_t%20%7C%20D_i%29%20%3D%20%5Cprod_k%20%28p_k%29%5E%7Bn_k%7D)

### Cosine Similarity 
To validate the association between mutational spectra and COSMIC signatures, cosine similarity between mutational spectra of the simulated samples and known signatures SBS5 and ID8 was calculated. 

### Signature Optimization with NNLS
Iterative NNLS decomposition was performed for each sample to identify the optimal combination of signatures that minimize reconstruction error. For WES and panel sequencing data the optimal number of signatures was between 2 and 5. The best solution was compared against alternative signature combinations using likelihood calculations to confirm its validity. For each sample, the likelihood of the optimal solution was compared to alternative solutions by calculating the likelihood using different combinations of signatures. For this calculation, the same likelihood formula as above was used, however instead of cluster probabilities (Di), the combined probability distribution (Ci) of each signature combination was calculated and used in place of Di. The combined probability distribution for every alternative solution was calculated using the following equation:

![Equation](https://latex.codecogs.com/png.latex?C_i%20%3D%20%5Ctext%7BexpA%7D%20%5Ccdot%20%5Ctext%7BsigA%7D%20%2B%20%5Ctext%7BexpB%7D%20%5Ccdot%20%5Ctext%7BsigB%7D%20%2B%20%5Ctext%7BexpC%7D%20%5Ccdot%20%5Ctext%7BsigC%7D)

The relative likelihood ratio was then calculated as the ratio of the likelihood of the best combination to that of the alternative solutions.

### Multivariate Analysis
A comprehensive feature set for each sample was constructed and was integrated into a gradient boosting model to derive a single predictive score for each sample and the relative importance of each feature. For each sample the following features were considered:

1. Cosine similarity with COSMIC signature SBS5
2. Cosine similarity with COSMIC signature ID8
3. Likelihood of the mutational spectrum belonging to SBS cluster (separate for all clusters)
4. Likelihood of the mutational spectrum belonging to INDEL cluster (separate for all clusters)
5. SBS5 exposure from the best NNLS decomposition
6. ID8 exposure from the best NNLS decomposition
7. Relative likelihood ratio for SBS
8. Relative likelihood ratio for INDEL
9. Total number of SNVs per sample
10. Total number of insertions per sample
11. Total number of deletions per sample

For the samples that did not have either mutation type (SBS or INDEL), the respective features were assigned a NA value. The weighting of the different features was handled by a Gradient Boosting Classifier (GBC), a machine learning technique based on decision trees. The GBC was implemented using the gbm3 package in R. A small bagging fraction (0.2) together with a  small learning rate (0.001-0.01)  and a large number of trees was used for more robust results. The best number of trees was determined using cross validation (CV) to minimize the overfitting. To further reduce overfitting, a 3-step training procedure was implemented to get the optimal classifier. A first model was trained (GBC-1) using a 5-fold CV on all the features generated from the artificial panel sequencing samples, to determine the best number of trees from the Bernoulli deviance. Next, a test and training set was selected through a stratified folding procedure on ERCC2-mutant vs ERCC2-wildtype and a second GBC (GBC-2) was trained on the train set, using the parameters determined by CV. The performance of GBC-2 ws tested on the test set. Testing the performance of GBC-2 on the test set avoids any artificial improvement in the performance due to overfitting, because the samples were not used in the training. Lastly, a final model, GBC-3 was trained using the whole dataset, without any test set, using the optimal setting determined by CV. As there is more data used for training GBC-3, it is expected to perform better. The performance of GBC-3 was evaluated based on specificity and accuracy measures of the predicted scores. Using the cutpoint R package, the optimal threshold for assigning the ERCC2-mutant vs ERCC2-wildtype labels was determined, as the value that maximizes the area under the receiver operator characteristic (ROC) curve using the whole feature dataset. 

## Detection of mutational signatures from real panel sequencing data using the developed computational model
The last step of the analysis was implementing the developed model on a real panel sequencing data. For this purpose, the MSK-IMPACT505 panel sequencing dataset, filtered only for bladder cancer samples was used. GBC-3 was implemented to predict the scores for each sample, and the previously determined threshold was used to classify the sample as either ERCC2-mutant or ERCC2-wildtype

# References
[1] L. B. Alexandrov et al., “Signatures of mutational processes in human cancer,” Nature, vol. 500, no. 7463, pp. 415–21, 2013, doi: https://doi.org/10.1038/nature12477.

[2] L. B. Alexandrov et al., “The repertoire of mutational signatures in human cancer,” Nature, vol. 578, no. 7793, pp. 94–101, Feb. 2020, doi: https://doi.org/10.1038/s41586-020-1943-3.

[3] J. G. Tate et al., “COSMIC: the Catalogue Of Somatic Mutations In Cancer,” Nucleic Acids Research, vol. 47, no. D1, pp. D941–D947, Oct. 2018, doi: https://doi.org/10.1093/nar/gky1015.

[4] S. Goodwin, J. D. McPherson, and W. R. McCombie, “Coming of age: ten years of next-generation sequencing technologies,” Nature reviews. Genetics, vol. 17, no. 6, pp. 333–51, 2016, doi: https://doi.org/10.1038/nrg.2016.49.

[5] R. Luthra, H. Chen, S. Roy-Chowdhuri, and R. Singh, “Next-Generation Sequencing in Clinical Molecular Diagnostics of Cancer: Advantages and Challenges,” Cancers, vol. 7, no. 4, pp. 2023–2036, Oct. 2015, doi: https://doi.org/10.3390/cancers7040874.

[6] J. A. Marteijn, H. Lans, W. Vermeulen, and J. H. J. Hoeijmakers, “Understanding nucleotide excision repair and its roles in cancer and ageing,” Nature Reviews Molecular Cell Biology, vol. 15, no. 7, pp. 465–481, Jun. 2014, doi: https://doi.org/10.1038/nrm3822.

[7] Nicolás Nieto Moreno, A. M. Olthof, and J. Q. Svejstrup, “Transcription-Coupled Nucleotide Excision Repair and the Transcriptional Response to UV-Induced DNA Damage,” Annual Review of Biochemistry, vol. 92, no. 1, pp. 81–113, Jun. 2023, doi: https://doi.org/10.1146/annurev-biochem-052621-091205.

[8] M. Jager et al., “Deficiency of nucleotide excision repair is associated with mutational signature observed in cancer,” Genome Research, vol. 29, no. 7, pp. 1067–1077, Jun. 2019, doi: https://doi.org/10.1101/gr.246223.118.

[9] J. A. Barbour et al., “ERCC2 mutations alter the genomic distribution pattern of somatic mutations and are independently prognostic in bladder cancer,” Cell Genomics, vol. 4, no. 8, 2024, doi: https://doi.org/10.1016/j.xgen.2024.100627. 

[10] V. Allen et al., “Somatic ERCC2 Mutations Correlate with Cisplatin Sensitivity in MuscleInvasive Urothelial Carcinoma,” Cancer Discov, vol. 4, no. 10, pp. 1140–1153, Sep. 2014, doi: https://doi.org/10.1158/21598290.CD140623.

[11] Q. Li et al., “ERCC2 Helicase Domain Mutations Confer Nucleotide Excision Repair Deficiency and Drive Cisplatin Sensitivity in Muscle-Invasive Bladder Cancer,” Clinical Cancer Research, vol. 25, no. 3, pp. 977–988, Feb. 2019, doi: https://doi.org/10.1158/1078-0432.CCR-18-1001.

[12] J. Börcsök et al., “Identification of a Synthetic Lethal Relationship between Nucleotide Excision Repair Deficiency and Irofulven Sensitivity in Urothelial Cancer,” Clinical Cancer Research: An Official Journal of the American Association for Cancer Research, vol. 27, no. 7, pp. 2011–2022, Apr. 2021, doi: https://doi.org/10.1158/1078-0432.CCR-20-3316.

[13] J. Kim et al., “Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors,” Nature Genetics, vol. 48, no. 6, pp. 600–606, Apr. 2016, doi: https://doi.org/10.1038/ng.3557.

[14] D. C. Gulhan, J. J.-K. Lee, G. E. M. Melloni, I. Cortés-Ciriano, and P. J. Park, “Detecting the mutational signature of homologous recombination deficiency in clinical samples,” Nature Genetics, vol. 51, no. 5, pp. 912–919, Apr. 2019, doi: https://doi.org/10.1038/s41588-019-0390-2.

[15] AACR Project GENIE Consortium, “AACR Project GENIE: Powering Precision Medicine through an International Consortium,” Cancer Discovery, vol. 7, no. 8, pp. 818–831, Jun. 2017, doi: https://doi.org/10.1158/2159-8290.cd-17-0151.
[16] F. Manders et al., “MutationalPatterns: the one stop shop for the analysis of mutational processes,” BMC Genomics, vol. 23, no. 1, Feb. 2022, doi: https://doi.org/10.1186/s12864-022-08357-3.

[17] R. Bonneville et al., “Landscape of Microsatellite Instability Across 39 Cancer Types,” JCO precision oncology, vol. 2017, 2017, doi: https://doi.org/10.1200/PO.17.00073.

[18] R. A. Gordon et al., “Comprehensive Molecular Characterization of MuscleInvasive Bladder Cancer,” Cell, vol. 171, no. 3, pp. 540-556.e25, 2017, doi: https://doi.org/10.1016/j.cell.2017.09.007.





