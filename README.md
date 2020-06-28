# DVA
**Computational methods for drug re-positioning identify potential anti-virals treatments against COVID-19**

Aanchal Mongia (IIIT Delhi), Emilie Chouzenoux (OPIS, Inria Saclay), Angshul Majumdar (IIIT Delhi)


TThis repository (DVA) contains:
*  DVA (Drug virus association database)
* Collection of matrix completion based computational techniques to predict anti-viral drug prediction for viruses


# Sources
* DrugBank: https://www.drugbank.ca/categories/DBCAT000066
* Antiviral drugs for viruses other than human immunodeficiency virus
* Approved antiviral drugs over the past 50 years
* Long-acting neuraminidase inhibitor laninamivir octanoate (cs-8958) versus oseltamivir as treatment for children with infuenza virus infection.
* Effectiveness of chloroquine and inflammatory cytokine response in patients with early persistent musculoskeletal pain and arthritis following chikungunya virus infection
* Heat shock protein 90 positively regulates chikungunya virus replication by stabilizing viral non-structural protein nsp2 during infection.
* Chikungunya virus: in vitro response to combination therapy with ribavirin and interferon alfa 2a.
* Structural basis for the inhibition of covid-19 virus main protease by carmofur, an antineoplastic drug
* Repurposing of the anti-malaria drug chloroquine for zika virus treatment and prophylaxis.
* Potential benefts of ibuprofen in the treatment of viral respiratory infections.
* ViPR: http://www.viprbrc.org/

The raw data can be found at: `./data_raw/database.xlsx`. The processed data has been created using the notebook `read_database.ipynb`. A schematic view of the DVA database curation and association prediction using it has been shown below.


![DVA-pipeline](./helper_functions/DVA.png)

The computational algorithms used to predict drug-virus association are available in: `helper_functions/alg_template`.
These are:
* Nuclear Norm Minimization based matrix completion [1]
* Matrix Facrorization based matrix completion [1]
* Deep matrix factorization [2]
* Graph regularized matrix factorization [3]
* Graph regularized matrix completion [4]
* Graph regularized binary matrix completion [5]

The results in the paper above can be reproduced by the following MATLAB scripts:

* `run.m`
* `./Experiments/novel_drugs_prediction.m`
* `./Experiments/coronavirus_pred.m`

# References
[1] Mongia, Aanchal, Debarka Sengupta, and Angshul Majumdar. "McImpute: Matrix completion based imputation for single cell RNA-seq data." Frontiers in genetics 10 (2019): 9.

[2] Mongia, Aanchal, Debarka Sengupta, and Angshul Majumdar. "deepMc: Deep Matrix Completion for Imputation of Single-Cell RNA-seq Data." Journal of Computational Biology (2019).

[3] Ezzat, Ali, et al. "Drug-target interaction prediction with graph regularized matrix factorization." IEEE/ACM transactions on computational biology and bioinformatics 14.3 (2016): 646-656.

[4] Mongia, Aanchal, and Angshul Majumdar. "Drug-target interaction prediction using Multi Graph Regularized Nuclear Norm Minimization." Plos one 15.1 (2020): e0226484.

[5] Mongia, Aanchal, Emilie Chouzenoux, and Angshul Majumdar. "Computational prediction of Drug-Disease association based on Graph-regularized one bit Matrix completion." bioRxiv (2020).
 
 

