---
permalink: /
title: "About"
excerpt: "About me"
author_profile: true
redirect_from: 
  - /about/
  - /about.html
---

Working as Director of AI & Drug Design. I am interested in applying machine learning, especially deep learning, techniques for solving challenging problems in chemistry, drug discovery. 

**Website is under development and template inspired from Xin"**

My primary research interests involve:  
- Quantitative Structure-Activity Relationships (QSAR) modeling.
- Development of data-driven and interpretable molecular representations via deep learning.
- Molecular informatics approach for visualizing and navigating chemical space.

## Recent Works & Projects

Reliable *in silico* approaches to replace animal testing for the evaluation of potential acute toxic effects are highly demanded by regulatory agencies. The first research project of my Ph.D. was developing a **hierarchical QSAR modeling** method that integrates binary, multiclass and regression models for predicting three [acute oral systemic toxicity endpoints](https://ntp.niehs.nih.gov/whatwestudy/niceatm/test-method-evaluations/acute-systemic-tox/models/index.html?utm_source=direct&utm_medium=prod&utm_campaign=ntpgolinks&utm_term=tox-models) used by a variety of regulatory bodies for toxicity test. This method represents (1) a promising alternative prediction method to animal testing; (2) a more efficient and powerful ensemble method compared to the model consensus (predictions averaging). This project is a collaboration with [Dr. Nicole Kleinstreuer](https://www.niehs.nih.gov/research/atniehs/dntp/assoc/niceatm/staff/kleinstreuer/index.cfm) at the NIEHS and published on [_Chemical Research in Toxicology_](https://pubs.acs.org/doi/10.1021/acs.chemrestox.9b00259). [Github](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling)

![Generative Models for Fragment-based drug design]
#(/images/HQSAR.png)

The SARS-CoV-2 virus contains a host of nonstructural proteins (NSPs) that contribute to its structure and viral function. Among them is the nonstructural protein 3 (NSP3), which contains a macrodomain (Mac1) that interferes with antiviral adenosine diphosphate (ADP)-ribosylation signaling. Catalytic mutations in Mac1 render viruses nonpathogenic, making this enzyme a promising target for antiviral development. For this reason, the third CACHE challenge focused on identifying binders of the Mac1 domain of NSP3 for the development of novel antivirals against SARS-CoV-2. To this end, we used available structural data of the NSP3 Mac1 domain in complex with known fragment binders as starting points for ligand discovery; our efforts were primarily focused on sub-sites of the ADP binding site in the NSP3 macrodomain. Then, using **Artificial intelligence (AI)-guided** and knowledge-based fragment merging and expansion approaches, we generated novel molecules that would serve as templates to identify highly similar compounds in the Enamine REAL database that would be commercially available. Our design yielded a library of 12,800 molecules, which was docked with our program FITTED to a representative crystal structure of NSP3. We ranked the predicted binding poses based on docking score, followed by visual pose analysis of the best 200 compounds. We finally selected and proposed 150 compounds for testing, followed by further shortlisting to yield a final list of 107 molecules. 91 compounds were purchased from Enamine and are being tested at the Structural Genomics Consortium (SGC). Our approach and findings will further contribute to our open science efforts, and we aim to continue to engage the scientific community.([_Chemrxiv_ ]([https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00430-x](https://chemrxiv.org/engage/chemrxiv/article-details/65c6e60b66c1381729521e8f))). 
#![Overview to AI-driven molecular_design](/images/LSTM_Fragment_based_drug_design.png)


**SMILES Pair Encoding** (SPE) is a data-driven substructure tokenization algorithm for deep learning. SPE learns a vocabulary of high frequency SMILES substrings from ChEMBL and then tokenizes new SMILES into a sequence of tokens for deep learning models. SPE splits SMILES into human-readable and chemically explainable substrings and shows superior performances on both generative and predictive tasks compared to the atom-level tokenization. [[ChemRxiv](https://doi.org/10.26434/chemrxiv.12339368.v1)] [[Github](https://github.com/XinhaoLi74/SmilesPE)]

#![SPE Overview](/images/SPE.PNG)




