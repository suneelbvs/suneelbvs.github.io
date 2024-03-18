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
**Project : Potency is All You need: Discovery of HDAC8 inhibitors**
Reliable *in silico* approaches to replace animal testing for the evaluation of potential acute toxic effects are highly demanded by regulatory agencies. The first research project of my Ph.D. was developing a **hierarchical QSAR modeling** method that integrates binary, multiclass and regression models for predicting three [acute oral systemic toxicity endpoints](https://ntp.niehs.nih.gov/whatwestudy/niceatm/test-method-evaluations/acute-systemic-tox/models/index.html?utm_source=direct&utm_medium=prod&utm_campaign=ntpgolinks&utm_term=tox-models) used by a variety of regulatory bodies for toxicity test. This method represents (1) a promising alternative prediction method to animal testing; (2) a more efficient and powerful ensemble method compared to the model consensus (predictions averaging). This project is a collaboration with [Dr. Nicole Kleinstreuer](https://www.niehs.nih.gov/research/atniehs/dntp/assoc/niceatm/staff/kleinstreuer/index.cfm) at the NIEHS and published on [_Chemical Research in Toxicology_](https://pubs.acs.org/doi/10.1021/acs.chemrestox.9b00259). [Github](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling)

**Project : Generative Models for Fragment-based drug design**


In the third CACHE challenge, we focused on developing new antivirals against SARS-CoV-2 by targeting the Mac1 domain of NSP3. We adopted an RNN-LSTM driven approach for molecular design, utilizing trained models for fragment expansion rather than designing whole molecules. This was based on structural data of the NSP3 Mac1 domain and its known fragment binders, particularly targeting the ADP binding site's sub-sites. hen, using **Artificial intelligence (AI)-guided** and knowledge-based fragment merging and expansion approaches, we generated novel molecules that would serve as templates to identify highly similar compounds in the Enamine REAL database that would be commercially available. We proposed 91 potential compounds, which are currently being tested. This work showcases our commitment to open science and engaging the scientific community in our research.([_Chemrxiv_ ]([https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00430-x](https://chemrxiv.org/engage/chemrxiv/article-details/65c6e60b66c1381729521e8f))). 
![Overview to AI-driven molecular_design](/images/LSTM_Fragment_based_drug_design.png)

**Project : Generative Models for Fragment-based drug design**

Treating acute myeloid leukemia (AML) by targeting FMS-like tyrosine kinase 3 (FLT-3) is considered an effective treatment strategy. By using AI-assisted hit optimization, we discovered a novel and highly selective compound with desired drug-like properties with which to target the FLT-3 (D835Y) mutant. In the current study, we applied an **AI-assisted de novo design approach** to identify a novel inhibitor of FLT-3 (D835Y). A recurrent neural network containing long short-term memory cells (LSTM) was implemented to generate potential candidates related to our in-house hit compound (PCW-1001). Approximately 10,416 hits were generated from 20 epochs, and the generated hits were further filtered using various toxicity and synthetic feasibility filters. Based on the docking and free energy ranking, the top compound was selected for synthesis and screening. Of these three compounds, PCW-A1001 proved to be highly selective for the FLT-3 (D835Y) mutant, with an IC50 of 764 nM, whereas the IC50 of FLT-3 WT was 2.54 μM. [[Frontiers in Molecular Biosciences]([https://doi.org/10.26434/chemrxiv.12339368.v1](https://www.frontiersin.org/articles/10.3389/fmolb.2022.1072028/full))] 

![PCW001](/images/AI_Driven_Molecular_Design_Workflow.png)




