---
title: Fragment-Based Generative Design
---

# Fragment-Based Generative Design (3rd CACHE Challenge)

In the third **CACHE** challenge we focused on developing new antivirals against
SARS-CoV-2 by targeting the **Mac1 domain of NSP3**. Rather than designing whole
molecules, we used an **RNN-LSTM** approach for *fragment expansion*.

![Fragment-based generative drug design](/images/LSTM_Fragment_based_drug_design.png)

## Method

- Started from structural data of the NSP3 Mac1 domain and its known fragment
  binders, focusing on sub-sites of the ADP binding site.
- Applied **AI-guided** and knowledge-based fragment merging and expansion to
  generate novel molecules.
- Used the generated molecules as templates to find highly similar,
  **commercially available** compounds in the Enamine REAL database.

## Outcome

We proposed **91 potential compounds**, which were taken forward for testing —
an effort grounded in open science and community engagement.

!!! quote "Preprint"
    *Generative models for fragment-based design against NSP3 Mac1* —
    [ChemRxiv (2024)](https://chemrxiv.org/engage/chemrxiv/article-details/65c6e60b66c1381729521e8f).
