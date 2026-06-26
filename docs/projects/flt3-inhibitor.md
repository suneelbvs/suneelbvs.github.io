---
title: Selective FLT-3 Inhibitor
---

# De novo Design of a Selective FLT-3 (D835Y) Inhibitor

Targeting FMS-like tyrosine kinase 3 (FLT-3) is an effective strategy for acute
myeloid leukemia (AML). Using **AI-assisted hit optimization**, we discovered a
novel, highly selective compound with drug-like properties against the
FLT-3 (D835Y) mutant.

![AI-driven molecular design workflow](/images/AI_Driven_Molecular_Design_Workflow.png)

## Approach

- A recurrent neural network with **long short-term memory (LSTM)** cells was
  trained to generate candidates related to our in-house hit compound
  (PCW-1001).
- Roughly **10,416 hits** were generated over 20 epochs.
- Candidates were filtered using toxicity and synthetic-feasibility filters,
  then ranked by **docking and free-energy** scoring.
- The top compounds were selected for synthesis and screening.

## Result

Of the synthesized compounds, **PCW-A1001** proved highly selective for the
FLT-3 (D835Y) mutant:

| Target            | IC~50~   |
| ----------------- | -------- |
| FLT-3 (D835Y)     | 764 nM   |
| FLT-3 (wild type) | 2.54 µM  |

!!! quote "Published"
    *AI-assisted de novo design of a selective FLT-3 (D835Y) inhibitor* —
    [Frontiers in Molecular Biosciences (2022)](https://www.frontiersin.org/articles/10.3389/fmolb.2022.1072028/full).
