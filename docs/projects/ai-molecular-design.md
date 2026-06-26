---
title: AI-Driven Molecular Design
---

# AI-Driven Molecular Design

Designing a new molecule with AI is rarely a single model call. In practice it
is a **loop** — generate, score, filter, learn — that gradually steers a vast
chemical space toward compounds that are potent, selective and actually
synthesizable.

![AI-driven molecular design workflow](/images/AI_Driven_Molecular_Design_Workflow.png)

## The workflow

1. **Frame the target.** Pin down the biology — the target (and its structure
   where available), known actives and fragments, and the property profile we
   want: potency, selectivity, solubility, synthetic accessibility.
2. **Generate.** A recurrent neural network with LSTM cells (or a transformer)
   learns the grammar of valid molecules and samples novel candidates, usually
   via fragment expansion or scaffold decoration rather than whole-molecule
   design.
3. **Filter aggressively.** Validity & novelty → drug-likeness & toxicity alerts
   → synthetic feasibility. Cheap filters protect the expensive ones.
4. **Score with structure.** Docking and, for the top tier, free-energy
   estimates rank the survivors down to a short, prioritized list.
5. **Close the loop.** Molecules that get made and tested become new training
   data, sharpening the generator and the scoring models each cycle.

!!! note "Why a pipeline beats a single clever model"
    Generation is maybe 20% of the work — filters and scoring carry the quality,
    and a human chemist stays in the loop at every ranking step.

## Outcome

This exact loop has produced real, testable chemistry — including a selective
FLT-3 inhibitor and 91 candidates proposed in the 3rd CACHE challenge. See the
[FLT-3 inhibitor](flt3-inhibitor.md) and
[fragment-based design](fragment-based-design.md) projects for details.
