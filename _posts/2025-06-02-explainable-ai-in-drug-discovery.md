---
title: 'Explainable AI in Drug Discovery: Why "Why" Matters'
date: 2025-06-02
permalink: /posts/2025/06/explainable-ai-in-drug-discovery/
excerpt: "A model that predicts activity is useful. A model that tells a chemist which substructure drives that activity is something you can act on."
tags:
  - Explainable AI
  - Drug Design
  - Machine Learning
---

A predictive model that says *"this molecule is active"* is helpful. A model that says *"this molecule is active **because of this substructure**"* changes how a chemist works. In drug discovery, where every synthesized compound costs real time and money, the **explanation is often worth more than the prediction**.

## The trust problem

Deep models are powerful precisely because they learn representations we did not hand them. That same flexibility makes them opaque. When a black box flags a compound, a medicinal chemist faces a hard question: *do I spend a week making this?* Without a reason, the rational answer is often "no" — and the model's value evaporates.

Explainability is what converts a probability into a **design hypothesis**.

## Three levels of explanation

Not all "explainable AI" means the same thing. It helps to separate:

1. **Global** — what has the model learned overall? Which descriptors or fragments dominate across the dataset?
2. **Local** — why *this* prediction for *this* molecule? Attribution methods (SHAP, integrated gradients, attention maps) highlight the atoms and bonds that pushed the score.
3. **Counterfactual** — what minimal change would flip the outcome? "Replace this nitro group and predicted toxicity drops" is a directly actionable suggestion.

Local and counterfactual explanations are the ones chemists reach for, because they map onto the next molecule to make.

## Mapping explanations back onto chemistry

The trick that makes this useful is grounding the explanation in the **molecular graph**. Atom-level attributions can be painted directly onto the 2D structure, turning a vector of numbers into a picture a chemist reads in seconds:

> Green where the model sees activity, red where it sees liability — overlaid on the structure itself.

This is where techniques like HQSAR-style fragment contributions and graph-attention weights earn their place: they speak the language of the bench.

## A few hard-won lessons

- **Faithfulness beats prettiness.** An explanation that looks convincing but doesn't reflect the model is worse than none.
- **Validate explanations like predictions.** If the model says a fragment drives activity, test whether removing it actually changes the assay.
- **Explainability is a design tool, not a compliance checkbox.** Used well, it shortens the generate–test loop.

---

The goal isn't to make AI *look* trustworthy — it's to make it genuinely **steerable**. When a chemist can see and challenge the reasoning, the model stops being an oracle and becomes a collaborator. That, more than any single accuracy number, is what gets AI adopted at the bench.
