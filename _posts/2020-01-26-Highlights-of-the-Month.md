---
title:  "Highlights of the Month: January 2020"
date: 2020-01-26
permalink: /posts/2020/01/Highlights-of-the-Month/
tags:
  - Highlights of the Month
---

Key words for this month: Functional Groups; Unsupervised/Self-supervised representation leaning; Generative Models; 

# 🎓 Research Papers 

### An algorithm to identify functional groups in organic molecules [[Journal of Cheminformatics]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z)
**Functional Groups (FGs)** are sets of connected atoms that determine the properties and reactivity of a molecule. Software systems to identify functional groups are mainly based on a predefined list of substructures. The author proposed an algorithm can automate extract the functional groups in a molecule without any pre-defined substructures. An analysis of the ChEMBL database was implemented and 3,080 functional groups were found. A following study of analyzing functional groups occurring in Natural Products can be found [here](https://pubs.acs.org/doi/10.1021/acs.jnatprod.8b01022).

### Unsupervised/Self-supervised representation learning for (1) [chemical](https://chemrxiv.org/articles/Inductive_Transfer_Learning_for_Molecular_Activity_Prediction_Next-Gen_QSAR_Models_with_MolPMoFiT/9978743/1) (2) [genomic](https://github.com/kheyer/Genomic-ULMFiT) and (3) [protein](https://arxiv.org/abs/1906.08230) data.

One of the main trends in machine learning recently is the rise of transfer learning in NLP. It refers to the idea of training a language model on a massive corpus (unlabeled) and then fine-tuning the trained language model to other specific tasks of interest. Transfer learning allows you to leverage the knowledge learned by the language models, which can give you a boost in performance and generalization while demanding much less labeled training data. Recently, researchers found this approach can also be applied to other sequence data such as molecules (SMILES), genome and protein.

### *De Novo* Molecular Design
(1) De novo generation of hit-like molecules from gene expression signatures using artificial intelligence [[Nature Communications]](https://www.nature.com/articles/s41467-019-13807-w) [[ChemRxiv]](https://chemrxiv.org/articles/De_Novo_Generation_of_Hit-like_Molecules_from_Gene_Expression_Signatures_Using_Artificial_Intelligence/7294388/1)

(2) Bidirectional Molecule Generation with Recurrent Neural Networks [[J. Chem. Inf. Model.]](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00943) [[Github]](https://github.com/ETHmodlab/BIMODAL)

(3) DeepScaffold: A Comprehensive Tool for Scaffold-Based De Novo Drug Discovery Using Deep Learning [[J. Chem. Inf. Model.]](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00727)  [[arXiv]](https://arxiv.org/abs/1908.07209) [[Github]](github.com/deep-scaffold)

### New Datasets
[The Natural Products Atlas: An Open Access Knowledge Base for Microbial Natural Products Discovery](https://pubs.acs.org/doi/10.1021/acscentsci.9b00806). A collection of **24,594 microbial natural products**  (The Natural Products Atlas, www.npatlas.org) contains referenced data for structure, compound names, source organisms, isolation references, total syntheses, and instances of structural reassignment. 

[Discovering the anticancer potential of non-oncology drugs by systematic viability profiling](https://www.nature.com/articles/s43018-019-0018-6). Screening of the entire collection of mostly non-cancer drugs for their anti-cancer capabilities.

### Rethinking drug design in the artificial intelligence era [[Nature Communications]](https://www.nature.com/articles/s41573-019-0050-3)
*Can artificial intelligence help us design better small-molecule drug candidates faster?* In this perspective, the authors discussed *five* **'grand challenges'** that need to be addressed in order to make the drug design with AI successful in the long run: 
- Obtaining appropriate datasets
    - understanding of the technical error and biological variability associated with the underlying data.
    - accidental misreporting of data
    - missing values
    - highly imbalanced data.
    - The questions of whether a compound is ‘active against a target’ or ‘toxic’ are much more complex and labeled with much greater difficulty and nuances
- Generating new hypotheses
    - what to make next? Mainly driven by human creativity and 'chemical intuition'.
- Optimizing in a multi- objective manner
- Reducing cycle times
- **Changing the research culture and creating an appropriate mindset.**

In addition to the five 'grand challenges', the authors also discussed some other areas where AI might be relevant to drug discovery:
- Data curation and the identification of potential mistakes in data reporting.
- How to represent molecules and proteins?
- low-data situations
- Identify areas in which AI can augment and support (rather than replace) chemists and drug designers to make their processes more productive
- Design AI systems in a way that allows them to observe and learn from human behaviour in feedback cycles that are deemed beneficial for both sides. 


# 📃 Articles and Blog Posts 

### How to (Not) Get a Job in Science

In this [blog post](https://practicalcheminformatics.blogspot.com/2020/01/how-to-not-get-job-in-science.html), **Patrick Walters** provides some really useful suggestions:
> - Be realistic about your experience.
> - The way that you approach a problem is far more critical than the specific problems you’ve worked on.
> - list your publications and provide hyperlinks: papers, preprints, posters, slides.
> - Your **cover letter and answers to website questions** are very important.  
> - Netrwork! Netrwork! Netrwork! **It isn’t easy, but it’s doable.**
    - Go to conferences and present a poster.
    - Get in torch with people who are doing work that you find interesting.
    - Publish
    - Make something useful. 
    - Use social media
    
### NLP Stories
[NLP Year in Review — 2019](https://medium.com/dair-ai/nlp-year-in-review-2019-fb8d523bcb19)

[2019: The Year of BERT](https://towardsdatascience.com/2019-the-year-of-bert-354e8106f7ba)

An amazing series of visualizations and illustrations by **Jay Alammar**
- [The Illustrated Transformer](https://jalammar.github.io/illustrated-transformer/)
- [The Illustrated BERT, ELMo, and co. (How NLP Cracked Transfer Learning)](https://jalammar.github.io/illustrated-bert/)

### Predicting Molecular Properties [[Kaggle Competition]](https://www.kaggle.com/c/champs-scalar-coupling/overview)
This is a competition (finished) that aims to predict magnetic interactions between atoms in a molecule. One interesting thing about this competition is all the top 3 teams used [Transformer](https://arxiv.org/abs/1706.03762) models in their solutions. ([#1 solution](https://www.kaggle.com/c/champs-scalar-coupling/discussion/106575), [#2 solution](https://www.kaggle.com/c/champs-scalar-coupling/discussion/106468), [#3 solution](https://www.kaggle.com/c/champs-scalar-coupling/discussion/106572)). More solutions and discussion can be found [here](https://www.kaggle.com/c/champs-scalar-coupling/discussion)

# ✨ Notable Mentions 

### Some cheminformatics/computational chemistry blogs I found useful and inspiring:
   - [RDKit Blog](https://rdkit.blogspot.com/) by **Greg Landrum** [Twitter](https://twitter.com/dr_greg_landrum?lang=en).
   - [IS LIFE WORTH LIVING?](https://iwatobipen.wordpress.com/) by **pen** [Twitter](https://twitter.com/iwatobipen).
   - [Practical Cheminformatics](http://practicalcheminformatics.blogspot.com/search?updated-max=2019-11-01T18:09:00-07:00&max-results=7) by **Patrick Walters** [Twitter](https://twitter.com/wpwalters).

### Conference
[AI Powered Drug Discovery and Manufacturing CONFERENCE 2020](https://www.aidm.mit.edu/) February 27 - 28, 2020, MIT, Cambridge, MA

