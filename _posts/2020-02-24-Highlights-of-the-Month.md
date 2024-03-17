---
title:  "Highlights of the Month: February 2020"
date: 2020-02-24
permalink: /posts/2020/02/Highlights-of-the-Month/
tags:
  - Highlights of the Month
---

**Key words** for this month: Reformer; TMAP; libmolgrid;

## Research Papers 🎓

### Cheminformatics
**[Visualization of very large high-dimensional data sets as minimum spanning trees](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0416-x)**. 
> [TMAP](http://tmap.gdb.tools/) is a new data visualization method designed for visualizing large, high-dimensional data sets containing chemical structures and associated properties while preserving both global and local features.

**Exploring chemical space using natural language processing methodologies for drug discovery** [Drug Discov Today](https://www.sciencedirect.com/science/article/pii/S1359644620300465?via%3Dihub) | [arXiv](https://arxiv.org/abs/2002.06053). 
> Text-based representations of chemicals and proteins can be thought of as unstructured languages codified by humans to describe domain-specific knowledge. Advances in natural language processing (NLP) methodologies in the processing of spoken languages accelerated the application of NLP to elucidate hidden knowledge in textual representations of these biochemical entities and then use it to construct models to predict molecular properties or to design novel molecules. This review outlines the impact made by these advances on drug discovery and aims to further the dialogue between medicinal chemists and computer scientists.

[Towards reproducible computational drug discovery](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0408-x)
> This review explores the following topics: (1) the current state-of-the-art on reproducible research, (2) research documentation (e.g. electronic laboratory notebook, Jupyter notebook, etc.), (3) science of reproducible research (i.e. comparison and contrast with related concepts as replicability, reusability and reliability), (4) model development in computational drug discovery, (5) computational issues on model development and deployment, (6) use case scenarios for streamlining the computational drug discovery protocol. 

[A Deep Learning Approach to Antibiotic Discovery](https://www.cell.com/cell/pdf/S0092-8674(20)30102-1.pdf)

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">A new machine-learning paper on antibiotic discovery has some pretty interesting results - a harbinger?<a href="https://t.co/v3nDF3D2Dl">https://t.co/v3nDF3D2Dl</a></p>&mdash; Derek Lowe (@Dereklowe) <a href="https://twitter.com/Dereklowe/status/1230523404094234624?ref_src=twsrc%5Etfw">February 20, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


[Discovery of Novel Chemical Reactions by Deep Generative Recurrent Neural Network](https://chemrxiv.org/articles/Discovery_of_Novel_Chemical_Reactions_by_Deep_Generative_Recurrent_Neural_Network/11635929/1)
> Condensed Graph of Reaction (CGR) encodes the structures of reactants and products into a single molecular graph (See Fig 1).

![png](/images/SMILES-CGR.png)

### NLP

[REALM : Retrieval-Augmented Language Model Pre-Training](https://kentonl.com/pub/gltpc.2020.pdf)
> To capture knowledge in a more modular and interpretable way, we augment language model pretraining with a latent knowledge retriever, which allows the model to retrieve and attend over documents from a large corpus such as Wikipedia, used during pre-training, fine-tuning and inference. For the first time, we show how to pre-train such a knowledge retriever in an unsupervised manner, using masked language modeling as the learning signal and backpropagating through a retrieval step that considers millions of documents.

[Compressive Transformers for Long-Range Sequence Modelling](https://arxiv.org/abs/1911.05507)
> A new model and dataset for long-range memory.

[ImageBERT: CROSS-MODAL PRE-TRAINING WITH LARGE-SCALE WEAK-SUPERVISED IMAGE-TEXT DATA](https://arxiv.org/abs/2001.07966). 
>ImageBERT is a new vision-language pre-trained model for image-text joint embedding from Microsoft. The model is pre-trained on four tasks simultaneously: Masked Language Modeling (MLM), Masked Object Classification (MOC), Masked Region Feature Regression (MRFR), and Image Text Matching (ITM). Along with the pre-trained model, the team also collected a Large-scale weAk-supervised Image-Text (LAIT) dataset from Web. One interesting finding of this paper is that **multi-stage pre-training strategy outperforms single-stage pre-training**.

[Exploring the Limits of Transfer Learning with a Unified Text-to-Text Transformer](https://arxiv.org/abs/1910.10683)
> Transfer learning, where a model is first pre-trained on a data-rich task before being fine-tuned on a downstream task, has emerged as a powerful technique in natural language processing (NLP). The effectiveness of transfer learning has given rise to a diversity of approaches, methodology, and practice. In this paper, we explore the landscape of transfer learning techniques for NLP by introducing a unified framework that converts every language problem into a text-to-text format.

### New dataset
[CrossDocked Dataset](https://chemrxiv.org/articles/3D_Convolutional_Neural_Networks_and_a_CrossDocked_Dataset_for_Structure-Based_Drug_Design/11833323)
> 22.5 million poses of ligands docked into multiple similar binding pockets across the Protein Data Bank and

### others

[A Scientist’s Guide to Social Media](https://pubs.acs.org/doi/10.1021/acscentsci.9b01273?utm_source=pubs_outreach_marketing&utm_medium=twitter&utm_campaign=0220_SAP_Scientists_Guide_to_Social_media&ref=pubs_outreach_marketing#.XknOLQj-kQA.twitter)

[Machine Learning in Python: Main developments and technology trends in data science, machine learning, and artificial intelligence](https://arxiv.org/abs/2002.04803)

## Software and Tools 💻 

[TMAP](http://tmap.gdb.tools/) is a new data visualization method designed for visualizing large, high-dimensional data sets containing chemical structures and associated properties while preserving both global and local features.

[Computer Vision Recipes](https://github.com/microsoft/computervision-recipes). 
>Best practices for computer vision from **Microsoft**. The goal for this repository is to build a comprehensive set of tools and examples that leverage recent advances in Computer Vision algorithms, neural architectures, and operationalizing such systems. This project is drawn from existing state-of-the-art libraries and build additional utility around loading image data, optimizing and evaluating models, and scaling up to the cloud. **The examples are provided as Jupyter notebooks and common utility functions**. All examples use **Pytorch** as the underlying deep learning library.

[ligmolgrid](https://github.com/gnina/libmolgrid)
> Represent 3D molecues using multidimensional arrays of voxelized molecular data for grid-based machine learning modeling.

## Articles and Blog Posts 📃

[Illustrating the Reformer](https://towardsdatascience.com/illustrating-the-reformer-393575ac6ba0). 

[Transformer - Illustration and code.ipynb](https://github.com/vinsis/math-and-ml-notes/blob/master/notebooks/Transformer%20-%20Illustration%20and%20code.ipynb)

[The Annotated GPT-2](https://amaarora.github.io/2020/02/18/annotatedGPT2.html)
> Re-implement OpenAI's GPT-2 in PyTorch using Hugging Face source code and try to explain all the magic that goes on inside the model.

[Train a language model from scratch using Hugging Face's transformers and tokenizers](https://huggingface.co/blog/how-to-train?utm_source=Deep+Learning+Weekly&utm_campaign=2a200317df-EMAIL_CAMPAIGN_2019_04_24_03_18_COPY_01&utm_medium=email&utm_term=0_384567b42d-2a200317df-103924405)

[FastHugs - Fastai-v2 + HuggingFace](http://www.ntentional.com/2020/02/18/fasthugs_demo.html)
> Fine-tuning a test classification model with HuggingFace transformer and the new fastai v2 library.

[Interactive filtering of predicted ligands](https://ljmartin.github.io/2020/01/31/Interative_filtering_of_predicted_ligands.html)
> Interesting application of ipywidgets, RDKit and jupyter notebook.

[A first look at JAX](https://www.pragmatic.ml/first-look-at-jax/). 
> **JAX** is Numpy on the CPU, GPU and TPU, with greate automatic differentiation for high performance machine learning research.

[Rethinking the Chemical Reaction as a Graph: Imaginary Transition Structures and Beyond](https://depth-first.com/articles/2020/02/24/rethinking-the-chemical-reaction-as-a-graph-imaginary-transition-structures-and-beyond/)

[Self-Supervised Learning with Image网](https://datasciencecastnet.home.blog/2020/02/22/self-supervised-learning-with-image%e7%bd%91/). How to get started with self-supervised learning, including examples of how to run and analyze experiments.

## Notable Mentions ✨

[The Art of the Algorithm: Machine Learning in Environmental Health Research](https://ehp.niehs.nih.gov/doi/10.1289/EHP6874). In this podcast, Dr. **Nicole Kleinstreuer** talked about the promise and potential pitfalls of artificial intelligence as it relates to environmental health research. She is the acting director of the *Interagency Center for the Evaluation of Alternative Toxicological Methods within the National Toxicology* Program at NIEHS. 

[Colab Pro]. 10 Dollars/month, it delivers better GPUs and longer runtimes.

[Build a Pro Deep Learning Workstation... for Half the Price](https://l7.curtisnorthcutt.com/build-pro-deep-learning-workstation). 



