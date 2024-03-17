---
title:  "Highlights of the Month: March 2020"
date: 2020-03-27
permalink: /posts/2020/03/Highlights-of-the-Month/
tags:
  - Highlights of the Month
---

**Key words** Sub-Word Tokenization; Ligand-Based Virtual Screening; Meta Learning; ELECTRA;

## Research Papers 🎓

### NLP:

[Pre-trained Models for Natural Language Processing: A Survey](https://arxiv.org/abs/2003.08271)
> A comprehensive review of Pre-trained Models for NLP

[A Primer in BERTology: What we know about how BERT works](https://arxiv.org/abs/2002.12327)

[ELECTRA: PRE-TRAINING TEXT ENCODERS AS DISCRIMINATORS RATHER THAN GENERATORS](https://openreview.net/pdf?id=r1xMH1BtvB). [[Github]](https://github.com/google-research/electra). [[Blog]](https://ai.googleblog.com/2020/03/more-efficient-nlp-model-pre-training.html)

### Sub-Word Tokenization:
- [Neural Machine Translation of Rare Words with Subword Units](https://www.aclweb.org/anthology/P16-1162/)

- [BPE-Dropout: Simple and Effective Subword Regularization](https://arxiv.org/abs/1910.13267)

- [Subword Regularization: Improving Neural Network Translation Models with Multiple Subword Candidates](https://arxiv.org/abs/1804.10959)

- [Neural Machine Translation with Byte-Level Subwords](https://arxiv.org/abs/1909.03341)


### Ligand-Based Virtual Screening:
> Ligand-based virtual screening utilizes information about molecules with known activity to predict the activity of new molecules. It assumes similar molecules tend to have similar activities/properties. Searching for active molecules can be viewed as similarity-based querying of a database storing molecules with unknown activity. Thus the use of molecular similarity measure is the cornerstone of the success of virtual screening. The following are papers of benchmarking different methods. 

- [Open-source platform to benchmark fingerprints for ligand-based virtual screening](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-26)

- [Heterogeneous Classifier Fusion for Ligand-Based Virtual Screening: Or, How Decision Making by Committee Can Be a Good Thing](https://pubs.acs.org/doi/abs/10.1021/ci400466r)

- [Benchmarking Platform for Ligand-Based Virtual Screening](https://ieeexplore.ieee.org/document/7822693)

- [Benchmarking Data Sets for the Evaluation of Virtual Ligand Screening Methods: Review and Perspectives](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5b00090)

- [Ligand-Based Virtual Screening Using Graph Edit Distance as Molecular Similarity Measure](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00820)

### Meta Learning:

[Meta-Learning: Learning to Learn Fast](https://lilianweng.github.io/lil-log/2018/11/30/meta-learning.html)

[Model-Agnostic Meta-Learning for Fast Adaptation of Deep Networks](https://arxiv.org/abs/1703.03400)

[Investigating Meta-Learning Algorithms for Low-Resource Natural Language Understanding Tasks](https://arxiv.org/abs/1908.10423)

[META-LEARNING INITIALIZATIONS FOR LOW-RESOURCE DRUG DISCOVERY](https://arxiv.org/abs/2003.05996)

### Others

[Molecule Attention Transformer](https://arxiv.org/abs/2002.08264). [[Code]](https://github.com/gmum/MAT/)
<p align='center'>
<img src="https://github.com/gmum/MAT/blob/master/assets/MAT.png" alt="architecture" width="600"/>
</p>


[Breaking Down Structural Diversity for Comprehensive Prediction of Ion-Neutral Collision Cross Sections](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.9b05772)

- A diverse set of 7405 molecues with CCS values.
- Machine learning models to predict CCS.

[A Deep Generative Model for Fragment-Based Molecule Generation](https://deepai.org/publication/a-deep-generative-model-for-fragment-based-molecule-generation)

## Software and Tools 💻 

[flair](https://github.com/flairNLP/flair)
> A very simple framework for state-of-the-art NLP. 

> **A text embedding library**. Flair has simple interfaces that allow you to use and combine different word and document embeddings, including our proposed Flair embeddings, BERT embeddings and ELMo embeddings.

[PyTorch implementation of SimCLR](https://github.com/Spijkervet/SimCLR)

[textualheatmap](https://github.com/AndreasMadsen/python-textualheatmap)
> Create interactive textual heatmaps for Jupiter notebooks.

[TensorFlow Quantum](https://ai.googleblog.com/2020/03/announcing-tensorflow-quantum-open.html)
> An Open Source Library for Quantum Machine Learning

<iframe width="560" height="315" src="https://www.youtube.com/embed/-o9AhIz1uvo" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

### Generative Teaching Networks:
[AI Generative Teaching Networks: Accelerating Neural Architecture Search by Learning to Generate Synthetic Training Data](https://eng.uber.com/generative-teaching-networks/)

[Implementation of Generative Teaching Networks for PyTorch](https://www.goodai.com/implementation-of-generative-teaching-networks-for-pytorch/).

## Articles and Blog Posts 📃

[Introducing BART](https://sshleifer.github.io/blog_v2/jupyter/2020/03/12/bart.html)

[Transformers are Graph Neural Networks](https://graphdeeplearning.github.io/post/transformers-are-gnns/)

[Graph Transformer tutorial](https://docs.dgl.ai/en/latest/tutorials/models/4_old_wines/7_transformer.html)

[From PyTorch to PyTorch Lightning — A gentle introduction](https://towardsdatascience.com/from-pytorch-to-pytorch-lightning-a-gentle-introduction-b371b7caaf09)
> PyTorch Lightning was created for professional researchers and PhD students working on AI research. PyTorch Lightning 

[WHAT IS TORCH.NN REALLY?](https://pytorch.org/tutorials/beginner/nn_tutorial.html)

[The Ultimate Guide to using the Python regex module](https://mlwhiz.com/blog/2019/09/01/regex/?utm_campaign=the-ultimate-guide-to-using-the-python-regex-module&utm_medium=social_link&utm_source=missinglettr-twitter)

[Tokenizers: How machines read](https://blog.floydhub.com/tokenization-nlp/)

[A Deep Dive into the Wonderful World of Preprocessing in NLP](http://mlexplained.com/2019/11/06/a-deep-dive-into-the-wonderful-world-of-preprocessing-in-nlp/)

[Machine Learning for Everyone](https://vas3k.com/blog/machine_learning/)
> Explain machine learning concepts with (1) simple words and (2) real-world examples.

[From PyTorch to JAX: towards neural net frameworks that purify stateful code](https://sjmielke.com/jax-purify.htm)

[Introducing DIET: state-of-the-art architecture that outperforms fine-tuning BERT and is 6X faster to train](https://blog.rasa.com/introducing-dual-intent-and-entity-transformer-diet-state-of-the-art-performance-on-a-lightweight-architecture/) 
> Dual Intent and Entity Transformer (DIET) is a multi-task transformer architecture that handles both intent classification and entity recognition together. 

[TRAINING ROBERTA FROM SCRATCH - THE MISSING GUIDE](https://zablo.net/blog/post/training-roberta-from-scratch-the-missing-guide-polish-language-model/)

[How to generate text: using different decoding methods for language generation with Transformers](https://huggingface.co/blog/how-to-generate)

[FROM Pre-trained Word Embeddings TO Pre-trained Language Models — Focus on BERT](https://towardsdatascience.com/from-pre-trained-word-embeddings-to-pre-trained-language-models-focus-on-bert-343815627598)

### Meta learning
> In meta-learning, there is a meta-learner and a learner. The meta-learner (or the agent) trains the learner (or the model) on a training set that contains a large number of different tasks. In this stage of meta-learning, the model will acquire a prior experience from training and will learn the common features representations of all the tasks. Then, whenever, there is a new task to learn, the model with its prior experience will be fine-tuned using the small amount of the new training data brought by that task.

- [From zero to research — An introduction to Meta-learning](https://medium.com/huggingface/from-zero-to-research-an-introduction-to-meta-learning-8e16e677f78a)

- [What is Model-Agnostic Meta-learning (MAML)?](https://towardsdatascience.com/model-agnostic-meta-learning-maml-8a245d9bc4ac)

## Notable Mentions ✨

[AI Curriculum](https://github.com/Machine-Learning-Tokyo/AI_Curriculum). Open Deep Learning and Reinforcement Learning lectures from top Universities like Stanford University, MIT, UC Berkeley.

[Huggingface's official notebook tutorials](https://github.com/huggingface/transformers/tree/master/notebooks)

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Today we&#39;re happy to release four new official notebook tutorials available in our documentation and in colab thanks to <a href="https://twitter.com/MorganFunto?ref_src=twsrc%5Etfw">@MorganFunto</a> to get started with tokenizers and transformer models in just seconds! <a href="https://t.co/zzBVWsEnef">https://t.co/zzBVWsEnef</a> (1/6) <a href="https://t.co/bpOYA6UTDk">pic.twitter.com/bpOYA6UTDk</a></p>&mdash; Hugging Face (@huggingface) <a href="https://twitter.com/huggingface/status/1235583660415488001?ref_src=twsrc%5Etfw">March 5, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

[Limitations of Graph Neural Networks (Stanford University)](https://www.youtube.com/watch?v=H6oOhElB3yE&feature=emb_logo)

<iframe width="560" height="315" src="https://www.youtube.com/embed/H6oOhElB3yE" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>



