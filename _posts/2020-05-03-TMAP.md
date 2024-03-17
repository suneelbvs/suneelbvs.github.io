---
title: 'A walkthrough of TMAP'
date: 2020-05-03
permalink: /posts/2020/05/TMAP/
tags:
  - TMAP

---

The notebook version is [here](https://github.com/XinhaoLi74/Cheminformatics-Notebooks/tree/master/TMAP).

## 1. What is TMAP?

TMAP is a really cool interactive visualization for big data. Examples can be found [here](http://tmap.gdb.tools/#ex-coil).

Tree MAP (TMAP) is an algorithm developed by Dr. [Probst](https://twitter.com/skepteis) for [visualization of very large high-dimensional data sets as minimum spanning trees](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0416-x).

> Visualizations based on TMAP are better suited than t‑SNE or UMAP for the exploration and interpretation of large data sets due to their tree‑like nature, increased local and global neighborhood and structure preservation, and the transparency of the methods the algorithm is based on.

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">A map of all approved, illicit, and experimental drugs.<br><br>A fully interactive Faerun visualization of <a href="https://twitter.com/DrugBankDB?ref_src=twsrc%5Etfw">@DrugBankDB</a> embedded using TMAP from 2²³ dimensional MHFP space. This combines 4 projects I have been working on so far during my PhD.<br><br>Check it out here: <a href="https://t.co/lSUhjIMfn1">https://t.co/lSUhjIMfn1</a> <a href="https://t.co/0wBbDPOzJO">pic.twitter.com/0wBbDPOzJO</a></p>&mdash; Daniel Probst (@skepteis) <a href="https://twitter.com/skepteis/status/1175064262869487617?ref_src=twsrc%5Etfw">September 20, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

The four projects mentioned are:
1. [A probabilistic molecular fingerprint for big data settings](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0321-8)
2. [SmilesDrawer: parsing and drawing SMILES‑encoded molecular structures using client‑side JavaScript](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425)
3. [FUn: a framework for interactive visualizations of large, high-dimensional datasets on the web](https://academic.oup.com/bioinformatics/article/34/8/1433/4657075)
4. [Visualization of very large high-dimensional data sets as minimum spanning trees](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-0416-x)

## 2. How it works?

TMAP consists of 4 phases:
- [LSH forest](http://infolab.stanford.edu/~bawa/Pub/similarity.pdf) indexing
- Construction of a c-approximate k-nearest neighbor graph
- Calculation of a minimum spanning tree ([MST](https://en.wikipedia.org/wiki/Minimum_spanning_tree)) of the c-approximate k-nearest neighbor graph
- Generation of a graph layout for the resulting MST

The first two steps are designed for big data setting and can be replaced by other k-nearest neighbor algorithms. In the resulting graph, the nodes are molecules (or other entities) and the edges are weighted by the distances between nodes. In the 3rd step, a tree is basicly a graph wihch doest not contrain cycle. A [MST](https://en.wikipedia.org/wiki/Minimum_spanning_tree) is the tree has the minimum sum of edge weights. The MST is ready for visualization. In the last step, we choose a layout to visualize the MST.

## 2. Code Walkthrough

### Packages

Installation of (1) RDKit (2) [TMAP](http://tmap.gdb.tools/) (3) [MHFP](https://github.com/reymond-group/mhfp) and (4) [Faerun](https://github.com/reymond-group/faerun-python).

```
conda install -c rdkit rdkit
conda install -c tmap tmap
pip install mhfp
pip install faerun
```


```python
import pandas as pd
import tmap
from faerun import Faerun
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import AllChem
```

### Data

Here we use ESOL dataset includes 1117 molecules.


```python
url = 'https://raw.githubusercontent.com/XinhaoLi74/molds/master/clean_data/ESOL.csv'

df = pd.read_csv(url)
df.shape
```




    (1117, 2)




```python
df.head(1)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>smiles</th>
      <th>logSolubility</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c...</td>
      <td>-0.77</td>
    </tr>
  </tbody>
</table>
</div>



### Step 1: Compute descriptors

[MHFP6](https://github.com/reymond-group/mhfp) 
> MHFP6 (MinHash fingerprint, up to six bonds) is a molecular fingerprint which encodes detailed substructures using the extended connectivity principle of ECFP in a fundamentally different manner, increasing the performance of exact nearest neighbor searches in benchmarking studies and enabling the application of locality sensitive hashing (LSH) approximate nearest neighbor search algorithms. 


```python
# The number of permutations used by the MinHashing algorithm
perm = 512

# Initializing the MHFP encoder with 512 permutations
enc = MHFPEncoder(perm)

# Create MHFP fingerprints from SMILES
# The fingerprint vectors have to be of the tm.VectorUint data type
fingerprints = [tmap.VectorUint(enc.encode(s)) for s in df["smiles"]]
```

Data types for TMAP: `VectorUnit` and `VectorFloat`

### Step 2: LSH indexing and coordinates generation


```python
# Initialize the LSH Forest
lf = tmap.LSHForest(perm)

# Add the Fingerprints to the LSH Forest and index
lf.batch_add(fingerprints)
lf.index()
```


```python
# Get the coordinates
x, y, s, t, _ = tmap.layout_from_lsh_forest(lf)
```

`x` and `y` are the coordinates of the nodes.  `s` and `t` store the indexes of start nodes and to nodes in the MST, respectively. 

### Step 3: Plotting


```python
# Now plot the data
faerun = Faerun(view="front", coords=False)
faerun.add_scatter(
    "ESOL_Basic",
    {   "x": x, 
        "y": y, 
        "c": list(df.logSolubility.values), 
        "labels": df["smiles"]},
    point_scale=5,
    colormap = ['rainbow'],
    has_legend=True,
    legend_title = ['ESOL (mol/L)'],
    categorical=[False],
    shader = 'smoothCircle'
)

faerun.add_tree("ESOL_Basic_tree", {"from": s, "to": t}, point_helper="ESOL_Basic")

# Choose the "smiles" template to display structure on hover
faerun.plot('ESOL_Basic', template="smiles", notebook_height=750)
```



<iframe
    width="100%"
    height="750"
    src="/images/TMAP/ESOL_Basic.html"
    frameborder="0"
    allowfullscreen
></iframe>




<a href='/images/TMAP/ESOL_Basic.html' target='_blank'>/images/TMAP/ESOL_Basic.html</a><br>


## 3. Advances

### 3.1. Add more legends

Sometime we want to plot multiple labels for a single dataset. There are two types of lables, continues and categorical.  The type of labels need to explicitly assign to the `categorical` in `faerun.add_scatter`.

We compute two categorical labels: (1) the `number of rings` (2) Linear molecules (`is_linear`): `0` if numrings > 1 and `1` if numrings = 0. 


```python
from rdkit.Chem import rdMolDescriptors
from rdkit import Chem
numrings = [rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(s)) for s in df["smiles"]]
set(numrings)
```




    {0, 1, 2, 3, 4, 5, 6, 7, 8}




```python
is_linear = [1 if r == 0 else 0 for r in numrings]
```

The molecules in ESOL datast contains rings from 0 to 8. Now we are going to plot three labels. We need to change setting in `faerun.add_scatter`, information of multiple labels are passed as lists.


```python
# Now plot the data
faerun = Faerun(view="front", coords=False)
faerun.add_scatter(
    "ESOL",
    {   "x": x, 
        "y": y, 
        "c": [list(df.logSolubility.values), numrings, is_linear], 
        "labels": df["smiles"]},
    point_scale=5,
    colormap = ['rainbow', 'Set1'],
    has_legend=True,
    categorical=[False, True, True],
    series_title = ['ESOL (mol/L)', 'Rings', 'is_linear'],
    legend_labels = [None, None, [(0, "No"), (1, "Yes")]],
    shader = 'smoothCircle'
)

faerun.add_tree("ESOL_tree", {"from": s, "to": t}, point_helper="ESOL")

# Choose the "smiles" template to display structure on hover
faerun.plot('ESOL', template="smiles", notebook_height=750)
```



<iframe
    width="100%"
    height="750"
    src="/images/TMAP/ESOL.html"
    frameborder="0"
    allowfullscreen
></iframe>




<a href='/images/TMAP/ESOL.html' target='_blank'>/images/TMAP/ESOL.html</a><br>


## 3.2. Use different descriptors/fingerprints

We can also use other descriptors/fingerprints. The descriptors/fingerprints need to be converted to [Minhash vectors](http://tmap.gdb.tools/#tmap.Minhash) first. It supports binary, indexed, string and also int and float weighted vectors as input and returns a list of Minhash vectors (List of `VectorUint`)

Methods for different types of input: 

`batch_from_binary_array`. MinHash vectors from **binary** arrays. The input vectors need to be a list of `VectorUchar`.

`batch_from_int_weight_array`. Create MinHash vectors from **integer** arrays (not only zeros and ones). The input vectors need to be a list of `VectorUint`.

`batch_from_sparse_binary_array`. Create MinHash vectors from **sparse binary** arrays. The input vectors need to be a list of `VectorUint`  – A list of vectors containing indices of ones in a binary array.

`batch_from_string_array`. Create MinHash vectors from **string** arrays. The input vector is a list of `list` or `string`.

`batch_from_weight_array`. Create MinHash vectors from **float** arrays. The input vectors need to be a list of `VectorFloat`. – A list of vectors containing float values. **Keyword Arguments**: **method** (str) – The weighted hashing method to use (`ICWS` (default) or `I2CWS`).


Let's try ECFP4 (binary).


```python
bits = 1024

mols = [Chem.MolFromSmiles(s) for s in df['smiles']]
ECFP4_fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,bits) for x in mols]
ecfp4_lists = [tmap.VectorUchar(list(fp)) for fp in ECFP4_fps]
```


```python
# Initialize the Minhash
enc = tmap.Minhash(bits)

# Initialize the LSH Forest
lf_ecfp4 = tmap.LSHForest(bits)

# Add the Fingerprints to the LSH Forest and index
lf_ecfp4.batch_add(enc.batch_from_binary_array(ecfp4_lists))
lf_ecfp4.index()
```


```python
x, y, s, t, _ = tmap.layout_from_lsh_forest(lf_ecfp4)
```


```python
# Now plot the data
faerun = Faerun(view="front", coords=False)
faerun.add_scatter(
    "ESOL_ECFP4",
    {   "x": x, 
        "y": y, 
        "c": [list(df.logSolubility.values), numrings, is_linear], 
        "labels": df["smiles"]},
    point_scale=5,
    colormap = ['rainbow', 'Set1'],
    has_legend=True,
    categorical=[False, True, True],
    series_title = ['ESOL (mol/L)', 'Rings', 'is_linear'],
    legend_labels = [None, None, [(0, "No"), (1, "Yes")]],
    shader = 'smoothCircle'
)

faerun.add_tree("ESOL_ECFP4_tree", {"from": s, "to": t}, point_helper="ESOL_ECFP4")

# Choose the "smiles" template to display structure on hover
faerun.plot("ESOL_ECFP4",template="smiles", notebook_height=750)
```



<iframe
    width="100%"
    height="750"
    src="/images/TMAP/ESOL_ECFP4.html"
    frameborder="0"
    allowfullscreen
></iframe>




<a href='/images/TMAP/ESOL_ECFP4.html' target='_blank'>/images/TMAP/ESOL_ECFP4.html</a><br>

