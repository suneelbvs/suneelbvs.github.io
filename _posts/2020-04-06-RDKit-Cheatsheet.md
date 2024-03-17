---
title: 'My RDKit Cheatsheet'
date: 2020-04-06
permalink: /posts/2020/04/RDKit-Cheatsheet/
tags:
  - Cheatsheet
  - RDKit

---

Cheatsheet for `RDKit` package in python: (1) Draw molecules in jupyter enviroment; (2) use with `Pandas Dataframe` (3) Descriptors/Fingerprints and (4) Similarity Search etc. 

## Installation

The `RDKit` pacakge only supports `conda` installation.

```
conda -c rdkit rdkit

```

## Setup

```python
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

rdkit.__version__
```

    '2020.03.1'



### Chem vs. AllChem

As mentioned in the [Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html#chem-vs-allchem):

> The majority of “basic” chemical functionality (e.g. reading/writing molecules, substructure searching, molecular cleanup, etc.) is in the `rdkit.Chem` module. More advanced, or less frequently used, functionality is in `rdkit.Chem.AllChem`. 

> If you find the Chem/AllChem thing annoying or confusing, you can use python’s “import … as …” syntax to remove the irritation:


   ```python
    from rdkit.Chem import AllChem as Chem
   ```

## Basic

Get a `RDKit molecule` from SMILES. `RDKit molecule` enable several features to handle molecules: drawing, computing fingerprints/properties, molecular curation etc.


```python
smiles = 'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'
mol = Chem.MolFromSmiles(smiles)
print(mol)
```

    <rdkit.Chem.rdchem.Mol object at 0x000001F84A4CEE90>
    

The RDKit molecules can be directly printed in jupyter enviroment.


```python
mol
```




![png](/images/rdkit_cheatsheet/output_9_0.png)



Convert a RDKit molecule to SMILES.


```python
smi = Chem.MolToSmiles(mol)
smi
```




    'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'



Convert a RDKit molecule to InchiKey.


```python
Chem.MolToInchiKey(mol)
```




    'VSIUFPQOEIKNCY-UHFFFAOYSA-N'



Convert a RDKit molecule to coordinative representation (which can be stored in `.sdf` file).


```python
mol_block = Chem.MolToMolBlock(mol)
print(mol_block)
```

    
         RDKit          2D
    
     23 24  0  0  0  0  0  0  0  0999 V2000
        5.2500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        3.7500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.7500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
       -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -5.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
       -6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -7.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -5.2500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -3.0000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -3.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
       -3.0000    5.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -3.7500    6.4952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.5000    5.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      2  3  1  0
      3  4  2  0
      3  5  1  0
      5  6  2  0
      6  7  1  0
      7  8  1  0
      8  9  2  0
      9 10  1  0
     10 11  1  0
     11 12  1  0
     12 13  1  0
     12 14  1  0
     10 15  2  0
     15 16  1  0
     16 17  1  0
     17 18  1  0
     17 19  1  0
     15 20  1  0
     20 21  2  0
     21 22  1  0
     22 23  2  0
     22  5  1  0
     21  8  1  0
    M  END
    
    

### Reading sets of molecules

Major types of molecular file formats:
1. `.csv` file that includes a column of `SMILES`. See `PandasTools` section.
2. `.smi/.txt` file that includes `SMILES`. Collect the SMILES as a list. The following code is an example to read a `.smi` file that contains one SMILES per line.

```python
file_name = 'somedata.smi'

with open(file_name, "r") as ins:
    smiles = []
    for line in ins:
        smiles.append(line.split('\n')[0])
print('# of SMILES:', len(smiles))
```
3. `.sdf` file that includes `atom coordinates`. Reading molecules from `.sdf` file. [Code Example](https://www.rdkit.org/docs/GettingStartedInPython.html#reading-sets-of-molecules)

### Draw molecules in Jupter environment

Print molecules in grid.


```python
smiles = [
    'N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c1ccccc1',
    'c1ccc2c(c1)ccc1c2ccc2c3ccccc3ccc21',
    'C=C(C)C1Cc2c(ccc3c2OC2COc4cc(OC)c(OC)cc4C2C3=O)O1',
    'ClC(Cl)=C(c1ccc(Cl)cc1)c1ccc(Cl)cc1'
]

mols = [Chem.MolFromSmiles(smi) for smi in smiles]
```


```python
Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(200, 200))
```




![png](/images/rdkit_cheatsheet/output_20_0.png)



## PandasTools

`PandasTools` enables using RDKit molecules as columns of a `Pandas Dataframe`. 


```python
import pandas as pd
from rdkit.Chem import PandasTools
```


```python
url = 'https://raw.githubusercontent.com/XinhaoLi74/molds/master/clean_data/ESOL.csv'

esol_data = pd.read_csv(url)
esol_data.head(1)
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



Add `ROMol` to Pandas Dataframe.


```python
PandasTools.AddMoleculeColumnToFrame(esol_data, smilesCol='smiles')
esol_data.head(1)
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
      <th>ROMol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c...</td>
      <td>-0.77</td>
      <td><img data-content="rdkit/molecule" src="data:image/png;base64,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" alt="Mol"/></td>
    </tr>
  </tbody>
</table>
</div>



`ROMol` column stores `rdchem.Mol` object.


```python
print(type(esol_data.ROMol[0]))
```

    <class 'rdkit.Chem.rdchem.Mol'>
    

Draw the structures in grid.


```python
PandasTools.FrameToGridImage(esol_data.head(8), legendsCol="logSolubility", molsPerRow=4)
```




![png](/images/rdkit_cheatsheet/output_29_0.png)



Adding new columns of properites use `Pandas` [map](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.map.html) method.


```python
esol_data["n_Atoms"] = esol_data['ROMol'].map(lambda x: x.GetNumAtoms())
esol_data.head(1)
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
      <th>ROMol</th>
      <th>n_Atoms</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c...</td>
      <td>-0.77</td>
      <td><img data-content="rdkit/molecule" src="data:image/png;base64,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" alt="Mol"/></td>
      <td>32</td>
    </tr>
  </tbody>
</table>
</div>



Before saving the dataframe as csv file, it is recommanded to drop the `ROMol` column.


```python
esol_data = esol_data.drop(['ROMol'], axis=1)
esol_data.head(1)
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
      <th>n_Atoms</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c...</td>
      <td>-0.77</td>
      <td>32</td>
    </tr>
  </tbody>
</table>
</div>



## Descriptors/Fingerprints

The RDKit has avariety of built-in functionality for generating molecular fingerprints/descriptors. A detialed description can be found [here](https://www.rdkit.org/docs/RDKit_Book.html#additional-information-about-the-fingerprints). 


```python
url = 'https://raw.githubusercontent.com/XinhaoLi74/molds/master/clean_data/ESOL.csv'

esol_data = pd.read_csv(url)
PandasTools.AddMoleculeColumnToFrame(esol_data, smilesCol='smiles')
esol_data.head(1)
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
      <th>ROMol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c...</td>
      <td>-0.77</td>
      <td><img data-content="rdkit/molecule" src="data:image/png;base64,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" alt="Mol"/></td>
    </tr>
  </tbody>
</table>
</div>



### Morgan Fingerprint (ECFPx)

`AllChem.GetMorganFingerprintAsBitVect` Parameters:
1. `radius`: no default value, usually set 2 for similarity search and 3 for machine learning.
2. `nBits`: number of bits, default is 2048. 1024 is also widely used.
3. [other parameterss](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) are ususlly left to default

More examples can be found in this [notebook](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/descriptors.ipynb) from my previous work.


```python
radius=3
nBits=1024

ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=radius, nBits=nBits) for x in esol_data['ROMol']]
```


```python
ECFP6[0]
```




    <rdkit.DataStructs.cDataStructs.ExplicitBitVect at 0x1f84e70f4e0>



ECFP6 fingerprint for each molecule has 1024 bits.


```python
len(ECFP6[0])
```




    1024



Save as a `.csv` file for futher use (e.g., machine learning). I usually save (1) SMILES as index and (2) each bit as a column to the csv file. 


```python
ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
ecfp6_bits = [list(l) for l in ECFP6]
df_morgan = pd.DataFrame(ecfp6_bits, index = esol_data.smiles, columns=ecfp6_name)
df_morgan.head(1)
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
      <th>Bit_0</th>
      <th>Bit_1</th>
      <th>Bit_2</th>
      <th>Bit_3</th>
      <th>Bit_4</th>
      <th>Bit_5</th>
      <th>Bit_6</th>
      <th>Bit_7</th>
      <th>Bit_8</th>
      <th>Bit_9</th>
      <th>...</th>
      <th>Bit_1014</th>
      <th>Bit_1015</th>
      <th>Bit_1016</th>
      <th>Bit_1017</th>
      <th>Bit_1018</th>
      <th>Bit_1019</th>
      <th>Bit_1020</th>
      <th>Bit_1021</th>
      <th>Bit_1022</th>
      <th>Bit_1023</th>
    </tr>
    <tr>
      <th>smiles</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c1ccccc1</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 1024 columns</p>
</div>



## Similarity Search

Compute the similarity of a reference molecule and a list of molecules. Here is an example of using ECFP4 fingerprint to compute the `Tanimoto Similarity` (the default metric of [DataStructs.FingerprintSimilarity](https://www.rdkit.org/docs/source/rdkit.DataStructs.html).


1. compute fingerprints


```python
ref_smiles = 'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'
ref_mol = Chem.MolFromSmiles(ref_smiles)
ref_ECFP4_fps = AllChem.GetMorganFingerprintAsBitVect(ref_mol,2)
```


```python
ref_mol
```




![png](/images/rdkit_cheatsheet/output_46_0.png)




```python
bulk_ECFP4_fps = [AllChem.GetMorganFingerprintAsBitVect(x,2) for x in esol_data['ROMol']]
```


```python
from rdkit import DataStructs

similarity_efcp4 = [DataStructs.FingerprintSimilarity(ref_ECFP4_fps,x) for x in bulk_ECFP4_fps]
```

We can also add the `similarity_efcp4` to the dataframe and visualize the structure and similarity.


```python
esol_data['Tanimoto_Similarity (ECFP4)'] = similarity_efcp4
PandasTools.FrameToGridImage(esol_data.head(8), legendsCol="Tanimoto_Similarity (ECFP4)", molsPerRow=4)
```




![png](/images/rdkit_cheatsheet/output_50_0.png)



Sort the result from highest to lowest.


```python
esol_data = esol_data.sort_values(['Tanimoto_Similarity (ECFP4)'], ascending=False)
PandasTools.FrameToGridImage(esol_data.head(8), legendsCol="Tanimoto_Similarity (ECFP4)", molsPerRow=4)
```




![png](/images/rdkit_cheatsheet/output_52_0.png)



## More Reading

1. [Offical documentation](https://www.rdkit.org/docs/).
2. [Getting Started with the RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
3. [The RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html)
4. [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
> This document provides example recipes of how to carry out particular tasks using the RDKit functionality from Python. The contents have been contributed by the RDKit community, tested with the latest RDKit release, and then compiled into this document.
