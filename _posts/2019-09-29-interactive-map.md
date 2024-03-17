---
title: 'Create an Interactive Visualization for Approved Small Molecule Drugs'
date: 2019-09-29
permalink: /posts/2019/09/interactive-map/
tags:

---

In the post, I want to generate an [interactive visualization](http://XinhaoLi74.github.io/images/samll_molecule_drug_ECFP4.html) of a chemical space. Each point in the map represents a molecule and close points have similar structures. When you move you mouse on a point, the name and structure of the moelcule will show up.


```python
from chembl_webresource_client.new_client import new_client
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

from sklearn.manifold import TSNE

from IPython.display import SVG
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.PandasTools import ChangeMoleculeRendering

#Bokeh library for plotting
import json
from bokeh.plotting import figure, show, output_notebook, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.transform import factor_cmap
from bokeh.plotting import figure, output_file, save
```



## Getting the small molecule drugs from ChEMBL

[ChEMBL](https://www.ebi.ac.uk/chembl/) is a manually curated database of bioactive molecules with drug-like properties. It brings together chemical, bioactivity and genomic data to aid the translation of genomic information into effective new drugs. Here we use the python client library of [ChEMBL API](https://github.com/chembl/chembl_webresource_client) to download the [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) for all of the approved small molecule drugs and put them into a Pandas dataframe. 


```python
molecule = new_client.molecule
approved_drugs = molecule.filter(max_phase=4)
small_molecule_drugs = [x for x in approved_drugs if x['molecule_type'] == 'Small molecule']
```

Extract information we need: (1) drug name (2) CHEMBL ID and (3) Canonical SMILES.


```python
struct_list = [(x['pref_name'], x['molecule_chembl_id'],x['molecule_structures'])for x in small_molecule_drugs if x]
smiles_list = [(a,b,c['canonical_smiles']) for (a,b,c) in struct_list if c]
```


```python
smiles_df = pd.DataFrame(smiles_list)
smiles_df.columns = ['Name','ChEMBL_ID','SMILES']
print(f'We downloaded {smiles_df.shape[0]} small molecule drugs from ChEMBL.')
```

    We downloaded 3147 small molecule drugs from ChEMBL.
    

Let check the data.


```python
smiles_df.head()
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
      <th>Name</th>
      <th>ChEMBL_ID</th>
      <th>SMILES</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PRAZOSIN</td>
      <td>CHEMBL2</td>
      <td>COc1cc2nc(nc(N)c2cc1OC)N3CCN(CC3)C(=O)c4occc4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NICOTINE</td>
      <td>CHEMBL3</td>
      <td>CN1CCC[C@H]1c2cccnc2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>OFLOXACIN</td>
      <td>CHEMBL4</td>
      <td>CC1COc2c(N3CCN(C)CC3)c(F)cc4C(=O)C(=CN1c24)C(=O)O</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NALIDIXIC ACID</td>
      <td>CHEMBL5</td>
      <td>CCN1C=C(C(=O)O)C(=O)c2ccc(C)nc12</td>
    </tr>
    <tr>
      <th>4</th>
      <td>INDOMETHACIN</td>
      <td>CHEMBL6</td>
      <td>COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c3ccc(Cl)cc3</td>
    </tr>
  </tbody>
</table>
</div>

## Molecular Fingerprints

Add RDKit Mol column to the dataframe.


```python
PandasTools.AddMoleculeColumnToFrame(smiles_df,smilesCol='SMILES')
```


```python
smiles_df.head(1)
```




<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Name</th>
      <th>ChEMBL_ID</th>
      <th>SMILES</th>
      <th>ROMol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PRAZOSIN</td>
      <td>CHEMBL2</td>
      <td>COc1cc2nc(nc(N)c2cc1OC)N3CCN(CC3)C(=O)c4occc4</td>
      <td><img src="data:image/png;base64,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" alt="Mol"/></td>
    </tr>
  </tbody>
</table>



Generare [ECFP4](https://docs.chemaxon.com/display/docs/Extended+Connectivity+Fingerprint+ECFP) fingerprint for each drug. 


```python
ECFP4_fps = [AllChem.GetMorganFingerprintAsBitVect(x,2) for x in smiles_df['ROMol']]
```

## Dimensionality Reduction of Features

Use t-SNE to reduce the dimension of fetures into 2 for visualization.


```python
tsne = TSNE(random_state=0).fit_transform(ECFP4_fps)
```

Define some functions for interactive visualization. Some of the codes are from [mol2vec_notebook](https://github.com/samoturk/mol2vec_notebooks/tree/master/Notebooks).


```python
def _prepareMol(mol,kekulize):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    return mc

def moltosvg(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    mc = _prepareMol(mol,kekulize)
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return SVG(svg.replace('svg:',''))
```


Get the image of molecules.
```python
svgs = [moltosvg(m).data for m in smiles_df.ROMol]
```

Generate the interactive visualization.
```python
ChangeMoleculeRendering(renderer='PNG')


source = ColumnDataSource(data=dict(x=tsne[:,0], y=tsne[:,1], desc= smiles_df.Name, 
                                    svgs=svgs))

hover = HoverTool(tooltips="""
    <div>
        <div>@svgs{safe}
        </div>
        <div>
            <span style="font-size: 17px; font-weight: bold;">@desc</span>
        </div>
    </div>
    """
)
interactive_map = figure(plot_width=1000, plot_height=1000, tools=['reset,box_zoom,wheel_zoom,zoom_in,zoom_out,pan',hover],
           title="Small Molecule Drug (ECFP4)")



interactive_map.circle('x', 'y', size=5, source=source, fill_alpha=0.2);


```

Save is as a ```html``` file. Please [click](http://XinhaoLi74.github.io/images/samll_molecule_drug_ECFP4.html) here to check the generated map. 
```python
output_file("interactive_map.html")
save(interactive_map)
```

All the code and data are available on my [github](https://github.com/XinhaoLi74/Bokeh-Interactive-Plots/blob/master/2019-09-29-interactive_map.ipynb).