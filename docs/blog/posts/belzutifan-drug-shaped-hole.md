---
date: 2026-08-16
title: "Belzutifan: Finding the Drug-Shaped Hole in an 'Undruggable' Protein"
categories:
  - Drug Discovery Case Studies
tags:
  - Computational Chemistry
  - Medicinal Chemistry
  - Structure-Based Drug Design
  - HIF-2α
  - Belzutifan
---

# Belzutifan: Finding the Drug-Shaped Hole in an “Undruggable” Protein

### A computational-chemistry + medicinal-chemistry design-history case study

Most of us learned the rule early: **transcription factors are undruggable**. Flat interfaces, no deep pocket, protein–protein interactions spread over hundreds of square angstroms. You do not put a small molecule on that and expect a drug.

Belzutifan (Welireg; formerly PT2977 / MK-6482) is a particularly useful counterexample — not because it proved every transcription factor druggable, but because it shows how the right structural observation can redefine the problem. The story moves from a buried water-filled cavity, through allosteric inhibition of HIF-2α, into a clinical pharmacokinetic failure caused by glucuronidation, and finally to a fluorine-enabled redesign.

<!-- more -->

<div class="bz-shell">
  <div class="bz-hero">
    <div>
      <span class="bz-kicker">INTERACTIVE DESIGN HISTORY</span>
      <h2>From cryptic pocket → clinical drug</h2>
      <p>Explore the HIF-2α–ARNT–belzutifan complex in 3D while following the decisions that transformed a supposedly undruggable transcription factor into an approved small-molecule target.</p>
    </div>
    <div class="bz-statgrid">
      <div class="bz-stat"><b>~290 Å³</b><span>buried PAS-B cavity</span></div>
      <div class="bz-stat"><b>7W80</b><span>HIF-2 complex + belzutifan</span></div>
      <div class="bz-stat"><b>72Q</b><span>PDB ligand ID</span></div>
      <div class="bz-stat"><b>2.75 Å</b><span>crystal resolution</span></div>
    </div>
  </div>

  <div class="bz-lab">
    <aside class="bz-steps" aria-label="Case-study stages">
      <button class="bz-step active" data-scene="complex"><span>01</span><b>The complex</b><small>HIF-2α + ARNT</small></button>
      <button class="bz-step" data-scene="hif"><span>02</span><b>The hidden pocket</b><small>HIF-2α PAS-B</small></button>
      <button class="bz-step" data-scene="ligand"><span>03</span><b>Belzutifan</b><small>Ligand 72Q</small></button>
      <button class="bz-step" data-scene="pocket"><span>04</span><b>Pocket residues</b><small>Shape + allostery</small></button>
      <button class="bz-step" data-scene="surface"><span>05</span><b>Buried cavity</b><small>Surface context</small></button>
      <button class="bz-step" data-scene="m252"><span>06</span><b>Allosteric relay</b><small>Focus M252</small></button>
    </aside>

    <div class="bz-viewwrap">
      <div id="bz-viewer" class="bz-viewer">
        <div id="bz-loading" class="bz-loading">Loading PDB 7W80…</div>
      </div>
      <div class="bz-viewbar">
        <div><span class="dot teal"></span> ARNT <span class="dot amber"></span> HIF-2α <span class="dot pink"></span> Belzutifan</div>
        <div class="bz-view-actions">
          <button id="bz-spin">Spin</button>
          <button id="bz-reset">Reset</button>
          <button id="bz-fullscreen">Fullscreen</button>
        </div>
      </div>
    </div>

    <div class="bz-insight" id="bz-insight">
      <span class="bz-kicker">01 · STRUCTURAL CONTEXT</span>
      <h3>HIF-2α does not act alone.</h3>
      <p>HIF-2α forms a transcriptionally active heterodimer with ARNT. The drug does not block the DNA-binding surface directly; it binds within the HIF-2α PAS-B domain and destabilizes the protein–protein interaction allosterically.</p>
      <div class="bz-callout"><b>Design question</b><br>Can a buried site control a distant protein–protein interface strongly enough to become pharmacology?</div>
    </div>
  </div>
</div>

## 1. Why HIF-2α was worth the fight

In clear-cell renal cell carcinoma, loss of functional VHL causes HIF-2α to accumulate rather than being efficiently degraded. HIF-2α then dimerizes with ARNT and drives a pro-tumor transcriptional program. The genetics made HIF-2α compelling; the modality made it difficult.

The central problem was obvious: **where do you put a small molecule on a transcription factor?**

## 2. The computational-chemistry lesson: read the water-filled hole

Structural work on the HIF-2α PAS-B domain revealed an unusually large, completely buried cavity in its hydrophobic core — approximately **290 Å³**. In the apo state, an ordered internal water network occupies that volume.

That observation changes the interpretation of “undruggable.” A buried hydrophobic cavity stabilized by structured water can be a **latent druggability signal** rather than an irrelevant crystallographic curiosity.

Three ideas matter computationally:

- **The static crystal structure can lie by omission.** An apparently sealed pocket may still be kinetically accessible because proteins breathe.
- **Desolvation and shape complementarity can dominate.** The ligand does not need an elaborate hydrogen-bond network if binding replaces energetically constrained water and packs efficiently.
- **Allostery must be modeled as part of the mechanism.** Occupancy of the PAS-B cavity perturbs the HIF-2α:ARNT interface rather than competing at DNA directly.

<div class="bz-decision">
  <span>COMPCHEM GATE</span>
  <h3>Do not ask only “Is there a visible pocket?”</h3>
  <p>Ask whether conformational sampling, internal waters, transient openings and cavity energetics reveal a pocket the static structure under-represents.</p>
</div>

## 3. Hit-to-lead: building into the cavity

Academic ligands established that the internal PAS-B cavity could be occupied. Structure-based optimization then advanced the series through compounds such as PT2399 and into **PT2385**, which clinically validated the target mechanism.

That was an enormous milestone: an orally administered small molecule could allosterically inhibit a transcription factor once considered outside conventional druggability rules.

But target validation was not the end of the design problem.

## 4. The medicinal-chemistry lesson: potency was not the fatal liability

PT2385 contained a benzylic secondary alcohol important for binding. Unfortunately, that same alcohol created a significant **phase-II metabolic liability**: glucuronidation, largely associated with UGT2B17, produced an inactive conjugate and variable exposure.

The naïve response would be to remove the hydroxyl. But if the hydroxyl is part of the pharmacophore, deleting it solves metabolism by breaking potency.

The better design question was:

> **Can we preserve the alcohol and make it less attractive to glucuronidation?**

That is a multi-parameter optimization problem, not a single-property optimization problem.

## 5. Fluorine as a multi-parameter design tool

The second-generation work changed the fluorination pattern and stereochemical arrangement around the indane alcohol, ultimately producing the all-cis **(1S,2S,3R)** configuration of belzutifan.

The redesign simultaneously addressed several objectives:

- slowed glucuronidation while preserving the potency-critical alcohol;
- improved potency through a more favorable bound geometry and local electronics;
- lowered lipophilicity and improved developability;
- retained an overall profile compatible with clinical advancement.

<div class="bz-compare">
  <div>
    <span>FIRST-GENERATION LESSON</span>
    <h3>PT2385</h3>
    <p>Target engagement worked. Clinical activity existed. Exposure variability and glucuronidation capped the molecule.</p>
  </div>
  <div class="arrow">→</div>
  <div>
    <span>REDESIGN OBJECTIVE</span>
    <h3>Preserve + detune</h3>
    <p>Keep the pharmacophoric OH, alter its electronic environment, and improve multiple properties at once.</p>
  </div>
  <div class="arrow">→</div>
  <div>
    <span>SECOND GENERATION</span>
    <h3>Belzutifan</h3>
    <p>A cleaner balance of potency, metabolism, lipophilicity and clinical exposure.</p>
  </div>
</div>

## 6. The allosteric mechanism in 3D

PDB **7W80** captures the heterodimeric HIF-2 complex with belzutifan. The structure shows belzutifan buried in the HIF-2α PAS-B pocket. Residue **M252** is particularly important in coupling pocket occupancy to destabilization of the HIF-2α–ARNT interface.

Use the interactive viewer above to switch between the full heterodimer, the HIF-2α chain, ligand-only view, selected pocket residues, surface context and the M252 allosteric relay.

## 7. What this case teaches drug designers

<div class="grid cards" markdown>

-   :material-water:{ .lg .middle } **Water can reveal druggability**

    ---
    A large buried cavity filled with structured water is not necessarily dead space. It can identify a latent ligand-binding opportunity.

-   :material-axis-arrow:{ .lg .middle } **Model motion, not only snapshots**

    ---
    A pocket can look sealed and still bind rapidly. Static crystallography should be complemented with conformational sampling and pocket-dynamics thinking.

-   :material-vector-intersection:{ .lg .middle } **Allostery changes the geometry of inhibition**

    ---
    You do not need to compete at the obvious functional interface if a remote cavity can destabilize that interface.

-   :material-flask-outline:{ .lg .middle } **Do not optimize potency in isolation**

    ---
    PT2385 already demonstrated biological activity. The decisive second-generation problem was metabolism and exposure.

-   :material-molecule:{ .lg .middle } **Preserve-and-detune can beat delete-and-replace**

    ---
    When a metabolically vulnerable group is also pharmacophoric, alter its local electronics or geometry before discarding it.

-   :material-human-male-board:{ .lg .middle } **Human judgment defines the gates**

    ---
    Computation can iterate tirelessly, but the highest-value decisions are often deciding what the real problem is and which compromise is acceptable.

</div>

## Final perspective

The Belzutifan story is valuable because the winning decisions happened at different layers of drug discovery. First came the structural judgment that a water-filled internal cavity was worth pursuing. Then came the mechanistic realization that filling it could allosterically weaken a transcription-factor complex. Finally came the medicinal-chemistry judgment that a clinical PK limitation should be solved by **preserving and electronically detuning** a necessary functional group rather than simply deleting it.

That is the broader lesson: automated loops perform iteration; **domain judgment defines what the loop should optimize toward**.

---

### Structure and primary reading

- **PDB 7W80** — crystal structure of the heterodimeric HIF-2 complex with belzutifan; ligand ID **72Q**.
- Key discovery literature includes early HIF-2α PAS-B cavity/ligand studies, PT2385 medicinal chemistry, the PT2977/belzutifan redesign, and subsequent structural work on HIF-2α inhibition and resistance.

<style>
.bz-shell{margin:1.2rem calc(50% - 50vw + 1rem) 2.4rem;max-width:calc(100vw - 2rem);border:1px solid var(--md-default-fg-color--lightest);border-radius:20px;overflow:hidden;background:var(--md-default-bg-color);box-shadow:0 20px 50px -28px rgba(14,116,144,.55)}
.bz-hero{display:grid;grid-template-columns:1.35fr 1fr;gap:1.2rem;padding:2rem;background:linear-gradient(125deg,#083344,#0e7490 50%,#0f766e);color:white}.bz-hero h2{color:white!important;border:0!important;margin:.2rem 0 .7rem!important;font-size:1.65rem}.bz-hero p{color:#cffafe;margin:0;max-width:50rem}.bz-kicker{font-size:.62rem;font-weight:800;letter-spacing:.14em;color:#fcd34d}.bz-statgrid{display:grid;grid-template-columns:1fr 1fr;gap:.6rem}.bz-stat{background:rgba(255,255,255,.09);border:1px solid rgba(255,255,255,.18);border-radius:12px;padding:.8rem}.bz-stat b{display:block;font-size:1.15rem;color:#fff}.bz-stat span{font-size:.68rem;color:#cffafe}.bz-lab{display:grid;grid-template-columns:180px minmax(420px,1.5fr) minmax(250px,.8fr);min-height:520px}.bz-steps{padding:.7rem;background:color-mix(in srgb,var(--md-default-bg-color) 94%,#0e7490 6%);border-right:1px solid var(--md-default-fg-color--lightest)}.bz-step{width:100%;display:grid;grid-template-columns:28px 1fr;text-align:left;column-gap:.45rem;padding:.65rem .55rem;margin:0 0 .4rem;border:1px solid transparent;background:transparent;border-radius:10px;color:var(--md-default-fg-color);cursor:pointer}.bz-step span{grid-row:1/3;font-size:.6rem;font-weight:800;color:#0891b2;padding-top:.12rem}.bz-step b{font-size:.74rem}.bz-step small{font-size:.61rem;color:var(--md-default-fg-color--light)}.bz-step:hover,.bz-step.active{background:color-mix(in srgb,var(--md-default-bg-color) 88%,#0891b2 12%);border-color:color-mix(in srgb,#0891b2 35%,transparent)}.bz-viewwrap{position:relative;min-width:0;background:#08141b}.bz-viewer{height:470px;position:relative}.bz-loading{position:absolute;inset:0;display:grid;place-items:center;color:#bae6fd;font-size:.8rem;letter-spacing:.04em}.bz-viewbar{min-height:50px;display:flex;align-items:center;justify-content:space-between;gap:.5rem;padding:.55rem .8rem;color:#dbeafe;background:#071017;font-size:.65rem;border-top:1px solid rgba(255,255,255,.08)}.dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin:0 .15rem 0 .5rem}.dot.teal{background:#22d3ee}.dot.amber{background:#fbbf24}.dot.pink{background:#f472b6}.bz-view-actions{display:flex;gap:.35rem}.bz-view-actions button{border:1px solid rgba(255,255,255,.18);background:#102630;color:#e0f2fe;border-radius:7px;padding:.28rem .5rem;cursor:pointer;font-size:.62rem}.bz-insight{padding:1.5rem;border-left:1px solid var(--md-default-fg-color--lightest);background:var(--md-default-bg-color)}.bz-insight h3{font-size:1.05rem;margin:.35rem 0 .75rem}.bz-insight p{font-size:.78rem;line-height:1.65}.bz-callout{margin-top:1.2rem;border-left:3px solid #f59e0b;padding:.75rem;background:color-mix(in srgb,var(--md-default-bg-color) 90%,#f59e0b 10%);font-size:.7rem}.bz-decision{margin:1.5rem 0;padding:1.1rem 1.2rem;border-radius:14px;background:linear-gradient(120deg,color-mix(in srgb,var(--md-default-bg-color) 88%,#0891b2 12%),color-mix(in srgb,var(--md-default-bg-color) 94%,#f59e0b 6%));border:1px solid var(--md-default-fg-color--lightest)}.bz-decision span,.bz-compare span{font-size:.6rem;font-weight:800;letter-spacing:.12em;color:#0e7490}.bz-decision h3{margin:.25rem 0 .35rem}.bz-decision p{margin:0}.bz-compare{display:grid;grid-template-columns:1fr auto 1fr auto 1fr;gap:.8rem;align-items:stretch;margin:1.6rem 0}.bz-compare>div:not(.arrow){padding:1rem;border:1px solid var(--md-default-fg-color--lightest);border-radius:12px}.bz-compare h3{margin:.25rem 0 .4rem}.bz-compare p{font-size:.72rem;margin:0}.bz-compare .arrow{display:grid;place-items:center;font-size:1.3rem;color:#0891b2}@media(max-width:1100px){.bz-lab{grid-template-columns:150px 1fr}.bz-insight{grid-column:1/-1;border-left:0;border-top:1px solid var(--md-default-fg-color--lightest)}.bz-viewer{height:430px}}@media(max-width:760px){.bz-shell{margin:1rem 0 2rem;max-width:none}.bz-hero{grid-template-columns:1fr;padding:1.25rem}.bz-lab{display:block}.bz-steps{display:grid;grid-template-columns:1fr 1fr;border-right:0;border-bottom:1px solid var(--md-default-fg-color--lightest)}.bz-step{margin:.15rem}.bz-viewer{height:360px}.bz-viewbar{align-items:flex-start;flex-direction:column}.bz-compare{grid-template-columns:1fr}.bz-compare .arrow{transform:rotate(90deg)}}
</style>

<script src="https://unpkg.com/ngl@2.0.0-dev.40/dist/ngl.js"></script>
<script>
(function(){
  const info={
    complex:{k:'01 · STRUCTURAL CONTEXT',t:'HIF-2α does not act alone.',p:'HIF-2α forms a transcriptionally active heterodimer with ARNT. Belzutifan acts from inside the HIF-2α PAS-B domain rather than blocking the DNA-binding surface directly.',q:'Can a buried site control a distant protein–protein interface strongly enough to become pharmacology?'},
    hif:{k:'02 · CRYPTIC DRUGGABILITY',t:'The useful pocket is buried inside HIF-2α.',p:'The PAS-B domain contains the internal cavity that made this target chemically tractable. The key computational insight is that an apparently closed protein can still sample ligand-accessible states.',q:'Treat structured internal water and transient pocket opening as design information, not noise.'},
    ligand:{k:'03 · LIGAND',t:'Belzutifan occupies the PAS-B cavity.',p:'Ligand 72Q is shown as ball-and-stick. Binding is dominated by fitting a hydrophobic internal cavity and replacing its pre-existing solvent network rather than by maximizing a dense polar interaction map.',q:'Shape complementarity and desolvation can be as important as hydrogen-bond counting.'},
    pocket:{k:'04 · POCKET ARCHITECTURE',t:'A compact hydrophobic environment surrounds the ligand.',p:'Representative pocket and relay residues are highlighted around belzutifan, including H248, M252, F254, F280, Y307, M309 and L319 in HIF-2α.',q:'The pocket should be read as a 3D energetic environment, not a list of contacts.'},
    surface:{k:'05 · BURIED CAVITY',t:'The ligand is hidden from bulk solvent.',p:'The translucent HIF-2α surface makes the central design paradox visible: the drug binds deeply inside a site that appears inaccessible in a single static structure.',q:'When a ligand reaches a “sealed” pocket quickly, the missing variable is protein dynamics.'},
    m252:{k:'06 · ALLOSTERIC RELAY',t:'M252 helps connect pocket occupancy to dimer disruption.',p:'Focusing on M252 illustrates why this is an allosteric drug: local packing changes inside PAS-B propagate toward the HIF-2α–ARNT interface and weaken heterodimerization.',q:'For allosteric projects, model the transmission pathway as carefully as the ligand pose.'}
  };
  const viewer=document.getElementById('bz-viewer'); if(!viewer||typeof NGL==='undefined') return;
  const stage=new NGL.Stage('bz-viewer',{backgroundColor:'#08141b'}); let comp=null, spinning=false;
  const resize=()=>stage.handleResize(); window.addEventListener('resize',resize);
  function rep(type,params){return comp.addRepresentation(type,params||{});}
  function render(scene){
    if(!comp)return; comp.removeAllRepresentations();
    if(scene==='complex'){
      rep('cartoon',{sele:':A',color:'#22d3ee',opacity:.9}); rep('cartoon',{sele:':B',color:'#fbbf24',opacity:.9}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.3}); comp.autoView();
    } else if(scene==='hif'){
      rep('cartoon',{sele:':B',color:'#fbbf24'}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.25}); comp.autoView(':B');
    } else if(scene==='ligand'){
      rep('cartoon',{sele:':B',color:'#64748b',opacity:.25}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.65}); comp.autoView('[72Q]',700);
    } else if(scene==='pocket'){
      rep('cartoon',{sele:':B',color:'#475569',opacity:.18}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.35}); rep('licorice',{sele:'(248:B or 252:B or 254:B or 280:B or 307:B or 309:B or 319:B)',color:'#fbbf24',radiusScale:.8}); comp.autoView('[72Q]',700);
    } else if(scene==='surface'){
      rep('surface',{sele:':B',color:'#fbbf24',opacity:.28,surfaceType:'av'}); rep('cartoon',{sele:':B',color:'#fbbf24',opacity:.22}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.45}); comp.autoView('[72Q]',700);
    } else if(scene==='m252'){
      rep('cartoon',{sele:':A',color:'#22d3ee',opacity:.3}); rep('cartoon',{sele:':B',color:'#fbbf24',opacity:.35}); rep('ball+stick',{sele:'[72Q]',color:'#f472b6',radiusScale:1.25}); rep('spacefill',{sele:'252:B',color:'#ef4444',radiusScale:.85}); comp.autoView('252:B',700);
    }
    const x=info[scene]; document.getElementById('bz-insight').innerHTML='<span class="bz-kicker">'+x.k+'</span><h3>'+x.t+'</h3><p>'+x.p+'</p><div class="bz-callout"><b>Design question</b><br>'+x.q+'</div>';
  }
  stage.loadFile('https://files.rcsb.org/download/7W80.pdb',{defaultRepresentation:false}).then(function(c){comp=c;document.getElementById('bz-loading').style.display='none';render('complex');});
  document.querySelectorAll('.bz-step').forEach(b=>b.addEventListener('click',function(){document.querySelectorAll('.bz-step').forEach(x=>x.classList.remove('active'));this.classList.add('active');render(this.dataset.scene);}));
  document.getElementById('bz-spin').addEventListener('click',function(){spinning=!spinning;stage.setSpin(spinning);this.textContent=spinning?'Stop':'Spin';});
  document.getElementById('bz-reset').addEventListener('click',()=>render(document.querySelector('.bz-step.active').dataset.scene));
  document.getElementById('bz-fullscreen').addEventListener('click',()=>viewer.requestFullscreen&&viewer.requestFullscreen());
})();
</script>
