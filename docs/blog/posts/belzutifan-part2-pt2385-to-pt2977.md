---
date: 2026-08-16
title: "Belzutifan Part 2: PT2385 → PT2977 — Fixing Metabolism Without Losing the Pharmacophore"
categories:
  - Drug Discovery Case Studies
tags:
  - Medicinal Chemistry
  - DMPK
  - Structure-Based Drug Design
  - HIF-2α
  - Belzutifan
---

# Belzutifan Part 2: PT2385 → PT2977

### Fixing glucuronidation without deleting the potency-critical alcohol

[← Part 1: Finding the Drug-Shaped Hole](belzutifan-drug-shaped-hole.md)

PT2385 proved that HIF-2α could be inhibited in patients. The next problem was not target biology. It was **exposure**. Extensive glucuronidation of a benzylic alcohol produced variable, dose-limited pharmacokinetics — yet that same alcohol was important for binding.

The second-generation design problem was therefore unusually clean:

> **Keep the pharmacophore. Change the chemistry around it enough to slow phase-II metabolism.**

<!-- more -->

<div class="mco-shell">
  <div class="mco-head">
    <div><span class="mco-kicker">INTERACTIVE MEDCHEM OPTIMIZATION</span><h2>Same pocket. Better molecule.</h2><p>Compare the experimental HIF-2 complexes for PT2385 and PT2977/belzutifan, then follow the design logic from clinical PK limitation to vicinal difluoro optimization.</p></div>
    <div class="mco-metrics"><div><b>6E3S</b><span>PT2385 · 79A</span></div><div><b>7W80</b><span>Belzutifan · 72Q</span></div><div><b>gem → vic</b><span>difluoro redesign</span></div><div><b>↓ UGT</b><span>phase-II metabolism</span></div></div>
  </div>

  <div class="mco-grid">
    <div class="mco-card"><div class="mco-title"><span>FIRST GENERATION</span><b>PT2385 · PDB 6E3S</b></div><div id="pt-view" class="mco-view"><div class="mco-load">Loading 6E3S…</div></div><div class="mco-controls"><button data-v="pt" data-mode="complex" class="active">Complex</button><button data-v="pt" data-mode="ligand">Ligand</button><button data-v="pt" data-mode="pocket">Pocket</button><button data-v="pt" data-mode="m252">M252</button></div></div>
    <div class="mco-card"><div class="mco-title"><span>SECOND GENERATION</span><b>PT2977 / Belzutifan · PDB 7W80</b></div><div id="bz2-view" class="mco-view"><div class="mco-load">Loading 7W80…</div></div><div class="mco-controls"><button data-v="bz" data-mode="complex" class="active">Complex</button><button data-v="bz" data-mode="ligand">Ligand</button><button data-v="bz" data-mode="pocket">Pocket</button><button data-v="bz" data-mode="m252">M252</button></div></div>
  </div>

  <div class="mco-switch">
    <button class="active" data-story="problem">1 · Clinical problem</button><button data-story="hypothesis">2 · Design hypothesis</button><button data-story="fluorine">3 · Fluorine move</button><button data-story="outcome">4 · Outcome</button>
  </div>
  <div id="mco-story" class="mco-story"><span>01 · CLINICAL PK</span><h3>The first molecule validated the target — but exposure was unreliable.</h3><p>PT2385 underwent extensive glucuronidation to an inactive metabolite. The key lesson is that excellent target biology does not rescue a molecule whose human exposure is capped by metabolism.</p><b>Optimization target:</b> reduce glucuronidation without sacrificing the benzylic alcohol required for potency.</div>
</div>

## 1. What actually failed?

The first-generation inhibitor **PT2385** showed clinical proof of concept, but its pharmacokinetics were variable and dose-limited because it was extensively converted to a glucuronide metabolite. The 2019 discovery paper identifies reduction of this phase-II metabolism as the central objective of the second-generation program.

This is an important distinction in lead optimization: **the target did not fail; the exposure profile did.**

## 2. Why not simply remove the alcohol?

Because a metabolically vulnerable group can simultaneously be a **binding pharmacophore**. Removing the benzylic OH would attack the glucuronidation handle, but could also remove a productive interaction and alter the preferred bound geometry.

So the medicinal-chemistry strategy became **preserve + detune**, rather than delete + replace.

<div class="mco-gate"><span>DESIGN GATE</span><h3>Can local stereoelectronics reduce conjugation while preserving binding?</h3><p>This is where fluorine becomes useful not as a generic metabolic blocker, but as a precisely positioned electronic and conformational perturbation.</p></div>

## 3. The key structural change: geminal → vicinal difluoro

PT2385 carries a **geminal difluoro** arrangement. In PT2977, the design moved to a **vicinal 2,3-difluoro** arrangement in an all-cis stereochemical series. The published optimization reports that this change produced the desired multi-parameter effect: **enhanced potency, lower lipophilicity, reduced phase-II metabolism, and substantially improved pharmacokinetics**.

<div class="mco-path">
  <div><span>PT2385</span><b>gem-difluoro + benzylic OH</b><small>potent, but high glucuronidation / variable exposure</small></div>
  <div class="mco-arrow">→</div>
  <div><span>DESIGN MOVE</span><b>vicinal F / F + defined stereochemistry</b><small>change electronics + conformation around the retained OH</small></div>
  <div class="mco-arrow">→</div>
  <div><span>PT2977</span><b>(1S,2S,3R) all-cis fluorohydrin</b><small>better potency / lipophilicity / PK balance</small></div>
</div>

## 4. What the enzyme kinetics added

The program did not infer glucuronidation only from microsomal disappearance. It synthesized glucuronide standards and quantified formation kinetics in **human intestinal microsomes** and **recombinant UGT2B17**, fitting substrate saturation data to Michaelis–Menten behavior.

That matters because the optimization target was mechanistic: **lower the rate of formation of the conjugated metabolite**.

For PT2385, the paper reports V~max~ values of **49 and 280 (pmol/min)/mg protein**, with K~m~ values of **20 and 11 µM** in human intestinal microsomes and recombinant UGT2B17, respectively.

## 5. The species translation was informative

The largest separation between PT2977 and PT2385 was not necessarily in rodents. The discovery paper reports **~9-fold higher dose-normalized AUC in dogs and ~20-fold higher in monkeys** for PT2977 versus PT2385, consistent with substantially less glucuronide formation in those higher species.

That supported the translational hypothesis: lowering glucuronidation should improve human exposure and reduce patient-to-patient variability.

## 6. Structural comparison: what should you look for?

The two viewers above are intentionally synchronized conceptually rather than geometrically. Use them to inspect:

- how both molecules occupy the same HIF-2α PAS-B cavity;
- retention of the alcohol-containing binding motif;
- the changed fluorination/stereochemical pattern around that motif;
- the proximity of **M252**, the allosteric relay residue implicated in antagonist-induced weakening of the HIF-2α–ARNT heterodimer.

The point is not that a dramatically different binding mode appeared. The medicinal-chemistry victory was achieving a **better systemic molecule while retaining the productive target-binding solution**.

## 7. The optimization lesson

<div class="grid cards" markdown>

-   :material-target:{ .lg .middle } **Name the real failure mode**

    ---
    PT2385 already had target engagement and clinical activity. The second-generation objective was pharmacokinetic reliability, not “more docking score.”

-   :material-flask:{ .lg .middle } **Measure the pathway you intend to fix**

    ---
    Direct glucuronide-formation kinetics provided a better mechanistic readout than relying only on generic metabolic-stability assays.

-   :material-molecule:{ .lg .middle } **Keep a necessary liability when needed**

    ---
    If the vulnerable functional group is pharmacophoric, redesign its local electronic environment before removing it.

-   :material-chart-bell-curve:{ .lg .middle } **Optimize the whole property vector**

    ---
    The successful move improved potency, lipophilicity, phase-II metabolism and PK simultaneously.

</div>

## Final perspective

PT2385 → PT2977 is a compact example of what senior medicinal chemistry looks like in practice. The key decision was not “which substituent improves potency?” It was **which molecular change preserves the binding solution while correcting the property that is limiting clinical exposure**.

That is why this case is so useful: the molecule did not need a new target, a new pocket, or a new mechanism. It needed a better **objective function**.

---

### Structural references

- **PDB 6E3S** — heterodimeric HIF-2 complex with PT2385; ligand **79A**.
- **PDB 7W80** — heterodimeric HIF-2 complex with belzutifan/PT2977; ligand **72Q**.
- Xu R. et al., *J. Med. Chem.* 2019, 62, 6876–6893 — discovery and optimization of PT2977.

<style>
.mco-shell{margin:1.2rem calc(50% - 50vw + 1rem) 2.2rem;max-width:calc(100vw - 2rem);border:1px solid var(--md-default-fg-color--lightest);border-radius:20px;overflow:hidden;background:var(--md-default-bg-color);box-shadow:0 20px 50px -30px rgba(14,116,144,.55)}.mco-head{display:grid;grid-template-columns:1.3fr 1fr;gap:1rem;padding:1.8rem 2rem;background:linear-gradient(125deg,#082f49,#0e7490 55%,#115e59);color:#fff}.mco-head h2{color:#fff!important;border:0!important;margin:.2rem 0 .6rem!important}.mco-head p{color:#cffafe;margin:0}.mco-kicker,.mco-title span,.mco-story>span,.mco-gate span,.mco-path span{font-size:.6rem;letter-spacing:.13em;font-weight:800;color:#fbbf24}.mco-metrics{display:grid;grid-template-columns:1fr 1fr;gap:.55rem}.mco-metrics div{padding:.75rem;border-radius:11px;background:rgba(255,255,255,.09);border:1px solid rgba(255,255,255,.18)}.mco-metrics b{display:block;font-size:1.05rem}.mco-metrics span{font-size:.63rem;color:#cffafe}.mco-grid{display:grid;grid-template-columns:1fr 1fr;background:#06131b;gap:1px}.mco-card{background:#091923}.mco-title{padding:.7rem .9rem;color:#fff}.mco-title b{display:block;font-size:.82rem}.mco-view{height:390px;position:relative;background:#061018}.mco-load{position:absolute;inset:0;display:grid;place-items:center;color:#a5f3fc;font-size:.75rem}.mco-controls{display:flex;gap:.3rem;padding:.6rem;background:#07141c}.mco-controls button,.mco-switch button{border:1px solid rgba(255,255,255,.12);background:#102630;color:#dff7ff;border-radius:8px;padding:.35rem .55rem;cursor:pointer;font-size:.63rem}.mco-controls button.active,.mco-controls button:hover,.mco-switch button.active{border-color:#22d3ee;background:#164e63}.mco-switch{display:flex;gap:.4rem;flex-wrap:wrap;padding:.8rem 1rem;background:#0b1c25}.mco-story{padding:1.1rem 1.2rem 1.3rem}.mco-story h3{margin:.25rem 0 .55rem}.mco-story p{margin:.3rem 0 .7rem}.mco-gate{margin:1.4rem 0;padding:1rem 1.1rem;border-left:4px solid #f59e0b;border-radius:10px;background:color-mix(in srgb,var(--md-default-bg-color) 91%,#f59e0b 9%)}.mco-gate h3{margin:.25rem 0}.mco-path{display:grid;grid-template-columns:1fr auto 1fr auto 1fr;gap:.7rem;align-items:stretch;margin:1.4rem 0}.mco-path>div:not(.mco-arrow){padding:1rem;border:1px solid var(--md-default-fg-color--lightest);border-radius:12px}.mco-path b,.mco-path small{display:block}.mco-path small{margin-top:.35rem;color:var(--md-default-fg-color--light)}.mco-arrow{display:grid;place-items:center;font-size:1.35rem;color:#0891b2}@media(max-width:850px){.mco-head,.mco-grid{grid-template-columns:1fr}.mco-view{height:340px}.mco-path{grid-template-columns:1fr}.mco-arrow{transform:rotate(90deg)}}
</style>
<script src="https://unpkg.com/ngl@2.0.0-dev.40/dist/ngl.js"></script>
<script>
(function(){
 const cfg={pt:{id:'6E3S',lig:'79A',el:'pt-view'},bz:{id:'7W80',lig:'72Q',el:'bz2-view'}}, stages={}, comps={};
 function show(key,mode){const c=comps[key];if(!c)return;c.removeAllRepresentations();const lig='['+cfg[key].lig+']';if(mode==='complex'){c.addRepresentation('cartoon',{sele:':A',color:'#22d3ee',opacity:.9});c.addRepresentation('cartoon',{sele:':B',color:'#fbbf24',opacity:.9});c.addRepresentation('ball+stick',{sele:lig,color:key==='pt'?'#fb923c':'#f472b6',radiusScale:1.25});c.autoView();}else if(mode==='ligand'){c.addRepresentation('cartoon',{sele:':B',color:'#64748b',opacity:.22});c.addRepresentation('ball+stick',{sele:lig,color:key==='pt'?'#fb923c':'#f472b6',radiusScale:1.65});c.autoView(lig,600);}else if(mode==='pocket'){c.addRepresentation('cartoon',{sele:':B',color:'#475569',opacity:.18});c.addRepresentation('ball+stick',{sele:lig,color:key==='pt'?'#fb923c':'#f472b6',radiusScale:1.3});c.addRepresentation('licorice',{sele:'(248:B or 252:B or 254:B or 280:B or 281:B or 293:B or 307:B or 309:B or 319:B)',color:'#fbbf24',radiusScale:.75});c.autoView(lig,600);}else{c.addRepresentation('cartoon',{sele:':A',color:'#22d3ee',opacity:.28});c.addRepresentation('cartoon',{sele:':B',color:'#fbbf24',opacity:.32});c.addRepresentation('ball+stick',{sele:lig,color:key==='pt'?'#fb923c':'#f472b6',radiusScale:1.2});c.addRepresentation('spacefill',{sele:'252:B',color:'#ef4444',radiusScale:.85});c.autoView('252:B',600);}}
 Object.keys(cfg).forEach(k=>{stages[k]=new NGL.Stage(cfg[k].el,{backgroundColor:'#061018'});stages[k].loadFile('https://files.rcsb.org/download/'+cfg[k].id+'.pdb',{defaultRepresentation:false}).then(c=>{comps[k]=c;document.querySelector('#'+cfg[k].el+' .mco-load').style.display='none';show(k,'complex');});});
 document.querySelectorAll('.mco-controls button').forEach(b=>b.addEventListener('click',function(){document.querySelectorAll('.mco-controls button[data-v="'+this.dataset.v+'"]').forEach(x=>x.classList.remove('active'));this.classList.add('active');show(this.dataset.v,this.dataset.mode);}));
 const stories={problem:['01 · CLINICAL PK','The first molecule validated the target — but exposure was unreliable.','PT2385 underwent extensive glucuronidation to an inactive metabolite. The key lesson is that excellent target biology does not rescue a molecule whose human exposure is capped by metabolism.','Optimization target: reduce glucuronidation without sacrificing the benzylic alcohol required for potency.'],hypothesis:['02 · MEDCHEM HYPOTHESIS','Preserve the OH. Detune its metabolic reactivity.','Instead of deleting a productive binding feature, the program changed the electronic and stereochemical environment around it.','The liability and the pharmacophore were the same atom-level region, so the solution had to be multiparameter.'],fluorine:['03 · STRUCTURAL EDIT','Move from geminal to vicinal difluoro substitution.','The vicinal fluorination pattern and defined all-cis stereochemistry changed local electronics, conformation and lipophilicity while retaining the alcohol.','Fluorine was used as a design tool, not merely as a metabolic blocker.'],outcome:['04 · PROPERTY VECTOR','One change improved several endpoints at once.','The reported result was enhanced potency, decreased lipophilicity, attenuated phase-II metabolism and substantially improved pharmacokinetics.','The success criterion was a better overall molecule — not optimization of one assay in isolation.']};
 document.querySelectorAll('.mco-switch button').forEach(b=>b.addEventListener('click',function(){document.querySelectorAll('.mco-switch button').forEach(x=>x.classList.remove('active'));this.classList.add('active');const s=stories[this.dataset.story];document.getElementById('mco-story').innerHTML='<span>'+s[0]+'</span><h3>'+s[1]+'</h3><p>'+s[2]+'</p><b>'+s[3]+'</b>';}));
 window.addEventListener('resize',()=>Object.values(stages).forEach(s=>s.handleResize()));
})();
</script>
