---
layout: home
title: "Suneel Kumar BVS"
permalink: /
excerpt: "Director of AI & Drug Design — building machine learning that accelerates drug discovery."
author_profile: false
redirect_from:
  - /about/
  - /about.html
---

<!-- ============================= HERO ============================= -->
<section class="hero">
  <div class="hero__inner">
    <div class="hero__copy">
      <span class="hero__eyebrow"><span class="dot"></span> AI &middot; Cheminformatics &middot; Drug Discovery</span>
      <h1 class="hero__title">Designing molecules with <span class="grad">artificial intelligence</span>.</h1>
      <p class="hero__subtitle">
        I'm Suneel Kumar BVS, Director of AI &amp; Drug Design at Molecular Forecaster.
        I build machine&nbsp;learning and deep&nbsp;learning systems that turn complex chemistry
        into faster, smarter decisions across the drug discovery pipeline.
      </p>
      <div class="hero__cta">
        <a class="btn-accent" href="/year-archive/">Read the blog</a>
        <a class="btn-ghost" href="/Resume/">View résumé</a>
      </div>
      <div class="hero__social">
        {% if site.author.github %}<a href="https://github.com/{{ site.author.github }}" aria-label="GitHub" title="GitHub"><i class="fab fa-github"></i></a>{% endif %}
        {% if site.author.linkedin %}<a href="https://www.linkedin.com/in/{{ site.author.linkedin }}" aria-label="LinkedIn" title="LinkedIn"><i class="fab fa-linkedin-in"></i></a>{% endif %}
        {% if site.author.twitter %}<a href="https://twitter.com/{{ site.author.twitter }}" aria-label="Twitter" title="Twitter"><i class="fab fa-twitter"></i></a>{% endif %}
        {% if site.author.googlescholar %}<a href="{{ site.author.googlescholar }}" aria-label="Google Scholar" title="Google Scholar"><i class="fas fa-graduation-cap"></i></a>{% endif %}
        {% if site.author.email %}<a href="mailto:{{ site.author.email }}" aria-label="Email" title="Email"><i class="fas fa-envelope"></i></a>{% endif %}
      </div>
    </div>
    <aside class="hero__card">
      <img class="hero__avatar" src="/images/{{ site.author.avatar }}" alt="{{ site.author.name }}">
      <h3>{{ site.author.name }}</h3>
      <p class="role">{{ site.author.bio }}</p>
      <ul class="hero__meta">
        <li><i class="fas fa-briefcase"></i> {{ site.author.employer }}</li>
        <li><i class="fas fa-map-marker-alt"></i> {{ site.author.location }}</li>
        <li><i class="fas fa-flask"></i> Generative &amp; explainable AI for chemistry</li>
      </ul>
    </aside>
  </div>
</section>

<!-- ============================= FOCUS ============================= -->
<section class="section">
  <div class="wrap">
    <div class="section__head">
      <span class="section__eyebrow">What I work on</span>
      <h2 class="section__title">Research &amp; focus areas</h2>
      <p class="section__lead">Bridging deep learning and medicinal chemistry to design better molecules, faster.</p>
    </div>
    <div class="card-grid">
      <div class="focus-card">
        <div class="focus-card__icon"><i class="fas fa-atom"></i></div>
        <h3>AI-driven molecular design</h3>
        <p>Generative models (RNN-LSTM, transformers) for de&nbsp;novo design, fragment expansion and hit optimization.</p>
      </div>
      <div class="focus-card">
        <div class="focus-card__icon"><i class="fas fa-lightbulb"></i></div>
        <h3>Explainable AI</h3>
        <p>Interpretable models that tell chemists <em>why</em> a prediction was made — trust you can act on.</p>
      </div>
      <div class="focus-card">
        <div class="focus-card__icon"><i class="fas fa-cube"></i></div>
        <h3>Structure-based discovery</h3>
        <p>Deep learning fused with docking and free-energy methods for structure-based drug discovery.</p>
      </div>
      <div class="focus-card">
        <div class="focus-card__icon"><i class="fas fa-shield-virus"></i></div>
        <h3>Toxicity forecasting</h3>
        <p>Data-driven AI/ML toxicity models — from theory to reproducible, practical pipelines.</p>
      </div>
    </div>
  </div>
</section>

<!-- ============================= STATS ============================= -->
<section class="section section--alt">
  <div class="wrap">
    <div class="stats">
      <div class="stat"><div class="stat__num">10+</div><div class="stat__label">Years in AI &amp; cheminformatics</div></div>
      <div class="stat"><div class="stat__num">91</div><div class="stat__label">AI-designed candidates proposed (CACHE)</div></div>
      <div class="stat"><div class="stat__num">764&nbsp;nM</div><div class="stat__label">IC<sub>50</sub> of lead FLT-3 inhibitor</div></div>
      <div class="stat"><div class="stat__num">1</div><div class="stat__label">Book chapter on AI toxicity models</div></div>
    </div>
  </div>
</section>

<!-- ============================= WORK ============================= -->
<section class="section">
  <div class="wrap">
    <div class="section__head">
      <span class="section__eyebrow">Selected projects</span>
      <h2 class="section__title">Featured work</h2>
      <p class="section__lead">A few recent projects at the intersection of machine learning and drug design.</p>
    </div>
    <div class="work-grid">
      <article class="work-card">
        <div class="work-card__media"><img src="/images/AI_ML%20Models_for_Toxicity_Forecasts.jpg" alt="Toxicity Forecasts book chapter"></div>
        <div class="work-card__body">
          <span class="work-card__tag">Book chapter</span>
          <h3>Toxicity Forecasts: Data-Driven AI/ML Models</h3>
          <p>A practical guide to AI/ML toxicity modeling with RDKit, DeepChem and scikit-learn — from theory to working code. Published by Apple Academic Press (Taylor &amp; Francis).</p>
          <a class="work-card__link" href="https://www.appleacademicpress.com/artificial-intelligence-for-chemical-sciences-concepts-models-and-applications-/9781774918326">Read more</a>
        </div>
      </article>
      <article class="work-card">
        <div class="work-card__media"><img src="/images/LSTM_Fragment_based_drug_design.png" alt="Generative fragment-based drug design"></div>
        <div class="work-card__body">
          <span class="work-card__tag">CACHE Challenge</span>
          <h3>Generative models for fragment-based design</h3>
          <p>An RNN-LSTM approach for fragment expansion against the SARS-CoV-2 NSP3 Mac1 domain — 91 novel candidates proposed and tested in an open-science effort.</p>
          <a class="work-card__link" href="https://chemrxiv.org/engage/chemrxiv/article-details/65c6e60b66c1381729521e8f">Read more</a>
        </div>
      </article>
      <article class="work-card">
        <div class="work-card__media"><img src="/images/AI_Driven_Molecular_Design_Workflow.png" alt="AI-driven molecular design workflow"></div>
        <div class="work-card__body">
          <span class="work-card__tag">Published</span>
          <h3>De novo design of a selective FLT-3 inhibitor</h3>
          <p>AI-assisted de&nbsp;novo design produced PCW-A1001, a highly selective FLT-3 (D835Y) inhibitor with an IC<sub>50</sub> of 764&nbsp;nM. Published in Frontiers in Molecular Biosciences.</p>
          <a class="work-card__link" href="https://www.frontiersin.org/articles/10.3389/fmolb.2022.1072028/full">Read more</a>
        </div>
      </article>
    </div>
  </div>
</section>

<!-- ============================= LATEST POSTS ============================= -->
<section class="section section--alt">
  <div class="wrap">
    <div class="section__head">
      <span class="section__eyebrow">From the blog</span>
      <h2 class="section__title">Latest writing</h2>
      <p class="section__lead">Notebooks, cheatsheets and notes on cheminformatics &amp; machine learning.</p>
    </div>
    <div class="posts-grid">
      {% for post in site.posts limit:3 %}
      <article class="post-card">
        {% if post.header.teaser %}
        <a href="{{ post.url }}" class="post-card__cover">
          <img src="{% if post.header.teaser contains '://' %}{{ post.header.teaser }}{% else %}/images/{{ post.header.teaser }}{% endif %}" alt="{{ post.title | escape }}">
        </a>
        {% else %}
        <a href="{{ post.url }}" class="post-card__cover post-card__cover--gen{% cycle '', ' alt' %}">
          <span class="ribbon">{{ post.tags | first | default: 'Notes' }}</span>
        </a>
        {% endif %}
        <div class="post-card__body">
          <div class="post-card__meta"><span>{{ post.date | date: "%b %-d, %Y" }}</span></div>
          <h3 class="post-card__title"><a href="{{ post.url }}">{{ post.title }}</a></h3>
          <div class="post-card__tags">
            {% for tag in post.tags limit:3 %}<span class="tag">{{ tag }}</span>{% endfor %}
          </div>
        </div>
      </article>
      {% endfor %}
    </div>
    <p style="margin-top:2rem;">
      <a class="btn-ghost--ink btn-ghost" href="/year-archive/">Browse all posts</a>
    </p>
  </div>
</section>

<!-- ============================= CTA ============================= -->
<section class="section">
  <div class="wrap">
    <div class="cta-band">
      <h2>Let's build something at the AI &times; chemistry frontier.</h2>
      <p>Open to collaborations on generative design, explainable AI and structure-based discovery.</p>
      <a class="btn-accent" href="mailto:{{ site.author.email }}">Get in touch</a>
    </div>
  </div>
</section>
