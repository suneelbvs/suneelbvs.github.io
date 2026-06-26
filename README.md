# suneelbvs.github.io

Personal site of **Suneel Kumar BVS** — Director, AI & Drug Design.
Built with [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).

## Local development

```bash
pip install -r requirements.txt
mkdocs serve        # live preview at http://127.0.0.1:8000
mkdocs build        # output to ./site
```

## Structure

```
docs/
  index.md            # home
  resume.md
  publications.md
  projects/           # project write-ups
  blog/posts/         # blog (Material blog plugin)
  images/             # all media
  stylesheets/extra.css   # teal + amber theme
mkdocs.yml            # site config, nav, palette, plugins
.github/workflows/    # build + deploy to GitHub Pages
```

## Deployment

Pushing to `main` triggers `.github/workflows/deploy.yml`, which builds the site
and publishes it to GitHub Pages.

> **Pages setting:** in repository **Settings → Pages**, set
> **Source → GitHub Actions** (the site is no longer built by the default
> Jekyll/branch deployment).
