project = 'Kraken 2'
copyright = '2026, Derrick Wood'
author = 'Derrick Wood'

release = '2.18.0'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output
html_title = 'Kraken 2 Documentation'
# html_theme = 'sphinx_rtd_theme'
html_theme = "alabaster"
html_show_sourcelink = False

# -- Options for EPUB output
# epub_show_urls = 'footnote'
