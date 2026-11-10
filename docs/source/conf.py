# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re
from pathlib import Path

DOCS_SOURCE = Path(__file__).resolve().parent
REPO_ROOT = DOCS_SOURCE.parents[1]

# -- Project information -----------------------------------------------------

project = 'SOS-ABS'
copyright = "2023, Centre National d'Etudes Spatiales (CNES)"
author = 'CNES, CS GROUP'
# version from info_Version.txt ("SOS-ABS version 5.1")
_m = re.search(r'([0-9]+(\.[0-9]+)+)', (REPO_ROOT / 'info_Version.txt').read_text())
release = _m.group(1) if _m else ''
version = release
today_fmt = '%Y-%m-%d'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx_copybutton',
    'sphinxcontrib.mermaid',
    # PDF/ePub version of the mermaid flowchart (algorithm.rst); needs graphviz
    'sphinx.ext.graphviz',
    # MyST markdown + Jupyter notebooks (myst_nb includes myst_parser)
    'myst_nb',
]

# the notebooks are stored with their outputs: SOS-ABS (compiled code and
# CKD data) is not available on Read the Docs
nb_execution_mode = 'off'

templates_path = ['_templates']
exclude_patterns = ['_build', '_readme.md', 'Thumbs.db', '.DS_Store',
                    '**.ipynb_checkpoints']

todo_include_todos = True

# the README included in index.rst starts at "##" once its title is skipped
suppress_warnings = ['myst.header']

# -- Math --------------------------------------------------------------------

# number labelled equations and refer to them as "Eq. (n)" with :eq:
math_eqref_format = 'Eq. ({number})'
math_numfig = True
numfig = True

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
pygments_style = 'sphinx'

html_theme_options = {
    'repository_url': 'https://github.com/CNES/RadiativeTransferCode-SOS',
    'repository_branch': 'master',
    'path_to_docs': 'docs/source',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_edit_page_button': True,
    'use_download_button': True,
    'navigation_with_keys': True,
    'show_toc_level': 2,
    'secondary_sidebar_items': ['page-toc', 'edit-this-page'],
}

html_title = 'SOS-ABS'
html_static_path = []

html_show_sourcelink = False
html_last_updated_fmt = today_fmt

htmlhelp_basename = 'sos_abs_doc'

# -- Options for LaTeX/PDF output --------------------------------------------

# xelatex handles the accented and Greek characters; the DejaVu fonts cover
# them and are available on most Linux systems (incl. Read the Docs)
latex_engine = 'xelatex'
latex_use_xindy = False
latex_elements = {
    'fontpkg': r'''
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
''',
}

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = [
    'amsmath',
    'colon_fence',
    'deflist',
    'dollarmath',
    'linkify',
    'replacements',
    'smartquotes',
]
