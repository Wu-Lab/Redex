# Configuration file for the Sphinx documentation builder.

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


project = "Redex"
copyright = "2026, Redex contributors"
author = "Redex contributors"
release = "v0.1.0"


extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_parser",
    "nbsphinx",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/yuanfanda/Redex",
    "use_repository_button": True,
}

html_static_path = []
html_show_sphinx = False

nbsphinx_execute = "never"
nbsphinx_allow_errors = True

autodoc_mock_imports = [
    "anndata",
    "dgl",
    "joblib",
    "matplotlib",
    "numpy",
    "ot",
    "pandas",
    "rpy2",
    "scanpy",
    "scipy",
    "sklearn",
    "torch",
    "tqdm",
]

