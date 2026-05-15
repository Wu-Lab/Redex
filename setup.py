from pathlib import Path

from setuptools import setup


ROOT = Path(__file__).parent


def read_requirements():
    requirements_file = ROOT / "requirements.txt"
    if not requirements_file.exists():
        return []
    return [
        line.strip()
        for line in requirements_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


setup(
    name="Redex",
    version="0.1.0",
    description="Spatially variable feature identification based on regional differential expression.",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8")
    if (ROOT / "README.md").exists()
    else "",
    long_description_content_type="text/markdown",
    author="Redex contributors",
    url="<YOUR_REPOSITORY_URL>",
    license="GPL-3.0",
    python_requires=">=3.8",
    package_dir={"": "src"},
    py_modules=["Redex", "method", "preprocess", "CCAmodel", "aug", "utils"],
    install_requires=read_requirements(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
