# Installation

## Create a virtual environment

Create and activate a clean conda environment:

```bash
git clone https://github.com/Wu-Lab/Redex.git
cd Redex
conda create -n redex python=3.8
conda activate redex
```



## Install Redex

Install the optional R dependencies used by the default `mclust` clustering workflow:

```bash
conda install -c conda-forge r=4.1.0
conda install -c conda-forge r-mclust
```

Install PyTorch and DGL. The following commands install CUDA 11.8 builds:

```bash
pip install torch==2.3.1+cu118 torchvision==0.18.1+cu118 torchaudio==2.3.1+cu118 --index-url https://download.pytorch.org/whl/cu118
pip install dgl==1.1.2+cu118 -f https://data.dgl.ai/wheels/cu118/repo.html
```

Install the remaining Python dependencies and the editable local package:

```bash
pip install -r requirements.txt
pip install -e .
```

## Jupyter kernel

Install a notebook kernel for the environment:

```bash
pip install ipykernel
python -m ipykernel install --user --name redex --display-name redex
```

