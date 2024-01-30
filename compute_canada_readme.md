# Python package installation

Run following commands to create virtual environment and activate it before installing python libraries:
```
module load python/3.10
virtualenv --no-download /project/[group_name]/[user_name]/nest_venv
source /project/[group_name]/[user_name]/fatema_venv/bin/activate
```
Please replace [group_name] and [user_name] with your assigned group and name.

Upgrade pip3:
```
pip3 install --no-index --upgrade pip
```

Install required libraries as follows:  
```
pip3 install ipython
pip3 install pandas
pip3 install altair vega_datasets
pip3 install flake8 pytest jinja2 sphinx m2r docutils
pip3 install collection
pip3 install stlearn
pip3 install scanpy
pip3 install qnorm
pip3 install csv pickle gzip matplotlib scipy sklearn 
```

Next, see which pytorch and CUDA versions are available:
```
avail_wheels "torch*"
module load cuda [press tab to see the available versions]
```

Next, based on the Pytorch and CUDA version, run with wheels links for installing. Here we are using pytorch 1.13.1 with CUDA 11.7:

```
pip3 install torch-scatter -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-sparse -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-cluster -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-geometric -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
```

# Submit GPU jobs

[to be]
