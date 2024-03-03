
## Mamba
I strongly advice you to use Mamba instead of Conda.

Link to [mamba](Mamba.readthedocs.io/en/latest/installation.html)

## Create the codna env.

`conda env create --file=requirements.yml` 
or 

`mamba env create --file=requirements.yml` 

if you just want to update a the environment use

`mamba env update -f requirements.yml`


## Run Tests
python -m unittest



## Run application with test model
python3 main.py --model1 tests/dat/e_coli_core.xml  --model2 tests/dat/e_coli_core.xml --output here.csv

