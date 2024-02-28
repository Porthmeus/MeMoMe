
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


## Execute bulk performance test
Executing 
`python3 -m tests.test_bulkPerformance  && google-chrome table.html.`
should generate a table called performace_table.html in the current folder.
