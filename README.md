
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

## Execute bulk performance test
Executing 
`python3 -m tests.test_bulkPerformance  && google-chrome performace_table.html.`
should generate a table called performace_table.html in the current folder.

# Development with docker
* You need to have docker installed

## Running the docker image for developers
cd into the MeMoMe folder
```bash
docker compose run dev # This will start the docker container interactively
```
You can now edit files locally with your favorite editor and them run MeMoMe in the container. In the container you can run all the development commands mentioned above.

## Running the docker image for production
Disclaimer: Currently the production image just runs the tets. The cli will be exposed later
cd into the MeMoMe folder
```bash
docker compose run prod # This will run the tests
```
