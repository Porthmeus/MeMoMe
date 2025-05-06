FROM continuumio/miniconda3:24.9.2-0 as conda

RUN conda install -n base -c conda-forge mamba

# Set working directory
WORKDIR /app

# Copy project files (optional)
COPY requirements.pinned.yml .

RUN mamba env create --file=requirements.pinned.yml
shell ["/bin/bash", "-c"]
run echo "source activate MeMoMe" >> /root/.bashrc && \
    source /root/.bashrc

#___________________________________________#

# Use a Miniconda base image
FROM conda as dev

# Install Mamba
RUN apt-get update && apt-get install vim -y    

#___________________________________________#

FROM conda  as prod

COPY . .
SHELL ["/bin/bash", "-c", "-l"]
ENTRYPOINT python -m unittest


