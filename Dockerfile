# Use a Miniconda base image
FROM continuumio/miniconda3:24.9.2-0 as dev

RUN apt-get update && apt-get install vim -y    

# Install Mamba
RUN conda install -n base -c conda-forge mamba

# Set working directory
WORKDIR /app

# Copy project files (optional)
COPY requirements.pinned.yml .

RUN mamba env create --file=requirements.pinned.yml
SHELL ["conda", "activate", "MeMoMe"]


FROM continuumio/miniconda3:24.9.2-0 as prod

RUN apt-get update && apt-get install vim -y    

# Install Mamba
RUN conda install -n base -c conda-forge mamba

# Set working directory
WORKDIR /app

# Copy project files (optional)
COPY . .

RUN mamba env create --file=requirements.pinned.yml
SHELL ["conda", "activate", "MeMoMe"]
RUN python -m unittest


