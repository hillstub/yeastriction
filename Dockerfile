# Use the official Python base image
FROM python:3.9-slim

# Set environment variables
ENV POETRY_VERSION=1.1.13
ENV PATH="/root/.poetry/bin:$PATH"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    unzip \
    && rm -rf /var/lib/apt/lists/*
# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python -

# Install Bowtie (version 1) and ViennaRNA
RUN wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip \
    && unzip bowtie-1.3.0-linux-x86_64.zip \
    && mv  bowtie-1.3.0-linux-x86_64 /opt/bowtie \
    && ln -s /opt/bowtie/bowtie /usr/local/bin/bowtie

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz \
    && tar -xzf ViennaRNA-2.4.17.tar.gz \
    && cd ViennaRNA-2.4.17 \
    && ./configure \
    && make \
    && make install

RUN wget https://github.com/primer3-org/primer3/archive/refs/tags/v2.5.0.tar.gz \
    && tar -xzf v2.5.0.tar.gz \
    && cd primer3-2.5.0/src \
    && make all \
    && cp primer3_core /usr/local/bin/primer3_core

# Set work directory
WORKDIR /app

RUN pip install poetry==1.4.2

# Copy the poetry.lock and pyproject.toml files
COPY poetry.lock pyproject.toml /app/

# Install project dependencies
RUN poetry install --no-root

# Copy the rest of the application code
COPY . /app

# Expose the port the app runs on
EXPOSE 8050

# Run the application
CMD ["poetry", "run", "python", "app.py"]
