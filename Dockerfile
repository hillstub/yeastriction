# Use the official Python base image
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install ViennaRNA, Primer3, and Bowtie (version 1)
RUN curl -LO https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz \
    && tar -xzf ViennaRNA-2.4.17.tar.gz \
    && cd ViennaRNA-2.4.17 \
    && ./configure --without-swig \
    && make \
    && make install

RUN curl -LO https://github.com/primer3-org/primer3/archive/refs/tags/v2.5.0.tar.gz \
    && tar -xzf v2.5.0.tar.gz \
    && cd primer3-2.5.0/src \
    && make all \
    && cp primer3_core /usr/local/bin/primer3_core

RUN curl -LO http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip \
    && unzip bowtie-1.3.0-linux-x86_64.zip \
    && mv  bowtie-1.3.0-linux-x86_64 /opt/bowtie \
    && ln -s /opt/bowtie/bowtie* /usr/local/bin/

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
CMD ["poetry", "run", "gunicorn", "-b", "0.0.0.0:8050", "app:server"]