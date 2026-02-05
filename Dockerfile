# Use Miniconda as base image
FROM continuumio/miniconda3:latest

# Metadata
LABEL maintainer="yhow" \
      description="Somatic Variant Calling Pipeline: From FastQ to Visualization"

# Set working directory
WORKDIR /pipeline

# Step 1: Install dependencies via Conda
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Set shell to execute subsequent commands within the gatk4 environment
SHELL ["conda", "run", "-n", "gatk4", "/bin/bash", "-c"]

# Step 2: Copy pipeline logic
COPY run_pipeline.sh .
COPY scr/ ./scr/

# Step 3: Grant execution permissions
RUN chmod +x run_pipeline.sh scr/*.sh

# Step 4: Create mount points for data persistence
RUN mkdir -p reads resources results aligned_reads plots tmp

# Step 5: Define container entrypoint
# Executes the master script within the gatk4 environment context
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "gatk4", "./run_pipeline.sh"]
