# Use the jupyter/minimal-notebook as the base image
FROM jupyter/minimal-notebook

# Copy the environment.yml file into the container
COPY environment.yml .

# Activate the notebook environment
RUN conda env create -f environment.yml && \
    conda clean --all -f -y && \
    echo "conda activate notebook" > ~/.bashrc

# Install IPython kernel for Jupyter Lab
RUN /bin/bash -c "source activate notebook && python -m ipykernel install --user --name=notebook"

# Set the working directory
WORKDIR /home/jovyan/work

# Copy the notebook file into the container
COPY wedock_win.ipynb ./

# Copy the Jupyter configuration file to set the default kernel
COPY jupyter_notebook_config.py /home/jovyan/.jupyter/

# Set the password to 'woshibianfuxia'
ENV JUPYTER_PASSWORD=argon2:$argon2id$v=19$m=10240,t=10,p=8$/hFy3/7uNm3MYlwJN/nclg$L+e1t+O8SIN3tBd3N09NTw

# Set the default notebook to 'wedock_win.ipynb'
CMD ["start-notebook.sh", "--NotebookApp.notebook_dir=/home/jovyan/work", "--notebook", "wedock_win.ipynb", "--NotebookApp.password=${JUPYTER_PASSWORD}"]

# Set GitHub repository label
LABEL repository="https://github.com/quantaosun/labodock_binder"

