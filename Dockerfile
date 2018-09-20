
# Start from an openSuSE image to mirror the development platform
FROM opensuse:latest

# Install (using the suse package manager) all system requirements to
# deploy the OS and python framework, including all its library dependencies.
RUN zypper install --no-confirm \
    python3 python3-devel \
    python3-pip \
    python3-setuptools \
    gsl gsl-devel \
    wget \
    gcc \
    git \
    gcc-fortran \
    freetype2-devel \
    libpng12-devel \
    && \
    \
    zypper install -d --no-confirm \
           python3-numpy python3-numpy-devel \
           python3-scipy \
           python3-matplotlib \
           python3-Django \
           python3-Pillow \
           python3-Cython



    # python3-numpy python3-numpy-devel \
    # python3-scipy \
    # python3-matplotlib \
    # python3-Django \
    # python3-Pillow \
    # python3-Cython \

# Define where the software lives, ...
WORKDIR /app

# ... and install the github-provided code into the docker container
ADD . /app

# Install all remaining python packages as specified in the requirements file
RUN pip3 install --upgrade pip && pip3 install --ignore-installed -r requirements.txt

# Compile the podi_cython module
RUN cd /app && python3 setup.py build_ext --inplace

# install astromatic software (Sextractor & Swarp)
RUN cd /tmp && \
    wget -nv https://www.astromatic.net/download/sextractor/sextractor-2.19.5-1.x86_64.rpm && \
    wget -nv https://www.astromatic.net/download/swarp/swarp-2.38.0-1.x86_64.rpm && \
    zypper --non-interactive --no-gpg-checks install sextractor-2.19.5-1.x86_64.rpm swarp-2.38.0-1.x86_64.rpm && \
    rm sextractor-2.19.5-1.x86_64.rpm swarp-2.38.0-1.x86_64.rpm

# Define the volumes we will need to run QuickReduce
VOLUME [\
    "/input", \
    "/output", \
    "/cals", \
    "/catalogs", \
    "/qr", \
    "/scratch", \
    "/mastercals" \
    ]



