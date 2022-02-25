#!/bin/bash -l

BASE=`pwd`
INSTALL_DIR=${BASE}
SOFTWARE_DIR=${BASE}/software
VIRTUAL_ENV=${BASE}/.virtualenvs/km

# Make virtual env
python -m venv ${BASE}/.virtualenvs/km

# Make software dir
mkdir -p ${SOFTWARE_DIR}

# Activate venv
source ${INSTALL_DIR}/.virtualenvs/km/bin/activate

## Download and install jellyfish
cd ${SOFTWARE_DIR}
curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz \
    --output ${SOFTWARE_DIR}/jellyfish-2.2.6.tar.gz
tar zxvf ${SOFTWARE_DIR}/jellyfish-2.2.6.tar.gz
cd ${SOFTWARE_DIR}/jellyfish-2.2.6

# General case
./configure --prefix=${INSTALL_DIR}/.virtualenvs/km --enable-python-binding
make -j 4 && make install && echo "==> Jellyfish installed"

## Download and install km
cd $INSTALL_DIR/software

# Clone km and setup venv
git clone https://github.com/iric-soft/km.git && cd km

# crucial for earlier versions of python (e.g. 3.5)
$VIRTUAL_ENV/bin/pip install pip setuptools --upgrade

# wheel package is necessary now that setuptools runs pip for dependecies
pip install wheel
python setup.py install && echo "==> km installed"

## Execute km on a small example
# Need to reload the virtual environment each time you open a new terminal
# with: source $INSTALL_DIR/.virtualenvs/km/bin/activate
echo "### Run km test ... ###"
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $INSTALL_DIR/software/km/data/jf/02H025_NPM1.jf | km find_report -t $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $INSTALL_DIR/software/km/data/jf/02H025_NPM1.jf -g
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa $INSTALL_DIR/software/km/data/jf/03H116_ITD.jf | km find_report -t $INSTALL_DIR/software/km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa
