#!/bin/bash -l

BASE=`pwd`
INSTALL_DIR=${BASE}
SOFTWARE_DIR=${BASE}/software
VIRTUAL_ENV=${BASE}/.virtualenvs/km

# install km
python3 -m venv .virtualenvs/km
source .virtualenvs/km/bin/activate
pip install --upgrade pip setuptools wheel
pip install km-walk && echo "==> km installed"

# install jellyfish with python binding
mkdir -p ${SOFTWARE_DIR}
cd ${SOFTWARE_DIR}
curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz \
	--output ${SOFTWARE_DIR}/jellyfish-2.2.6.tar.gz
tar zxvf ${SOFTWARE_DIR}/jellyfish-2.2.6.tar.gz
cd ${SOFTWARE_DIR}/jellyfish-2.2.6
./configure --prefix=${INSTALL_DIR}/.virtualenvs/km --enable-python-binding
make -j 4 && make install && echo "==> Jellyfish installed"
