#!/bin/bash

sudo apt-get install libarpack2-dev libsuperlu3-dev gfortran libatlas-dev libblas-dev libboost-system-dev libboost-date-time-dev libboost-filesystem-dev

sudo apt-get install gcc g++ g++-4.7
sudo apt-get install texlive dvipng libfreetype6-dev libpng-dev
sudo apt-get install python-pip
sudo pip install matplotlib numpy scikit-learn statsmodels

cd ./src
make


