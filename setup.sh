#!/bin/bash

echo -e "Starting installation of required c++ libraries...\n\n"
sudo apt-get install libarpack2-dev libsuperlu3-dev gfortran libatlas-dev libblas-dev libboost-system-dev libboost-date-time-dev libboost-filesystem-dev
echo -e "Starting installation required gcc...\n\n"
sudo apt-get install gcc g++ g++-4.7
echo -e "Starting installation required LaTex libraries(for report generating)...\n\n"
sudo apt-get install texlive-latex-base texlive-latex-extra dvipng libfreetype6-dev libpng-dev
echo -e "Starting installation python...\n\n"
sudo apt-get install python python-pip && sudo easy_install distribute & sudo apt-get install python-dev
echo -e "Starting installation required python modules...\n\n"
sudo pip install matplotlib numpy scikit-learn statsmodels
echo -e "Finished installation!\n\n\n"

cd ./
sudo chmod 777  -R ./

cd ./src
make


