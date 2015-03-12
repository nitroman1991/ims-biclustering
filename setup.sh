#!/bin/bash

echo -e "Starting installation of required c++ libraries...\n\n"
sudo apt-get -y install libarpack2-dev libsuperlu3-dev gfortran libatlas-dev libblas-dev libboost-system-dev libboost-date-time-dev libboost-filesystem-dev
if [ $? -ne 0 ]; then
	echo -e "Error with installing c++ libraries, aborting..."
	return 1
fi
echo -e "Starting installation required gcc...\n\n"
sudo apt-get -y install gcc g++ g++-4.7
if [ $? -ne 0 ]; then
	echo -e "Error with installing gcc, aborting..."
	return 1
fi
echo -e "Starting installation required LaTex libraries(for report generating)...\n\n"
sudo apt-get -y install texlive-latex-base texlive-latex-extra dvipng libfreetype6-dev libpng-dev
if [ $? -ne 0 ]; then
	echo -e "Error with installing LaTex libraries, aborting..."
	return 1
fi
echo -e "Starting installation python...\n\n"
sudo apt-get -y install python python-pip && sudo easy_install distribute & sudo apt-get install python-dev
if [ $? -ne 0 ]; then
	echo -e "Error with installing c++ libraries, aborting..."
	return 1
fi
echo -e "Starting installation required python modules...\n\n"
sudo pip install matplotlib numpy scikit-learn statsmodels
if [ $? -ne 0 ]; then
	echo -e "Error with installing python modules, aborting..."
	return 1
fi
echo -e "Finished installation!\n\n\n"

cd ./
sudo mkdir ./data
sudo chmod 777  -R ./

cd ./src
make


