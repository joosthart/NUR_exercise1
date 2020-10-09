#!/bin/bash

echo "Running code for Exercise 1."

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist create it!"
  mkdir plots
fi

echo "Creating the data directory if it does not exist"
if [ ! -d "data" ]; then
  echo "Directory does not exist create it!"
  mkdir data
fi

echo "Creating the output directory if it does not exist"
if [ ! -d "output" ]; then
  echo "Directory does not exist create it!"
  mkdir output
fi

echo "Download data for problem 1."
if [ ! -d "data/CoolingTables" ]; then
  wget https://www.strw.leidenuniv.nl/WSS08/coolingtables_highres.tar.gz
  tar -xzf coolingtables_highres.tar.gz
  rm coolingtables_highres.tar.gz
  mv CoolingTables data/CoolingTables
fi


echo "Creating virtual environment and install dependencies."
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt

echo "Run the first problem ..."
python3 code/problem1.py

# Video
if [ ! -f "coolingrate.mp4" ]; then
  echo "Combining images into .mp4"
  ffmpeg -framerate 25 -pattern_type glob -i "plots/coolingrate_z*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 coolingrate.mp4 -y
fi

echo "Download data for problem 2."
if [ ! -f "data/wgs.dat" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/wgs.dat
  mv wgs.dat data/wgs.dat
fi

if [ ! -f "data/wss.dat" ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/wss.dat
  mv wss.dat data/wss.dat
fi

echo "Run the second problem ..."
python3 code/problem2.py

echo "Run the third problem ..."
python3 code/problem3.py

echo "Generating the pdf"
pdflatex -interaction=batchmode -jobname=solutions main.tex

echo "All done! Results saved to solutions.pdf."



