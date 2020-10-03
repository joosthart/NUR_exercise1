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

echo "Download data for problem 1."
if [ ! -d "data/coolingtables_highres" ]; then
  # wget https://www.strw.leidenuniv.nl/WSS08/coolingtables_highres.tar.gz
  tar -xzf coolingtables_highres.tar.gz
  mv CoolingTables data/CoolingTables
fi

# Script that returns a plot
echo "Run the first script ..."
python3 problem1.py

# # Script that pipes output to a file
# echo "Run the second script ..."
# python3 helloworld.py > helloworld.txt

# # Script that saves data to a file
# echo "Run the third script ..."
# python3 cos.py

# # Script that generates movie frames
# echo "Run the fourth script ..."
# python3 sinemovie.py

# # code that makes a movie of the movie frames
# ffmpeg -framerate 25 -pattern_type glob -i "plots/snap*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 sinemovie.mp4

# echo "Generating the pdf"

# pdflatex template.tex
# bibtex template.aux
# pdflatex template.tex
# pdflatex template.tex


