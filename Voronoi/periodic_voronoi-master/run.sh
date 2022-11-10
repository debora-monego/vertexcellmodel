#!/usr/bin/bash


HOME=/Users/deboramonego/Documents/Voronoi/periodic_voronoi-master

N=$1
echo "N = "$N

folder="voronoi"
# mkdir $folder
cd $folder

# run voronoi diagram without periodic boundaries
# python $HOME/voronoi.py $N

# plot voronoi diagram
# python $HOME/plot.py "voronoi.jpg"


# run voronoi diagram with periodic boundaries
python $HOME/periodic_voronoi.py $N

# plot voronoi diagram
python $HOME/periodic_plot.py "periodic_voronoi_poly_${N}cells.jpg" 



