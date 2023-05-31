#!/usr/bin/bash

HOME=/Users/deboramonego/Documents/VertexCellModel/Voronoi

N=$1
lx=$2
ly=$3

echo "N = "$N "box size = ["$lx","$ly"]"

folder="initial_conf"
# mkdir $folder
cd $folder

# run voronoi diagram without periodic boundaries
# python $HOME/voronoi.py $N

# plot voronoi diagram
# python $HOME/plot.py "voronoi.jpg"


# run voronoi diagram with periodic boundaries
python3 $HOME/periodic_voronoi.py $N $lx $ly

# plot voronoi diagram
python3 $HOME/periodic_plot.py "periodic_voronoi_poly_${N}cells.jpg" $lx $ly
