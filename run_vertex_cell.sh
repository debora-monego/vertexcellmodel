#! /bin/bash

HOME=/Users/deboramonego/Documents/VertexCellModel-main
cd $HOME
g++ -std=c++11 -g src/main.cpp src/parser.cpp src/polygon.cpp src/geometry.cpp src/vector.cpp src/energy.cpp src/force.cpp src/transitionT1.cpp -o vertexcell

N=100

# Parameters
lx=10
ly=10

ka=1                   # area force coefficient
A0=1                    # current preferred area for polygon
P0=1                   # preferred perimeter for polygon
J=1                     # J value

gamma_f=0.04           # multiplication factor for gamma (gamma = gamma_f * ka * A0)
Lambda_f=0.12          # multiplication factor for Lambda (Lambda = Lambda_f * ka * sqrt(pow(A0, 3)))

lmin=0.04              # edge rearangment treshold
ksep=1.5                # new separation after transition (* original separation)

xi=0.2                  # motility coefficient
eta=0.1                 # noise scalling coefficient

delta_t=0.01            # timestep
T=0.05                    # total simulation time

T1_enabled=true         # enable T1 transitions

#for ka in $(seq 0.1 0.1 1) $(seq 1 1 10)
for ka in 1
do
    echo "Running ka = ${ka}..."
    folder="results/${N}cells_${ka}ka_T1"
    mkdir -p "$folder"
    ./vertexcell periodic_vertices.txt periodic_edges.txt periodic_polygons.txt "$lx" "$ly" "$ka" "$A0" "$P0" "$J" "$lmin" "$ksep" "$xi" "$eta" "$delta_t" "$T" "$T1_enabled" "$gamma_f" "$Lambda_f" "$folder"
done