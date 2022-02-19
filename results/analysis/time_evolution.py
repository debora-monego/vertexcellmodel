# Written by Lucas Kampman
# 2019/07/24

# gif_generator.py: This Python 3.4 script reads in the lattice files produced by Alex Devanny's
# implementation of a 2D 2Type Potts simulation ata given interval (usually 1). It then
# processes them into numpy arrays, assigning cells to either red or green based on cell
# type, where the exact shade is determined by the cell's ID number. An animation is made
# from images generated from these arrays.

# Depending on the values of save_animation and show_animation, the animation can be saved to
# an .mp4 file and/or shown to the user, respectively.

# Code snippets have been stolen from various places, mostly the MatPlotLib documentation.

# Command line execution format:

# python[3] gif_generator.py OUTPUT_DIRECTORY_PATH NUMBER_OF_TIMESTEPS START_TIMESTEP

# Dependencies
# Python 3.7.4
# numpy==1.16.4
# Pillow==6.1.0
# matplotlib==2.2.4

import os
import time
from os import listdir
from os.path import isfile, join
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import numpy as np
from PIL import Image


# The index of the first output lattice file used to generate a movie
start_step = 0

# The index of the final output lattice file used to generate a movie
end_step = 1000

# The length of the edge of the square lattice
lattice_edge = 150

# The number of cells in the simulation
# (only used to scale color values so there's a good range)
num_cells = 30

# The interval between timesteps that are printed in lattice files (e.g. if printing every 10 timesteps)
step_interval = 25

# The number of types of cells in the simulation
# Currently only works with 0 or 1, but could be extended by modifying the
# lattice_to_frame() function appropriately.
num_cell_types = 1

# Set to true to save a .mp4 file of the animation
save_animation = True

# Set to true to show the animation in a matplotlib window
show_animation = False

# Set to true to print progress messages to the terminal
verbose = True

# Array to store data read in from simulation output lattice files
sim_array = []

# Number of frames to skip when generating an animation
# (e.g. skip = 10 will skip 9 frames for every 1 that's animated
skip = 1000

sim_array = []

print(os.path.dirname('.'))
data_path = '../.'

def main():
    read_args()
    read_lattice_files()
    # lattices_to_animation()
    lattices_to_FuncAnimation()

def read_args():
    if (len(sys.argv) > 0):
        data_path = sys.argv[0]
        if (len(sys.argv) > 1):
            end_step = sys.argv[1]
            if (len(sys.argv) > 2):
                start_step = sys.argv[2]
                end_step = sys.argv[1] + sys.argv[2]
                if (len(sys.argv) > 3):
                    print("Too many arguments")

def read_lattice_files():
    filesInDirectory = [f for f in listdir(data_path) if isfile(join(data_path, f)) and str(f).__contains__('lattice_') and str(f).__contains__('.txt')]
    ints_in_filenames = [[int(s) for s in str(f).split("_") if s.isdigit()] for f in filesInDirectory]
    step_to_file = [[ints_in_filenames[i], filesInDirectory[i]] for i in range(len(filesInDirectory))]

    
    # filesInDirectory.sort()
    step_to_file.sort()

    ordered_trimmed_files = [step_to_file[i][1] for i in range(1, len(step_to_file))]
    if (verbose):
        for arr in ordered_trimmed_files:
            print(arr)

    global end_step

    if (len(ordered_trimmed_files) > end_step):
        del ordered_trimmed_files[end_step:]
        if (verbose):
            print("Truncated to ",end_step,"output files.")
    else:
        end_step = len(ordered_trimmed_files)

    for file in ordered_trimmed_files:
        with open(data_path+'/'+file, mode='r') as curr_file:
            if (verbose):
                print("reading in file: {}".format(str(curr_file)))

            curr_step_arr = []
            line = curr_file.readline()
            line_count = 0
            while(line):
                curr_step_arr.append(list(map(int, line.split())))
                line = curr_file.readline()
            sim_array.append(curr_step_arr)
    if (verbose):
        print('Lattice files read in successfully.')

# Image generation

def lattice_to_frame(step_index, _ax):
    # Data for plotting
    if (verbose):
        print('Generating RGB array for timestep {0}\r'.format(step_index)),

    ID = 0
    frame = np.zeros((lattice_edge, lattice_edge, 3))
    for y in range(lattice_edge):
        for x in range(lattice_edge):
            ID = sim_array[step_index][x][y]
            if ID==0:
                frame[x,y] = [0.9, 0.9, 0.9]
            elif ID % 2 == 1:
                frame[x,y] = [0.1, 0.1, 0.5 + ID/(4*num_cells)]
            else:
                frame[x,y] = [0.5 + ID/(4*num_cells), 0.1, 0.1]
    return frame

# After data is/are read in from individual lattice files into the sim_array that stores all of this information,
# timepoints are processed sequentially into RGB arrays and then animated using FuncAnimation.
#
# After that, depending on settings, the plot can be saved and/or displayed.
def lattices_to_FuncAnimation():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_ylim(0, lattice_edge)

    frames = [lattice_to_frame(i, ax) for i in range(start_step, end_step, skip)]

    # Function used by the FuncAnimation class to generate the first frame of the video
    def init():
        return plt.imshow(frames[0], animated=True)

    # Function used by the FuncAnimation class to generate each subsequent frame of the video
    def animate(nframe):
        plt.cla()
        if (verbose):
            print("Generating artist frame number ", nframe, "of ", len(frames))
        im = plt.imshow(frames[nframe], animated=True)
        plt.title('timestep = {}'.format((nframe + start_step) * step_interval * skip))
        plt.gca().invert_yaxis()

    if (verbose):
        print("Generating animation.")
    ani = animation.FuncAnimation(fig, animate, frames=(int((end_step-start_step)/skip)), interval=40)  # , repeat_delay=1000

    if (save_animation):
        if (verbose):
            print("Saving animation to .mp4")
        # Writer options: imagemagick, ffmpeg
        # Imagemagick mp4 can only be opened in certain programs (not ppt), use ffmpeg as writer for more reliable saving
        ani.save('{0}_simulation_{1}_steps.mp4'.format(time.strftime("%Y%m%d-%H%M%S"), (end_step - start_step)),writer="ffmpeg",fps=10)

    if (show_animation):
        if (verbose):
            print("Displaying animation.")
        plt.show()

# Run the main method.
main()
