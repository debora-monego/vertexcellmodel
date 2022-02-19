#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

token = 'timestep'
chunks = []
current_chunk = []

for line in open('outvertex.txt'):
   if line.startswith(token) and current_chunk: 
      # if line starts with token and the current chunk is not empty
      chunks.append(current_chunk[:]) #  add not empty chunk to chunks
      current_chunk = [] #  make current chunk blank
   # just append a line to the current chunk on each iteration
   current_chunk.append(line)

chunks.append(current_chunk) 
print chunks