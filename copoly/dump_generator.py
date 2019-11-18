import numpy as np
import subprocess as sb
from itertools import islice
from math import ceil
import re

def get_number_of_timesteps(filename,header_size=9):
    print("Getting number of timesteps")
    with open('tmp.out','w') as tmp:
        first_header = sb.check_output(['head','-n',str(header_size),filename]).decode('utf-8')
        num_entries = int(re.findall(r'[0-9]+',' '.join(first_header.split('\n')))[1])
    num_lines = int(sb.check_output(['wc','-l'],stdin=open(filename,'r')).decode('utf-8'))
    return(int(num_lines/(num_entries+header_size)))

def read_dump(filename,header_size=9,scale=False):
    cur_line_number=0
    with open(filename,'r') as datafile:
        while True:
            try:
                header = ''.join(datafile.readline() for i in range(header_size))
                if header=='':
                    break
                header_data = re.findall(r'[0-9]+',header)
                timestep=header_data[0]
                num_entries = int(header_data[1])
                box_bounds = np.array([list(map(float,line.split(' '))) for line in header.split('\n')[-5:-2]])
                lines = [datafile.readline() for i in range(num_entries)]
                data1D = np.fromstring(' '.join(lines).replace('\n',' '),sep = ' ')
                try:
                    data = data1D.reshape((num_entries,ceil(len(data1D)/num_entries)))
                except ZeroDivisionError:
                    print("No entries for timestep returning empty data frame")
                    data=[]
                if scale:
                    data[:,-3:] = data[:,-3:]*(box_bounds[0,1]-box_bounds[0,0])+box_bounds[0,0]
                yield(timestep,data)
            except EOFError as e:
                break
