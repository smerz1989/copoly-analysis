import numpy as np
import subprocess as sb
from itertools import islice
from math import ceil
import re
import bz2
import os

def get_number_of_timesteps(filename,header_size=9):
    print("Getting number of timesteps")
    _ ,file_extension = os.path.splitext(filename)
    compress = file_extension=='.bz2'
    with open('tmp.out','w') as tmp:
        if compress:
            unzip_process = sb.Popen(['bzcat',filename],stdout=sb.PIPE)
            head_process = sb.Popen(['head','-n',str(header_size)],stdin=unzip_process.stdout,stdout=sb.PIPE)
            unzip_process.stdout.close()
            first_header_bytes, err =head_process.communicate()
            first_header = first_header_bytes.decode('utf-8')
        else:
            first_header = sb.check_output(['head','-n',str(header_size),filename]).decode('utf-8')
        num_entries = int(re.findall(r'[0-9]+',' '.join(first_header.split('\n')))[1])
    if compress:
        unzip_process = sb.Popen(['bzcat',filename],stdout=sb.PIPE)
        wc_process = sb.Popen(['wc','-l'],stdin=unzip_process.stdout,stdout=sb.PIPE)
        unzip_process.stdout.close()
        num_lines, err = wc_process.communicate()
        num_lines = int(num_lines)
    else:
        num_lines = int(sb.check_output(['wc','-l'],stdin=open(filename,'r')).decode('utf-8'))
    return(int(num_lines/(num_entries+header_size)))


def read_dump_compressed(filename,header_size=9,scale=False):
    cur_line_number=0
    with bz2.open(filename,'r') as datafile:
        while True:
            try:
                header = ''.join(datafile.readline().decode('utf-8') for i in range(header_size))
                if header=='':
                    break
                header_data = re.findall(r'[0-9]+',header)
                timestep=header_data[0]
                num_entries = int(header_data[1])
                box_bounds = np.array([list(map(float,line.split(' '))) for line in header.split('\n')[-5:-2]])
                lines = [datafile.readline().decode('utf-8') for i in range(num_entries)]
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

def read_dump(filename,header_size=9,scale=False):
    cur_line_number=0
    _ ,file_extension = os.path.splitext(filename)
    if file_extension=='.lammpstrj':
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
    else:
       with bz2.open(filename,'r') as datafile:
            while True:
                try:
                    header = ''.join(datafile.readline().decode('utf-8') for i in range(header_size))
                    if header=='':
                        break
                    header_data = re.findall(r'[0-9]+',header)
                    timestep=header_data[0]
                    num_entries = int(header_data[1])
                    box_bounds = np.array([list(map(float,line.split(' '))) for line in header.split('\n')[-5:-2]])
                    lines = [datafile.readline().decode('utf-8') for i in range(num_entries)]
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

