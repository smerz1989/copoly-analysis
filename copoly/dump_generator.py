import numpy as np
import subprocess as sb
from itertools import islice
from math import ceil
import re

#def read_file_bash(filename,start,finish):
#    sb.check_output(['head'])


def read_dump(filename,header_size=9):
    cur_line_number=0
    with open(filename,'r') as datafile:
        while True:
            #header = ''.join(islice(datafile,cur_line_number,cur_line_number+9))
            header = ''.join(datafile.readline() for i in range(header_size))
            #print(header)
            header_data = re.findall(r'[0-9]+',header)
            num_entries = int(header_data[1])
            #print(num_entries)
            lines = [datafile.readline() for i in range(num_entries)]
            #lines = islice(datafile,cur_line_number+num_entries)
            #cur_line_number +=num_entries+9
            data1D = np.fromstring(' '.join(lines).replace('\n',' '),sep = ' ')
            data = data1D.reshape((num_entries,ceil(len(data1D)/num_entries)))
            yield(data)
