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
            try:
                header = ''.join(datafile.readline() for i in range(header_size))
                if header=='':
                    break
                header_data = re.findall(r'[0-9]+',header)
                num_entries = int(header_data[1])
                lines = [datafile.readline() for i in range(num_entries)]
                data1D = np.fromstring(' '.join(lines).replace('\n',' '),sep = ' ')
                data = data1D.reshape((num_entries,ceil(len(data1D)/num_entries)))
                yield(data)
            except EOFError as e:
                break