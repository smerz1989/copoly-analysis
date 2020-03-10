from dump_generator import read_dump_w_header
import sys
import numpy as np
from io import StringIO

filename = sys.argv[1]
isbonddump = sys.argv[2]
with open('{}.filtered'.format(filename),'w') as trjfile:
	trjfile.write('')

d = read_dump_w_header(filename)
for snapshot in d:
    if int(snapshot[0])%20000==0:
        print("On timestep {}".format(snapshot[0]))
        with open('{}.filtered'.format(filename),'a') as filtered_file, StringIO() as s:
            filtered_file.write(snapshot[1])
            if isbonddump=="yes":
                np.savetxt(s,snapshot[2],fmt='%d %d %d')
            else:
                np.savetxt(s,snapshot[2],fmt='%d %d %f %f %f')
            filtered_file.write(s.getvalue())
