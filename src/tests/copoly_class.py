import numpy as np
import subprocess as sb
import os, sys

class Copolymer(object):    
    def __init__(self,sequence,sep=1.1,lt_dir="./lts/",name="Copoly",outfile='poly.lt'):
        self.sequence = sequence
        self.sep = sep   
        self.lt_dir = os.path.abspath(lt_dir)
        self.name = name
        self.outfile = outfile

    def write_coords(self):
        coords = np.array([[0,0+self.sep*i,0] for i in range(len(self.sequence))])
        np.savetxt('coords.raw',coords,'%0.4f %0.4f %0.4f')
 
    def write_sequence(self):
        seq_array = np.array([seq for seq in self.sequence])
        np.savetxt('sequence.txt',seq_array,fmt='%s')

    def create_ltfile(self):
        self.write_coords()
        self.write_sequence()
        sb.call(["mv","sequence.txt","coords.raw",self.lt_dir])
        os.chdir(self.lt_dir)
        sb.call(["./genpoly_lt.py","-bond","SASA","SA2","SA1",
                                   "-header","import \"a_bead.lt\"\nimport \"b_bead.lt\"\nimport \"graft_bead.lt\"",
                                   "-polymer-name",self.name,
                                   "-inherits","COPOLYFF",
                                   "-sequence","sequence.txt"],stdin=open('coords.raw','r'),stdout=open(self.outfile,'w'))

