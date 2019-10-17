"""This module contains the molecule class and all classes and functions
related to molecules.  This includes Bond, Angle, and Dihedral classes as well as helper functions that take in atom and bond lists and generate molecules.

"""
from itertools import groupby, permutations
import itertools
import networkx as ntwkx
import numpy as np
from math import *
from copoly.dump import dump
from copoly import dump_generator
from subprocess import call, check_output
import re
import igraph
import subprocess as sb
import pandas as pd

class Molecule(object):
    """This class is used to represent a molecule in a simulation and holds all
        the objects related to a molecule including the related Atom, Bond, Angle, and Dihedral Objects.

        Parameters
        ----------
        molID : int
            The unique integer identifier of the molecule as set in the LAMMPS input file
        atoms : list of type Atom
            The Atom objects associated with this molecule
        bonds : list of type Bond
            The Bond objects associated with this molecule
        angles : list of type Angle
            The Angle objects associated with this molecule
        dihedrals : list of type Dihedral
            The Dihedral objects associated with this molecule

    """
    def __init__(self,molID,atoms,bonds,angles,dihedrals):
        self.molID = molID
        self.atoms = atoms
        self.bonds = bonds
        self.graph = self.atomsAsGraph() if self.bonds is not None else None

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return all([self.atoms==other.atoms,self.bonds==other.bonds])
        else:
            return False
    def __neq__(self,other):
        return not self.__eq__(other)

    def get_com(self):
        """Calculates the center of mass of the molecule.

        Returns
        -------
        com : array of floats
            The center of mass of the molecule.
        """
        positions = np.array([atom.position for atom in self.atoms])
        com = np.mean(positions,axis=0)
        return(com)

    def setAnchorAtom(self,atomID):
        """Sets the anchor atom of the molecule this is the atom of the molecule that is anchored to the nanoparticle.  Setting this important when using the CBMC regrowth move.
        
        Parameters
        ----------
        atomID : int
            The unique atom ID identifier of the Atom object that is anchored to the nanoparticle.
        """
        atom = self.getAtomByID(atomID)
        if(atom!=None):
            self.anchorAtom = atom
            return True
        else:
            return False
    
    def atomsAsGraph(self):
        """Returns the graph data structure of the atoms in the molecule defined by the connectivity defined by the Bond objects of the molecule.  
        The nodes of the graph are the atom ID's of the atoms.

        Returns
        -------
        Networkx Graph Object
            A Networkx Graph object with the nodes being the atom ID's of the atoms in the molecule and the connections defined by the Molecule's Bonds
        """
        return molecule2graph(self.atoms,self.bonds)
    def getAtomByID(self,atomID):
        """Returns the Atom object associated with the given atom ID as long as the atom is associated with the molecule.

        Parameters
        ----------
        atomID : int
            The unique atom ID identifier of the Atom object which you hope to retrieve

        Returns
        -------
            Atom
                The Atom object associated with the atom ID or None if the Atom is not associated with this molecule.
        """
        for atom in self.atoms:
            if(atom.atomID==atomID):
                return atom
        return None

    def getAtomByMolIndex(self,index):
        """Returns the Atom by it's index in the molecule where the index is defined as the number of bonds away from the anchor atom 
        (i.e. the atom at index 1 is the atom directly connected to the anchor atom)

        Parameters
        ----------
        index : int
            The molecular index of the Atom object one wishes to retrieve where the index is defined as the number of bonds away from the anchor atom.

        Returns
        -------
        Atom Object
            The Atom Object located at the specified index if the index is greater than the number of atoms in the molecule minus one then 
            the function returns None as this is out of the range of the atom list.
        """
        if(index>(len(self.atoms)-1)):
            return None
        successor_dict= ntwkx.dfs_successors(self.graph,source=self.anchorAtom.atomID)
        currentID = self.anchorAtom.atomID
        for i in range(index):
            currentID = successor_dict[currentID][0]
        return self.getAtomByID(currentID)

    def align_to_vector(self,vector):
        molecule_vector = self.get_com()-self.anchorAtom.position
        if np.linalg.norm(vector)==0. or np.linalg.norm(molecule_vector)==0.:
            raise ValueError("Alignment vector passed in must have a non-zero magnitude.")
        anchor_position = self.anchorAtom.position
        axis_rotation = np.cross(molecule_vector,vector)
        angle = acos(np.dot(molecule_vector/np.linalg.norm(molecule_vector),vector/np.linalg.norm(vector)))
        rotate_atoms = [atom for atom in self.atoms if not (atom.atomID==self.anchorAtom.atomID)]
        for atom in rotate_atoms:
            atom.position = rot_quat((atom.position-anchor_position),angle,axis_rotation)+anchor_position
        

    def move_atoms(self,move):
        for atom in self.atoms:
            atom.position+=move

    def move_atoms_by_index(self,move,index):
        for i in range(index,len(self.atoms)):
            self.getAtomByMolIndex(i).position+=move


class Bond(object):
    """The Bond object represents a bond between two atoms. The format is similar to a LAMMPS bond, therefore a 
    Bond object consists of a Bond ID which uniquely defines the bond, a bond type, and the atom ID's of the two atoms involved in the bond.
    
    Parameters
    ----------
    bondID : int
        The unique integer identifying the bonds same as the one defined in the LAMMPS input file
    bondType : int
        An integer which represents the type of bond this is, the same as the bond type number defined in LAMMPS input file.
    atom1 : int
        The atom ID associated with the first atom in the bond
    atom2 : int
        The atom ID associated with the second atom in the bond
    """
    def __init__(self,bondID,bondType,atom1,atom2):
        self.bondID = int(bondID)
        self.bondType = int(bondType)
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return self.bondID == other.bondID
        else:
            return False
    def __neq__(self,other):
        return not self.__eq__(other)

class Angle(object):
    """The Angle object represents the angle between three connected atoms.  The format is the same as in the Angles section of a LAMMPS input file.
    
    Parameters
    ----------
    angleID : int
        The unique integer identifier of the angle, the same as the one defined in the LAMMPS inpit file
    angleType : int
        The integer identifier of the angle type which is the same as the angle type number defined in the LAMMPS input file
    atom1 : int
        The atom ID of the first atom associated with the angle, same as the one defined in LAMMPS input file.
    atom2 : int
        The atom ID of the second atom associated with the angle, same as the one defined in LAMMPS input file.
    atom3 : int
        The atom ID of the third atom associated with the angle, same as the one defined in LAMMPS input file.
    """
    def __init__(self, angleID,angleType,atom1,atom2,atom3):
        self.angleID = int(angleID)
        self.angleType =int(angleType)
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return other.angleID == self.angleID
        else:
            return False
    def __neq__(self,other):
        return not self.__eq__(other)


def load_bond_trajectory_old(filename):
    d = dump(filename)
    times = d.time()
    bond_snapshots = {}
    for timestep in times:
        b_type,batom1,batom2 = d.vecs(timestep,"c_bcomp[1]","c_bcomp[2]","c_bcomp[3]")
        bonds = {}
        for bid,(btype,batom1,batom2) in enumerate(zip(b_type,batom1,batom2)):
            bonds[bid] = Bond(bid,btype,batom1,batom2)
        bond_snapshots[timestep]=bonds
    return(bond_snapshots)

def load_bond_trajectory(bdump):
    times = bdump.time()
    for timestep in times:
        b_type,batom1,batom2 = bdump.vecs(timestep,"c_bcomp[1]","c_bcomp[2]","c_bcomp[3]")
        bonds = np.concatenate((np.array(b_type)[:,np.newaxis],np.array(batom1)[:,np.newaxis],np.array(batom2)[:,np.newaxis]),axis=1).astype(int)
        yield((timestep,bonds))

def load_atom_trajectory(adump):
    times = adump.time()
    for timestep in times:
        #adump.unscale(timestep)
        adump.sort(timestep)
        a_id,a_type,xs,ys,zs = adump.vecs(timestep,"id","type","x","y","z")
        atom_coords = np.concatenate((np.array(a_id)[:,np.newaxis],
                               np.array(a_type)[:,np.newaxis],
                               np.array(xs)[:,np.newaxis],
                               np.array(ys)[:,np.newaxis],
                               np.array(zs)[:,np.newaxis]),axis=1).astype(float)
        yield((timestep,atom_coords))

def set_anchor_atoms(molecules,anchortype):
    """For every molecule in molecules set the anchot atom to the atom with atom type anchortype.

    Parameters
    ----------
    molecules : list of Molecule
    A list of Molecule objects that will have their anchor atoms set to the atom with type anchortype.
    anchortype : int
    The atom type of the atom to set as the anchor for each Molecule in molecules.
    """
    for key, molecule in molecules.items():
        anchorIDs = [atom.atomID for atom in molecule.atoms if atom.atomType==anchortype]
        if len(anchorIDs)>0:
            molecule.setAnchorAtom(anchorIDs[0])

def construct_molecule_trajectory(datafile,bondtrajectory,atomtrajectory=None):
    """From a LAMMPS input file construct a list of Molecule objects based on the molecules in the LAMMPS input file.

    Parameters
    ----------
    filename : str
        The name of the LAMMPS input file that contains the molecules

    Returns
    -------
    molecules : Molecule List
        A list of Molecule objects with the data specified by the LAMMPS input file passed in.
    """
    atoms,atoms_array = loadAtoms(datafile)
    bond_snapshots = load_bond_trajectory(bondtrajectory)
    atom_snapshots = None if atomtrajectory is None else load_atom_trajectory(atomtrajectory)
    for timestep, bonds in bond_snapshots:
        if not atom_snapshots is None:
            timestep,atom_coords = next(atom_snapshots)
        else:
            atom_coords=None
        yield((timestep,SimulationSnapshot(atoms,bonds,atoms_array,atom_coords=atom_coords)))

def construct_molecule_trajectory_from_generators(datafile,bondfile,atomfile):
    """From a LAMMPS input file construct a list of Molecule objects based on the molecules in the LAMMPS input file.

    Parameters
    ----------
    filename : str
        The name of the LAMMPS input file that contains the molecules

    Returns
    -------
    molecules : Molecule List
        A list of Molecule objects with the data specified by the LAMMPS input file passed in.
    """
    atoms,atoms_array = loadAtoms(datafile)
    bond_snapshots = dump_generator.read_dump(bondfile)
    atom_snapshots = dump_generator.read_dump(atomfile,scale=False)
    for (timestep,traj_array), (timestep,bonds) in zip(atom_snapshots,bond_snapshots):
        yield(SimulationSnapshot(atoms,bonds,atoms_array,atom_coords=traj_array))



def rot_quat(vector,theta,rot_axis):
    """Rotates a vector about a specified axis a specified angle theta using the quaternion method

    Parameters
    ----------
    vector : float vector
        A vector of three elements that represents the X, Y, Z coordinates of the vector that one wishes to rotate
    theta : float
        The angle in radians by which the vector rotates.
    rot_axis : float vector
        A vector of three elements which represents the X,Y,Z elements of the rotation axis.

    Returns
    -------
    new_vector : float vector
        A vector of three elements representing the X,Y,Z coordinates of the old vector after rotation
    """
    if np.linalg.norm(rot_axis)==0.:
        raise ValueError("The rotation axis must have a non-zero magnitude in order for rotation about the axis to make sense.")
    rot_axis = rot_axis/np.linalg.norm(rot_axis)
    vector_mag = np.linalg.norm(vector)
    quat = np.array([cos(theta/2),sin(theta/2)*rot_axis[0],sin(theta/2)*rot_axis[1],sin(theta/2)*rot_axis[2]])
    quat_inverse = np.array([cos(theta/2),-sin(theta/2)*rot_axis[0],-sin(theta/2)*rot_axis[1],-sin(theta/2)*rot_axis[2]])

    vect_quat = np.array([0,vector[0],vector[1],vector[2]])/vector_mag
    new_vector = quat_mult(quat_mult(quat,vect_quat),quat_inverse)
    return new_vector[1:]*vector_mag

def quat_mult(q1,q2):
    w = w1*w2-x1*x2-y1*y2-z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 + y1*w2 + z1*x2 - x1*z2
    z = w1*z2 + z1*w2 + x1*y2 - y1*x2
    return np.array([w,x,y,z])


class SimulationSnapshot(object):
    """A class encapsulating the atom position and bond topology of a snapshot of atoms in a NP copoly simulation.

    Parameters
    ----------
    atoms : list of type Atom
        A list of all atoms within the simulation.
    bonds : list of type Bond
        A list of all bonds within the simulation snapshot.
    """ 
    def __init__(self,atoms,bonds,atoms_array=None,atom_coords=None,anchor_atom_type=8):
        self.atoms = {atom.atomID: atom for atom in atoms if not atom.atomType==1}
        self.bonds = bonds
        self.atoms_array = atoms_array
        if not atom_coords is None:
            self.atoms_array[:,3:6]=atom_coords[:,2:5]
        self.anchor_atoms = atoms_array[atoms_array[:,2].astype(int)==anchor_atom_type]
        self.topology = self.create_topology_network(atoms_array,self.bonds)
        self.molecules = list(self.topology.components())
        self.relabel_molecules()
        self.monomers = [molecule for molecule in self.molecules if len(molecule)==3]
        self.chains = [molecule for molecule in self.molecules if len(molecule)>3]
        self.chain_lengths = [int(len(molecule)/3) for molecule in self.molecules]
    
    def guess_simulation_bounds(self):
        xmin,xmax = (np.amin(self.atoms_array[:,3]),np.amax(self.atoms_array[:,3]))
        ymin,ymax = (np.amin(self.atoms_array[:,4]),np.amax(self.atoms_array[:,4]))
        zmin,zmax = (np.amin(self.atoms_array[:,5]),np.amax(self.atoms_array[:,5]))
        return(np.array([[xmin,xmax],[ymin,ymax],[zmin,zmax]]))

    def relabel_molecules(self):
        for i,molecule in enumerate(self.molecules):
            self.atoms_array[molecule+self.min_node-1,1]=i+self.min_node

    def unwrap_molecules(self,bounds):
        for molecule in self.molecules:
            coords = self.atoms_array[molecule+self.min_node-1,3:]
            for i,coord in enumerate(coords):
                if not i==0:
                    disps = coords[i]-coords[i-1]
                    correction = [0 if disp<bound else -disp for disp,bound in zip(disps,bounds)]
                    self.atoms_array[molecule[i],3:]=self.atoms_array[molecule[i],3:]+correction

    def to_LAMMPS_datafile(self):
        sim_bounds = self.guess_simulation_bounds()
        #self.unwrap_molecules([sim_bounds[0,1]-sim_bounds[0,0],
        #                       sim_bounds[1,1]-sim_bounds[1,0],
        #                       sim_bounds[2,1]-sim_bounds[2,0]])
        print("Guessed sim_bounds are x: {}-{}, y: {}-{}, z: {}-{}".format(sim_bounds[0,0],sim_bounds[0,1],
                                                                           sim_bounds[1,0],sim_bounds[1,1],
                                                                           sim_bounds[2,0],sim_bounds[2,1]))
        np.savetxt('atom_tmp.txt',self.atoms_array,fmt='%d\t%d\t%d\t%0.4f\t%0.4f\t%0.4f')
        index_bonds = np.concatenate(((np.arange(len(self.bonds))+1)[:,np.newaxis],self.bonds),axis=1)
        np.savetxt('bonds_tmp.txt',index_bonds,fmt='%d\t%d\t%d\t%d')
        unique_atoms = int(np.amax(np.unique(self.atoms_array[:,2])))
        unique_bonds = int(np.amax(np.unique(self.bonds[:,0])))
        header = ("LAMMPS Description\n\n"
                  "{} atoms\n"
                  "{} bonds\n\n"
                  "{} atom types\n"
                  "{} bond types\n\n"
                  "{} {} xlo xhi\n"
                  "{} {} ylo yhi\n"
                  "{} {} zlo zhi\n\n"
                  "Masses\n\n").format(len(self.atoms_array),len(self.bonds),
                                       unique_atoms,unique_bonds,
                                       sim_bounds[0,0],sim_bounds[0,1],
                                        sim_bounds[1,0],sim_bounds[1,1],
                                        sim_bounds[2,0],sim_bounds[2,1])
        masses = "\n".join(["{} {}".format(i+1,mass) for i,mass in enumerate(range(unique_atoms))])
        with open('tmpdata.data','w') as datafile:
            datafile.write(header)
            datafile.write(masses)
            datafile.write("\n\nAtoms\n\n")
            datafile.write(open('atom_tmp.txt','r').read())
            datafile.write("\nBonds\n\n")
            datafile.write(open('bonds_tmp.txt','r').read())

    def visualize_snapshot(self,imagefilename,destfolder=''):
        self.to_LAMMPS_datafile()
        monomer_molnumbers = [str(molnumber) for molnumber in self.atoms_array[np.array(self.monomers)[:,0]+self.min_node-1,1].astype(int)]
        sb.call(["ovitos","ovito_template.py","-f","tmpdata.data",'-d',destfolder,"-m",",".join(monomer_molnumbers),'-i',imagefilename]) 
        #sb.call(["ovitos","ovito_template.py","-f","tmpdata.data","-m","0"])  
 
    def create_topology_network(self,atoms, bonds, anchor_atom_type=8):
        G = igraph.Graph()
        nodes = self.atoms_array[:,0].astype(int)
        #nodes = list(self.atoms.keys())
        #atom_attributes = list(self.atoms.values())
        atom_attributes = self.atoms_array
        G.add_vertices(nodes)
        #G.vs["Atom"] = atoms
        G.vs["Atom"] = self.atoms_array
        self.min_node = np.min(G.vs['name'])
        #self.min_node = int(np.amin(self.atoms_array[:,0]))
        G.add_edges(zip(bonds[:,1].astype(int)-self.min_node,self.bonds[:,2].astype(int)-self.min_node))
        return(G)
    
    def get_number_chains(self):
        return(len([length for length in self.chain_lengths if length>1]))

    def get_number_monomers(self):
        return(len([length for length in self.chain_lengths if length==1]))

    def get_monomer_type(self,monomer):
        node_types = [node['Atom'][2] for node in self.topology.vs]
        unique_types, type_counts = np.unique(np.array(node_types),return_counts=True)
        return(unique_types[type_counts==1])

    def unravel_chain_generator(self,chain_generator):
        [chainnode for node in chain_generator]

    def get_monomer_type_fraction(self,atom_type=3):
        atoms = np.array(self.monomers).flatten()
        atom_types = [self.topology.vs[atom]['Atom'][2] for atom in atoms]
        types, counts = np.unique(atom_types,return_counts=True)
        return(counts[types==atom_type][0]/len(self.monomers))

    def get_chain_type_fraction(self,atom_type=3):
        atoms = itertools.chain(self.chains)
        atom_types = [self.topology.vs[atom]['Atom'][2] for atom in atoms]
        types, counts = np.unique(atom_types,return_counts=True)
        return(counts[types==atom_type])

    def get_dop(self):
        return(np.mean([length for length in self.chain_lengths if length>0.5]))

    def get_pdi(self):
        length_length2 = np.array([[length,length**2] for length in self.chain_lengths if length>1])
        if len(length_length2):
            n_ave, m_ave = (np.mean(length_length2,axis=0)[0], np.mean(length_length2,axis=0)[1])
        else:
            n_ave, m_ave = (1,1)
        return(m_ave/n_ave)
    
    def get_chain_sequence(self,chain):
        sequence = ' '.join([str(int(atom['Atom'][2])) for atom in chain])
        return(sequence)

    def get_chain_sequence_filtered(self,chain,filter_atom_types=(2,5,6,7,8),a_type_id=3,b_type_id=4):
        sequence = ''.join([str(int(atom['Atom'][2])) for atom in chain if int(atom['Atom'][2]) not in filter_atom_types])
        return(sequence)

    def get_sequence_probs(self,sequence,filter_atom_type=(5,6,7),a_type_id=3,b_type_id=4):
        filter_str = r'['+','.join([str(atom_type) for atom_type in filter_atom_type])+']\ ?'
        str_sequence = self.get_chain_sequence(sequence)
        filtered_seq = re.sub(filter_str,"",str_sequence).strip()
        filtered_seq = re.sub(r'\ ',"",filtered_seq)
        filtered_seq_one_shift = filtered_seq[-1]+filtered_seq[:-1]
        numAA = len(re.findall(r'(?='+str(a_type_id)*2+')',filtered_seq))
        numBB = len(re.findall(r'(?='+str(b_type_id)*2+')',filtered_seq)) 
        numAB = len(re.findall(r'(?='+str(a_type_id)+str(b_type_id)+')',filtered_seq))
        numBA = len(re.findall(r'(?='+str(b_type_id)+str(a_type_id)+')',filtered_seq))
        #numAA_oneshift = len(re.findall(r'(?='+str(a_type_id)*2+')',filtered_seq_one_shift))
        #numBB_oneshift = len(re.findall(r'(?='+str(b_type_id)*2+')',filtered_seq_one_shift)) 
        #numAB_oneshift = len(re.findall(r'(?='+str(a_type_id)+str(b_type_id)+')',filtered_seq_one_shift))
        #numBA_oneshift = len(re.findall(r'(?='+str(b_type_id)+str(a_type_id)+')',filtered_seq_one_shift))
        return((numAA,numBB,numAB,numBA))        

    def get_all_probs(self):
        sequences = self.get_sequences()
        chain_counts = zip(*[self.get_sequence_probs(seq) for seq in sequences])
        chain_sums = [sum(count) for count in chain_counts]
        chain_probs = np.array(chain_sums)/sum(chain_sums) if sum(chain_sums)>0 else (0,0,0,0) 
        return(chain_probs) 

    def get_conditional_probs(self,a_type_id=3,b_type_id=4):
        sequences = self.get_sequences()
        str_sequences = [self.get_chain_sequence_filtered(seq) for seq in sequences]
        numAA,numBB,numAB,numBA,numA,numB = (0,0,0,0,0,0)
        for seq in str_sequences:
            numAA += len(re.findall(r'(?='+str(a_type_id)*2+')',seq))
            numBB += len(re.findall(r'(?='+str(b_type_id)*2+')',seq)) 
            numAB += len(re.findall(r'(?='+str(a_type_id)+str(b_type_id)+')',seq))
            numBA += len(re.findall(r'(?='+str(b_type_id)+str(a_type_id)+')',seq))
            numA += len(re.findall(r'(?='+str(a_type_id)+')',seq[:-1]))
            numB += len(re.findall(r'(?='+str(b_type_id)+')',seq[:-1]))
        probs = ((numAA/numA),(numBB/numB),(numAB/numA),(numBA/numB)) if numA!=0 and numB!=0 else (0,0,0,0)
        return(probs)


    def get_sequences(self):
        sequences = [self.topology.bfsiter(self.topology.vs[int(anchor[0])-self.min_node]) for anchor in self.anchor_atoms]
        return(sequences) 


class Atom(object):
    """The Atom class represents an atom described by the LAMMPS full atom style.

    Parameters
    ----------
    atomID : int
        The unique integer identifier of the atom as specified in the LAMMPS input file.
    molID : int
        The unique integer identifier of the molecule associated with this atom as specified in the LAMMPS input file.
    atomType : int
        The integer that identifies the atom type as specified in the LAMMPS input file.
    charge : float, optional
        The charge on the atom defaults to 0.
    position : float vector, optional
        A three element vector which represent the X,Y,Z coordinates of the atom.  It defaults to a vector of [X=0,Y=0,Z=0].
    """
    def __init__(self,atomID,molID,atomType,charge=0.,position=[0.,0.,0.]):
        self.atomID = int(atomID)
        self.molID = int(molID)
        self.atomType = int(atomType)
        self.charge = float(charge)
        self.position = position
    
    def __eq__(self,other):
        if isinstance(other, self.__class__):
            return self.atomID == other.atomID
        else:
            return False
    def __neq__(self,other):
        return not self.__eq__(other)
    def get_pos(self):
        return self.position
    def get_charge(self):
        return self.charge
    def get_type(self):
        return self.atomType
    def get_mol_ID(self):
        return self.molID
    def get_atom_ID(self):
        return self.atomID

def loadAtoms(filename,style="angle"):
    """Loads the atoms from a LAMMPS  input file and returns a list of Atom object which represent those atoms.
    
    Parameters
    ----------
    filename : str
        The name of the LAMMPS input file which contain the atoms

    Returns
    -------
    Atom List
        A list of Atom objects which contains the atom info in the given LAMMPS input file.
    """
    with open('tmp.out','w') as temp_file:
        call(["awk",'/Atoms/{flag=1;next}/Bonds/{flag=0}flag',filename],stdout=temp_file)
    atoms_array = np.loadtxt("tmp.out",skiprows=1)
    call(["rm","tmp.out"])
    atom_list = [Atom(atom[0],atom[1],atom[2],0.,atom[3:6]) for atom in atoms_array]
    return (atom_list,atoms_array)


def loadBonds(filename,style="angle"):
    """Loads the atoms from a LAMMPS  input file and returns a list of Atom object which represent those atoms.
    
    Parameters
    ----------
    filename : str
        The name of the LAMMPS input file which contain the atoms

    Returns
    -------
    Atom List
        A list of Atom objects which contains the atom info in the given LAMMPS input file.
    """
    with open('tmp.out','w') as temp_file:
        call(["awk",'/Bonds/{flag=1;next}/Angles/{flag=0}flag',filename],stdout=temp_file)
    bonds_array = np.loadtxt("tmp.out",skiprows=1,dtype=int)
    call(["rm","tmp.out"])
    #atom_list = [Atom(atom[0],atom[1],atom[2],0.,atom[3:6]) for atom in atoms_array]
    return(bonds_array)















