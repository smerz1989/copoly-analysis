import simulation_class as smc
import argparse
import os

parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-N','--total_monomers',dest='total_monomers',type=int,default=3000, help="Total number of monomers in the system")
parser.add_argument('-fA','--fractionA',dest='fA',type=float,default=0.5, help="Fraction of monomers of species A")
parser.add_argument('--epsAA',dest='eAA',type=float,default=1.0,help="Attraction between monomers AA in kbT")
parser.add_argument('--epsBB',dest='eBB',type=float,default=1.0,help="Attraction between monomers BB in kbT")
parser.add_argument('--epsAB',dest='eAB',type=float,default=1.0,help="Attraction between monomers AB in kbT")
parser.add_argument('-p',dest='p',type=float,default=0.9,help="Monomer conversion where simulation will end")
parser.add_argument('-f','--folder',dest='dest_folder')
parser.add_argument('--servername',dest='servername',default=None)
parser.add_argument('--send_to_cluster',dest='send_to_cluster',type=bool,default=False)

args = parser.parse_args()

epsilons = (args.eAA,args.eBB,args.eAB)

project_path = os.path.abspath(args.dest_folder)

sim = smc.Simulation(total_monomers=args.total_monomers,
                     monomer_attractions=epsilons,
                     p = args.p,monomer_A_fraction=args.fA,
                     send_to_cluster=args.send_to_cluster,
                     servername=args.servername)
sim.compile_simulation()
sim.move_simulation_files(project_path)
sim.start_simulation()
