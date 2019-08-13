
class Simulation(object):
    def __init__(self,total_monomers=3000,monomer_A_fraction,p=0.9,dump_frequency):
        self.total_monomers = total_monomers
		self.monomer_A_fraction = monomer_A_fraction
		self.p = p
		self.dump_frequency = dump_frequency
