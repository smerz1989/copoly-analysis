import ovito.io as io
import ovito.modifiers as mods
import ovito.vis as vis
import math
import sys
import argparse

parser = argparse.ArgumentParser(description='Flexible Python script template for use with Ovito visualization tool.')
parser.add_argument('-f','--filename',dest='filename', help="Parent directory where simulations are kept")
parser.add_argument('-d','--dest_folder',dest='dest_folder',help="Location where results are placed.")
parser.add_argument('-m','--monlist', dest='mon_list',default='0',help="Comma separated list of monomers to be set to transparent")
parser.add_argument('-i','--imagename',dest='imagename',default='testfile.png',help='Name to give rendered image.')
args = parser.parse_args()

filename = args.filename
#filename = sys.argv[1]
#mol1 = int(sys.argv[2])

def render_overlay(painter, **args):
    painter.drawText(10, 10, "Hello world")

monomer_bead_radius = "0.3"
peripheral_radius = "0.15"
anchor_bead_radius = "0.35"
np_bead_radius = "0.5"
transparency="0.6"

print(args.mon_list)
#import pdb;pdb.set_trace()
node = io.import_file(filename,atom_style="angle")

node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==1"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[np_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(160/255,160/255,160/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==2"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(204/255,153/255,255/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==3"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(102/255,102/255,255/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==4"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(255/255,0/255,127/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==5"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[peripheral_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(255/255,252/255,242/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==6"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(255/255,0/255,255/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==7"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(127/255,0/255,255/255)))
node.modifiers.append(mods.SelectExpressionModifier(expression="ParticleType==8"))
node.modifiers.append(mods.ComputePropertyModifier(expressions=[anchor_bead_radius],output_property="Radius",only_selected=True))
node.modifiers.append(mods.AssignColorModifier(color=(204/255,255/255,229/255)))
#node.modifiers.append(mods.ComputePropertyModifier(expressions=[monomer_bead_radius],output_property="Transparency",only_selected=False))
monomer_expression = " || ".join(["MoleculeIdentifier == {}".format(mol) for mol in args.mon_list.split(",")]) 
print(monomer_expression)
node.modifiers.append(mods.SelectExpressionModifier(expression=monomer_expression))
#node.modifiers.append(mods.ComputePropertyModifier(expressions=[transparency],output_property="Transparency",only_selected=True))
#node.modifiers.append(mods.InvertSelectionModifier())
node.modifiers.append(mods.DeleteSelectedParticlesModifier())

io.export_file(node,'output.data',format="lammps_data",atom_style="angle")

vp = vis.Viewport()
vp.type = vis.Viewport.Type.PERSPECTIVE
vp.camera_pos = (-20,20,-20)
vp.camera_dir = (3,-2,3)
vp.fov = math.radians(60.0)

settings = vis.RenderSettings()
settings.renderer = vis.TachyonRenderer()
settings.renderer.shadows = True
settings.filename = args.dest_folder+'/'+args.imagename
settings.size = (800,600)

node.add_to_scene()
vp.render(settings)
