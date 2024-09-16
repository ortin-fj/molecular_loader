from ase.io import read
import bpy
import bmesh
import mathutils
import numpy as np
import math
from ase.neighborlist import NeighborList, natural_cutoffs





#####################  SETTINGS OF THE SCRIPT ####################

## Path of the file
file_input='C:/Users/trabajo/Desktop/code/proyectos/molecular_loader/test.xyz'


### In case two atoms were to be bonded but the program could not recognize the bond, you may provide a list of lists containing the ID of the atoms that
### should be linked. An example: 
#explicit = [[74,52],[74,53],[74,54],[74,55],[74,56],[74,52],[74,58],[74,59],[74,60],[74,61],[39,40],[39,43],[39,46],[39,71]]
### Remember that these are Python indexes, i.e., the index of the atoms starts at 0! Your atom Nº 1 in VMD/Pymol/Molden/Jmol... is 0 in this script.

## By default its empty.
explicit=[[]]
if explicit[0]:
    for i in range(len(explicit)):
        j,k=explicit[i]
        explicit.append([k,j])

## Slow factor. Default=1. Increase if the trajectory looks stumbling
slow_factor=1

#################################################################







### Default settings, customize at will
default_colour = (0.0129016, 0.0129016, 0.0129016, 1)
default_radius = 0.25
default_bond_cutoff = 1.7


### Default properties assigned to each atom
### Copied from: https://github.com/nanotech-empa/blender-atomic-loader
default_properties = {
    "H": {"colour": (1.0, 1.0, 1.0, 1.0), "radius": 0.2652},
"He": {"colour": (0.8509803921568627, 1.0, 1.0, 1.0), "radius": 0.2048},
"Li": {"colour": (0.8, 0.5019607843137255, 1.0, 1.0), "radius": 0.9580},
"Be": {"colour": (0.7607843137254902, 1.0, 0.0, 1.0), "radius": 0.6938},
"B": {"colour": (1.0, 0.7098039215686275, 0.7098039215686275, 1.0), "radius": 0.5616},
"C": {"colour": (0.367, 0.367, 0.367, 1.0), "radius": 0.4625},
"N": {"colour": (0.18823529411764706, 0.3137254901960784, 0.9725490196078431, 1.0), "radius": 0.4625},
"O": {"colour": (1.0, 0.050980392156862744, 0.050980392156862744, 1.0), "radius": 0.4625},
"F": {"colour": (0.5647058823529412, 0.8784313725490196, 0.3137254901960784, 1.0), "radius": 0.3304},
"Ne": {"colour": (0.7019607843137254, 0.8901960784313725, 0.9607843137254902, 1.0), "radius": 0.2511},
"Na": {"colour": (0.6705882352941176, 0.3607843137254902, 0.9490196078431372, 1.0), "radius": 1.1893},
"Mg": {"colour": (0.5411764705882353, 1.0, 0.0, 1.0), "radius": 0.9911},
"Al": {"colour": (0.7490196078431373, 0.6509803921568628, 0.6509803921568628, 1.0), "radius": 0.8259},
"Si": {"colour": (0.9411764705882353, 0.7843137254901961, 0.6274509803921569, 1.0), "radius": 0.7268},
"P": {"colour": (1.0, 0.5019607843137255, 0.0, 1.0), "radius": 0.6607},
"S": {"colour": (1.0, 1.0, 0.18823529411764706, 1.0), "radius": 0.6607},
"Cl": {"colour": (0.12156862745098039, 0.9411764705882353, 0.12156862745098039, 1.0), "radius": 0.6607},
"Ar": {"colour": (0.5019607843137255, 0.8196078431372549, 0.8901960784313725, 1.0), "radius": 0.4691},
"K": {"colour": (0.5607843137254902, 0.25098039215686274, 0.8313725490196079, 1.0), "radius": 1.4536},
"Ca": {"colour": (0.23921568627450981, 1.0, 0.0, 1.0), "radius": 1.1893},
"Sc": {"colour": (0.9019607843137255, 0.9019607843137255, 0.9019607843137255, 1.0), "radius": 1.0571},
"Ti": {"colour": (0.7490196078431373, 0.7607843137254902, 0.7803921568627451, 1.0), "radius": 0.9250},
"V": {"colour": (0.6509803921568628, 0.6509803921568628, 0.6705882352941176, 1.0), "radius": 0.8920},
"Cr": {"colour": (0.5411764705882353, 0.6, 0.7803921568627451, 1.0), "radius": 0.9250},
"Mn": {"colour": (0.611764705882353, 0.47843137254901963, 0.7803921568627451, 1.0), "radius": 0.9250},
"Fe": {"colour": (0.8784313725490196, 0.4, 0.2, 1.0), "radius": 0.9250},
"Co": {"colour": (0.9411764705882353, 0.5647058823529412, 0.6274509803921569, 1.0), "radius": 0.8920},
"Ni": {"colour": (0.3137254901960784, 0.8156862745098039, 0.3137254901960784, 1.0), "radius": 0.8920},
"Cu": {"colour": (0.7843137254901961, 0.5019607843137255, 0.2, 1.0), "radius": 0.8920},
"Zn": {"colour": (0.49019607843137253, 0.5019607843137255, 0.6901960784313725, 1.0), "radius": 0.8920},
"Ga": {"colour": (0.7607843137254902, 0.5607843137254902, 0.5607843137254902, 1.0), "radius": 0.8589},
"Ge": {"colour": (0.4, 0.5607843137254902, 0.5607843137254902, 1.0), "radius": 0.8259},
"As": {"colour": (0.7411764705882353, 0.5019607843137255, 0.8901960784313725, 1.0), "radius": 0.7598},
"Se": {"colour": (1.0, 0.6313725490196078, 0.0, 1.0), "radius": 0.7598},
"Br": {"colour": (0.6509803921568628, 0.1607843137254902, 0.1607843137254902, 1.0), "radius": 0.7598},
"Kr": {"colour": (0.3607843137254902, 0.7215686274509804, 0.8196078431372549, 1.0), "radius": 0.5814},
"Rb": {"colour": (0.4392156862745098, 0.1803921568627451, 0.6901960784313725, 1.0), "radius": 1.5527},
"Sr": {"colour": (0.0, 1.0, 0.0, 1.0), "radius": 1.3214},
"Y": {"colour": (0.5803921568627451, 1.0, 1.0, 1.0), "radius": 1.2223},
"Zr": {"colour": (0.5803921568627451, 0.8784313725490196, 0.8784313725490196, 1.0), "radius": 1.0241},
"Nb": {"colour": (0.45098039215686275, 0.7607843137254902, 0.788235294117647, 1.0), "radius": 0.9580},
"Mo": {"colour": (0.32941176470588235, 0.7098039215686275, 0.7098039215686275, 1.0), "radius": 0.9580},
"Tc": {"colour": (0.23137254901960785, 0.6196078431372549, 0.6196078431372549, 1.0), "radius": 0.8920},
"Ru": {"colour": (0.1411764705882353, 0.5607843137254902, 0.5607843137254902, 1.0), "radius": 0.8589},
"Rh": {"colour": (0.0392156862745098, 0.49019607843137253, 0.5490196078431373, 1.0), "radius": 0.8920},
"Pd": {"colour": (0.0, 0.4117647058823529, 0.5215686274509804, 1.0), "radius": 0.9250},
"Ag": {"colour": (0.7529411764705882, 0.7529411764705882, 0.7529411764705882, 1.0), "radius": 1.0571},
"Cd": {"colour": (1.0, 0.8509803921568627, 0.5607843137254902, 1.0), "radius": 1.0241},
"In": {"colour": (0.6509803921568628, 0.4588235294117647, 0.45098039215686275, 1.0), "radius": 1.0241},
"Sn": {"colour": (0.4, 0.5019607843137255, 0.5019607843137255, 1.0), "radius": 0.9580},
"Sb": {"colour": (0.6196078431372549, 0.38823529411764707, 0.7098039215686275, 1.0), "radius": 0.9580},
"Te": {"colour": (0.8313725490196079, 0.47843137254901963, 0.0, 1.0), "radius": 0.9250},
"I": {"colour": (0.5803921568627451, 0.0, 0.5803921568627451, 1.0), "radius": 0.9250},
"Xe": {"colour": (0.25882352941176473, 0.6196078431372549, 0.6901960784313725, 1.0), "radius": 0.7136},
"Cs": {"colour": (0.3411764705882353, 0.09019607843137255, 0.5607843137254902, 1.0), "radius": 1.7179},
"Ba": {"colour": (0.0, 0.788235294117647, 0.0, 1.0), "radius": 1.4205},
"La": {"colour": (0.4392156862745098, 0.8313725490196079, 1.0, 1.0), "radius": 1.2884},
"Ce": {"colour": (1.0, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Pr": {"colour": (0.8509803921568627, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Nd": {"colour": (0.7803921568627451, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Pm": {"colour": (0.6392156862745098, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Sm": {"colour": (0.5607843137254902, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Eu": {"colour": (0.3803921568627451, 1.0, 0.7803921568627451, 1.0), "radius": 1.2223},
"Gd": {"colour": (0.27058823529411763, 1.0, 0.7803921568627451, 1.0), "radius": 1.1893},
"Tb": {"colour": (0.18823529411764706, 1.0, 0.7803921568627451, 1.0), "radius": 1.1562},
"Dy": {"colour": (0.12156862745098039, 1.0, 0.7803921568627451, 1.0), "radius": 1.1562},
"Ho": {"colour": (0.0, 1.0, 0.611764705882353, 1.0), "radius": 1.1562},
"Er": {"colour": (0.0, 0.9019607843137255, 0.4588235294117647, 1.0), "radius": 1.1562},
"Tm": {"colour": (0.0, 0.8313725490196079, 0.3215686274509804, 1.0), "radius": 1.1562},
"Yb": {"colour": (0.0, 0.7490196078431373, 0.2196078431372549, 1.0), "radius": 1.1562},
"Lu": {"colour": (0.0, 0.6705882352941176, 0.1411764705882353, 1.0), "radius": 1.1562},
"Hf": {"colour": (0.30196078431372547, 0.7607843137254902, 1.0, 1.0), "radius": 1.0241},
"Ta": {"colour": (0.30196078431372547, 0.6509803921568628, 1.0, 1.0), "radius": 0.9580},
"W": {"colour": (0.12941176470588237, 0.5803921568627451, 0.8392156862745098, 1.0), "radius": 0.8920},
"Re": {"colour": (0.14901960784313725, 0.49019607843137253, 0.6705882352941176, 1.0), "radius": 0.8920},
"Os": {"colour": (0.14901960784313725, 0.4, 0.5882352941176471, 1.0), "radius": 0.8589},
"Ir": {"colour": (0.09019607843137255, 0.32941176470588235, 0.5294117647058824, 1.0), "radius": 0.8920},
"Pt": {"colour": (0.8156862745098039, 0.8156862745098039, 0.8784313725490196, 1.0), "radius": 0.8920},
"Au": {"colour": (1.0, 0.8196078431372549, 0.13725490196078433, 1.0), "radius": 0.8920},
"Hg": {"colour": (0.7215686274509804, 0.7215686274509804, 0.8156862745098039, 1.0), "radius": 0.9911},
"Tl": {"colour": (0.6509803921568628, 0.32941176470588235, 0.30196078431372547, 1.0), "radius": 1.2554},
"Pb": {"colour": (0.3411764705882353, 0.34901960784313724, 0.3803921568627451, 1.0), "radius": 1.1893},
"Bi": {"colour": (0.6196078431372549, 0.30980392156862746, 0.7098039215686275, 1.0), "radius": 1.0571},
"Po": {"colour": (0.6705882352941176, 0.3607843137254902, 0.0, 1.0), "radius": 1.2554},
"At": {"colour": (0.4588235294117647, 0.30980392156862746, 0.27058823529411763, 1.0), "radius": 0.8391},
"Rn": {"colour": (0.25882352941176473, 0.5098039215686274, 0.5882352941176471, 1.0), "radius": 0.7929},
"Fr": {"colour": (0.25882352941176473, 0.0, 0.4, 1.0), "radius": 1.1562},
"Ra": {"colour": (0.0, 0.49019607843137253, 0.0, 1.0), "radius": 1.1562},
"Ac": {"colour": (0.4392156862745098, 0.6705882352941176, 0.9803921568627451, 1.0), "radius": 1.2884},
"Th": {"colour": (0.0, 0.7294117647058823, 1.0, 1.0), "radius": 1.1893},
"Pa": {"colour": (0.0, 0.6313725490196078, 1.0, 1.0), "radius": 1.1893},
"U": {"colour": (0.0, 0.5607843137254902, 1.0, 1.0), "radius": 1.1562},
"Np": {"colour": (0.0, 0.5019607843137255, 1.0, 1.0), "radius": 1.1562},
"Pu": {"colour": (0.0, 0.4196078431372549, 1.0, 1.0), "radius": 1.1562},
"Am": {"colour": (0.32941176470588235, 0.3607843137254902, 0.9490196078431372, 1.0), "radius": 1.1562},
"Cm": {"colour": (0.47058823529411764, 0.3607843137254902, 0.8901960784313725, 1.0), "radius": 1.1562},
"Bk": {"colour": (0.5411764705882353, 0.30980392156862746, 0.8901960784313725, 1.0), "radius": 1.1562},
"Cf": {"colour": (0.6313725490196078, 0.21176470588235294, 0.8313725490196079, 1.0), "radius": 1.1562},
"Es": {"colour": (0.7019607843137254, 0.12156862745098039, 0.8313725490196079, 1.0), "radius": 1.1562},
"Fm": {"colour": (0.7019607843137254, 0.12156862745098039, 0.7294117647058823, 1.0), "radius": 1.1562},
"Md": {"colour": (0.7019607843137254, 0.050980392156862744, 0.6509803921568628, 1.0), "radius": 1.1562},
"No": {"colour": (0.7411764705882353, 0.050980392156862744, 0.5294117647058824, 1.0), "radius": 1.1562},
"Lr": {"colour": (0.7803921568627451, 0.0, 0.4, 1.0), "radius": 1.1562},
"Rf": {"colour": (0.8, 0.0, 0.34901960784313724, 1.0), "radius": 1.1562},
"Db": {"colour": (0.8196078431372549, 0.0, 0.30980392156862746, 1.0), "radius": 1.1562},
"Sg": {"colour": (0.8509803921568627, 0.0, 0.27058823529411763, 1.0), "radius": 1.1562},
"Bh": {"colour": (0.8784313725490196, 0.0, 0.2196078431372549, 1.0), "radius": 1.1562},
"Hs": {"colour": (0.9019607843137255, 0.0, 0.1803921568627451, 1.0), "radius": 1.1562},
"Mt": {"colour": (0.9215686274509803, 0.0, 0.14901960784313725, 1.0), "radius": 1.1562}
}

### This definition of a function to create a material has also been copied from :
###  https://github.com/nanotech-empa/blender-atomic-loader
def create_basic_material(mat_name='C', colour=default_properties['C']["colour"]):      
    # Switch to Cycles 
    bpy.context.scene.render.engine = 'CYCLES'
    
    # check whether the material already exists
    if bpy.data.materials.get(mat_name):
        mat = bpy.data.materials[mat_name]
    else:
        # create the material
        mat = bpy.data.materials.new(mat_name)
        mat.diffuse_color = colour # viewport color RGBA
        mat.use_nodes = True
    
        # get the material nodes
        nodes = mat.node_tree.nodes
        
        # clear all nodes to start clean
        for node in nodes:
            nodes.remove(node)
        
        # create glossy node
        node_glossy = nodes.new(type='ShaderNodeBsdfGlossy')
        node_glossy.inputs[0].default_value = (1,1,1,1)  # RGBA
        node_glossy.inputs[1].default_value = 0.223 # roughness
        node_glossy.location = -200,100
        
        # create diffuse node
        node_diffuse = nodes.new(type='ShaderNodeBsdfDiffuse')
        node_diffuse.inputs[0].default_value = colour  # RGBA
        node_diffuse.inputs[1].default_value = 0.202 # roughness
        node_diffuse.location = -200,-100
        
        # create mix shader node
        node_mix = nodes.new(type='ShaderNodeMixShader')
        node_mix.inputs[0].default_value = 0.925
        node_mix.location = 0,0
        
        # create output node
        node_output = nodes.new(type='ShaderNodeOutputMaterial')   
        node_output.location = 200,0
        
        # link nodes
        links = mat.node_tree.links
        link_diff_mix = links.new(node_diffuse.outputs[0], node_mix.inputs[2])
        link_gloss_mix = links.new(node_glossy.outputs[0], node_mix.inputs[1])
        
        # link mix to output node
        link_mix_out = links.new(node_mix.outputs[0], node_output.inputs[0])
    
    return mat



### Reads atomic information

mol = read(file_input, index=0)
at_name = mol.get_chemical_symbols()
at_loc = mol.get_positions()
first_lengths = {}

for element in at_name:
    locals()[str(element) + '_mat'] = create_basic_material(mat_name='{}_mat'.format(element), colour=default_properties['{}'.format(element)]["colour"])


################ ATOMS #####################
for i,j in enumerate(at_name):
    atom_name='Atom{}'.format(i+1)
    mesh = bpy.data.meshes.new(atom_name)
    basic_sphere = bpy.data.objects.new(atom_name, mesh)

    # Add the object into the scene.
    bpy.context.collection.objects.link(basic_sphere)
    
    # Select the newly created object
    bpy.context.view_layer.objects.active = basic_sphere
    basic_sphere.select_set(True)
    
    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()

    # create a location matrix
    mat_loc = mathutils.Matrix.Translation(at_loc[i])
    bmesh.ops.create_uvsphere(bm, u_segments=18, v_segments=16, radius=default_properties[at_name[i]]['radius'], matrix=mat_loc)
    bm.to_mesh(mesh)
    bm.free()


# smooth the surface
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.ops.object.shade_smooth()

    # put the origin to the geometry
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')

    # Assign it to the context object
    obj = bpy.context.active_object

#    obj.keyframe_insert(data_path='location', frame=0)
    
    if at_name[i] in default_properties:
        obj.data.materials.append(locals()[str(at_name[i]) + "_mat"])

############# BONDS #######################

ntrl_ctff = natural_cutoffs(mol)
nl = NeighborList(cutoffs=ntrl_ctff)
nl.update(mol)
con_mat=nl.get_connectivity_matrix().asformat('array')

pos1 = list(zip(*np.where(con_mat == 1)))
pos=pos1.copy()

for i,j in pos1:
    pos.append([j,i])
    


for i,j in pos:
    if i == j:
        continue
    ### Atom positions
    x1=at_loc[i][0]
    y1=at_loc[i][1]
    z1=at_loc[i][2]
    x2=at_loc[j][0]
    y2=at_loc[j][1]
    z2=at_loc[j][2]


    ### Bond positions
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1    
    dist = math.sqrt(dx**2 + dy**2 + dz**2)
    
    ### If the bonded atoms are equal (i.e, two carbon atoms) it draws a cylinder of n length, otherwise two cylinders of n/2 length
    if default_properties[at_name[i]]["radius"] < default_properties[at_name[j]]["radius"]:
        curr_rad = default_properties[at_name[i]]["radius"]/2
    else:
        curr_rad = default_properties[at_name[j]]["radius"]/2
        
    ### Centers the cylinder
    parameter = dist - default_properties[at_name[i]]["radius"] - default_properties[at_name[j]]["radius"]
    if at_name[i] == at_name[j]:
        first_length = dist
    else:
        first_length = default_properties[at_name[i]]["radius"] + parameter/2

    
    ### Set the cylinder
    bpy.ops.mesh.primitive_cylinder_add(
    vertices = 32, 
    radius = curr_rad,
    depth = first_length,
    location = (x1+first_length*(dx/dist)/2,  y1+first_length*(dy/dist)/2, z1+first_length*(dz/dist)/2)  
    ) 
       
    entry = first_length
    if 'bond{}_{}'.format(i+1,j+1) not in first_lengths:
        first_lengths['bond{}_{}'.format(i+1,j+1)] = entry
    
    
    ### It lines up the cylinder with the atoms 
    phi = math.atan2(dy, dx) 
    theta = math.acos(dz/dist)    
    bpy.context.object.rotation_euler[1] = theta 
    bpy.context.object.rotation_euler[2] = phi 
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')
    bpy.ops.object.shade_smooth()
    obj = bpy.context.active_object
    obj.name = 'bond{}_{}'.format(i+1,j+1)
    obj.data

    ### Appends the material
    if at_name[i] in default_properties:
        obj.data.materials.append(locals()[str(at_name[i]) + "_mat"])
        


############## KEYFRAMES #################



mol = read(file_input, index="1:")


# Atoms
for i,image in enumerate(mol):
    for j in range(image.get_number_of_atoms()):
        obj = bpy.context.scene.objects['Atom{}'.format(j+1)]
        obj.location = image.get_positions()[j]
        obj.keyframe_insert(data_path = 'location', frame = i*slow_factor)



# Bonds
for i,image in enumerate(mol):
    seen = set()
    at_loc = image.get_positions()
    for j,k in pos:
        pair = frozenset((j,k))
        if pair in seen and at_name[j] == at_name[k]:
            continue
        seen.add(pair)
    
        x1=at_loc[j][0]
        y1=at_loc[j][1]
        z1=at_loc[j][2]
        x2=at_loc[k][0]
        y2=at_loc[k][1]
        z2=at_loc[k][2]
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1    
        dist = math.sqrt(dx**2 + dy**2 + dz**2)
    
    
        parameter = dist - default_properties[at_name[j]]["radius"] - default_properties[at_name[k]]["radius"]
        
        
        if at_name[j] == at_name[k]:
            length = dist
        else:
            length = default_properties[at_name[j]]["radius"] + parameter/2
    
        # Rotation
        phi = math.atan2(dy, dx) 
        theta = math.acos(dz/dist) 
        obj = bpy.context.scene.objects['bond{}_{}'.format(j+1,k+1)]   
        obj.rotation_euler[1] = theta 
        obj.rotation_euler[2] = phi
        obj.keyframe_insert(data_path='rotation_euler', frame = i*slow_factor)
        
        #Location
        obj.location = (x1+length*(dx/dist)/2,  y1+length*(dy/dist)/2, z1+length*(dz/dist)/2)
        obj.keyframe_insert(data_path='location', frame = i*slow_factor)
        
        #Scale
        
        obj.scale[2] = (length/first_lengths['bond{}_{}'.format(j+1,k+1)])
        obj.keyframe_insert(data_path='scale', frame = i*slow_factor)

        
    
    
        