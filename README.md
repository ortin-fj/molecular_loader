# molecular_loader
Python script that allows you to easily upload single geometries or molecular dynamics trajectories in Blender.


Welcome to Molecular Loader for Blender! Here you will find a simple code that you may use to load molecules into Blender. Blender already has built-in function to upload molecules, but this code can also upload molecular dynamics trajectories! You can use this script to make beautiful movies for your scientific conferences. You can upload any type of structure file that ASE supports (XYZ,PDB...). 

The package "Molecular Nodes" from Brady Johnston (https://github.com/BradyAJohnston/MolecularNodes) is far superior in this aspect, but this code allows to easily have full control of every single bond in the simulation. 

I would like to acknowledge "Blender Atomic Loader" (https://github.com/nanotech-empa/blender-atomic-loader). I copied the function that creates the material from this code. This code was designed for Blender v.2 and unfortunately it doesn't work anymore (at least in my PC).


##########    WHAT YOU WILL NEED   ##########

1. Blender, versions 3 or 4.
2. Python packages: ase, mathutils, numpy, math. The packages bpy and bmesh should already be included in Blender's API.


###########   HOW TO USE THIS CODE ############

1. This code is designed to be used from the Python API of Blender (tested on Blender 3.6 and Blender 4.1). You need to open the file molecular_loader.py in the text editor of Blender.

2. You need to change the variable 'file_input' to the path of the molecule you want to upload. 

3. You can execute the code!



###########   EXTRA DETAILS ############

1. The list of bonds in the system is created considering the covalent radii of each atom as implemented in ase. However, for some special cases (high-energy structures or species containing metals) not all the bonds may be drawn. In those cases, the variable 'explicit' should be substituted by a list of lists containing the index of the pairs of atoms whose bond you would like to add. For example, if I want to add a bond between my first and third atoms, and another one between my first and fourth atoms:

explicit = [[0,2],[0,3]]

Note that in Python the indexes start at 0!

2. The code is prepared to assign one snapshot of the trajectory to each frame of the movie. If you want to spread the snapshots you can change the value of 'slow_factor'. If slow_factor=n, the code will insert one snapshot every n frames. 

3. Once the keyframes have been assigned for the position of the atoms in the movie, it is  
annoying to change the orientation of the species. Make sure to rotate the species in your structure file prior to creating the final scene.

4. This code is based on a brute-force approach. While large files may be uploaded using this code (I got to upload 700 atoms in 500 snapshots, i.e., thousands of objects in Blender), these tasks may take a while. Loading medium-sized files shouldn't take longer than a minute or two. The package includes an .xyz file with a tioproline molecule (test.xyz) as a test case. In my PC, the loading is instantaneous. 

5. I hope this code is useful for you. Feel free to modify it and share it with more people. If you have any suggestion I will be happy to hear it!
