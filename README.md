# Discrete orthogonal geodesic net (DOG) editor
A simple DOG editor, implementing methods from the following papers:

"Discrete geodesic nets for modeling developable surfaces"  
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung
  

"The shape space of discrete orthogonal geodesic nets"  
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung

The content is as follows  :

ext/          -- External dependencies (libigl, Pardiso solver)  
src/          -- Source code  

The implementation needs to solve a sparse linear system, and requires PARDISO.

Due to licensing restrictions, we unfortunately cannot include PARDISO in this archive. To use the PARDISO solver please save the .dylib/.so/.dll file of the latest release (6.0.0) in the directory 'ext/pardiso' and compile using CMake.

To load a planar mesh invoke the binary as follows (e.g. on OSX):  
$ ./dog_edtior planar

To specify a different resolution (height width), load it by:  
$ ./dog_edtior planar 40 30

After editing, one can save the mesh in the editor. The resulting mesh could later be loaded by calling:  
$ ./dog_edtior MESH_PATH.obj