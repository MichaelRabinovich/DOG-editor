# Discrete orthogonal geodesic net (DOG) editor
This contains source code for DOG exploration implementing methods from the following ACM TOG and Siggraph Asia papers:

"Discrete geodesic nets for modeling developable surfaces"
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung

"The shape space of discrete orthogonal geodesic nets"
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung

It currently only supports mesh parametrization minimizing the symmetric Dirichlet isometric energy.

The content is as follows

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