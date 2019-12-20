# Discrete orthogonal geodesic net (DOG) editor
A simple DOG editor, implementing methods from the following papers:

"Discrete geodesic nets for modeling developable surfaces"  
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung
  

"The shape space of discrete orthogonal geodesic nets"  
Michael Rabinovich, Tim Hoffmann and Olga Sorkine-Hornung

The editor currently does not support editing of curved folded surfaces, but I will update another version that does soon. Please contact me via email (michael.rabinovich.27@gmail.com) in case you have any questions or issues.

The content is as follows:  

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

This simple editor supports vertex based editing, similar to the one demonstrated here: https://www.youtube.com/watch?v=rd5mg6VsfnA&t=78s  

The user can select vertex handles by first moving to selection mode (by rolling the corresponding 'Edit mode' button or simply pressing 's') and then marking vertices.  

![Alt text](data/images/selection.jpg?raw=true "Selecting vertex handles")  
After vertices have been selected, the user can translate them by moving to translate mode (by rolling the corresponding 'Edit mode' button or simply pressing 'd').  
  
![Alt text](data/images/translation.jpg?raw=true "Translating handles")  
The handles can be cleared/resetted by pressing 'c'. Optimization can be turned off/on by unticking/ticking the "is_optimizing" box or by pressing the space button. This can be useful if one wants to set new constraints after a short deformation, as otherwise the sheet will often quickly flatten.


The default editing scheme use a bending objective and an isometry objective with a high weight. The objective weights can be changed in the menu. Note that one can remove the isometry objective entirely, but in that case it is important to add some weight for the edge regularizing objective.  

To visualize the Gauss Map or the surface rulings, simple change the ViewMode.
![Alt text](data/images/GaussAndRulings.jpg?raw=true "Visualizing the Gauss map of the surface and the DOG discrete rulings")  

