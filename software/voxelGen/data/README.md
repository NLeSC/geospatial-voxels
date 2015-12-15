Files description:
==================
* **AHN_Noordereiland.txt** is a point cloud data set taken from Actueel Hoogtebestand Nederland [(AHN)](http://www.ahn.nl/index.html).
* **Around_Noordereiland.rar** covers the larger area around Noordereiland as shown in the paper (GIS datasets).
* **Noordereiland_BatchVoxels.txt** is a sample output of our algorithms that contains voxels as tuples of [X,Y,Z,R,G,B].
* **RoofSurface, WallSurface, and GroundSurface** are OBJ files generated from the cityGML files of the municipality of Rotterdam by the FME (Feature Manipulation Engine) software application. Note that some data and attributes have been lost through the process of conversion. We have attributed distinict colours to the objects contained within each file in our voxelization process. The cityGML models are at the Level of Detail 2 (LOD2), which means they represent the roof shape in detail but the walls are only extrusions towards the ground. Read more about [cityGML and semantic level of detail](http://www.citygml.org/?id=1539). Check also the [original cityGML files](http://www.rotterdam.nl/links_rotterdam_3d).
