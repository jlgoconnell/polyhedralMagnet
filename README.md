# polyhedralMagnet
Code to calculate interactions between polyhedral permanent magnets including magnetic fields, forces, and torques.

This code requires the use of the third party Geom2d and Geom3d libraries.

## polyhedronForce.m

polyhedronForce.m semi-analytically calculates the forces and torques on one magnet due to the field of the other.

#### The input arguments are:
* verticesA: An (n x 3) matrix of the vertices of the magnet producing the field (magnet A).
* verticesB: An (n x 3) matrix of the vertices of the magnet we are calculating the force and torque on (magnet B).
* magA: The (1 x 3) magnetisation vector of magnet A.
* magB: The (1 x 3) magnetisation vector of magnet B.
* meshnum: A positive integer describing how much the mesh should be refined. A larger meshnum results in a more accurate calculation but takes longer to compute.
* torquept: The (1 x 3) point about which the torque on magnet B should be calculated.
* d: An (n x 3) matrix specifying the range of displacements of magnet B. If this argument is omitted, the function assumes no additional displacements.

#### The output arguments are:
* F: An (n x 3) matrix describing the calculated force(s) on magnet B. This has identical dimensions to the input argument d.
* T: An (n x 3) matrix describing the calculated torque(s) on magnet B. This has identical dimensions to the input argument d.
* t: A vector of length n describing the time it took for each force and torque calculation.

## polyhedronField.m

polyhedronField.m analytically calculates the magnetic field at any point(s) due to a polyhedral permanent magnet.

#### The input arguments are:
* vertices: An (n x 3) matrix of the vertices of the magnet.
* mag: A vector of length 3 describing the magnetisation vector of the magnet.
* obspt: An (n x 3) matrix of points at which the magnetic field should be calculated.
* Fac: Geometric information about the magnet which can be used to speed up calculation time if numerous calculations on the same magnet are performed. Generally this argument should be omitted, in which case the function will generate the geometric information when necessary.

#### The output arguments are:
* B: An (n x 3) matrix describing the magnetic field strength in Teslas at each point in obspt.

## Other functions

trapField.m analytically calculates the magnetic field due to a magnetically charged trapezium. This function is used in polyhedronField.m and should not be used by itself.

trapDecomp.m decomposes a polygon into a series of trapezia to be used by trapField.m. It should not be used by itself.

myDot.m is simply a fast dot product calculator of two 3-length vectors. It is used instead of Matlab's dot() function since it is faster for this specific case.
