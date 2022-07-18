# RandomPlacementOnSurfaceV3

This script was originally written to allow me to randomly place oxygen atoms atop a Au(221) surface. I reworked the script to allow for more options, but there were still several issues that arose with version 2. The latest version has the following upgrades:
- It now places not only atoms but whole molecules reandomly.
- Pre-adsorbed molecules/atoms on the surface can be taken into consideration.
- Molecules/atoms are placed better now.
- The final structure with all original atoms is outputted and not just the surface in the previous version.

To use the script, at least 3 files are needed.
The one file contains the structure.
The second file contains the cell vectors, list of surface atoms and list of adsorbates.
The last set of files are the files containing the atoms/molecules that should be added. These files need to end with the extension .atom.

The layout of the cell vectors file (and must be named "Vectors.txt") should look like:

```
1.7702580000000001E+01    0.0000000000000000E+00    0.0000000000000000E+00
0.0000000000000000E+00    1.7702580000000001E+01    0.0000000000000000E+00
0.0000000000000000E+00    0.0000000000000000E+00    1.7867090000000008E+01
99-144
None
```

The first 3 lines are the cell vectors, the fourth line is the list or range of surface atoms and the last line contains the list or range of adsorbates. If no adsorbates are on the structure, use the word None.

For the .atom files, the layout should look like:

```
N
A x1 y1 z1
B x2 y2 z2
C x3 y3 z3
...
Z xk yk zk
```

where ```N``` is the number of copies of the molecule in this file that should be placed on the surface, ```A-Z``` is the atoms that make up the molecule and next to it, is its coordinates. Make sure that the atom that will be mostly interacting with the surface is zeroed, i.e., it is at the origin, with no other atom having a lower z-value. Also make sure that the x and y coordinates are also zeroed. 

When using the script, it is a general form of:

```
python -i RandomPlacement.py -in <STRUCTURE_FILE_WITH_SURFACE> -out <NAME_OF_OUTPUT_FILE> -plot <Y or N> -mdist <REAL_VALUE_MIN_DISTANCE_BETWEEN_SURFACE_ATOMS_AND_ADSORBED_ATOMS> -xdist <REAL_VALUE_MAX_DISTANCE_BETWEEN_SURFACE_ATOM_AND_ADSORBED_ATOM> -adist <REAL_VALUE_MIN_DISTANCE_BETWEEN_ADSORBED_ATOMS>
```

with ```-plot``` you can decide if you want to make a simple plot of the surface (Y) or not (N). ```-mdist``` is purely the minimum distance between the adsorbed atom and the surface atoms, while ```-xdist``` is the maximum distance, and ```-adist``` is the minimum distance between adsorbed atoms (these being the central atoms). 

The script works by reading in the cell vectors, the list of surface atoms and the list of adsorbates. It then separates the surface atoms and does a convex hull to find the outer most atoms that make up the surface. These atoms/points are used to create a polygon, with this polygon split into a bunch of triangles. Points are then placed ranomly using the equation:

<x4,y4,z4> = a1<x1,y,z1> + a2<x2,y2,z2>

<x3,y3,z3> = <0,0,0>  and a1,a2 ∈ [0,1], and a1+a2<=1.

Once the atoms/molecules have been randomly placed, the script checks that it is within a certain distance to the surface atoms and other adsorbates. It does this on a 3x3x1 unit cell. If it meets the criteria, the molecule is placed and the script conintues.
