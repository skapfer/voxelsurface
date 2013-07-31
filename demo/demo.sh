#!/bin/bash
# voxelize some surfaces.

# the demo surface is about [-0.05:1.05]^3 cut from Koch's C(Y) surface,
# where the unit cell is [0:1]^3.
# file format is .poly.

# inflated-output = inflate the surface.
# inflated-volfrac = set the volume fraction.
# discret = number of voxels to produce, here 64x64x64
../voxelsurface --surface=koch_CY.tuc.plus.poly --inflated-output=sheet.bin --inflation-volfrac=0.1 --discret=64

# You can use more than one surface, then the union is used.
# Here, the 4srs branched surface [not included]
#../voxelsurface --multi-surface 4YSD.rescaled.surf0.poly,4YSD.rescaled.surf1.poly,4YSD.rescaled.surf2.poly,4YSD.rescaled.surf3.poly --inflated-output=out.bin --discret=64

# filled-output = fill one of the labyrinths.
../voxelsurface --surface=koch_CY.tuc.plus.poly --filled-output=labyrinths.bin --discret=64

