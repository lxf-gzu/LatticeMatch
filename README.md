# LatticeMatch


## Licence info:
Copyright (C) 2012 Andreas Grois

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## Description
This is a small script that calculates possible lattice matches as discussed in:

Fritz, Torsten: Molecular Architecture in Heteroepitaxially Grown Organic Thin Films
Dresden: sfps - Wissenschaftlicher Fachverlag, 1999
ISBN: 3-934264-50-6


The advantage of this script compared to the original work of T. Fritz is, that this
program does not use a brute-force approach, but instead determines ranges of possible
values for the angle between the first lattice vectors of the interface unit cell of
substrate and adlayer analytically.
To do this, the user has to enter the substrate interface unit cell, giving the length
of the lattice vectors a1 and a2, and also the angle between them, alpha.
For the interface unit cell of the adlayer the user has to input ranges in which the
magnitude of the lattice vectors b1 and b2 may lie, as well as a range for the angle
between them, beta. The program then uses these inequalities together with the
underdefined set of equations describing a coincident lattice match to determine ranges
in which the angle theta between a1 and b1 may lie.

According to the work quoted above, a coincident lattice match occurs when both elements
in one of the columns of the epitaxy matrix are integer, while in case all entries are
integer, one is dealing with a commensurate lattice match.

The epitaxy matrix, as given in the work quoted above, reads:

```math
( px, qy )  
( qx, py )
```
with:
```math
px=b1*sin(alpha-theta)/(a1*sin(alpha))
qx=b2*sin(alpha-theta-beta)/(a1*sin(alpha))
qy=b1*sin(theta)/(a2*sin(alpha))
py=b2*sin(theta+beta)/(a2*sin(alpha))
```
There are two reasons for the length of this program:
o) asin isn't unique.
o) ranges of angles are a pain to deal with.
Both points together made writing the correct conditions for upper and lower bounds
of the result ranges a pain, and I'm pretty certain there might still be the one or the other
bug hidden.

##Usage
Using this program is quite easy: Just supply the input as command line parameters in this order:
a1, a2, alpha, b1min, b1max, b2min, b2max, betamin, betamax
Use degrees for the angles.

The output consists of several ranges of possible values for theta, one range per line.

##Contact
To contact the author either use electronic mail: Andreas Grois <andreas.grois@jku.at>
or write to:
Andreas Grois, Institute for Semiconductor and SolidState physics, Johannes Kepler University,
Altenbergerstraße 69, 4040 Linz, AUSTRIA
