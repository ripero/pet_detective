# pet_detective
This project provides a software that solves the Pet Detective game,
which is part of the Lumosity suite.

## Compilation
The source code is written in Fortran 2008.
You can compile it with your favourite Fortran 2008 compiler,
or with your favourite Fortran 95 compiler
with adequate partial support for Fortran 2003 and 2008.
In particular, it can be compiled with
[gfortran](https://gcc.gnu.org/wiki/GFortran) v4.8.4 and above
(earlier versions may also be fine).
```bash
$ cd src
$ gfortran pd.F90 -o pd
```
## Input files
In order to run the software, you will need 5 input files.
Some samples are provided in the `samples` directory.

### `lengths.dat`
Pet Detective takes place in a subset
of an embedding rectangular grid.

File `lengths.dat` must contain a single line
with two integer numbers separated by a space.
* The first number should be the number of points
along the horizontal direction of the embedding rectangular grid.
* The second number should be the number of points
along the vertical direction of the embedding rectangular grid.

Most of the higher levels take place in a 4x6 grid;
in that case, the content of file `lengths.dat` should be:
```
4 6
```

### `missing_points.dat`
In the first levels,
some of the nodes of the embedding rectangular grid
may not be displayed.

This information should be encoded in file `missing_points.dat`
with the following syntax:
* First line: a single integer number, the number of missing points.
If all points are present, this number should be a zero.
* One additional line for every missing point.
Each line should contain two integer numbers separated by a space:
   * The first number should be the x coordinate of the missing point.
Horizontal coordinates increase from left to right,
start at one, and increase by unit steps.
   * The second number should be the y coordinate of the missing point.
Vertical coordinates increase from top to bottom
(attention, this is not standard!),
start at one, and increase by unit steps.

For example, if in the second row starting from the top,
the first three points starting from the left are missing,
the content of file `missing_points.dat` should be:
```
3
1 2
2 2
3 2
```

### `missing_links.dat`
Most of the nodes of the grid will have links to all their neighbours -
2 links if they are in a corner, 3 links if they are on a side,
and 4 links if they are in the central part of the grid.
However, some links are usually missing.

This information should be encoded in file `missing_links.dat`
with the following syntax:
* First line: a single integer number, the number of missing links.
If all links between existing nodes are present,
this number should be a zero.
* One additional line for every missing link.
For each link, the coordinates of the two points
between which the link is missing should be written
as a sequence of four integer numbers separated by spaces:
   * The first two numbers should be the x and y coordinates
   of one of the points in the (missing) link.
   * The third and fourth numbers should be the x and y coordinates
   of the other point in the (missing) link.

For example, if the only missing link is
the one between the top left point and the point to its right,
then the contents of file `missing_links.dat` should be:
```
1
1 1   2 1
```

### `car.dat`
This file should always have two lines:
* The first line must consist of a single integer number:
the total distance that the car can travel.
* The second line should be the x andd y coordinates
of the initial position of the car.

For example, if the car can travel a distance of 14 intersections
and starts in the node immediately to the right
of the top left one, the contents of file `car.dat` should be:
```
14
2 1
```

### `pets.dat`
This file contains all the information
related to the pets and the homes.

Its syntax is the following one:
* First line: a single integer number, the total number of pets.
* One additional line for every pet,
each of them consisting of four integer numbers:
   * The first two numbers should be the x and y coordinates
   of the position where the pet is waiting.
   * The third and fourth numbers should be the x and y coordinates
   of the position of the home of the pet.

For example, if the grid is 3x3,
there are two pets in the rightmost and leftmost nodes of the top line,
and their homes are the diagonally opposite nodes of the lattice,
the contents of file `pets.dat` should be:
```
2
1 1   3 3
3 1   1 3
```

## Running the software and understanding the output
Just run the binary from the directory
that contains the five files described above.
```
$ ls
car.dat  lengths.dat  missing_links.dat  missing_points.dat  pets.dat
$ ../src/pd
 1--- 2--- 3--- 4
 |    |         |
 5--- 6--- 7--- 8
 |    |    |     
 9---10   11---12
 |    |         |
13---14---15---16
 |    |    |    |
17   18---19---20
 |    |         |
21   22---23---24
Computing distances between nodes...... done
Number of calls to DFS: 54655
SOLUTION FOUND:
14  15  16  24  23  22  13  17  9  10  18  19  20  11  7  6  5  1  2  4  3
```
The output starts with a sketch of the grid,
with all their nodes and existing links.
Each node is given a number.
If a solution is found, the sequence of nodes
that forms the solution will be provided in the last line.
