# BoundaryLayer
Efficient and robust two-dimensional boundary-layer code for incinerator modeling. 
The governing equations are solved using a Von Mises transformation in a nondimensional setting. 

## Installation
BoundaryLayer model can be run using the pre-existing executables or compiling the code on
user’s local machine from the source. All operating parameters needed to run the model are specified in the 'input.dat' file. 

### Pre-existing compiled versions
1. Download the OS-specific compiled versions of the BoundaryLayer model from the `run' folder. 
2. Copy the ‘input.dat’ file in the same folder.
3. To run the model, doubleclick or type `./boundaryLayer' in the terminal.

### Compiling from the source code
#### Installation requirements
1. Compiled Cantera code (Version 3.0 or more)
2. Python and Boost
3. Sundials
4. Microsoft Visual Studio (Required on Windows)

#### Compilation on Windows
1. Open the file 'SConstruct' from the `compile/windows' folder.
2. Change the directory paths to the local directory paths where Cantera and Boost suite are installed.
3. Open MSVC command prompt in the same folder and run `scons’. The code should compile without any errors.
4. Copy the ‘input.dat’ file in the same folder.
5. Run the command `boundaryLayer.exe’ to run the program.
6. The code can also be run by double-clicking on ‘boundaryLayer.exe’.

#### Compilation on Linux or Mac
1. Open the file 'Makefile' from the `compile/linux' folder.
2. Change the directory paths to the local directory paths where Cantera and Boost suite are installed.
3. Copy the ‘input.dat’ file in the same folder.
4. Open the terminal in the same folder and run `make all’. The code should compile without any errors.
5. Doubleclick or type the command `./boundaryLayer’ to run the program.

## Citation
Please cite this package along with the research article, when used in a scholarly work. 

_`Computationally efficient and robust boundary-layer code for incinerator modeling', G. Kogekar, CF. Goldsmith,  In preparation, 2023_

## License
BoundaryLayer is released under the MIT license; see LICENSE for details.


