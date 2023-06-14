#include "main.h"

int main()
{
	/* This is the main subroutine which initilizes, solves and 
	then post-processes the boundary-layer model				*/

	cout<<"\nThis code simulates a steady-state, two dimensional boundary-layer model\n \n";

	//Create a BoundaryLayer class instance
	/* The BoundaryLayer constructor is called, where all input values are assigned 
	to BoundaryLayer global variables	*/
	BoundaryLayer BL_instance("input.dat");
	
	//Sundials solve function
	int err = solveBoundaryLayer(BL_instance);

	//Post-processing
	return(0);
}

