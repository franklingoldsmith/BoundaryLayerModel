/**
*  @file PackBed.h
*     Header file for PackBed module
*/
using namespace std;

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>			// Geometry 
#include <stdio.h>      
#include <math.h>     
#include <cmath>
#include <memory>

// Include Boost libraries

#include <boost/config.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/lu.hpp> 
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;

//Include Cantera files
#include <cantera/thermo.h>
#include <cantera/kinetics.h>
#include <cantera/transport.h>
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include <cantera/kinetics/GasKinetics.h>
//#include <cantera/kinetics/ImplicitSurfChem.h>
//#include "cantera/base/fmt.h"
#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"

using namespace Cantera;
typedef std::vector<double> vector_double;

class BoundaryLayer
{
public:
	// Constructor and destructor functions
	BoundaryLayer(std::string filename);
	~BoundaryLayer();

	// Miscellaneous functions
	void readInputFile(std::string);
	void resizeArrays();
	void updateEosType(std::string);
	int withinRange(double x, double lowerBound, double upperBound);
	int withinRange(double x, double lowerBound);
    
    // Solver functions
	void residual(double);
	void getInitialParams();
	
	// Update phase states
	void initializeState();

	// Variable functions
	double vel(int j)
	{
		if(isnan(m_y[m_nv*j + offset_u]))
		{
			cout<<"\n vel is NaN";
		}
		return m_y[m_nv*j + offset_u];
	}

	double T(int j)
	{
		if(isnan(m_y[m_nv*j + offset_T]))
		{
			cout<<"\n j = "<<j<<"\t T is NaN"; 
		}
		return m_y[m_nv*j + offset_T];
	}

	double p(int j)
	{
		if(isnan(m_y[m_nv*j + offset_p]))
		{
			cout<<"\n p is NaN";
		}
		return m_y[m_nv*j + offset_p];
	}

	double Y(int j, int k)
	{
		if(isnan(m_y[m_nv*j + offset_Ygas + k]))
		{
			cout<<"\n yk is NaN";
		}
		return m_y[m_nv*j + offset_Ygas + k];
	}

	double X(int j, int k)
	{
		if(isnan(m_mmw[j]*Y(j,k)/m_Wk[k]))
		{
			cout<<"\n xk is NaN";
		}
		return m_mmw[j]*Y(j,k)/m_Wk[k];
	}

	double rSqr(int j)
	{
		if(isnan(m_y[m_nv*j + offset_rSqr]))
		{
			cout<<"\n r2 is NaN";
		}
		return m_y[m_nv*j + offset_rSqr];
	}

	// Returns rho*u*r2 value at the plus interface
	double rhou_r2(int j)
	{
		double r2, u;
		r2 = 0.5 *(rSqr(j) + rSqr(j + 1));
		u = 0.5 *(vel(j) + vel(j + 1));
		return m_rho_plus[j]*u*r2/rho0;
	}

	// Returns rho*u*r value at the plus interface
	double rhou_r(int j)
	{
		double r, u;
		r = 0.5 *(m_r[j] + m_r[j + 1]);
		u = 0.5 *(vel(j) + vel(j + 1));
		return m_rho_plus[j]*u*r/rho0;
	}

	double du_dpsi(int j)
	{
		if(isnan((vel(j+1) - vel(j))/(m_psi[j+1] - m_psi[j])))
		{
			cout<<"\n du_dpsi is NaN";
		}
		return (vel(j+1) - vel(j))/(m_psi[j+1] - m_psi[j]);
	}

	double dp_dpsi(int j)
	{
		if(isnan((p(j+1) - p(j))/(m_psi[j+1] - m_psi[j])))
		{
			cout<<"\n dp_dpsi is NaN";
		}
		return (p(j+1) - p(j))/(m_psi[j+1] - m_psi[j]);
	}

	double dT_dpsi(int j)
	{
		if(isnan((T(j+1) - T(j))/(m_psi[j+1] - m_psi[j])))
		{
			cout<<"\n dT_dpsi is NaN";
		}
		return (T(j+1) - T(j))/(m_psi[j+1] - m_psi[j]);
	}

	void setGas(int j);
	void setGasAtMidpoint(int j);
	void updateThermo();
	void updateFluxes();

	// Component names
	string componentName(int);

	// Destruction efficiencies
	double calculateDEs();

	// Wall temperature profile
	void readCSV(std::string filename){
        // The data is stored in m_csvData array
        std::ifstream f(filename);
        // Check if file is open
        if(!f.is_open()) throw std::runtime_error("Could not open file");

        std::string line, col;
        double val;

        // Read the first line with column names
        if(f.good())
        {
            // Read the first line in the file
            std::getline(f, line);
            std::stringstream ss(line);

            // Extract column name
            while(std::getline(ss, col, ',')) {
            
                // Initialize and add <colname, double vector> pairs to result
                m_csvData.push_back({col, std::vector<double> {}});
            }
        }

        // Read remaining lines
        while(std::getline(f, line)) {
        
            int ind = 0;
            std::stringstream ss(line);
            
            // Extract data
            while(ss >> val){
            
                // Add the current value to the 'ind' column's values vector
                m_csvData.at(ind).second.push_back(val);
            
                // If the next token is a comma, ignore it and move on
                if(ss.peek() == ',') ss.ignore();
            
                 // Increment the column index
                ind++;
            }
        }

		ncsvData = m_csvData.at(0).second.size();

		// Close file
		f.close();
	}
	int getMeshIndex(double);
	double getT_wall(double);

	/* Variables*/
public:

	/* EosType: returns type of EoS
	0: Ideal gas EoS
	*/

	int EosType = -1;
    
	//Filename:
	string m_filenameGas;								// Gas-phase Mechanism file name
	string m_gasPhase;									// Name of the gas phase
	
	//Geometry:
	double m_Lbed = 0.0;							// Length of the reactor [m]
	double m_Rin = 0;								// Inner radius of the reactor tube
	double m_Rout;									// Outer radius of the reactor tube
	double m_Dh;									// Hydraulic diameter of the reactor tube
	int solve_annulus;                              // Annular geometry
	int m_psiPoints = 0;							// Number of points in radial direction
	double m_stretchFactor = 1.1;						// Stretch factor
	vector_double m_psi_plus, m_delta_psi;
	vector_double multiplier;
	
	// Species related properties
	int m_nsp = 0;					// Number of gas-phase species
	vector_double m_Wk = {};		// Molecular weights of gas-phase species
	Array2D m_xMole = {};		// Mole fractions of gas-phase species
	
	// Cantera objects
	ThermoPhase* m_gas = NULL;
	GasKinetics* m_kin = NULL;
	Transport* m_tran = NULL;

	/* Define variable pointers
	u = 0, r = 1, p = 2, T = 3 and Ygas = 1 to nsp */
	
	int offset_u = 0, offset_rSqr = 1, offset_p = 2, offset_T = 3;
	int offset_Ygas = 4;
	int m_Neq = 0;									// Number of equations
	int m_nv = 0;									// Number of variables for each grid point
	vector_double m_rSqr = {};						// Vector for r^2
	vector_double m_r = {};							// Vector for r
    vector_double m_y = {};						    // Solution vector y
    vector_double m_yPrime = {};					// Vector for yprime
    vector_double m_res = {};					    // Vector for residual
	vector_double m_y0 = {};						// Initial solver vector
    vector_double m_yPrime0 = {};				    // Initial solver vector
    vector_double m_ID = {};                        // Vector for id  (= 1 for differential, 0 for algebraic)
    vector_double m_xk = {};						// Mole fractions
	vector_double m_diff = {};                      // mixture diffusion coefficients
	vector_double m_rho = {};						// Density at the grid point
	vector_double m_rho_plus = {};					// Density at the plus interface
	vector_double m_visc = {};						// Viscosity at the grid point
	vector_double m_visc_plus = {};					// Viscosity at the plus interface
	vector_double m_lamda_plus = {};				// Thermal conductivity at the plus interface
	vector_double m_mmw = {};						// Mean molecular weight at the grid point
	vector_double m_mmw_plus = {};					// Mean molecular weight at the plus interface
	vector_double m_ybar = {};						// Temporary vector 
	Array2D m_flux = {};							// fluxes at the interface

	// Boundary conditions
	double m_Twall_0 = 1;						// Initial guess for the wall temperature
	string m_wallBC = {};					// Heat BC type (Isothermal of flux)
	double m_Tenv = 0; 							// Temperature of the environment
	double m_Twall = 0.0; 						// Wall temperature
	double m_hEnv = 0;							// Convective heat transfer coefficient to the surroundings
	double m_kWall = 0; 						// Conductivity of the wall material
	double m_th_Wall = 0;						// Thickness of the wall 
	double m_thermalR = 0;						// Thermal resistance 
	vector_double m_yInlet = {};			    // Inlet mass fractions
    vector_double m_xInlet = {};			    // Inlet mole fractions
	double m_inletFlux = 0;			            // Total Inlet flux
	vector_double m_outFlow = {};				// Outlet mass flow rate vector 
	vector_double m_yEnd = {};				    // Outlet mass fraction at the centerline 
	vector_double m_xEnd = {};				    // Outlet mole fraction at the centerline 
	double m_inletPres = 0.0;					// Inlet pressure 
	double m_inletVel = 0.0;				    // Inlet velocity
	double m_inletTemp = 0.0;				    // Inlet temperature
	double flowRate_sccm = 0;

	// Flag to specify inlet mass flow rate or velocity or Reynolds number
	int flowRate = -1;							// 0 for velocity, 1 for sccm, 2 for Re
    
    // Add flag for inlet mass or mole fraction
    int is_moleFrac = -1;                       // 0 for mole fractions, 1 for mass fractions

	// Index of the species used for DE calculation
	string m_DE_speciesName;

	// Solver inputs
	double solver_ATOL = 0.0, solver_RTOL = 0.0; 	// Absolute and relative tolerances for the solver
	double solver_tEnd = -1;
	double solver_initTimeStep = -1;
	double solver_maxDt = 1e-6;
	int solver_print = 1;
    int solver_Method = -1;
    int solver_nPass = -1;
    vector_double m_atol = {};                      // absolute tolerance vector
    double m_tCurrent = 0;                          // Current time
	vector_double m_psi = {};                       // Radial coordinate of the grid point
	double pRef, rhou_D, cp_dT;

	// Inlet values
	double mu0, cp0, lambda0, rho0, u0, deltaT;
	double Re0, Pr0, Ec0;
	vector_double m_diff0 = {};
	vector_double Pe0 = {};

	// Temperature profile
	string m_file_TempProfile;
	int is_TempProfile = 0;
	std::vector<std::pair<std::string, std::vector<double>>> m_csvData;
	int ncsvData;
};
