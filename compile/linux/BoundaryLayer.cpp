#include "BoundaryLayer.h"

/* 
    This is the constructor for the BoundaryLayer class.
    It reads the input file and initiazes the BoundaryLayer class object. 
*/
BoundaryLayer::BoundaryLayer(std::string filename):
	m_gas(0),
	m_kin(0),
	m_tran(0)
{
	// Read the input file
    readInputFile(filename);
        
	/* Define variable pointers
	u = 0, r2 = 1, p = 2, T = 3 and Ygas = 1 to nsp */
	
	offset_u = 0;
    offset_rSqr = 1;
    offset_p = 2;
    offset_T = 3;
    offset_Ygas = 4;				// offset_Ygas to (offset_Ygas + m_nsp-1)

    // Get total number of variables in the solution vector
    m_nv = offset_Ygas + m_nsp;

    // Get total number of equations
    m_Neq = m_nv * m_psiPoints;

	// Resize function
	resizeArrays();

    // Initialize the thermo state (called only once at the start)
    initializeState();

    // Calculate dimensionless constants at the inlet
    Re0 = rho0*u0*m_Dh/mu0; 
    Pr0 = cp0*mu0/lambda0;
    Ec0 = u0*u0/cp0/deltaT;
    Pe0.resize(m_nsp,0.0);
    for(int i = 0; i< m_nsp; i++)
    {
        Pe0[i] = m_Dh*u0/m_diff0[i];
    }
    // Calculate Schmidt number
    double Sc = mu0/rho0/m_Dh;
    cout<<"\n Re = "<<Re0<<"\t Pr = "<<Pr0<<"\t Ec = "<<Ec0<<"\t Sc = "<<Sc;

    // Calculate thermal resistance
    if(m_wallBC == "FLUX")
    {
        double Rplust = m_Rout + m_th_Wall;
        m_thermalR = log(Rplust/m_Rout)/m_kWall + 1/(Rplust*m_hEnv);
        m_thermalR *= lambda0;
    }

    // Setup radial mesh (in terms of psi)
    
    // Non-dimensional value
    //double m_psi0 = 0.125*(m_Rout-m_Rin)/(m_Rout+m_Rin);
    double m_psi0 = 0.5*(m_Rout/m_Dh)*(m_Rout/m_Dh) - 0.5*(m_Rin/m_Dh)*(m_Rin/m_Dh);
    double ratio;
    m_psi[0] = 0;//0.5*(m_Rin/m_Dh)*(m_Rin/m_Dh);
    for(int  i = 1; i < m_psiPoints; i++)
    {
        ratio = (m_psiPoints - 1 - i)/double((m_psiPoints - 1));
        m_psi[i] = m_psi0 * (1 - pow(ratio, m_stretchFactor));
    }
    /*// Linear interpolation
    double slope_m = (m_Rout*m_Rout - m_Rin*m_Rin)/m_psiPoints;
    double c = m_Rin*m_Rin;
    for(int  i = 0; i < m_psiPoints; i++)
    {
        m_psi[i] = (slope_m*i + c)/2;
    }*/

    for(int  i = 0; i < (m_psiPoints-1); i++)
    {
        m_psi_plus[i] = 0.5 * (m_psi[i] + m_psi[i + 1]);
        if(i > 0)
        {
            m_delta_psi[i] = m_psi_plus[i] - m_psi_plus[i-1];
        } else
        {
            m_delta_psi[0] = 2 * (m_psi_plus[0] - m_psi[0]);
        }
    }

    // Set length in dimensionless form
    solver_tEnd = solver_tEnd/m_Dh;

    // Populate multiplier vector (Required while saving dimensionless solution) 
    multiplier.resize(m_nv, 1.0);
    multiplier[offset_u] = u0;
    multiplier[offset_rSqr] = m_Dh*m_Dh;
    multiplier[offset_p] = pRef;
    multiplier[offset_T] = deltaT;

    // Get initial parameters
    getInitialParams();
}

/* 
    This is the destructor for the BoundaryLayer class.
*/
BoundaryLayer::~BoundaryLayer()
{
	// BoundaryLayer destructor
}

/* 
    This function initializes the arrays and matrices used in calculations.
    It also sets thermodyanamic states of the gas phase based on the initial conditions.
*/
void BoundaryLayer::initializeState()
{
	//Inlet

    if (is_moleFrac == 0)
    {
        m_gas->setMoleFractions_NoNorm(m_xInlet.data());
        m_gas->setState_TP(m_inletTemp, m_inletPres);
    }
    else
    {
        m_gas->setMassFractions_NoNorm(m_yInlet.data());
        m_gas->setState_TP(m_inletTemp, m_inletPres);
    }
    m_gas->getMassFractions(m_yInlet.data());

    // Save molecular weights
	m_Wk = m_gas->molecularWeights();

    // Get thermodynamic properties at the inlet
    rho0 = m_gas->density();
    cp0 = m_gas->cp_mass();
    lambda0 = m_tran->thermalConductivity();
    mu0 = m_tran->viscosity();
    if(flowRate ==2)
    {
        u0 = Re0*mu0/(rho0*m_Dh);
    }
    else {
        u0 = m_inletVel;
    }    
    cout<<"\n Inlet velocity = "<<u0;

    //Calculate inlet flux
    m_inletFlux = u0 * (m_gas->density());
    cout<<"\n m_inletFlux =  "<<m_inletFlux;

    if(m_wallBC == "FLUX")
    {
        deltaT = m_Tenv - m_inletTemp;
    } else 
    {
        deltaT = m_Twall - m_inletTemp;
    }

    pRef = rho0*u0*u0;
    rhou_D = rho0*u0/m_Dh;
    cp_dT = cp0 * deltaT;

    m_diff0.resize(m_nsp);
    m_tran->getMixDiffCoeffs(&m_diff0[0]);
    cout<<"\n The setup has been initialized. \n";
}

/* 
    Residual function:
    This function calculates residual functions used to solve BoundaryLayer DAEs.
	y is the solution vector and yprime is the derivative.
	resid() is passed to a IDA function call.
*/
void BoundaryLayer::residual(double t)
{
	int iSpecies, j, ind;
	double pres, Temp, rho, u, r_sqr;
	double rho_u;
    vector_double wdot(m_nsp, 0.0);
    vector_double hk(m_nsp, 0.0), cpk(m_nsp, 0.0);
    double grad_fluxVel, grad_cond, sumEnthalpyFlux, sumFluxEnergy;
    double Phi = 0.0, cp_mass = 0.0;
    vector_double grad_fluxk(m_nsp, 0.0);

    // Update thermodynamics at each grid point and calculate corresponding fluxes
    updateThermo();

    // Calculate outflow of all species
    for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
    {
        m_outFlow[iSpecies] = 0.0;
        for (j = 0; j < (m_psiPoints-1); j++)
        {
            m_outFlow[iSpecies] += vel(j)*m_rho[j]*m_r[j]*(m_r[j+1] - m_r[j])*Y(j,iSpecies);
        }
        // Converting to dimensional form
        m_outFlow[iSpecies] *= 2* M_PI* u0* (m_Dh*m_Dh);

        // Save outlet mass fraction at r = 0
        m_yEnd[iSpecies] = Y(0,iSpecies);
        m_xEnd[iSpecies] = X(0,iSpecies);
    }    

    // Centerline (or inner surface of annulus)

    m_res[offset_rSqr] = m_rSqr[0] - (m_Rin*m_Rin)/(m_Dh*m_Dh);
    m_res[offset_u] = vel(1) - vel(0);      // Zero velocity gradient
    m_res[offset_p] = (p(1) - p(0));        // Zero pressure gradient
    m_res[offset_T] = (T(1) - T(0));        // Zero temperature gradient

    for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
    {
        m_res[offset_Ygas + iSpecies] = Y(1,iSpecies) - Y(0,iSpecies);   // Zero gradient
    }
    // Interior mesh points
    for (j = 1; j < (m_psiPoints-1); j++)
    {
        ind = j * m_nv;

        // Convert to dimensional parameters
        Temp = T(j) * deltaT + m_inletTemp;
        pres = p(j) * pRef;

        //Species mass fractions
        std::vector<double> yMass(m_y.begin() + ind + offset_Ygas, m_y.begin() + ind + offset_Ygas + m_nsp);
    
        //Set state
        m_gas->setMassFractions_NoNorm(yMass.data());
        m_gas->setState_TP(Temp, pres);

        //Get density using given EoS
        rho = (m_gas->density())/rho0;              // dimensionless
        rho_u = vel(j) * rho;                       // dimensionless

        // Get net progress rates
        m_kin->getNetProductionRates(wdot.data());
        
        // Get partial molar enthalpies
        m_gas->getPartialMolarEnthalpies(hk.data());   //m_gas->getEnthalpy_RT

        // Get partial specific heat
        m_gas->getPartialMolarCp(cpk.data());

        //Get total specific heat
        cp_mass = m_gas->cp_mass();

        // Calculate the gradient of the velocity term
        grad_fluxVel = (m_visc_plus[j]*rhou_r2(j)*du_dpsi(j) 
                        - m_visc_plus[j-1]*rhou_r2(j-1)*du_dpsi(j-1))
                        / (m_delta_psi[j] * mu0);

        // Calculate the gradient of the species flux terms
        for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
        {
            double r1 = 0.5 * (m_r[j] + m_r[j-1]);
            double r2 = 0.5 * (m_r[j] + m_r[j+1]);
            grad_fluxk[iSpecies] = (r2 * m_flux(j,iSpecies) 
                                    - r1 * m_flux(j-1,iSpecies))
                                    / m_delta_psi[j];
        }

        // Calculate the gradient of the conduction term
        grad_cond = (m_lamda_plus[j]*rhou_r2(j)*dT_dpsi(j)
                     - m_lamda_plus[j-1]*rhou_r2(j-1)*dT_dpsi(j-1))
                     / (m_delta_psi[j] * lambda0);
        
        // Calculate dTdpsi at point j
        double dTdpsi = 0.5*(T(j+1) - T(j-1)) / m_delta_psi[j];
        double dudpsi = 0.5*(vel(j+1) - vel(j-1)) / m_delta_psi[j];

        // Calculate energy flux terms
        sumEnthalpyFlux = 0;
        sumFluxEnergy = 0;
        for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
        {
            // Calculate \sum (cp_k * j_k * rho * u * r * dTdpsi)
            sumFluxEnergy += 0.5 * (m_flux(j,iSpecies) + m_flux(j-1,iSpecies))
                             * (cpk[iSpecies]/cp0) * rho_u * m_r[j] * dTdpsi / Pe0[iSpecies];
            // Calculate \sum (h_k * W_k * wdot_k)
            sumEnthalpyFlux += hk[iSpecies] * m_Wk[iSpecies] * wdot[iSpecies] / (rhou_D*cp_dT);
        }

        /*// Calculate viscous dissipation term
        Phi = rho_u * m_r[j] * dudpsi;
        Phi = m_visc[j]/mu0 * Phi * Phi;
        Phi *= (mu0*u0*u0)/(m_D*m_D); */
        
        // Residual functions ->
        
        // Velocity u
        m_res[ind + offset_u] = m_yPrime[ind + offset_u] + m_yPrime[ind + offset_p]/rho_u
                        - grad_fluxVel/Re0;

        // Radial r^2
        m_res[ind + offset_rSqr] = - (m_rSqr[j] - m_rSqr[j-1])/(m_psi[j] - m_psi[j-1]) 
                                    + 4/(rho_u + vel(j-1)*m_rho[j-1]/rho0);
        
        // Pressure
        m_res[ind + offset_p] = p(j+1) - p(j);

        // Mass fractions
        double sumY = 0.0;
        for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
        {
            sumY += Y(j,iSpecies);
            m_res[ind + offset_Ygas + iSpecies] = m_yPrime[ind + offset_Ygas + iSpecies] 
                                                  + (grad_fluxk[iSpecies] * rho_u /Pe0[iSpecies] 
                                                  - m_Wk[iSpecies]* wdot[iSpecies]/rhou_D)
                                                  / rho_u;
        }

        // Last species
        /*int k = m_nsp -1;
        double a = 1 - sumY;
        m_res[ind + offset_Ygas + m_nsp -1] =  a;//(1 - sumY);*/
    
        // Temperature
        m_res[ind + offset_T] = m_yPrime[ind + offset_T] - (Ec0 * vel(j) * m_yPrime[ind + offset_p]
                                + rho_u*grad_cond / (Re0*Pr0) - sumFluxEnergy
                                - sumEnthalpyFlux  + Phi)
                                /(rho_u * cp_mass/cp0);
    }

    // Reactor wall
    j = m_psiPoints -1;
    ind = j * m_nv;
    m_res[ind + offset_u] = vel(j);                         // No slip BC
    m_res[ind + offset_rSqr] = m_rSqr[j] - (m_Rout/m_Dh)*(m_Rout/m_Dh);            // r^2 = Rout^2
    m_res[ind + offset_p] = - (m_rSqr[j] - m_rSqr[j-1])/(m_psi[j] - m_psi[j-1])
                            + 4/(rho_u + vel(j-1)*m_rho[j-1]/rho0);  // implicit BC

    // Flux boundary condition
    if(m_wallBC == "FLUX")
    {
        // Convert to dimensional parameters
        Temp = T(j) * deltaT + m_inletTemp;
        pres = p(j) * pRef;

        //Species mass fractions
        std::vector<double> yMass(m_y.begin() + ind + offset_Ygas, m_y.begin() + ind + offset_Ygas + m_nsp);
    
        //Set state
        m_gas->setMassFractions_NoNorm(yMass.data());
        m_gas->setState_TP(Temp, pres);

        // Get net progress rates
        m_kin->getNetProductionRates(wdot.data());
            
        // Get partial molar enthalpies
        m_gas->getPartialMolarEnthalpies(hk.data());   //m_gas->getEnthalpy_RT

        // Calculate energy flux terms
        sumEnthalpyFlux = 0;
        for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
        {
            // Calculate \sum (h_k * W_k * wdot_k)
            sumEnthalpyFlux += hk[iSpecies] * m_Wk[iSpecies] * wdot[iSpecies] / (rhou_D*cp_dT);
        }
        double m_Rout2 = m_Rout*m_Rout; 
        sumEnthalpyFlux *= Re0*Pr0*((m_Rout+m_th_Wall)*(m_Rout+m_th_Wall) - m_Rout2)/(m_Dh*m_Dh);
        grad_cond = - (m_lamda_plus[j-1]/lambda0)*rho_u*(m_Rout2/(m_Dh*m_Dh))
                    *(T(j)-T(j-1))/(m_psi[j] - m_psi[j-1]);
        m_res[ind + offset_T] = (T(j) - 1)/m_thermalR - grad_cond - sumEnthalpyFlux;
        if(t==0)
        {
            m_Twall_0 = 5 * (m_psi[j] - m_psi[j-1])/((m_lamda_plus[j-1]/lambda0)*rho_u*(m_Rout2/(m_Dh*m_Dh))) + T(j-1);
        }
    } else {
        m_res[ind + offset_T] = T(j) - 1;                 // No slip BC
    }
    // Temperature profile
    if(is_TempProfile)
    {
        Temp = (getT_wall(t) - m_inletTemp)/ deltaT;
        m_res[ind + offset_T] = T(j) - Temp;
    }
    for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
    {
        m_res[ind + offset_Ygas + iSpecies] = Y(j,iSpecies) - Y(j-1,iSpecies);   // Zero gradient
    }
}

/* 
    Set the initial solution vector y0
*/
void BoundaryLayer::getInitialParams()
{
    int iGrid, ind, i, k;
    
    // Calculate r^2 based on the initial condition
    m_rSqr[0] = (m_Rin/m_Dh)*(m_Rin/m_Dh);
    for (iGrid = 1; iGrid < m_psiPoints; iGrid++)
    {
        // Non-dimensional
        m_rSqr[iGrid] = 2 * (m_psi[iGrid] - m_psi[iGrid-1]) + m_rSqr[iGrid-1];
    }
    
    for (iGrid = 0; iGrid < m_psiPoints; iGrid++)
    {
        ind = iGrid*m_nv;

        // Give constant velocity profile
        m_y0[ind + offset_u] = 1; //m_inletVel;

        // Give parabolic profile of initial velocity
        //m_y0[ind + offset_u] = m_inletVel * (1 - m_rSqr[iGrid]/(m_Rout*m_Rout));

        // r^2
        m_y0[ind + offset_rSqr] = m_rSqr[iGrid];

        //species mass fractions
        std::copy(m_yInlet.begin(), m_yInlet.end(), &m_y0[ind + offset_Ygas]);

        //pressure
        m_y0[ind + offset_p] = m_inletPres/pRef;

        // Give linear profile for initial temperature
        //m_y0[ind + offset_T] = m_Twall + (m_Rout - sqrt(m_rSqr[iGrid]))/m_Rout * (m_inletTemp - m_Twall);
        m_y0[ind + offset_T] = 0;//m_inletTemp;
    } 

    // set no-slip boundary at the wall
    m_y0[m_nv*(m_psiPoints-1) + offset_u] = 0.0;

    // Pressure
    m_y0[m_nv*(m_psiPoints-1) + offset_p] = m_inletPres/pRef;
    
    // set T = Twall at the wall
    m_y0[m_nv*(m_psiPoints-1) + offset_T] = 1;// m_Twall;
    // set r2 = R2 at the wall
    m_y0[m_nv*(m_psiPoints-1) + offset_rSqr] = (m_Rout/m_Dh)*(m_Rout/m_Dh); //m_Rout*m_Rout;
    
    //Set m_y = m_y0
    for (i = 0; i < m_Neq; i++)
    {
        m_y[i] = m_y0[i];

        // Set ID to 1 and m_yPrime to 0 as a default
        m_ID[i] = 1;
        m_yPrime0[i] = 0;
        m_yPrime[i] = 0;

        // Set tolerances
        m_atol[i] = solver_ATOL;
    }
    
    //Calculate yPrime0
    residual(0);

    // Set wall temperature based on the residual calculation
    m_y0[m_nv*(m_psiPoints-1) + offset_T] = m_Twall_0;
    //m_y[m_nv*(m_psiPoints-1) + offset_T] = m_y0[m_nv*(m_psiPoints-1) + offset_T];

    for (iGrid = 0; iGrid < m_psiPoints; iGrid++)
    {
        ind = m_nv * iGrid;
        m_atol[ind + offset_rSqr] = 1e-6;
        m_atol[ind + offset_u] = 1e-4;
        for (k = 0; k < m_nv; k++)
        {
            //Set ID and absolute tolerances
            if(iGrid == 0 || iGrid == (m_psiPoints-1))
            {
                //m_atol[ind + 0] = 1e-2*m_inletVel;
                m_ID[ind + k] = 0;
                if(iGrid == (m_psiPoints-1))
                {
                    m_atol[ind + offset_T] = 1e-3;
                }
                /*if(k == offset_p && iGrid == (m_psiPoints-1))
                {
                    m_atol[ind + k] = 1e-6;
                }*/
            }
            else {
                //m_atol[ind + 0] = 1e-2;
                m_ID[ind + offset_p] = 0;
                m_ID[ind + offset_rSqr] = 0;
                if(k != offset_p && k != offset_rSqr)
                {
                    m_yPrime0[ind + k] = -m_res[ind + k];
                }
            }
        }
    }
}

double BoundaryLayer::calculateDEs()
{
    int k = m_gas->speciesIndex(m_DE_speciesName);
    cout<<"\n DE species = "<<m_DE_speciesName;
    double inFlow = m_inletFlux*m_yInlet[k]*M_PI*(m_Rout*m_Rout - m_Rin*m_Rin);  // dimensional form  
    double outFlow =  m_outFlow[k]; // dimensional form
    cout<<"\n Outflow = "<<outFlow<<"\n inflow = "<<inFlow;
    double DE_value = 100*(inFlow-outFlow)/inFlow;
    double DE_centerY = 100*(m_yInlet[k] - m_yEnd[k])/m_yInlet[k];
    //double DE_centerY = 100*(m_inletFlux*m_yInlet[k] - m_yEnd[k]*vel(0)*m_rho[0])/m_inletFlux*m_yInlet[k];
    double DE_centerX = 100*(m_xInlet[k] - m_xEnd[k])/m_xInlet[k];
    cout<<"\n DE_value = "<<DE_value<<"\n DE_centerY = "<<DE_centerY<<"\n DE_centerX = "<<DE_centerX;

    // Write DEs to the file
    std::ofstream outfile;
    outfile.open("DEs_PFAS.dat", ios::app);
    outfile<<"\n";
    outfile<<m_inletTemp<<"\t"<<m_Twall<<"\t"<<Re0<<"\t"<<u0<<"\t"<<DE_value<<"\t"<<DE_centerX<<"\t"<<DE_centerY<<"\t"<<inFlow<<"\t"<<outFlow;
    outfile.close();
    return DE_value;
}

/* Evaluation of required thermodynamic, kinetics properties  */

void BoundaryLayer::updateThermo()
{    
    // Update thermo and required properties at the mesh point
    for (int j = 0; j < m_psiPoints; j++)
    {
        // Update r^2 and r values
        m_rSqr[j] = m_y[j*m_nv + offset_rSqr];
        m_r[j] = sqrt(m_rSqr[j]);
        if(m_rSqr[j] < 0.0)
        {
            m_r[j] = 0.0;
        }
        setGas(j);                  // m_rho is updated after this call.
        if(j < (m_psiPoints-1))
        {
            setGasAtMidpoint(j);    // Update thermo at the interface
        }
    }

    // Update fluxes 
    updateFluxes();
}

void BoundaryLayer::setGas(int j)
{
    int ind = j * m_nv;

    // Convert temperature and pressure to dimensional values
    double Temp = T(j)*deltaT + m_inletTemp;
    double pres = p(j)*pRef;

    //Species mass fractions
    std::vector<double> yMass(m_y.begin() + ind + offset_Ygas, m_y.begin() + ind + offset_Ygas + m_nsp);
    
    //Set state
    for(int l = 0; l < m_nsp; l++)
    {
        if(j == 0)
        {
            //cout<<"\t l = "<<l<<"\t"<<yMass[l];
        }
        if(yMass[l] < 0.0)
        {
            //yMass[l] = 0;
        }
    }

    m_gas->setMassFractions_NoNorm(yMass.data());   
    m_gas->setState_TP(Temp, pres);
    
    //Get density using given EoS
    m_rho[j] = (m_gas->density());

    // Mean mol. wt.
    m_mmw[j] = m_gas->meanMolecularWeight();

    // Get viscosity
    m_visc[j] = m_tran->viscosity();

    // Save mole fractions
    vector_double xCurr(m_nsp, 0.0);
    m_gas->getMoleFractions(xCurr.data());
    for(int k = 0; k < m_nsp; k++)
    {
        m_xMole(k,j) = xCurr[k];
    }
}

void BoundaryLayer::setGasAtMidpoint(int j)
{
    double Temp = 0.5*(T(j)+T(j+1));
    double pres = 0.5*(p(j)+p(j+1));

    // Convert to dimensional values
    Temp = deltaT*Temp + m_inletTemp;
    pres = pRef * pres;
    
    // Mole fractions at the j point
    int ind = j*m_nv + offset_Ygas;
    std::vector<double> yMass1(m_y.begin() + ind, m_y.begin() + ind + m_nsp);

    // Mole fractions at the (j+1) point
    ind = (j+1)*m_nv + offset_Ygas;
    std::vector<double> yMass2(m_y.begin() + ind, m_y.begin() + ind + m_nsp);
    
    // Get average mole fractions
    for (int k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yMass1[k] + yMass2[k]);
    }
    m_gas->setMassFractions_NoNorm(m_ybar.data());
    m_gas->setState_TP(Temp, pres);
    m_mmw_plus[j] = m_gas->meanMolecularWeight();
    m_rho_plus[j] = m_gas->density();
    m_visc_plus[j] = m_tran->viscosity();
    m_lamda_plus[j] = m_tran->thermalConductivity();

    // Calculate diffusion coefficients with current temperature
    m_tran->getMixDiffCoeffs(&m_diff[j*m_nsp]);
}

void BoundaryLayer::updateFluxes()
{
    double const1, xAvg, yAvg, pAvg, sum;
    double rho;
    for (int j = 0; j < (m_psiPoints-1); j++) 
    {
        sum = 0.0;
        rho = m_rho_plus[j]/rho0;
        const1 = - rho*rhou_r(j)/m_mmw_plus[j];
        pAvg = 0.5 * (p(j) + p(j+1));
        for (int k = 0; k < m_nsp; k++) 
        {
            xAvg = 0.5 * (X(j,k) + X(j+1,k));
            yAvg = 0.5 * (Y(j,k) + Y(j+1,k));
            m_flux(j,k) = (X(j+1,k) - X(j,k))/ (m_psi[j+1] - m_psi[j]) 
                          + (xAvg - yAvg)/pAvg * dp_dpsi(j);
            m_flux(j,k) *= const1*m_Wk[k]*m_diff[k+m_nsp*j]/m_diff0[k];
            sum += m_flux(j,k);
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (int k = 0; k < m_nsp; k++) {
            //m_flux(j,k) += sum*Y(j,k);
        }
    }
    for (int k = 0; k < m_nsp; k++) 
    {
        m_flux(m_psiPoints-1,k) = 0.0;
    }
}
/*
	This function reads the mechanism and creates gas object based on the given equation of state
*/
void BoundaryLayer::updateEosType(std::string eos_type)
{
	// Sets type of Eos and read phase
	EosType = 0;

	// Read Gas-Phase
    if(m_filenameGas.empty())
    {
        cout<<"\n Error: No gas phase mechanism file is specified. Program will exit now...\n";
        exit(-1);
    }
    ifstream file_temp;
    file_temp.open(m_filenameGas);
    if(!file_temp)
    {
        cout<<"\n Error: Specified file doesn't exist in current folder. Program will exit now...\n";
        exit(-1);
    }

    cout << "Gas-phase name = "<< m_gasPhase <<endl; 
    if(m_gasPhase.empty())
    {
        cout<<"\n Error: No gas phase is specified. Program will exit now... \n";
        exit(-1);
    }
    if (boost::algorithm::ends_with(m_filenameGas, "xml"))
    {
        XML_Node* xgas = get_XML_File(m_filenameGas); 
		XML_Node* const xg = xgas->findNameID("phase", m_gasPhase);
		m_gas = newPhase(*xg);	
		std::vector<ThermoPhase*> phases1{m_gas};
		m_kin = new GasKinetics();
		importKinetics(*xg, phases1, m_kin);
		m_tran = newDefaultTransportMgr(m_gas);		
    } else if (boost::algorithm::ends_with(m_filenameGas, "yaml"))
    {
        auto sol = newSolution(m_filenameGas, m_gasPhase, "");
        m_gas = newPhase(m_filenameGas, m_gasPhase);
        auto gasTP = sol->thermo();
        auto kin = newKinetics({m_gas}, m_filenameGas, m_gasPhase);
        m_kin = dynamic_cast<GasKinetics*>(kin.release());
        m_tran = newDefaultTransportMgr(m_gas);
    } else 
    {
       cout<<"\n Mechanism file format is not supported."; 
    }
    m_nsp = m_gas->nSpecies();
	cout << "\n Number of gas species = " << m_nsp << endl;
	cout << "Number of gas reactions = " << m_kin->nReactions() << endl;

    // Read the gas phase
	if (eos_type == "IDEAL") {
		EosType = 0;
        //Ideal gas EoS			
		cout << "Using Ideal gas equation of state.... \n";
        
	} /*else if (eos_type == "REDLICHKWONG") {

		EosType = 1;
        cout<<"\n Error: Redlich-Kwong Equation of state is not implemented. Program will exit now..";
        exit(-1);
	} else if (eos_type == "PENGROBINSON") {

		EosType = 1;
        cout<<"\n Error: Peng-Robinson Equation of state is not implemented. Program will exit now..";
        exit(-1);
    */
    else {
        cout<<"\n Error: Equation of state is not specified. Program will exit now..";
        cout<<"\n Valid options are IDEAL or REDLICHKWONG or PENGROBINSON \n";
        exit(-1);
	}
}

/*
	Resize all arrays to their appropriate dimensions. 
*/
void BoundaryLayer::resizeArrays()
{
	//Solution arrays
	m_y.resize(m_Neq,0.0);
	m_yPrime.resize(m_Neq,0.0);
	m_y0.resize(m_Neq,0.0);	
    m_yPrime0.resize(m_Neq, 0.0);
    m_res.resize(m_Neq, 0.0);
    m_ID.resize(m_Neq, 0.0);
    m_atol.resize(m_Neq, 0.0);
    m_flux.resize(m_psiPoints, m_nsp, 0.0);
    m_diff.resize((m_psiPoints-1) * m_nsp, 0.0);
    m_psi.resize(m_psiPoints, 0.0);
    m_psi_plus.resize(m_psiPoints-1, 0.0);
    m_delta_psi.resize(m_psiPoints-1, 0.0);
    m_rSqr.resize(m_psiPoints,0.0);
    m_r.resize(m_psiPoints,0.0);
    m_xk.resize(m_nsp, 0.0);
    m_visc.resize(m_psiPoints, 0.0);
    m_visc_plus.resize(m_psiPoints-1, 0.0);
    m_lamda_plus.resize(m_psiPoints-1, 0.0);
    m_rho.resize(m_psiPoints, 0.0);
    m_rho_plus.resize(m_psiPoints-1, 0.0);
    m_ybar.resize(m_nsp, 0.0);
    m_mmw.resize(m_psiPoints, 0.0);
    m_mmw_plus.resize(m_psiPoints-1, 0.0);
    m_Wk.resize(m_nsp,0.0);
    m_xMole.resize(m_nsp, m_psiPoints, 0.0);
    m_outFlow.resize(m_nsp,0.0);
    m_yEnd.resize(m_nsp,0.0);
    m_xEnd.resize(m_nsp,0.0);
}

/*
	Read the input file
*/
void BoundaryLayer::readInputFile(std::string filename)
{
	std::ifstream input(filename);
	std::map<std::string, double> vars;
	std::string line;
	std::string varName;
	double varValue;
    int iSpecies;
	std::string memSpeciesName;
    std::string eosType;
    std::string solmethod;
	
	cout<<"\nReading from the input file \n";

	// Open file
	while(input.good())
	{
		getline(input, line);
		istringstream ss(line);
        ss >> varName;
        if (varName == "GASFILENAME")
        {
            ss >> m_filenameGas;
        } 
        else if (varName == "GASPHASE")
        {
            ss >> m_gasPhase;
        } 
        else if (varName == "EOSTYPE")
		{
			ss >> eosType;
		}
        else if (varName == "SPECIESFRAC")
        {
            std::string spFrac;
            ss >> spFrac;
            if(spFrac == "MOLE")
            {
                is_moleFrac = 0;
            } 
            else if (spFrac == "MASS")
            {
                is_moleFrac = 1;
            } 
            else {
                cout<<"\n Error: Invalid option to specify input composition";
                cout<<"\n Valid options are MASS or MOLE \n";
                exit(-1);
            }            
        }
        else if (varName == "DE_SPECIES")
        {
            ss >> m_DE_speciesName; 
        }
        else if (varName == "VEL")
        {
            flowRate = 0;
        }
        else if (varName == "SCCM")
        {
            flowRate = 1;      
        }
        else if(varName == "RE")
        {
            flowRate = 2;
        }
        else if (varName == "SOLVEMETHOD")
        {
            ss >> solmethod;
            if(solmethod == "BAND")
            {
                solver_Method = 0;
            } 
            else if (solmethod == "DENSE")
            {
                solver_Method = 1;
            } 
            else {
                cout<<"\n Error: Invalid option to specify solver method";
                cout<<"\n Valid options are BAND or DENSE \n";
                exit(-1);
            }            
        } else if (varName == "TEMPPROFILE_FILE")
        {
            ss >> m_file_TempProfile;
        } else if (varName == "WALL_BC") {
            ss >> m_wallBC;
        }
        ss >> varValue;
        // Read the parameter value
		vars[varName] = varValue;
	}

    //Geometry
    m_Rout = vars["CHANNELRAD"];
    m_Rin = vars["INNER_RAD"];
    m_Dh = 2 * (m_Rout-m_Rin);
    if(m_Rin > 0)
    {
        solve_annulus = 1;
        m_Rin = 0.5*(m_Rout+m_Rin); // Change the reference axis
    }
    if(!withinRange(m_Rout,0))
    {
        cout<<"\n Error: Channel radius can not be negative. \n";
        exit(-1);
    }
    m_psiPoints = vars["NUMPOINTS_PSI"];
    m_stretchFactor = vars["STRETCH_FACTOR"];

	//EoS type
	updateEosType(eosType);

    //Wall boundary condition
    m_Twall = vars["TWALL"];
    m_hEnv = vars["HENV"];
    m_kWall = vars["KWALL"];
    m_th_Wall = vars["THICKNESS_WALL"];
    m_Tenv = vars["TENV"];

	//Boundary condition
	m_yInlet.resize(m_nsp,0.0);
    m_xInlet.resize(m_nsp, 0.0);
	m_inletTemp = vars["INLETTEMP"];
    if(!withinRange(m_inletTemp, 273))
    {
        cout<<"\n Error: Inlet temperature can not be negative. \n";
        exit(-1);
    }
	m_inletPres = vars["INLETPRES"];

    // Temperature profile
    is_TempProfile = vars["TEMP_PROFILE"];
    if(is_TempProfile)
    {
        readCSV(m_file_TempProfile);
    }

    if(!withinRange(m_inletPres,0))
    {
        cout<<"\n Error: Inlet pressure can not be negative. \n";
        exit(-1);
    }
	
    if (flowRate == 0)
	{
		m_inletVel = vars["VEL"];
        cout << "\n Inlet velocity = " << m_inletVel;
	}
	else if (flowRate == 1)
	{
		flowRate_sccm = vars["SCCM"];
		double area_out = M_PI*(m_Rout*m_Rout - m_Rin*m_Rin);

		//Convert to standard velocity
		m_inletVel = flowRate_sccm * 1e-6 / 60 / area_out;		
		//Add pressure and temperature correction
		m_inletVel = m_inletVel * m_inletTemp / 273 * 1e5/m_inletPres;
        cout << "\n Inlet velocity = " << m_inletVel;
	} 
    else if (flowRate == 2) {
        Re0 = vars["RE"];
    } 
    else {
        cout<<"\n Error: Invalid option to specify inlet flow rate";
        cout<<"\n Valid options are SCCM or VEL or RE \n";
        exit(-1);
    }

	// Apply initial conditions    
    for (iSpecies = 0; iSpecies < m_nsp; iSpecies++)
    {
        if (is_moleFrac == 0)
        {
            m_xInlet[iSpecies] = 0.0; 
            m_xInlet[iSpecies] = vars[m_gas->speciesName(iSpecies)];
        }
        else if (is_moleFrac == 1)
        {
            m_yInlet[iSpecies] = 0.0; 
            m_yInlet[iSpecies] = vars[m_gas->speciesName(iSpecies)];
        }
        else
        {
            cout<<"\n Error: No composition is specified";
            cout<<"\n Valid options are MOLE or MASS \n";
            exit(-1);
        }
    }

	// Solver inputs
	solver_ATOL = vars["ATOL"]; 
	solver_RTOL = vars["RTOL"]; 
	solver_tEnd = vars["LENGTH"];  
    solver_initTimeStep = vars["INITDT"];
    solver_maxDt = vars["MAXDT"];
    solver_print = vars["PRINTSTEP"];
    
    cout << "\nFinished reading the input file \n";
	// close file
	input.close();
}

string BoundaryLayer::componentName (int n)
{
    switch (n) {
        case 0:
            return "Velocity";
        case 1:
            return "r2";
        case 2:
            return "Pressure";
        case 3:
            return "Temperature";
        default:
            if (n >= offset_Ygas) {
                return (m_gas->speciesName(n - offset_Ygas));
            } else {
                return "<unknown>";
            }
        }
}
 
/*
	Function to check the input is within the allowed range
*/
int BoundaryLayer::withinRange(double x, double lowerBound, double upperBound)
{
    return(x>=lowerBound && x<=upperBound);
}

/*
	Function to check the input is within the allowed range
*/
int BoundaryLayer::withinRange(double x, double lowerBound)
{
    return(x>=lowerBound);
}

// Temperature profile
double BoundaryLayer::getT_wall(double z)
{
    double Twall;
    double zm, slope, zp, ym, yp;

    // Find nearest index to 't' from m_csvData
    int j = getMeshIndex(z*m_Dh); // Conversion to dimensional value
    //Interpolate to get the Temperature at the location 't'
    zm = m_csvData.at(0).second.at(j);
    ym = m_csvData.at(1).second.at(j);
    if(j != (ncsvData-1))
    {
        zp = m_csvData.at(0).second.at(j+1);
        slope = (z - zm) / (zp - zm);

        yp = m_csvData.at(1).second.at(j+1);
        Twall = ym + slope * (yp - ym);
        return Twall;
    } else {
        return ym;
    }
}

int BoundaryLayer::getMeshIndex(double z)
{
	int i, ind;
	if (z == 0)
	{
		ind = 0;
	}
	else {
		for (i = 0; i < (ncsvData - 1); i++)
		{
			if (m_csvData.at(0).second.at(i+1) >= z)
			{
				ind = i;
				break;
			}
            if (z >= m_csvData.at(0).second.at(ncsvData-1))
            {
                ind = ncsvData-1;
            }
		}
	}
	return ind;
}