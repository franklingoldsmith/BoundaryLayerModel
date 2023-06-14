#include "solver.h"
#include<ctime>
#include <chrono>

int solveBoundaryLayer(BoundaryLayer BL_instance)
{
    using namespace std::chrono;
    
    int nerr = 0, maxIt = 2;
    void* ida_mem;
    N_Vector y, yPrime, absTol, id, constraints;
    int retval, iout = 0;
    realtype tret,tout;
    sunindextype mu, ml;
    SUNMatrix A;
	SUNLinearSolver LS;				//Linear solver used
    SUNContext ctx;

    y = NULL;
    yPrime = NULL;
    absTol = NULL;
    constraints = NULL;
    A = NULL;
    LS = NULL;
    id = NULL;

    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    BL_ptr = &BL_instance;
    N_eqs = BL_ptr->m_Neq;
    nsp = BL_ptr->m_nsp;
    nv = BL_ptr->m_nv;
    points = BL_ptr->m_psiPoints;
    
    sol_ATOL = 1e-18;
    sol_RTOL = 1e-6;
    sol_Nout = 1e10;
    sol_tStart = 0.0;
    sol_tEnd = BL_ptr->solver_tEnd;
    sol_Method = BL_ptr->solver_Method;
    sol_initTimeStep = BL_ptr->solver_initTimeStep;
    sol_maxTimestep = BL_ptr->solver_maxDt;
    sol_print = BL_ptr->solver_print;

    /* Allocate N-vectors. */
    y = N_VNew_Serial(N_eqs, ctx);  //N_VNew_Serial(N_eqs, ctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) return(1);
    yPrime = N_VClone(y);//N_VNew_Serial(N_eqs); //N_VClone(y);
    if (check_retval((void*)yPrime, "N_VNew_Serial", 0)) return(1);
    absTol = N_VClone(y);//N_VNew_Serial(N_eqs); //N_VClone(y);
    if (check_retval((void*)absTol, "N_VNew_Serial", 0)) return(1);
    id = N_VClone(y);//N_VNew_Serial(N_eqs); //N_VClone(y);
    if (check_retval((void*)id, "N_VNew_Serial", 0)) return(1);
    constraints = N_VClone(y); //N_VNew_Serial(N_eqs); //N_VClone(y);
    if (check_retval((void*)constraints, "N_VNew_Serial", 0)) return(1);
    
    /* Set constraints to all 1's for nonnegative solution values. */
    N_VConst(ONE, constraints);

    realtype* cData = N_VGetArrayPointer(constraints);
    for (int i = 0 ; i < points; i++)
    {
        cData[nv * i] = 0;          // No constraint for velocity
    }

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    SetInitialProfile(y, yPrime, absTol, id);

    /* Call IDACreate and IDAMalloc to initialize IDA memory */
    ida_mem = IDACreate(ctx);   //IDACreate(ctx);
    if (check_retval((void*)ida_mem, "IDACreate", 0)) return(1);

    retval = IDASetId(ida_mem, id);
    if (check_retval(&retval, "IDASetId", 1)) return(1);

    retval = IDASetConstraints(ida_mem, constraints);
    if (check_retval(&retval, "IDASetConstraints", 1)) return(1);
    N_VDestroy(constraints);

    retval = IDAInit(ida_mem, res_fun, sol_tStart, y, yPrime);
    if (check_retval(&retval, "IDAInit", 1)) return(1);

    /* Call IDADense and set up the linear solver. */
    //A = SUNDenseMatrix(N_eqs, N_eqs);
    A = SUNDenseMatrix(N_eqs, N_eqs, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(y,A, ctx);
    //LS = SUNDenseLinearSolver(y, A);
    if (check_retval((void*)LS, "SUNDenseLinearSolver", 0)) return(1);

    cout << "\n Dense linear solver is used. \n";

    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(ida_mem, sol_RTOL, absTol);
    if (check_retval(&retval, "IDASStolerances", 1)) return(1);

    /* Attach the matrix and linear solver */
    //retval = IDADlsSetLinearSolver(ida_mem, LS, A);
    retval = IDASetLinearSolver(ida_mem, LS, A);
    if (check_retval(&retval, "IDADlsSetLinearSolver", 1)) return(1);

    /* Call IDACalcIC (with default options) to correct the initial values. */
    /*tout = RCONST(0.001);
    retval = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, tout);
    if(check_retval(&retval, "IDACalcIC", 1)) return(1);*/
  
    // set options

    /*retval = IDASetMaxOrd(ida_mem, 5);
    if (check_retval(&retval, "IDASetMaxOrd", 1)) return(retval);

    retval = IDASetMaxNumSteps(ida_mem, 200000);
    if (check_retval(&retval, "IDASetMaxNumSteps", 1)) return(retval);

    retval = IDASetMaxStep(ida_mem, sol_maxTimestep);
    if (check_retval(&retval, "IDASetMaxStep", 1)) return(retval);

    retval = IDASetInitStep(ida_mem, sol_initTimeStep);
    if (check_retval(&retval, "IDASetInitStep", 1)) return(retval);*/

    retval = IDASetStopTime(ida_mem, sol_tEnd);
    if (check_retval(&retval, "IDASetStopTime", 1)) return(retval);

    /*retval = IDASetMaxErrTestFails(ida_mem, 50);
    if (check_retval(&retval, "IDASetMaxErrTestFails", 1)) return(retval);

    retval = IDASetMaxNonlinIters(ida_mem, 50);
    if (check_retval(&retval, "IDASetMaxNonlinIters", 1)) return(retval);

    retval = IDASetMaxConvFails(ida_mem, 20);
    if (check_retval(&retval, "IDASetMaxConvFails", 1)) return(retval);*/

    tret = 0.0;
    //cout << "\n" << tret << "\t";

    // Write solution at the radial location 'r' at grid point i
    int gridIndex[] = { 0, points/5, 2*points/5, 3*points/5, 4*points/5, points-2};
    std::vector<std::string> file_names(6);
    std::vector<std::string> file_namesMole(6);
    for(int i = 0; i < 6; i++)
    {
        file_names[i] = "i_" + to_string(int(gridIndex[i])) + ".csv";
        file_namesMole[i] = "i_" + to_string(int(gridIndex[i])) + "_mole.csv";
        
        std::ofstream f(file_names[i]);
        // Print titles first
        f << "z" << ", ";
        for(int n = 0; n < nv; n++)
        {
            f << BL_ptr->componentName(n) << ", ";
        }
        f<<"\n";
        f.close();

        std::ofstream f1(file_namesMole[i]);
        // Print titles first
        f1 << "z" << ", ";
        for(int n = 0; n < nv; n++)
        {
            f1 << BL_ptr->componentName(n) << ", ";
        }
        f1<<"\n";
        f1.close();
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (iout = 0, tout = sol_tEnd; iout <= sol_Nout; iout++) 
    {
        retval = IDASolve(ida_mem, tout, &tret, y, yPrime, IDA_ONE_STEP);
        if (check_retval(&retval, "IDASolve", 1)) return(retval);
        cout<<"\n i = "<<iout <<"\t z_nondim = "<<tret;    
            
        if (iout % sol_print == 0)
        {
            std::string filename = "sol_" + to_string(double(tret*(BL_ptr->m_Dh))) + ".csv";
            std::string filename2 = "nondim_sol_" + to_string(double(tret)) + ".csv";
            writeCSV(ida_mem, y, filename);
            writeCSV_nonDim(ida_mem, y, filename2);
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double, std::milli> time_span = t2 - t1;
        // Write solution at the radial location 'r' at grid point i
        for(int  i = 0; i<5; i++)
        {
            writeCSV_axial(ida_mem, y, file_names[i], tret, gridIndex[i]);
            writeCSV_axialMole(ida_mem, y, file_namesMole[i], tret, gridIndex[i]);
        }
        if (retval == IDA_SUCCESS && iout < sol_Nout) {
            //Save output every nth time-step
            if (tret >= sol_tEnd && retval == IDA_SUCCESS)
            {
                double T = BL_ptr-> m_inletTemp;
                double p = BL_ptr-> m_inletPres/1e5;
                double v = BL_ptr-> m_inletVel;
                double Re = BL_ptr->Re0;
                std::string filename = "BLsolution_" + to_string(int(T)) + "K_" + to_string(int(p)) + "Bar_" + to_string(int(Re)) + "vel.csv";
                writeCSV(ida_mem, y, filename);
                std::string filename2 = "nondim_BLsolution_" + to_string(int(T)) + "K_" + to_string(int(p)) + "Bar_" + to_string(int(Re)) + "vel.csv";
                writeCSV_nonDim(ida_mem, y, filename2);

                // Calculate destruction efficiencies (DEs)
                //cout<<"\n Destruction efficiencies: ";
                //cout<<"\n CF4_DE = "<<BL_ptr->calculateDEs("CF4");
                cout<<"\n Simulation time = "<<time_span.count();
                cout<<"\n The program exited successfully.";
                break;
            }
        }
    }

    IDAFree(&ida_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    N_VDestroy(y);
    N_VDestroy(yPrime);
    N_VDestroy(absTol);
    N_VDestroy(id);

    return(0);
}

/* 
    Residual function called by IDA
*/
static int res_fun(realtype t, N_Vector y, N_Vector yPrime, N_Vector res, void *user_ptr)
{
    //cout<<"\t y_vec address"<<&y;
    realtype* yData = NULL,* yPrimeData = NULL,* resData = NULL;
    yData = N_VGetArrayPointer(y);
    yPrimeData = N_VGetArrayPointer(yPrime);
    resData = N_VGetArrayPointer(res);
    
    res_count++;
    int tempCount = 0;
    // Reshape the solution vector
	for (int i = 0; i < N_eqs; i++)
    {
        BL_ptr->m_y[i] = yData[i];
        BL_ptr->m_yPrime[i] = yPrimeData[i];
        BL_ptr->m_res[i] = resData[i];
    }
    
    // Call residual function from the PackBed class
    BL_ptr->residual(t);

	// Reshape the solution vector
    for (int i = 0; i < N_eqs; i++)
    {
        yData[i] = BL_ptr->m_y[i];
        yPrimeData[i] = BL_ptr->m_yPrime[i];
        resData[i] = BL_ptr->m_res[i];
    }
    return(0);
}

/* 
    Set initial profile (i.e. set y0)
*/
static int SetInitialProfile(N_Vector y, N_Vector yp, N_Vector absTol, N_Vector id)
{
    
    int iVar = 0;
    realtype* yData = NULL, *yPrimeData = NULL, *absTolData = NULL, *idData = NULL;
    yData = N_VGetArrayPointer(y);
    yPrimeData = N_VGetArrayPointer(yp);
    absTolData = N_VGetArrayPointer(absTol);
    idData = N_VGetArrayPointer(id);
    
    //Initialize y
    for (int i = 0; i < N_eqs; i++)
    {
        yData[i] = BL_ptr->m_y0[i];
        yPrimeData[i] = BL_ptr->m_yPrime0[i];
        absTolData[i] = BL_ptr->m_atol[i];
        idData[i] = BL_ptr->m_ID[i];
    }
    return(0);
}

// Write final solution to the file				

static void WriteOutput(void* mem, N_Vector ySol, realtype t)
{
	int iVar;
    std::ofstream outfile;
    realtype* yData = N_VGetArrayPointer(ySol);
    outfile.open("BLSolutionSteady.dat");
    int counter = 0;
    for(int i = 0; i < points; i++)
    {
        for(int n = 0; n < nv; n++)
        {
            outfile << yData[counter] << "\t";
            counter++;
        }
        outfile<<"\n";
    }
    outfile.close();
}

static int check_retval(void *returnvalue, const char *funcname, int opt)
{	
	int *retval;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && returnvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	/* Check if retval < 0 */
	else if (opt == 1) {
		retval = (int *)returnvalue;
		if (*retval < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
				funcname, *retval);
			return(1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && returnvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	return(0);
}

static void writeCSV_nonDim(void* mem, N_Vector ySol, std::string filename) 
{
    std::ofstream f(filename);
    realtype* yData = N_VGetArrayPointer(ySol);
    int iVar;

    // Print titles first
    f << "r" << ", ";
    for(int n = 0; n < nv; n++)
    {
        f << BL_ptr->componentName(n) << ", ";
    }
    f<<"\n";
    int counter = 0;
    int offset = BL_ptr->offset_rSqr;
    
    for(int i = 0; i < points; i++)
    {
        // Print radial distance first
        f << sqrt(yData[i * nv + offset]) << ", ";

        for(int n = 0; n < nv; n++)
        {
            f << yData[counter] << ",";
            counter++;
        }
        f<<"\n";
    }
    f.close();
}

static void writeCSV(void* mem, N_Vector ySol, std::string filename) 
{
    std::ofstream f(filename);
    realtype* yData = N_VGetArrayPointer(ySol);
    int iVar;
    double r, val;

    // Print titles first
    f << "r" << ", ";
    for(int n = 0; n < nv; n++)
    {
        f << BL_ptr->componentName(n) << ", ";
    }
    f<<"\n";
    int counter = 0;
    int offset = BL_ptr->offset_rSqr;
    
    for(int i = 0; i < points; i++)
    {
        // Print radial distance first
        r = (BL_ptr->multiplier[offset])* (yData[i * nv + offset]);
        r = sqrt(r);
        f << r << ", ";

        for(int n = 0; n < nv; n++)
        {
            val = yData[counter]*(BL_ptr->multiplier[n]);
            if(n == (BL_ptr->offset_T))
            {
                val += BL_ptr->m_inletTemp;
            }
            f << val << ",";
            counter++;
        }
        f<<"\n";
    }
    f.close();
}

static void writeCSV_axial(void* mem, N_Vector ySol, std::string filename, double t, int i) 
{
    std::ofstream f;
    f.open(filename, std::ios_base::app);
    realtype* yData = N_VGetArrayPointer(ySol);
    double val;
    f << t*(BL_ptr->m_Dh) << ", ";

    for(int n = 0; n < nv; n++)
    {
        val = yData[i*nv + n]*(BL_ptr->multiplier[n]);
        if(n == (BL_ptr->offset_T))
        {
            val += BL_ptr->m_inletTemp;
        }
        f << val << ",";
    }
    f<<"\n";
    f.close();
}

static void writeCSV_axialMole(void* mem, N_Vector ySol, std::string filename, double t, int i) 
{
    std::ofstream f;
    f.open(filename, std::ios_base::app);
    realtype* yData = N_VGetArrayPointer(ySol);
    double val;
    f << t*(BL_ptr->m_Dh) << ", ";

    for(int n = 0; n < nv; n++)
    {
        val = yData[i*nv + n]*(BL_ptr->multiplier[n]);
        if(n == (BL_ptr->offset_T))
        {
            val += BL_ptr->m_inletTemp;
        }
        // Print mole fractions
        int indY = BL_ptr->offset_Ygas;
        if(n >= indY && n < (BL_ptr->m_nsp + indY))
        {
            val = BL_ptr->m_xMole(n - indY,i);
        }
        f << val << ",";
    }
    f<<"\n";
    f.close();
}