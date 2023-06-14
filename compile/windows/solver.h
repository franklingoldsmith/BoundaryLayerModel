using namespace std;
#include "BoundaryLayer.h"

//Sundial
#include <ida/ida.h>                              /* prototypes for IDA subroutine                */
//#include <ida/ida_dense.h>                        /* access to dense IDA matrix                   */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunmatrix/sunmatrix_dense.h>            /* access to dense SUNMatrix                    */
#include <sunlinsol/sunlinsol_dense.h>            /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h>             /* access to band SUNMatrix                     */
#include <sunlinsol/sunlinsol_band.h>             /* access to band SUNLinearSolver               */
#include <sundials/sundials_types.h>              /* definition of realtype                       */
#include <sundials/sundials_math.h>               /* contains the macros ABS, SUNSQR, and EXP     */ 
#include <sunlinsol/sunlinsol_spgmr.h>            /* access to SPGMR SUNLinearSolver              */
#include <sunlinsol/sunlinsol_spfgmr.h>           /* access to SPGMR SUNLinearSolver              */
#include <sunlinsol/sunlinsol_spbcgs.h>           /* access to SPBCGS SUNLinearSolver             */
#include <sunlinsol/sunlinsol_sptfqmr.h>          /* access to SPTFQMR SUNLinearSolver            */
#include <ida/ida_direct.h>
#include <ida/ida_spils.h>

#define ZERO            RCONST(0.)
#define ONE             RCONST(1.0)

#define BAND            RCONST(0)
#define DENSE           RCONST(1)
#define SPGMR           RCONST(2)
#define SPFGMR          RCONST(3)
#define SPBCGS          RCONST(4)
#define SPTFQMR         RCONST(5)
#define KLU             RCONST(6)

/* Shared Problem Constants */

static double sol_ATOL, sol_RTOL, sol_Nout, sol_tEnd, sol_tStart, sol_initTimeStep, sol_maxTimestep;
static int sol_Method, sol_nPass, sol_nGridSave, sol_print;
static double sol_maxNonlinIters, sol_maxNonlinConvFails, sol_setSuppressAlg, sol_maxord, sol_maxsteps, sol_h0, sol_tstop, sol_maxErrTestFails;

/* Private Helper Functions */

int solveBoundaryLayer(BoundaryLayer);                      // Solve function, calls residual function in PackBed class
void PrintIntro(void);
void PrintHeader(void);
static void WriteOutput(void* mem, N_Vector c, realtype t);
static void writeCSV(void* mem, N_Vector c, std::string f);
static void writeCSV_nonDim(void* mem, N_Vector c, std::string f);
static void writeCSV_axial(void* mem, N_Vector c, std::string f, double t, int i);
static void writeCSV_axialMole(void* mem, N_Vector c, std::string f, double t, int i);
static void PrintErrInfo(int nerr);
static void PrintFinalStats(void* mem, realtype t, N_Vector y, int linsolver);

/* Functions Called by the Solver */

static int res_fun(realtype t, N_Vector y, N_Vector ydot, N_Vector res, void *user_ptr);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

static sunindextype N_eqs;

// Define a PackBed class instance with a local scope
static BoundaryLayer *BL_ptr;
static int nsp, nspSurf;
static int nv, points;

//Initialize solution vector
static int SetInitialProfile(N_Vector y, N_Vector yp, N_Vector tol, N_Vector id);

static int res_count = 0;