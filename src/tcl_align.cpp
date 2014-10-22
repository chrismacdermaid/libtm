/**
 * @file   tcl_align.cpp
 * @author  <macdercm@mtdoom>
 * @date   Wed Oct  22 17:32:54 2014
 *
 * @brief  routines for aligning two sets of coordinates
 *
 *  See align.h for algorithm details
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <tcl.h>

#include "math_extra.h"
#include "tcl_align.h"
#include "tcl_math.h"

/**
 * @def SIGN
 *
 * @brief returns 0 if x = 0, -1 if x < 0; 1 if x > 0
 */

#define SIGN(x) (( x > 0 ) - ( x < 0 ))

/**
 * @def MMAX
 *
 * Number of rows in matrix to perform SVD on.
 */

#define MMAX 3

/**
 * @def NMAX
 *
 * Number of columns in matrix to perform SVD on.
 */

#define NMAX 3

/**
 * @def MAXITER
 *
 * Maximum number of iterations when doing the alignment
 */

#define MAXITER 1


/**
 * @def R_SMALL
 *
 * A small number that's not 0.
 *
 */

#define R_SMALL 0.000000001

static void normalize3d(double *v) {
    double vlen;
    vlen = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (vlen > R_SMALL) {
        v[0] /= vlen;
        v[1] /= vlen;
        v[2] /= vlen;
    } else {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
    }
}

using namespace MathExtra;

void MatrixFitRMS(int n, double **p, double **q, double *w, double m[4][4]) {

    // Defines
    double mm[3][3],aa[3][3];
    double tol, sig, gam;
    double sg, bb, cc, tmp;
    int a, b, maxiter, iters, ix, iy, iz;

    // Initialize
    identity3(mm);
    zero3(aa);

    /* RMS fit tolerance */
    tol = 1e-15;

    /* maximum number of fitting iterations */
    maxiter = 1000;

    // Calculate COMs
    double com1[3] = { 0.0 };
    double com2[3] = { 0.0 };
    double mass = 0, tm = 0;

    for (int i = 0; i < n; i++) {

        mass = w[i];
        tm += mass;

        // COM of group 1
        com1[0] += mass * p[i][0];
        com1[1] += mass * p[i][1];
        com1[2] += mass * p[i][2];

        // COM of group 2
        com2[0] += mass * q[i][0];
        com2[1] += mass * q[i][1];
        com2[2] += mass * q[i][2];

    }

    if (tm == 0) {
        fprintf(stdout, "Align: Bad mass in COM calculation");
    }

    com1[0] /= tm;
    com1[1] /= tm;
    com1[2] /= tm;

    com2[0] /= tm;
    com2[1] /= tm;
    com2[2] /= tm;

    // Move the coordinates to the COM frame and calculate
    // Correlation
    double cor[MMAX*NMAX] = { 0.0 };
    for (int i = 0; i < n; i++) {
        double x1 = w[i] * (p[i][0] - com1[0]);
        double y1 = w[i] * (p[i][1] - com1[1]);
        double z1 = w[i] * (p[i][2] - com1[2]);

        double x2 = q[i][0] - com2[0];
        double y2 = q[i][1] - com2[1];
        double z2 = q[i][2] - com2[2];

        cor[0] += x2 * x1;
        cor[1] += x2 * y1;
        cor[2] += x2 * z1;
        cor[3] += y2 * x1;
        cor[4] += y2 * y1;
        cor[5] += y2 * z1;
        cor[6] += z2 * x1;
        cor[7] += z2 * y1;
        cor[8] += z2 * z1;
    }

    for (a=0; a<3; a++) for (b=0; b<3; b++) aa[a][b] = cor[3*a+b];

    if(n>1) {
        /* Primary iteration scheme to determine rotation matrix for molecule 2 */
        iters = 0;
        while(1) {
            /* IX, IY, and IZ rotate 1-2-3, 2-3-1, 3-1-2, etc.*/
            iz = (iters+1) % 3;
            iy = (iz+1) % 3;
            ix = (iy+1) % 3;
            sig = aa[iz][iy] - aa[iy][iz];
            gam = aa[iy][iy] + aa[iz][iz];

            if(iters>=maxiter) {
                fprintf(stdout, "Matrix: Warning: no convergence (%1.8f<%1.8f after %d iterations).\n",(float)tol,(float)gam,iters);
            }

            /* Determine size of off-diagonal element.  If off-diagonals exceed the
               diagonal elements * tolerance, perform Jacobi rotation. */
            tmp = sig*sig + gam*gam;
            sg = sqrt(tmp);
            if((sg!=0.0F) &&(fabs(sig)>(tol*fabs(gam)))) {
                sg = 1.0F / sg;
                for(a=0;a<3;a++)
                {
                    bb = gam*aa[iy][a] + sig*aa[iz][a];
                    cc = gam*aa[iz][a] - sig*aa[iy][a];
                    aa[iy][a] = bb*sg;
                    aa[iz][a] = cc*sg;

                    bb = gam*mm[iy][a] + sig*mm[iz][a];
                    cc = gam*mm[iz][a] - sig*mm[iy][a];
                    mm[iy][a] = bb*sg;
                    mm[iz][a] = cc*sg;
                }
            } else {
                break;
            }
            iters++;
        }
    }
    /* At this point, we should have a converged rotation matrix (mm).  Calculate
       the weighted RMS error. */

    // Normalize (is this necessary?)
    normalize3d(mm[0]);
    normalize3d(mm[1]);
    normalize3d(mm[2]);

    // generate the pre/post translation matrix
    double post[4][4];
    double pre[4][4];

    moveby(com2,post);
    moveto(com1,pre);

    // Generate the 4x4 transformation matrix that contains the pre-translation
    double mm4[4][4];
    mm4[0][0] = mm[0][0]; mm4[0][1] = mm[1][0]; mm4[0][2] = mm[2][0]; mm4[0][3] = 0.0;
    mm4[1][0] = mm[0][1]; mm4[1][1] = mm[1][1]; mm4[1][2] = mm[2][1]; mm4[1][3] = 0.0;
    mm4[2][0] = mm[0][2]; mm4[2][1] = mm[1][2]; mm4[2][2] = mm[2][2]; mm4[2][3] = 0.0;
    mm4[3][0] = 0.0;      mm4[3][1] = 0.0;      mm4[3][2] = 0.0;      mm4[3][3] = 1.0;

    // Multiply the stuff out m = [post]*[temp]*[pre]
    double temp4[4][4];
    times4(mm4,pre,temp4);
    times4(post,temp4,m);

    // Print out iterative statistics
    /**
        fprintf(stdout, "\n----------------------- VMD/PYMOL Alignment Statistics ----------------------\n");

        fprintf(stdout, "Number of Iterations = %d\n", iters);

        fprintf(stdout, "\nCOM1: %10.4f %10.4f %10.4f\t COM2: %10.4f %10.4f %10.4f\n",
                com1[0], com1[1], com1[2], com2[0], com2[1], com2[2]);

        fprintf(stdout, "\nCorrelation matrix =\n");
        write3(aa);

        fprintf(stdout, "\nOptimized 4x4 Transformation Matrix = \n");
        write4(m);

        fprintf(stdout, "\n-----------------------------------------------------------------------------\n");
    **/
}

static int tcl_align(ClientData /*clientdata*/, Tcl_Interp *interp,
                     int objc, Tcl_Obj * const objv[])
{

    if (objc != 4) {
        Tcl_SetResult(interp, (char *) "Align expects three arguments: v1 v2 w", TCL_STATIC );
        return TCL_ERROR;
    }

    int N = 0; int M = 0; int NW=0;
    Tcl_Obj **pp = NULL; Tcl_Obj **qq = NULL; Tcl_Obj **ww = NULL;

    if (Tcl_ListObjGetElements(interp, objv[1], &N,  &pp) != TCL_OK || N == 0 ||
        Tcl_ListObjGetElements(interp, objv[2], &M,  &qq) != TCL_OK || M == 0 ||
        Tcl_ListObjGetElements(interp, objv[3], &NW, &ww) != TCL_OK || NW == 0 ||
        N != M || N != NW)
    {

        Tcl_SetResult(interp, (char *) "Badly formed array or dimension mismatch", TCL_STATIC );
        return TCL_ERROR;
    }

    //  Memory for coordinates and weights
    double **p = NULL;
    double **q = NULL;
    double  *w = NULL;

    p = new double*[N];
    for(int i = 0; i < N; ++i)
        p[i] = new double[3];

    q = new double*[M];
    for(int i = 0; i < M; ++i)
        q[i] = new double[3];

    w = new double[N];

    // Get coordinates from TCL lists
    tcl_get_array_obj(objv[1], interp, p);
    tcl_get_array_obj(objv[2], interp, q);
    tcl_get_vector_obj(objv[3], w, interp);

    double m[4][4];
    zero4(m);

    MatrixFitRMS(N,p,q,w,m);

    // Produce TCL matrix
    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][3])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    // Free memory
    for(int i = 0; i < N; ++i)
        delete [] p[i];
    delete [] p;

    for(int i = 0; i < M; ++i)
        delete [] q[i];
    delete [] q;

    delete [] w;

    return TCL_OK;
}


// +---------------------------------+
// | Add tcl commands to interpreter |
// +---------------------------------+

int align_init(Tcl_Interp *interp)
{

    /** Each command is checked against commands currently loaded into the interpretor.
     * By default, existing commands are NOT overwritten to preserve VMDs commands
     * when necessary
     */

    Tcl_CmdInfo cmdinfo;

    // align
    if (Tcl_GetCommandInfo(interp, (char *) "align", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "align", tcl_align,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    return TCL_OK;

}

int align_destroy(Tcl_Interp *interp) {
    Tcl_DeleteCommand(interp, (char *) "align");
}
