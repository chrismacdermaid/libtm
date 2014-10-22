/**
 * @file   tcl_math.cpp
 * @author macdercm <macdercm@mtdoom>
 * @date   Fri Jan 24 15:29:11 2014
 *
 * @brief  Add some rudimentary vector/matrix functions to tcl
 *
 * This is shamelessly ripped and improved on from VMDs implementations
 */

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: TclVec.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.42 $      $Date: 2010/12/16 04:08:43 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A C-based implementation of some performance-critical Tcl callable
 *   routines in VMD.  The C implementation outperforms a raw Tcl version
 *   by a factor of three or so.  The performance advantage helps
 *   significantly when doing analysis in VMD.
 ***************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "tcl.h"
#include "main.h"
#include "tcl_math.h"
#include "math_extra.h"

using namespace TCLMATH_NS;
using namespace MathExtra;

/**
 * @brief calculate the next power of two for an integer
 *
 * 1000 returns 1024, 255 returns 256 etc...
 *
 * @return TCL_OK/TCL_ERROR
 */

static int tcl_nextpow2(ClientData /*clientdata*/, Tcl_Interp *interp,
                        int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "nextpow2 expects a number", "\n", NULL);
        return TCL_ERROR;
    }

    double d1 = 0.0;
    if (Tcl_GetDoubleFromObj(interp, objv[1], &d1) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "nextpow2: non-numeric in first argument", TCL_STATIC );
        return TCL_ERROR;
    }

    double d2 = 0.0;
    d2 = nextpow2(d1);

    Tcl_SetObjResult(interp, Tcl_NewLongObj(d2));
    return TCL_OK;
}


/**
 * @brief returns vector 2 subtracted from vector 1
 *
 * @param clientdata null
 * @param interp tcl interpreter
 * @param argc number of arguments passed == 3
 * @param objv array of arguments pased
 *
 * @return TCL_OK/TCL_ERROR
 */
static int tcl_vecsub(ClientData /*clientdata*/, Tcl_Interp *interp,
                      int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "vecsub expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1] and objv[2]
    int num1 = 0, num2 = 0;
    Tcl_Obj **data1, **data2;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
        return TCL_ERROR;

    // Check to make sure the vectors are equal in length
    if (num1 != num2) {
        Tcl_AppendResult(interp, "vecsub expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
    for (int i = 0; i < num1; i++) {
        double d1 = 0, d2 = 0;
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecsub: non-numeric in first argument", TCL_STATIC );
            return TCL_ERROR;
        }
        if (Tcl_GetDoubleFromObj(interp, data2[i], &d2) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecsub: non-numeric in second argument", TCL_STATIC );
            return TCL_ERROR;
        }
        Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(d1 - d2));
    }
    Tcl_SetObjResult(interp, tcl_result);
    return TCL_OK;
}

/**
 * Calculates the vector sum of two vectors of aribitrary length
 *
 * @param ClientData null
 * @param Tcl_Interp tcl interpreter
 * @param argc number of arguments passed
 * @param argv array of passed arguments
 *
 * @return TCL_OK/TCL_ERROR
 */
// Function:  vecadd v1 v2 {v3 ...}
//  Returns: the sum of vectors; must all be the same length
//  The increase in speed from Tcl to C++ is 4561 / 255 == 18 fold
static int tcl_vecadd(ClientData, Tcl_Interp *interp, int argc,
                      Tcl_Obj * const objv[]) {
    if (argc < 3) {
        Tcl_WrongNumArgs(interp, 1, objv, (char *)"vec1 vec2 ?vec3? ?vec4? ...");
        return TCL_ERROR;
    }
    int num;
    Tcl_Obj **data;
    if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) {
        return TCL_ERROR;
    }
    double *sum = new double[num];
    int i;
    for (i=0; i<num; i++) {
        if (Tcl_GetDoubleFromObj(interp, data[i], sum+i) != TCL_OK) {
            delete [] sum;
            return TCL_ERROR;
        }
    }
    // do the sums on the rest
    int num2;
    for (int term=2; term < argc; term++) {
        if (Tcl_ListObjGetElements(interp, objv[term], &num2, &data) != TCL_OK) {
            delete [] sum;
            return TCL_ERROR;
        }
        if (num != num2) {
            Tcl_SetResult(interp, (char *) "vecadd: two vectors don't have the same size", TCL_STATIC);
            delete [] sum;
            return TCL_ERROR;
        }
        for (i=0; i<num; i++) {
            double df;
            if (Tcl_GetDoubleFromObj(interp, data[i], &df) != TCL_OK) {
                delete [] sum;
                return TCL_ERROR;
            }
            sum[i] += df;
        }
    }

    // and return the result
    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
    for (i=0; i<num; i++) {
        Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(sum[i]));
    }
    Tcl_SetObjResult(interp, tcl_result);
    delete [] sum;
    return TCL_OK;
}

/**
 * Calculates the vector (tensor) product of two vectors of aribitrary length
 * Term-by-term multiply
 *
 * @param ClientData null
 * @param Tcl_Interp tcl interpreter
 * @param argc number of arguments passed
 * @param argv array of passed arguments
 *
 * @return TCL_OK/TCL_ERROR
 */
static int tcl_vecmul(ClientData /*clientdata*/, Tcl_Interp *interp,
                      int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "vecmul expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1] and objv[2]
    int num1 = 0, num2 = 0;
    Tcl_Obj **data1, **data2;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
        return TCL_ERROR;

    // Check to make sure the vectors are equal in length
    if (num1 != num2) {
        Tcl_AppendResult(interp, "vecmul expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
    for (int i = 0; i < num1; i++) {
        double d1 = 0, d2 = 0;
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecmul: non-numeric in first argument", TCL_STATIC );
            return TCL_ERROR;
        }
        if (Tcl_GetDoubleFromObj(interp, data2[i], &d2) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecmul: non-numeric in second argument", TCL_STATIC );
            return TCL_ERROR;
        }
        Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(d1 * d2));
    }
    Tcl_SetObjResult(interp, tcl_result);
    return TCL_OK;
}

static int tcl_vecsum(ClientData /*clientdata*/, Tcl_Interp *interp,
                      int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "vecsum expects a vector\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    double sum = 0.0;
    for (int i = 0; i < num1; i++) {
        double d1 = 0;
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_AppendResult(interp, "vecsum non-numeric vector\n", NULL);
            return TCL_ERROR;

        }

        sum += d1;
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(sum));
    return TCL_OK;
}

/**
 * Calculates the dot product of two vectors of aribitrary length
 *
 * @param ClientData null
 * @param Tcl_Interp tcl interpreter
 * @param argc number of arguments passed
 * @param argv array of passed arguments
 *
 * @return TCL_OK/TCL_ERROR
 */

static int tcl_vecdot(ClientData /*clientdata*/, Tcl_Interp *interp,
                      int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "vecdot expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1] and objv[2]
    int num1 = 0, num2 = 0;
    Tcl_Obj **data1, **data2;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
        return TCL_ERROR;

    // Check to make sure the vectors are equal in length
    if (num1 != num2) {
        Tcl_AppendResult(interp, "vecdot expects two vectors of equal length\n", NULL);
        return TCL_ERROR;
    }

    double ans = 0.0;
    for (int i = 0; i < num1; i++) {
        double d1 = 0, d2 = 0;
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecdot: non-numeric in first argument\n", TCL_STATIC );
            return TCL_ERROR;
        }
        if (Tcl_GetDoubleFromObj(interp, data2[i], &d2) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecdot: non-numeric in second argument\n", TCL_STATIC );
            return TCL_ERROR;
        }
        ans += (d1 * d2);
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(ans));
    return TCL_OK;
}

/**
 * Calculates the cross product of two vectors of length 3
 *
 * @param ClientData null
 * @param Tcl_Interp tcl interpreter
 * @param argc number of arguments passed
 * @param argv array of passed arguments
 *
 * @return TCL_OK/TCL_ERROR
 */

static int tcl_veccross(ClientData /*clientdata*/, Tcl_Interp *interp,
                        int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "veccross expects two vectors of length 3\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1] and objv[2]
    int num1 = 0, num2 = 0;
    Tcl_Obj **data1, **data2;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
        return TCL_ERROR;

    // Check to make sure the vectors are equal in length and have length 3
    if (num1 != num2 || num1 != 3) {
        Tcl_AppendResult(interp, "veccross expects two vectors of length 3\n", NULL);
        return TCL_ERROR;
    }

    double v1[3] = { 0.0 }, v2[3] = { 0.0 }, v3[3] = { 0.0 };

    // Get the vectors
    tcl_get_vector_obj(objv[1], v1, interp);
    tcl_get_vector_obj(objv[2], v2, interp);

    // Take cross product
    cross3(v1, v2, v3);

    // Set the result
    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
    if (Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(v3[0])) != TCL_OK
        || Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(v3[1])) != TCL_OK
        || Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(v3[2])) != TCL_OK)
        return TCL_ERROR;

    Tcl_SetObjResult(interp, tcl_result);

    return TCL_OK;
}

static int tcl_vecinvert(ClientData /*clientdata*/, Tcl_Interp *interp,
                         int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "vecinvert expects two arguments\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    // Set the result
    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);

    double d1 = 0.0;
    for (int i = 0; i < num1; i++)
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK || Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(-d1)) != TCL_OK)
            return TCL_ERROR;

    Tcl_SetObjResult(interp, tcl_result);

    return TCL_OK;
}

static int tcl_vecnorm(ClientData /*clientdata*/, Tcl_Interp *interp,
                       int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "vecinvert expects one vector\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    // Calculate the sum
    double scale = 0.0, d1 = 0.0;
    for (int i = 0; i < num1; i++) {
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_AppendResult(interp, "vecinvert non-numeric vector\n", NULL);
            return TCL_ERROR;

        }

        scale += d1 * d1;
    }

    scale = 1.0 / sqrt(scale);

    // Make an object to store the result
    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);

    // Apply the scaling to the vector elements
    for (int i = 0; i < num1; i++)
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK || Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(d1 * scale)) != TCL_OK)
            return TCL_ERROR;

    // Set the result
    Tcl_SetObjResult(interp, tcl_result);

    return TCL_OK;
}

static int tcl_vecscale(ClientData /*clientdata*/, Tcl_Interp *interp,
                        int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "vecscale expects two arguments\n", NULL);
        return TCL_ERROR;
    }

    double scale = 0.0;
    // Get scale factor
    if (Tcl_GetDoubleFromObj(interp, objv[1], &scale) != TCL_OK)
        return TCL_ERROR;

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[2], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    // Set the result
    Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);

    double d1 = 0.0;
    for (int i = 0; i < num1; i++)
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK || Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(d1 * scale)) != TCL_OK)
            return TCL_ERROR;

    Tcl_SetObjResult(interp, tcl_result);

    return TCL_OK;
}

static int tcl_veclength(ClientData /*clientdata*/, Tcl_Interp *interp,
                         int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "veclength expects one vector\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    // Calculate the sum
    double length = 0.0, d1 = 0.0;
    for (int i = 0; i < num1; i++) {
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_AppendResult(interp, "veclength non-numeric vector\n", NULL);
            return TCL_ERROR;

        }

        length += d1 * d1;
    }

    length = sqrt(length);

    // Set the result
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(length));

    return TCL_OK;
}

static int tcl_veclength2(ClientData /*clientdata*/, Tcl_Interp *interp,
                          int argc, Tcl_Obj * const objv[])
{

    // Check number of arguments
    if (argc != 2) {
        Tcl_AppendResult(interp, "veclength2 expects one vector\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1]
    int num1 = 0;
    Tcl_Obj **data1 = NULL;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;

    // Calculate the sum
    double lengthsq = 0.0, d1 = 0.0;
    for (int i = 0; i < num1; i++) {
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_AppendResult(interp, "veclength2 non-numeric vector\n", NULL);
            return TCL_ERROR;
        }

        lengthsq += d1 * d1;
    }

    // Set the result
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(lengthsq));

    return TCL_OK;
}

static int tcl_vecdist(ClientData /*clientdata*/, Tcl_Interp *interp,
                       int argc, Tcl_Obj * const objv[]) {

    // Check number of arguments
    if (argc != 3) {
        Tcl_AppendResult(interp, "vecdist expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    // Get the individual lists from objv[1] and objv[2]
    int num1 = 0, num2 = 0;
    Tcl_Obj **data1, **data2;
    if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
        return TCL_ERROR;

    // Check to make sure the vectors are equal in length
    if (num1 != num2) {
        Tcl_AppendResult(interp, "vecdist expects two vectors of equal length", "\n", NULL);
        return TCL_ERROR;
    }

    double distsq = 0.0;
    for (int i = 0; i < num1; i++) {
        double d1 = 0, d2 = 0;
        if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecdist: non-numeric in first argument", TCL_STATIC );
            return TCL_ERROR;
        }
        if (Tcl_GetDoubleFromObj(interp, data2[i], &d2) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "vecdist: non-numeric in second argument", TCL_STATIC );
            return TCL_ERROR;
        }

        distsq += (d1 - d2) * (d1 - d2);
    }

    double dist = sqrt(distsq);

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(dist));
    return TCL_OK;
}

static double* obj_getdoublearray(Tcl_Interp *interp, Tcl_Obj *const objv[], int *len)
{
    int num;

    Tcl_Obj **data;
    if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK)
        return NULL;

    double *list = (double*) malloc(num*sizeof(double));
    if (list == NULL)
        return NULL;

    for (int i=0; i<num; i++) {
        double tmp;
        if (Tcl_GetDoubleFromObj(interp, data[i], &tmp) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "veclength: non-numeric in vector", TCL_STATIC);
            free(list);
            return NULL;
        }
        list[i] = tmp;
    }

    *len = num;

    return list;
}

/**
 * @brief Calculate the average of the values
 * @brief in the list
 * @param ClientData null
 * @param interp Interpreter instance
 * @param argc number of passed arguments
 * @param argv Passed arguments
 **/

static int obj_vecmean(ClientData, Tcl_Interp *interp,
                       int argc, Tcl_Obj *const objv[])
{
    if (argc != 2) {
        Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
        return TCL_ERROR;
    }

    int num;
    double *list = obj_getdoublearray(interp, objv, &num);
    if (list == NULL)
        return TCL_ERROR;

    double sum = 0.;
    for (int i=0; i<num; i++) {
        sum += list[i];
    }
    sum /= (double) num;
    free(list);

    Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
    Tcl_SetDoubleObj(tcl_result, sum);
    return TCL_OK;
}

static int obj_vecstddev(ClientData, Tcl_Interp *interp,
                         int argc, Tcl_Obj *const objv[])
{
    if (argc != 2) {
        Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
        return TCL_ERROR;
    }

    int i, num;
    double* list = obj_getdoublearray(interp, objv, &num);
    if (list == NULL)
        return TCL_ERROR;

    double mean = 0.;
    for (i=0; i<num; i++) {
        mean += list[i];
    }
    mean /= (double) num;

    double stddev = 0.;
    for (i=0; i<num; i++) {
        double tmp = list[i] - mean;
        stddev += tmp * tmp;
    }
    stddev /= (double) num;
    stddev = sqrt(stddev);
    free(list);

    Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
    Tcl_SetDoubleObj(tcl_result, stddev);
    return TCL_OK;
}

/**
 * @brief calculate transpose of 4x4 matrix
 *
 * @param ClientData
 * @param interp tcl interpreter
 * @param argc number of arguments passed
 * @param objv tcl obj vector
 *
 * @return TCL_OK, TCL_ERROR
 */

static int tcl_transtranspose(ClientData, Tcl_Interp *interp,
                              int argc, Tcl_Obj * const objv[])
{

    /**
     * objv[0] = transtranspose
     * objv[1] = 4x4 matrix
     */

    if (argc < 2) {
        Tcl_AppendResult(interp, "transtranspose: expected a 4x4 matrix", "\n", NULL);
        return TCL_ERROR;
    }

    double m1[4][4];
    double m2[4][4];
    zero4(m1);
    zero4(m2);

    // get matrix
    if (tcl_get_matrix("transtranspose:", interp, objv[1], m1) != TCL_OK) {
        return TCL_ERROR;
    }

    transpose4(m1, m2);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m2[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m2[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m2[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m2[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

/**
 * @brief produce 4x4 rotation matrix for rotation about axis a given magnitude
 *
 * @param ClientData null
 * @param interp tcl interpreter
 * @param argc number of arguments passed
 * @param argv array of arguments passed
 *
 * @return TCL_OK/TCL_ERROR
 */

static int tcl_transaxis(ClientData, Tcl_Interp *interp,
                         int argc, Tcl_Obj * const objv[])
{

    /**
     * objv[0] = transaxis
     * objv[1] = <x,y,z>
     * objv[2] = amount
     * objv[3] = unit <deg rad>
     */

    if (argc < 3) {
        Tcl_AppendResult(interp, "transaxis: not enough arguments", "\n", NULL);
        return TCL_ERROR;
    }

    // Check requested axis and set vector
    double u[3] = { 0.0 };
    const char *axis = Tcl_GetStringFromObj(objv[1], NULL);

    if (strcmp(axis, "x") == 0) {
        u[0] = 1.0;
    } else if (strcmp(axis, "y") == 0) {
        u[1] = 1.0;
    } else if (strcmp(axis, "z") == 0) {
        u[2] = 1.0;
    } else {
        Tcl_AppendResult(interp, "transaxis: axis must be x, y, or z", "\n", NULL);
        return TCL_ERROR;
    }

    // Magnitude of rotation
    double theta = 0.0;
    if (Tcl_GetDoubleFromObj(interp, objv[2], &theta) != TCL_OK)
        return TCL_ERROR;

    //check if the user specified units, if not assume input is
    //in degrees
    if (argc == 4) {
        const char *units = Tcl_GetStringFromObj(objv[3], NULL);
        if (strcmp(units, "deg") == 0) {
            theta *= DEG2RAD;
        } else if (strcmp(units, "rad") == 0) {
            // Do nothing
        } else {
            Tcl_AppendResult(interp, "transaxis: rotation units should be deg or rad", "\n", NULL);
            return TCL_ERROR;
        }
    } else {
        // Assume Degrees
        theta *= DEG2RAD;
    }

    // 4x4 vector, zero it
    double m[4][4];
    zero4(m);

    // Translation vector
    double v[3] = { 0.0 };

    //Get the rotation matrix
    axis_angle_to_mat_trans4(theta, u, v, m);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

static int tcl_transabout(ClientData, Tcl_Interp *interp,
                          int argc, Tcl_Obj * const objv[])
{

    /**
     * objv[0] = transabout
     * objv[1] = u {0.0 0.0 0.0}
     * objv[2] = amount
     * objv[3] = quantity <deg rad>
     */

    if (argc < 4) {
        Tcl_AppendResult(interp, "transabout: not enough arguments", "\n", NULL);
        return TCL_ERROR;
    }

    // Convert rotation vector from TCL to double vector
    double u[3] = { 0.0 };
    if (tcl_get_vector_obj(objv[1], u, interp) != TCL_OK) {
        Tcl_AppendResult(interp, "transabout: misformed rotational vector", "\n", NULL);
        return TCL_ERROR;
    }

    // Magnitude of rotation
    double theta = 0.0;
    if (Tcl_GetDoubleFromObj(interp, objv[2], &theta) != TCL_OK)
        return TCL_ERROR;

    //check if the user specified units, if not assume input is
    //in degrees
    if (argc == 4) {
        const char *units = Tcl_GetStringFromObj(objv[3], NULL);
        if (strcmp(units, "deg") == 0) {
            theta *= DEG2RAD;
        } else if (strcmp(units, "rad") == 0) {
            // Do nothing
        } else {
            Tcl_AppendResult(interp, "transabout: rotation units should be deg or rad", "\n", NULL);
            return TCL_ERROR;
        }
    } else {
        // Assume Degrees
        theta *= DEG2RAD;
    }

    // 4x4 vector, zero it
    double m[4][4];
    zero4(m);

    // Translation vector
    double v[3] = { 0.0 };

    //Get the rotation matrix
    axis_angle_to_mat_trans4(theta, u, v, m);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

static int tcl_transvec(ClientData /*clientdata*/, Tcl_Interp *interp,
                        int argc, Tcl_Obj * const objv[])
{

    /**
     * objv[0] = transvec
     * objv[1] = 'x', 'y', or 'z'
     * objv[2] = 3-vector, <x,y,z>
     */

    if (argc < 3) {
        Tcl_AppendResult(interp, "transvec: not enough arguments", "\n", NULL);
        return TCL_ERROR;
    }

    const char *axis = Tcl_GetStringFromObj(objv[1], NULL);

    // Convert vector from TCL to double vector
    double v[3] = { 0.0 };
    if (tcl_get_vector_obj(objv[2], v, interp) != TCL_OK) {
        Tcl_AppendResult(interp, "transvec: misformed rotational vector", "\n", NULL);
        return TCL_ERROR;
    }

    double m[4][4];

    axistovec(axis, v, m);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

static int tcl_transvecinv(ClientData /*clientdata*/, Tcl_Interp *interp,
                           int argc, Tcl_Obj * const objv[])
{

    /**
     * objv[0] = transvecinv
     * objv[1] = 'x', 'y', or 'z'
     * objv[2] = 3-vector, <x,y,z>
     */

    if (argc < 3) {
        Tcl_AppendResult(interp, "transvecinv: not enough arguments", "\n", NULL);
        return TCL_ERROR;
    }

    const char *axis = Tcl_GetStringFromObj(objv[1], NULL);

    // Convert vector from TCL to double vector
    double v[3] = { 0.0 };
    if (tcl_get_vector_obj(objv[2], v, interp) != TCL_OK) {
        Tcl_AppendResult(interp, "transvecinv: misformed rotational vector", "\n", NULL);
        return TCL_ERROR;
    }

    double m[4][4];

    vectoaxis(axis, v, m);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

static int tcl_matmult(ClientData /*clientdata*/, Tcl_Interp *interp,
                       int argc, Tcl_Obj * const objv[])
{

    if (argc < 3) {
        Tcl_AppendResult(interp, "matmult: expected two 4x4 matrices", "\n", NULL);
        return TCL_ERROR;
    }

    double m1[4][4];
    double m2[4][4];
    double m3[4][4];
    zero4(m1);
    zero4(m2);
    zero4(m3);

    // get matrix 1
    if (tcl_get_matrix("matmult:", interp, objv[1], m1) != TCL_OK) {
        return TCL_ERROR;
    }

    // get matrix 2
    if (tcl_get_matrix("matmult:", interp, objv[2], m2) != TCL_OK) {
        return TCL_ERROR;
    }

    // Multiply m1 on m2 return in m3
    times4(m1, m2, m3);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    for (int i = 0; i < 4; i++) {
        Tcl_Obj *innerlist = Tcl_NewListObj(0, NULL);
        if (Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m3[i][0])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m3[i][1])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m3[i][2])) != TCL_OK
            || Tcl_ListObjAppendElement(interp, innerlist, Tcl_NewDoubleObj(m3[i][3])) != TCL_OK || Tcl_ListObjAppendElement(interp, outerlist, innerlist) != TCL_OK)
            return TCL_ERROR;
    }

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

static int tcl_matvec(ClientData /*clientdata*/, Tcl_Interp *interp,
                      int argc, Tcl_Obj * const objv[])
{

    if (argc < 3) {
        Tcl_AppendResult(interp, "matmult: expected a 4x4 matrix and a vector length 4\n", NULL);
        return TCL_ERROR;
    }

    double m1[4][4];
    double v1[4] = { 0.0 }, v2[4] = { 0.0 };
    zero4(m1);

    v1[3] = v2[3] = 1.0;

    // get matrix
    if (tcl_get_matrix("matmult:", interp, objv[1], m1) != TCL_OK) {
        return TCL_ERROR;
    }

    // get vector
    if (tcl_get_vector_obj(objv[2], v1, interp) != TCL_OK) {
        return TCL_ERROR;
    }

    // Multiply m1 on m2 return in m3
    matvec4(m1, v1, v2);

    // Outer list to hold 4 inner lists
    Tcl_Obj *outerlist = Tcl_NewListObj(0, NULL);

    if (Tcl_ListObjAppendElement(interp, outerlist, Tcl_NewDoubleObj(v2[0])) != TCL_OK
        || Tcl_ListObjAppendElement(interp, outerlist, Tcl_NewDoubleObj(v2[1])) != TCL_OK
        || Tcl_ListObjAppendElement(interp, outerlist, Tcl_NewDoubleObj(v2[2])) != TCL_OK)
        return TCL_ERROR;

    // Set result to constructed list
    Tcl_SetObjResult(interp, outerlist);

    return TCL_OK;
}

// +------------------+
// |                  |
// | HELPER FUNCTIONS |
// |                  |
// +------------------+

/**
 * @brief parses tcl string containg a list to a double vector
 *
 * @param s the string containg the tcl list
 * @param val parsed vector in double vector (returned)
 * @param interp tcl interpreter
 *
 * @return TCL_OK/TCL_ERROR
 */

int tcl_get_vector(const char *s, double *val, Tcl_Interp *interp)
{

    int num;
    const char **pos;
    if (Tcl_SplitList(interp, s, &num, &pos) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "need three data elements for a vector", TCL_STATIC );
        return TCL_ERROR;
    }
    if (num != 3) {
        Tcl_SetResult(interp, (char *) "need three numbers for a vector", TCL_STATIC );
        return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, pos[0], val + 0) != TCL_OK || Tcl_GetDouble(interp, pos[1], val + 1) != TCL_OK || Tcl_GetDouble(interp, pos[2], val + 2) != TCL_OK) {
        ckfree((char *) pos);
        // free of tcl data
        return TCL_ERROR;
    }

    ckfree((char *) pos);
    // free of tcl data
    return TCL_OK;
}

/**
 * @brief Convert a vector of arbitrary length in TCL_Obj to double[]
 *
 * @param s tcl object containg a list
 * @param val the double vector to return into
 * @param interp the tcl interpreter
 *
 * @return TCL_OK/TCL_ERROR
 */
int tcl_get_vector_obj(Tcl_Obj *s, double *val, Tcl_Interp *interp)
{

    int num = 0;
    Tcl_Obj **pos = NULL;
    if (Tcl_ListObjGetElements(interp, s, &num, &pos) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "vector expected", TCL_STATIC );
        return TCL_ERROR;
    }

    for (int i = 0; i < num; i++)
        if (Tcl_GetDoubleFromObj(interp, pos[i], val + i) != TCL_OK)
            return TCL_ERROR;

    return TCL_OK;
}

/**
 * @brief Convert an array of arbitrary length in TCL_Obj to double[][]
 *
 * @param s tcl object containg a list
 * @param val the double vector to return into
 * @param interp the tcl interpreter
 *
 * @return TCL_OK/TCL_ERROR
 */
int tcl_get_array_obj(Tcl_Obj *s, Tcl_Interp *interp, double **a)
{

    int num = 0;
    Tcl_Obj **pos = NULL;
    if (Tcl_ListObjGetElements(interp, s, &num, &pos) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Badly formed array or dimension mismatch", TCL_STATIC );
        return TCL_ERROR;
    }

    // Loop over list of list, return double into a
    int num2 = 0;
    Tcl_Obj **pos2 = NULL;
    for (int i = 0; i < num; i++) {

        if (Tcl_ListObjGetElements(interp, pos[i], &num2, &pos2) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Badly formed array or dimension mismatch", TCL_STATIC );
            return TCL_ERROR;
        }

        double *b = a[i];
        for (int j = 0; j < num2; j++)
            if (Tcl_GetDoubleFromObj(interp, pos2[j], b + j) != TCL_OK) {
                Tcl_SetResult(interp, (char *) "Non-Numeric in array", TCL_STATIC );
                return TCL_ERROR;
            }
    }

    return TCL_OK;
}


/**
 * @brief parse a tcl object containg a 4x4 matrix to a double [4][4] array
 *
 * @param fctn a string describing what function this is being called from, used in error
 * @param interp tcl interpreter
 * @param s tcl object containg the string representation of the matrix
 *
 * @return TCL_OK/TCL_ERROR
 */

int tcl_get_matrix(const char *fctn, Tcl_Interp *interp, Tcl_Obj *s, double mat[4][4])
{
    int num_rows = 0;
    Tcl_Obj **data_rows = NULL;
    if (Tcl_ListObjGetElements(interp, s, &num_rows, &data_rows) != TCL_OK) {
        char tmpstring[256];
        sprintf(tmpstring, "%s: badly formed matrix", fctn);
        Tcl_SetResult(interp, tmpstring, TCL_VOLATILE );
        return TCL_ERROR;
    }
    if (num_rows != 4) {
        char tmpstring[256];
        sprintf(tmpstring, "%s: need a 4x4 matrix", fctn);
        Tcl_SetResult(interp, tmpstring, TCL_VOLATILE );
        return TCL_ERROR;
    }
    int num_row[4] = { 0 };
    Tcl_Obj **data_row[4];
    if (Tcl_ListObjGetElements(interp, data_rows[0], num_row + 0, data_row + 0) != TCL_OK || num_row[0] != 4
        || Tcl_ListObjGetElements(interp, data_rows[1], num_row + 1, data_row + 1) != TCL_OK || num_row[1] != 4
        || Tcl_ListObjGetElements(interp, data_rows[2], num_row + 2, data_row + 2) != TCL_OK || num_row[2] != 4
        || Tcl_ListObjGetElements(interp, data_rows[3], num_row + 3, data_row + 3) != TCL_OK || num_row[3] != 4) {
        Tcl_AppendResult(interp, fctn, ": poorly formed matrix", NULL);
        return TCL_ERROR;
    }

    // now get the numbers
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double tmp = 0.0;
            if (Tcl_GetDoubleFromObj(interp, data_row[i][j], &tmp) != TCL_OK) {
                char tmpstring[256];
                sprintf(tmpstring, "%s: non-numeric in matrix", fctn);
                Tcl_SetResult(interp, tmpstring, TCL_VOLATILE );
                return TCL_ERROR;
            } else {
                mat[i][j] = (double) tmp;  // Matrix4 is transpose of Tcl's matrix
            }
        }
    }
    return TCL_OK;
}


/**
 * @brief Load commands into the interpretor
 *
 * @param interp tcl interpreter
 * @return TCL_OK/TCL_ERROR
 */

int matvec_init(Tcl_Interp *interp)
{

    /** Each command is checked against commands currently loaded into the interpretor.
     * By default, existing commands are NOT overwritten to preserve VMDs commands
     * when necessary
     */

    Tcl_CmdInfo cmdinfo;

    // Next power of two
    if (Tcl_GetCommandInfo(interp, (char *) "nextpow2", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "nextpow2", tcl_nextpow2,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Transaxis command
    if (Tcl_GetCommandInfo(interp, (char *) "transaxis", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "transaxis", tcl_transaxis,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Transabout command
    if (Tcl_GetCommandInfo(interp, (char *) "transabout", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "transabout", tcl_transabout,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecsub command
    if (Tcl_GetCommandInfo(interp, (char *) "vecsub", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecsub", tcl_vecsub,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecadd command
    if (Tcl_GetCommandInfo(interp, (char *) "vecadd", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecadd", tcl_vecadd,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecmul command
    if (Tcl_GetCommandInfo(interp, (char *) "vecmul", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecmul", tcl_vecmul,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecsum command
    if (Tcl_GetCommandInfo(interp, (char *) "vecsum", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecsum", tcl_vecsum,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecdot command
    if (Tcl_GetCommandInfo(interp, (char *) "vecdot", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecdot", tcl_vecdot,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Veccross command
    if (Tcl_GetCommandInfo(interp, (char *) "veccross", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "veccross", tcl_veccross,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecinvert command
    if (Tcl_GetCommandInfo(interp, (char *) "vecinvert", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecinvert", tcl_vecinvert,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecnorm command
    if (Tcl_GetCommandInfo(interp, (char *) "vecnorm", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecnorm", tcl_vecnorm,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecscale command
    if (Tcl_GetCommandInfo(interp, (char *) "vecscale", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecscale", tcl_vecscale,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Veclength command
    if (Tcl_GetCommandInfo(interp, (char *) "veclength", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "veclength", tcl_veclength,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Veclength2 command
    if (Tcl_GetCommandInfo(interp, (char *) "veclength2", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "veclength2", tcl_veclength2,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecdist command
    if (Tcl_GetCommandInfo(interp, (char *) "vecdist", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "vecdist", tcl_vecdist,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Matmult command
    if (Tcl_GetCommandInfo(interp, (char *) "matmult", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "matmult", tcl_matmult,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Matvec command
    if (Tcl_GetCommandInfo(interp, (char *) "matvec", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "matvec", tcl_matvec,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Transvec command
    if (Tcl_GetCommandInfo(interp, (char *) "transvec", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "transvec", tcl_transvec,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Transvecinv command
    if (Tcl_GetCommandInfo(interp, (char *) "transvecinv", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "transvecinv", tcl_transvecinv,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Transtranspose command
    if (Tcl_GetCommandInfo(interp, (char *) "transtranspose", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, (char *) "transtranspose", tcl_transtranspose,
                             (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    // Vecmean
    if (Tcl_GetCommandInfo(interp, (char *) "vecmean", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, "vecmean", obj_vecmean,
                             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

    ///Vecstdev
    if (Tcl_GetCommandInfo(interp, (char *) "vecstdev", &cmdinfo) == 0)
        Tcl_CreateObjCommand(interp, "vecstddev", obj_vecstddev,
                             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

    return TCL_OK;
}


/**
 * @brief Delete commands from the interpretor
 *
 * @param interp tcl interpreter
 * @return TCL_OK/TCL_ERROR
 */
int matvec_destroy(Tcl_Interp *interp)
{
    Tcl_DeleteCommand(interp, (char *) "nextpow2");
    Tcl_DeleteCommand(interp, (char *) "transaxis");
    Tcl_DeleteCommand(interp, (char *) "transabout");
    Tcl_DeleteCommand(interp, (char *) "vecsub");
    Tcl_DeleteCommand(interp, (char *) "vecadd");
    Tcl_DeleteCommand(interp, (char *) "vecmul");
    Tcl_DeleteCommand(interp, (char *) "vecsum");
    Tcl_DeleteCommand(interp, (char *) "vecdot");
    Tcl_DeleteCommand(interp, (char *) "veccross");
    Tcl_DeleteCommand(interp, (char *) "nextpow2");
    Tcl_DeleteCommand(interp, (char *) "vecinvert");
    Tcl_DeleteCommand(interp, (char *) "vecnorm");
    Tcl_DeleteCommand(interp, (char *) "vecscale");
    Tcl_DeleteCommand(interp, (char *) "veclength");
    Tcl_DeleteCommand(interp, (char *) "veclength2");
    Tcl_DeleteCommand(interp, (char *) "vecdist");
    Tcl_DeleteCommand(interp, (char *) "matmult");
    Tcl_DeleteCommand(interp, (char *) "matvec");
    Tcl_DeleteCommand(interp, (char *) "transvec");
    Tcl_DeleteCommand(interp, (char *) "transvecinv");
    Tcl_DeleteCommand(interp, (char *) "transtranspose");
    Tcl_DeleteCommand(interp, (char *) "vecmean");
    Tcl_DeleteCommand(interp, (char *) "vecstdev");

    return TCL_OK;
}
