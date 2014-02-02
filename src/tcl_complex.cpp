/**
 * @file   tcl_complex.cpp
 * @author macdercm <macdercm@mtdoom>
 * @date   Fri Jan 31 14:14:42 2014
 *
 * @brief  Complex math routines for TCL
 *
 *
 */


#include <iostream>
#include <complex>
#include <vector>

#include "main.h"
#include "tcl.h"

#include "tcl_math.h"
#include "tcl_complex.h"
#include "math_extra.h"

using namespace TCLMATH_NS;
using namespace MathExtra;

/**
 * @brief Calculate the complex conjugate of the provided list
 *
 *
 * @return list of complex conjugates
 */
static int tcl_complex_conj(ClientData /*clientdata*/, Tcl_Interp *interp,
                            int objc, Tcl_Obj * const objv[]) {

    std::vector<std::complex<double> > v;

    if ( objc != 2 ) {
        Tcl_WrongNumArgs(interp, objc, objv, "complex_conj expects a list of complex numbers");
        return TCL_ERROR;
    }
    Tcl_Obj* list = objv[1];

    int count = 0;
    if ( Tcl_ListObjLength(interp, list, &count) != TCL_OK ) {
        Tcl_AppendResult(interp, "Couldn't get complex list length");
        return TCL_ERROR;
    }

    if (tcl_get_complex_vector(NULL, interp, list, v) != TCL_OK ) {
        Tcl_AppendResult(interp, "Couldn't retrieve complex list, malformed list?");
        return TCL_ERROR;
    }

    // Take the conjugate of each of the values in the vector
    for(std::vector<std::complex<double> >::iterator
                it = v.begin(); it != v.end(); ++it) {
        *it = std::conj(*it);
    }

    // Construct a list of complex objects
    list = Tcl_NewListObj(0, NULL);
    if (tcl_put_complex_list(NULL, interp, list, v) != TCL_OK ) {
        Tcl_AppendResult(interp, "Couldn't construct complex list");
        return TCL_ERROR;
    }

    // Return the list to the user
    Tcl_SetObjResult(interp, list);

    return TCL_OK;
}

static int tcl_veclength2_complex(ClientData /*clientdata*/, Tcl_Interp *interp,
                            int objc, Tcl_Obj * const objv[]) {

int tcl_get_complex_vector(ClientData /*clientdata*/, Tcl_Interp *interp,
                           Tcl_Obj * const list, std::vector<std::complex<double> > &v) {

    int count = 0;
    if ( Tcl_ListObjLength(interp, list, &count) != TCL_OK ) {
        Tcl_AppendResult(interp, "Couldn't get complex list length");
        return TCL_ERROR;
    }

    /// Get complex values from the list
    for ( size_t i=0; i<count; ++i ) {

        Tcl_Obj *x;
        if ( Tcl_ListObjIndex(interp, list, i, &x) != TCL_OK ) {
            Tcl_AppendResult(interp, "Invalid index");
            return TCL_ERROR;
        }

        Tcl_Obj *xRe, *xIm;
        if ( Tcl_ListObjIndex(interp, x, 0, &xRe) != TCL_OK ) {
            Tcl_AppendResult(interp, "Can't get real component");
            return TCL_ERROR;
        }

        if ( Tcl_ListObjIndex(interp, x, 1, &xIm) != TCL_OK ) {
            Tcl_AppendResult(interp, "Can't get imag component");
            return TCL_ERROR;
        }

        double Re = 0.0, Im = 0.0;
        if ( Tcl_GetDoubleFromObj(interp, xRe, &Re) != TCL_OK ) {
            Tcl_AppendResult(interp, "Real part not a double");
            return TCL_ERROR;
        }

        if ( Tcl_GetDoubleFromObj(interp, xIm, &Im) != TCL_OK ) {
            Tcl_AppendResult(interp, "Imag part not a double");
            return TCL_ERROR;
        }

        // store vector of complex values
        v.push_back(std::complex<double> (Re, Im));
    }

    return TCL_OK;
}



int tcl_put_complex_list(ClientData /*clientdata*/, Tcl_Interp *interp,
                         Tcl_Obj *list, std::vector<std::complex<double> > &v) {

    Tcl_Obj* Re;
    Tcl_Obj* Im;
    Tcl_Obj* C;

    for(std::vector<std::complex<double> >::iterator
                it = v.begin(); it != v.end(); ++it) {

        Re = Tcl_NewDoubleObj(std::real(*it));
        Im = Tcl_NewDoubleObj(std::imag(*it));
        C = Tcl_NewListObj(0, NULL);

        if ( Tcl_ListObjAppendElement(interp, C, Re) != TCL_OK ) {
            Tcl_AppendResult(interp, "Couldn't append real part");
            return TCL_ERROR;
        }
        if ( Tcl_ListObjAppendElement(interp, C, Im) != TCL_OK ) {
            Tcl_AppendResult(interp, "Couldn't append imaginary part");
            return TCL_ERROR;
        }
        if ( Tcl_ListObjAppendElement(interp, list, C) != TCL_OK ) {
            Tcl_AppendResult(interp, "Couldn't append complex result");
            return TCL_ERROR;
        }
    }

    Tcl_SetObjResult(interp, list);

    return TCL_OK;
}

int complex_init(Tcl_Interp *interp) {

    // complex_conj
    Tcl_CreateObjCommand(interp, (char *) "complex_conj", tcl_complex_conj,
                         (ClientData)NULL, (Tcl_CmdDeleteProc *) NULL);

    return TCL_OK;
}
