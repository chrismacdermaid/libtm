/**
 * @file   main.h
 * @author macdercm <macdercm@mtdoom>
 * @date   Fri Jan 24 11:12:59 2014
 *
 * @brief  main header
 *
 *
 */

#ifndef TM_MAIN_H
#define TM_MAIN_H

#include "tcl.h"

namespace TCLMATH_NS {

//! Defines a macro for an unused variable. At compilation time, if the user uses
//! the variable, an error is thrown, and compilation fails.
//! Example: int reinit_tcl_atomsel_obj(int UNUSED(argc), Tcl_//! Obj * const UNUSED(objv[]))
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

}

#endif
