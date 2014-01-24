/**
 * @file   tcl_math.h
 * @author macdercm <macdercm@mtdoom>
 * @date   Fri Jan 24 11:18:17 2014
 * 
 * @brief  Vector/Matrix routines for tcl excised from VMD
 * 
 * 
 */

#ifndef TM_TCLMATH_H
#define TM_TCLMATH_H

/* +---------------------------------+  */
/* | Vector and Matrix Commands   | */
/* +---------------------------------+  */

int matvec_init(Tcl_Interp *interp);

/* +---------+   */
/* | HELPERS |   */
/* +---------+   */

int tcl_get_vector(const char *s, double *val, Tcl_Interp *interp); /**< Get vector from TCL String */

int tcl_get_vector_obj(Tcl_Obj *s, double *val, Tcl_Interp *interp); /**< Get vector from TCL Obj */

int tcl_get_matrix(const char *fctn, Tcl_Interp *interp, Tcl_Obj *s, double mat[4][4]);

#endif
