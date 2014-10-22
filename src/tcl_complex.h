/**
 * @file   tcl_complex.h
 * @author Chris MacDermaid <macdercm@mtdoomVM.localdomain>
 * @date   Sun Feb  2 12:55:19 2014
 *
 * @brief  Header for tcl_complex routines
 *
 *
 */


#ifndef TM_TCLCOMPLEX_H
#define TM_TCLCOMPLEX_H

#include <vector>
#include <complex>

int complex_init(Tcl_Interp *interp);
int complex_destroy(Tcl_Interp *interp);

// +---------+
// | HELPERS |
// +---------+

int tcl_get_complex_vector(ClientData /*clientdata*/, Tcl_Interp *,
                           Tcl_Obj * const, std::vector<std::complex<double> > &);

int tcl_put_complex_list(ClientData /*clientdata*/, Tcl_Interp *,
                         Tcl_Obj *, std::vector<std::complex<double> > &);

int tcl_put_real_list(ClientData /*clientdata*/, Tcl_Interp *,
                         Tcl_Obj *, std::vector<std::complex<double> > &);
#endif
