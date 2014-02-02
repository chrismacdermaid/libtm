// -*-c++-*-

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <tcl.h>
#include "main.h"
#include "tcl_math.h"
#include "tcl_complex.h"

using namespace TCLMATH_NS;

/* Tcl plugin wrapper. it parses the arguments, calls the
   real code and then passes the result to the interpreter. */
int tcl_tm(ClientData UNUSED(clientdata), Tcl_Interp *interp,
           int UNUSED(objc), Tcl_Obj *const UNUSED(objv[]))
{
    // Load the matrix/vector commands
    matvec_init(interp);

    // Load complex commands
    complex_init(interp);

    return TCL_OK;
}

/**
 * Register the plugin with the TCL interpreter
 *
 */

extern "C" {

    /* register the plugin with the tcl interpreters */
#if defined(TMTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

    BOOL APIENTRY DllMain( HANDLE hModule,
                           DWORD  ul_reason_for_call,
                           LPVOID lpReserved )
    {
        return TRUE;
    }

    EXTERN int Tm_Init(Tcl_Interp *interp)

#else

            int Tm_Init(Tcl_Interp *interp)

#endif
    {

#if defined(USE_TCL_STUBS)
        if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
            return TCL_ERROR;
        if (Tcl_PkgRequire(interp, "Tcl", TCL_VERSION, 0) == NULL)
            return TCL_ERROR;
#endif

        if (Tcl_PkgProvide(interp, PACKAGE_NAME, PACKAGE_VERSION) != TCL_OK)
            return TCL_ERROR;

        // Create the initializer command
        Tcl_CreateObjCommand(interp,"tm",tcl_tm,
                             (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

        // Call the initializer command so tm is loaded automatically
        Tcl_Eval(interp,"tm");

        return TCL_OK;
    }

#if defined(TMTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

    BOOL APIENTRY DllMain( HANDLE hModule,
                           DWORD  ul_reason_for_call,
                           LPVOID lpReserved )
    {
        return TRUE;
    }

    EXTERN int Tm_SafeInit(Tcl_Interp *)

#else

            int Tm_SafeInit(Tcl_Interp *)

#endif
    {
        return TCL_OK;
    }

#if defined(TMTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

    BOOL APIENTRY DllMain( HANDLE hModule,
                           DWORD  ul_reason_for_call,
                           LPVOID lpReserved )
    {
        return TRUE;
    }

    EXTERN int Tm_SafeUnload(Tcl_Interp *)

#else

            int Tm_SafeUnload(Tcl_Interp *)

#endif
    {
        return TCL_OK;
    }

#if defined(TMTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

    BOOL APIENTRY DllMain( HANDLE hModule,
                           DWORD  ul_reason_for_call,
                           LPVOID lpReserved )
    {
        return TRUE;
    }

    EXTERN int Tm_Unload(Tcl_Interp *)

#else

            int Tm_Unload(Tcl_Interp *)

#endif
    {
        return TCL_OK;
    }
}
