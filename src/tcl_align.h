/**
 * @file   align.h
 * @author  <chris@cabal.nfs>
 * @date   Mon May  7 16:10:08 2012
 *
 * @brief  Compute the optimal alignment matrix for two coordinate systems
 *
 * Uses the Kabsch algorithm and the ACML singular value decomposition
 * routines to determine the optimal rotation matrix that superimposes
 * coordinates A with coordinates B.
 *
 * Acta Cryst. (1976). A32, 922-923    [ doi:10.1107/S0567739476001873 ]
 *
 * A solution for the best rotation to relate two sets of vectors. W. Kabsch
 * Abstract: A simple procedure is derived which determines a best rotation
 * of a given vector set into a second vector set by minimizing the weighted
 * sum of squared deviations. The method is generalized for any given metric
 * constraint on the transformation.
 *
 * This is the same approach used by VMD and Pymol.
 *
 * http://en.wikipedia.org/wiki/Kabsch_algorithm
 *
 * Coutsias EA, Seok C, Dill KA (2004). "Using quaternions to calculate RMSD".
 * Comput Chem 25 (15): 1849–1857. doi:10.1002/jcc.20110. PMID 15376254.
 *
 * Kabsch W (1976). "A solution for the best rotation to relate two sets of vectors".
 * Acta Crystallographica 32 (5): 922–923. doi:10.1107/S0567739476001873.
 *
 * For SVD, see the ACML. dgesdd_c_example.c
 *
 * The algorithm from VMD/PYMOL is also provided, and gives identical results,
 * however, the vmd/pymol method is:
 *
 * 0.00058 seconds / 0.00016 seconds = 3.625x faster...
 *
 * So we use it instead...
 */

#ifndef TM_TCLALIGN_H
#define TM_TCLALIGN_H

int align_init(Tcl_Interp *interp);
int align_destroy (Tcl_Interp *interp);

void MatrixFitRMS(int n, double **p, double **q, double *w, double m[4][4]); /**< Alignment algorithm from VMD/PYMOL */

#endif
