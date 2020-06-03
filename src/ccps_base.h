/* $Id: ccps_base.h,v 1.19 2009-05-06 15:51:52 Alex_Pletzer Exp $ */

/* for size_t */
#include <stdlib.h> 

#include "fpreproc/f77name.h"

#ifndef CCPS_BASE_H
#define CCPS_BASE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Get number of components
 * @param number (out)
 */
#define psGetNumComponents F77NAME(ccps_get_num_components)
  void psGetNumComponents(int* numb);

/**
 * Get component index number 
 * @param name capitilized component name (e.g. NBI, FUS, ...)
 * @param number index (out)
 */
#define psGetComponentNumber F77NAME(ccps_get_component_number)
  void psGetComponentNumber(const char* name, int* number, size_t);

/**
 * Get component string descriptor 
 * @param name capitilized component name (e.g. NBI, FUS, ...)
 * @param descr   (out)
 */
#define psGetComponentDescriptor F77NAME(ccps_get_component_descriptor)
  void psGetComponentDescriptor(const char* name, char* descr, size_t, size_t);

/**
 * Global (static) method to set debug level
 * @param debug debug level (0, 1, 2, or 3)
 * @param ierr error code (0=OK)
 */
#define psSetDebugLevel F77NAME(ccps_set_debug_level)
void psSetDebugLevel(const int* debug, int* ierr);

/**
 * Global (static) method to get version id
 * @param version_id version id
 * @param size length of Fortran character version_id 
 * @param ierr error code (0=OK)
 */
#define psGetVersionId F77NAME(ccps_get_version_id)
void psGetVersionId(char* version_id, int* size, int* ierr, size_t);

/** 
 * Global (static) method to convert integer species data to MKSA units
 * @param zAtom atomic number (or number of protons), -1 for electrons
 * @param zCharge  1 <= ionic charge <= zAtom
 * @param amu atomic weight, nearest integer (ignored unless 1 <= zAtom <=2)
 * @param qAtom charge of fully stripped ion (out)
 * @param qCharge actual charge (out)
 * @param mass mass in kg (out)
 * @param ierr error code (0=OK)
 */
#define psSpeciesConvert F77NAME(ccps_species_convert)
void psSpeciesConvert(const int* zAtom, const int* zCharge, const int* amu,
		      double* aAtom, double* qCharge, double* mass, int* ierr);
/**
 * Generate labels for each allocated species list (thermal, beam, RF minorities, 
 * fusion)
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psLabelSpecies F77NAME(ccps_label_species)
void psLabelSpecies(int* iobj, int* ierr);
  
/**
 * Create reduced list of species out of the union of thermal, beam, RF minority,
 * and fusion species lists
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psMergeSpeciesLists F77NAME(ccps_merge_species_lists)
  void psMergeSpeciesLists(int* iobj, int* ierr);

/**
 * Construct neutral species corresponding to each atomic species
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psNeutralSpecies F77NAME(ccps_neutral_species)
void psNeutralSpecies(int* iobj, int* ierr);

/**
 * Update hash table; lock (or unlock) initialization sections of state
 * @param iobj reference to opaque plasma state object
 * @param mdescr_iflag if !=0 => lock machine description section
 * @param sconfig_iflag if !=0 => lock shot description section
 * @param simInit_iflag if !=0 => lock simulation description section
 * @param ierr error code (0=OK)
 */
#define psUpdateHashCodes F77NAME(ccps_update_hashcodes)
void psUpdateHashCodes(int* iobj, const int* mdescr_iflag, const int* sconfig_iflag, const int* simInit_iflag, int* ierr);

/**
 * Restore state's hash codes from a previously saved copy
 * @param iobj reference to opaque plasma state object 
 * @param ierr error code (0=OK)
 */
#define psRestoreHashCodes F77NAME(ccps_restore_hash_codes)
  void psRestoreHashCodes(int* iobj, int* ierr);

/**
 * Save a copy of the state's current hash codes in the state object's
 * @param iobj reference to opaque plasma state object 
 * @param ierr error code (0=OK)
 */
#define psSaveHashCodes F77NAME(ccps_save_hash_codes)
  void psSaveHashCodes(int* iobj, int* ierr);

/**
 * Set selected profiles to zero
 * @param iobj reference to opaque plasma state object 
 * @param cclist component list, set elements to 1 to activate
 * @param save_hash set to 1 to save hash codes
 * @param ierr error code (0=OK)
 */
#define psClearProfs F77NAME(ccps_clear_profs)
  void psClearProfs(int* iobj, const int cclist[], const int* save_hash, 
		    int* ierr);

/**
 * Write file containing updates
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param hashflag .NE.0 to update I/O hash table
 * @param ierr error code (0=OK)
 */
#define psWriteUpdateFile F77NAME(ccps_write_update_file)
void psWriteUpdateFile(int* iobj, const char* filepath, const int* hashflag, int* ierr, size_t);

/**
 * Invert psi -> rho map
 * @param iobj reference to opaque plasma state object
 * @param tol tolerance for Psi->rho map
 * @param npsi number of psi values in & rhopsi values out
 * @param psi psi values in
 * @param rhopsi rho values out
 * @param n_outside number of values not in range
 * @param ierr error code (0=OK)
 */
#define psPsiMap F77NAME(ccps_psimap)
void psPsiMap(int* iobj, const double* tol, const int* npsi, const double* psi, double* rhopsi, int* n_outside, int* ierr);

/**
 * Allocate a portion of a state object (those elements for which the size is nonzero)
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psAlloc F77NAME(ccps_alloc)
void psAlloc(int* iobj, int* ierr);

/**
 * Copy state objects
 * @param iobj_from reference to opaque plasma state object
 * @param iobj_to reference to opaque plasma state object
 * @param cclist list of elements to be copied (0=no copied, 1=copied)
 * @param inodims 1 if output should contain only copied data
 * @param iwant_1deq 1 to copy (nrho_eq) dim, even if cclist(ps_cnum_EQ)==0
 * @param ierr error code (0=OK)
 */
#define psCopyPlasmaState F77NAME(ccps_copy_plasma_state)
  void psCopyPlasmaState(int* iobj_from, int* iobj_to, const int cclist[],
			 const int* inodims, const int* iwant_1deq, int* ierr);

/**
 * Free a state object
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psFree F77NAME(ccps_free)
void psFree(int* iobj, int* ierr);

/**
 * Read state object contents from file
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psGetPlasmaState F77NAME(ccps_get_plasma_state)
void psGetPlasmaState(int* iobj, const char* filepath, int* ierr, size_t);

/**
 * Retrieve label information for specified profile
 * @param iobj reference to opaque plasma state object
 * @param id profile ID
 * @param label item label (up to 120 characters)
 * @param units item physical units (up to 32 characters)
 * @param component component owning the item (up to 32 characters)
 * @param ierr error code (0=OK)
 */
#define psGetProfLabel F77NAME(ccps_get_prof_label)
void psGetProfLabel(int* iobj, const int* id, char* label, char* units, char* component, int* ierr, size_t, size_t, size_t);

/**
 * Construct a state object
 * @param iobj reference to opaque plasma state object
 * @param bytename null-terminated tag name byte string (up to 32 characters)
 * @param ierr error code (0=OK)
 */
#define psInit F77NAME(ccps_init)
void psInit(int* iobj, const char* bytename, int* ierr, size_t);

/**
 * Interpolate a 1d state element profile at a single point
 * @param iobj reference to opaque plasma state object
 * @param id profile ID
 * @param ideriv  0 for value, 1 for d/dx, 2 for d2/dx2
 * @param xval value at which to interpolate
 * @param result interpolated value returned
 * @param ierr error code (0=OK)
 */
#define psIntrp1dSingle F77NAME(ccps_intrp_1d_single)
void psIntrp1dSingle(int* iobj, const int* id, const int* ideriv, const double* xval, double* result, int* ierr);

/**
 * 
 * @param iobj reference to opaque plasma state object
 * @param id profile ID
 * @param ideriv 0 for value, 1 for d/dx, 2 for d2/dx2
 * @param nvals size of target & results vectors
 * @param xvals values at which to interpolate
 * @param results interpolated values
 * @param ierr error code (0=OK)
 */
#define psIntrp1dVector F77NAME(ccps_intrp_1d_vector)
void psIntrp1dVector(int* iobj, const int* id, const int* ideriv, const int* nvals, const double* xvals, double* results, int* ierr);

/**
 * Interpolate a 2d f(x1,x2) state element profile at a single point
 * @param iobj reference to opaque plasma state object
 * @param id profile ID
 * @param ideriv1 0 for value, 1 for d/dx1, 2 for d2/d[x1]2
 * @param ideriv2 0 for value, 1 for d/dx2, 2 for d2/d[x2]2
 * @param x1val value at which to interpolate
 * @param x2val value at which to interpolate
 * @param result interpolation result
 * @param ierr error code (0=OK)
 */
#define psIntrp2dSingle F77NAME(ccps_intrp_2d_single)
void psIntrp2dSingle(int* iobj, const int* id, const int* ideriv1, const int* ideriv2, const double* x1val, const double* x2val, double* result, int* ierr);

/**
 * Interpolate a 2d state element profile at a vector of target values
 * @param iobj reference to opaque plasma state object
 * @param id profile ID
 * @param ideriv1 0 for value, 1 for d/dx1, 2 for d2/d[x1]2
 * @param ideriv2 0 for value, 1 for d/dx2, 2 for d2/d[x2]2
 * @param nvals size of target & results vectors
 * @param x1vals values at which to interpolate
 * @param x2vals values at which to interpolate
 * @param results interpolation results
 * @param ierr error code (0=OK)
 */
#define psIntrp2dVector F77NAME(ccps_intrp_2d_vector)
void psIntrp2dVector(int* iobj, const int* id, const int* ideriv1, const int* ideriv2, const int* nvals, const double* x1vals, const double* x2vals, double* results, int* ierr);

/**
 * Read file containing machine description
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param action null terminated string, if 'NEW' or 'INIT' then erase plasma state
 * @param g_filepath null terminated file path string to override box size and limiter contour (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psMDescrRead F77NAME(ccps_mdescr_read)
  void psMDescrRead(int* iobj, const char* filepath, const char* action, const char* g_filepath, int* ierr, size_t, size_t, size_t);

/**
 * Write commented machine description namelist from contents of state
 * @param iobj reference to opaque plasma state object
 * @param filename null terminated file name (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psMDescrWrite F77NAME(ccps_mdescr_write)
void psMDescrWrite(int* iobj, const char* filename, const int* ierr, size_t);

/**
 * Write geqdsk file based on contents of plasma state
 * @param iobj reference to opaque plasma state object
 * @param filename null terminated file name (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psGeqdskWrite F77NAME(ccps_wr_geqdsk)
void psGeqdskWrite(int* iobj, const char* filename, const int* ierr, size_t);

/**
 * Read file containing updates-- i.e. state element value changes
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param finalflag  ! .NE.0 to update state internal tables
 * @param ierr error code (0=OK)
 */
#define psReadUpdateFile F77NAME(ccps_read_update_file)
void psReadUpdateFile(int* iobj, const char* filepath, const int* finalflag, int* ierr, size_t);

/**
 * Compute and return flux surface averaged or integrated metrics over a specified rho grid
 * @param iobj reference to opaque plasma state object
 * @param metric_name null terminated string specifying the desired metric... (case insensitive)
 * @param ngrid size of RHO grid provided
 * @param rho_grid the RHO grid (strict ascending order in range 0 to 1)
 * @param nresult size of result expected (=ngrid or =ngrid-1)
 * @param result results returned
 * @param ierr error code (0=OK)
 */
#define psRhoMetric F77NAME(ccps_rho_metric)
void psRhoMetric(int* iobj, const char* metric_name, const int* ngrid, const int* rho_grid, const int* nresult, double* result, int* ierr, size_t);

/**
 * Rezone a 1d profile f(rho), rho=sqrt(tor.flux) flux label
 * @param iobj reference to opaque plasma state object
 * @param ID profile to be rezoned. If quantity is specified as weighted then the rezoning will be weighted
 * @param ngrid size of RHO grid provided
 * @param rho_grid the RHO grid (strict ascending order in range 0 to 1)
 * @param result[ngrid-1] results returned ***(ngrid-1) values***
 * @param curflag set .NE.0 for current density rezone
 * @param normflag set .NE.0 for area or volume norm
 * @param smooflag set .NE.0 for smoothing
 * @param ierr error code (0=OK)
 */
#define psRhoRezone F77NAME(ccps_rho_rezone)
void psRhoRezone(int* iobj, const int* id, const int* ngrid, const double* rho_grid, double* result, const int* curflag, const int* normflag, const int* smooflag, int* ierr);

/**
 * Read file containing shot configuration data
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psSConfigRead F77NAME(ccps_sconfig_read)
void psSConfigRead(int* iobj, const char* filepath, int* ierr, size_t);

/**
 * Write file containing shot configuration data
 * @param iobj reference to opaque plasma state object
 * @param filepath null terminated file path string (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psSConfigWrite F77NAME(ccps_sconfig_write)
void psSConfigWrite(int* iobj, const char* filepath, int* ierr, size_t);

/**
 * Update hash table and interpolation data associated with state
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psStateMemoryUpdate F77NAME(ccps_state_memory_update)
void psStateMemoryUpdate(int* iobj, int* ierr);

/**
 * Utility routine for setting beam energy fractions
 * @param iobj reference to opaque plasma state object
 * @param ichk if 1 set error code is table has not been initialized
 * @param ierr error code (0=OK)
 */
#define psChkNbfrac F77NAME(ccps_chk_nbfrac)
  void psChkNbfrac(int* iobj, const int* ichk, int* ierr);

/**
 * Write state object contents to file
 * @param iobj reference to opaque plasma state object
 * @param filepath  null terminated file path string (max 256 chars)
 * @param ierr error code (0=OK)
 */
#define psStorePlasmaState F77NAME(ccps_store_plasma_state)
void psStorePlasmaState(int* iobj, const char* filepath, int* ierr, size_t);

/**
 * Update MHD equilibrium stored in state object, from G-eqdsk file
 * @param iobj reference to opaque plasma state object
 * @param g_filepath null terminated G-eqdsk file (max 256 chars)
 * @param bdy_crat curvature ratio limit for boundary
 * @param kcur_option method to obtain q and I, values must be in 
                      {-1 auto-choose, 0=take efit values, 1=compute from psirz}
 * @param rho_curbrk value from which current profile is interpolated to edge (eg 0.9)
 * @param npprime number of values for pressure profile smoothing (0 for no smoothing)
 * @param rho_psmoo rho values around which profile should be smoothed
 * @param delrho_psmoo widths of smoothing kernels
 * @param ierr error code (0=OK)
 */
#define psUpdateEquilibrium F77NAME(ccps_update_equilibrium)
  void psUpdateEquilibrium(int* iobj, const char* g_filepath, 
			   const double* bdy_crat, const int* kcur_option, 
			   const double* rho_curbrk, const int* npprime, 
			   const double* rho_psmoo, const double* delrho_psmoo, 
			   int* ierr, size_t);

/**
 * Verify that profile of enclosed toroidal flux is consistent with rho coordinate
 * @param iobj reference to opaque plasma state object
 * @param reltol demanded relative accuracy (~0.02)
 * @param nrhowk size of working grid for <1/R^2> integration (can be zero)
 * @param rhowk working grid (see above)
 * @param update_phit set to zero if toroidal flux should not be updated
 * @param update_psi set to zero if poloidal flux should not be updated
 * @param update_q set to zero if safety factor should not be updated (update_psi and update_q cannot be both != 0)
 * @param check_phit set to one if toroidal flux should be checked
 * @param check_q set to one is safety factor should be checked
 * @param relerr_max max relative error in phi or q tests
 * @param ierr error code (0=OK)
 */
#define psMhdeqVerify F77NAME(ccps_mhdeq_verify)
void psMhdeqVerify(int* iobj, const double* reltol, 
         const int* nrhowk, const double* rhowk, 
         const int* update_phit, const int* update_psi, const int* update_q, 
	 const int* check_phit, const int* check_q, 
         double* relerr_max, int* ierr);

/**
 * Derive flux surface averages and profiles from equilibrium geometry and fields
 * @param iobj reference to opaque plasma state object
 * @param action specify action to be taken (eg "Everything")
 * @param nrhowk size of working grid (can be zero)
 * @param rhowk working grid
 * @param ierr error code (0=OK)
 */
#define psMhdeqDerive F77NAME(ccps_mhdeq_derive)
  void psMhdeqDerive(int* iobj, const char* action, 
          const int* nrhowk, const double* rhowk, int* ierr, size_t);

/**
 * Summarize difference between two states, based on hash code comparison.
 * @param iobj1 reference to first opaque plasma state object
 * @param lbl1 label for 1st state
 * @param iobj2 reference to second opaque plasma state object
 * @param lbl2 label for 2nd state
 * @param sstr show all diffs for section, profile, dimension containg this string
 * @param icomp number of hash comparisons done
 * @param idiff number of DIFFERENCES found
 */
#define psHashDiff F77NAME(ccps_hash_diff)
  void psHashDiff(int* iobj1, const char* lbl1, int* iobj2, const char* lbl2,
		  const char* sstr, int* icomp, int* idiff, 
		  size_t, size_t, size_t);

/**
 * Store plasma state data in xplasma and enable interpolation services
 * scalar data is written if section hashes indicate change, or if tests
 * indicate true. Grid data items are written if hashes indicate change, 
 * or if grid IDs cannot be found in xplasma
 * @param iobj reference to opaque plasma state object
 * @param ierr error code (0=OK)
 */
#define psXload F77NAME(ccps_xload)
  void psXload(int* iobj, int* ierr);

/**
 * Load xplasma from state structure; only load items that changed
 * or are not found in the xplasma object (small items may be loaded
 * regardless).  Check for rebuild of equilibrium in xplasma.
 * @param iobj reference to opaque plasma state object
 * @param iout output Fortran unit number
 * @param icheck set to 0 to skip equilibrium check
 * @param ierr error code (0=OK)
 */
#define psXloadCheck F77NAME(ccps_xload_check)
  void psXloadCheck(int* iobj, const int* iout, const int* icheck, int* ierr);  

#ifdef __cplusplus
}
#endif

#endif /* CCPS_BASE_H */
