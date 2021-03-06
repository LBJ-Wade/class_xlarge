/** @file perturbations.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include "thermodynamics.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"

#define _scalars_ ((ppt->has_scalars == _TRUE_) && (index_md == ppt->index_md_scalars))
#define _vectors_ ((ppt->has_vectors == _TRUE_) && (index_md == ppt->index_md_vectors))
#define _tensors_ ((ppt->has_tensors == _TRUE_) && (index_md == ppt->index_md_tensors))

#define _set_source_(index) ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index][index_tau * ppt->k_size[index_md] + index_k]

/**
 * flags for various approximation schemes
 * (tca = tight-coupling approximation,
 *  rsa = radiation streaming approximation,
 *  ufa = massless neutrinos / ultra-relativistic relics fluid approximation)
 *
 * CAUTION: must be listed below in chronological order, and cannot be
 * reversible. When integrating equations for a given mode, it is only
 * possible to switch from left to right in the lists below.
 */

//@{

enum tca_flags {tca_on, tca_off};
enum rsa_flags {rsa_off, rsa_on};
enum ufa_flags {ufa_off, ufa_on};
enum ncdmfa_flags {ncdmfa_off, ncdmfa_on};

//@}

/**
 * labels for the way in which each approximation scheme is implemented
 */

//@{

enum tca_method {first_order_MB,first_order_CAMB,first_order_CLASS,second_order_CRS,second_order_CLASS,compromise_CLASS};
enum rsa_method {rsa_null,rsa_MD,rsa_MD_with_reio,rsa_none};
enum ufa_method {ufa_mb,ufa_hu,ufa_CLASS,ufa_none};
enum ncdmfa_method {ncdmfa_mb,ncdmfa_hu,ncdmfa_CLASS,ncdmfa_none};
enum tensor_methods {tm_photons_only,tm_massless_approximation,tm_exact};

//@}

/**
 * List of coded gauges. More gauges can in principle be defined.
 */

//@{

enum possible_gauges {
  newtonian, /**< newtonian (or longitudinal) gauge */
  synchronous /**< synchronous gauge with \f$ \theta_{cdm} = 0 \f$ by convention */
};

//@}


//@{

/**
 * maximum number of k-values for perturbation output
 */
#define _MAX_NUMBER_OF_K_FILES_ 30

//@}


/**
 * Structure containing everything about perturbations that other
 * modules need to know, in particular tabuled values of the source
 * functions \f$ S(k, \tau) \f$ for all requested modes
 * (scalar/vector/tensor), initial conditions, types (temperature,
 * E-polarization, B-polarisation, lensing potential, etc), multipole
 * l and wavenumber k.
 *
 */

struct perturbs
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantitites are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  short has_perturbations; /**< do we need to compute perturbations at all ? */

  short has_cls; /**< do we need any harmonic space spectrum C_l (and hence Bessel functions, transfer functions, ...)? */

  short has_scalars; /**< do we need scalars? */
  short has_vectors; /**< do we need vectors? */
  short has_tensors; /**< do we need tensors? */

  short has_ad;      /**< do we need adiabatic mode? */
  short has_bi;      /**< do we need isocurvature bi mode? */
  short has_cdi;     /**< do we need isocurvature cdi mode? */
  short has_nid;     /**< do we need isocurvature nid mode? */
  short has_niv;     /**< do we need isocurvature niv mode? */

  /* perturbed recombination */
  /* Do we want to consider perturbed temperature and ionization fraction? */
  short has_perturbed_recombination;
  /** Neutrino contribution to tensors */
  enum tensor_methods tensor_method;  /**< way to treat neutrinos in tensor perturbations(neglect, approximate as massless, take exact equations) */

  short evolve_tensor_ur;             /**< will we evolve ur tensor perturbations (either becasue we have ur species, or we have ncdm species with massless approximation) ? */
  short evolve_tensor_ncdm;             /**< will we evolve ncdm tensor perturbations (if we have ncdm species and we use the exact method) ? */

  short has_cl_cmb_temperature;       /**< do we need Cl's for CMB temperature? */
  short has_cl_cmb_polarization;      /**< do we need Cl's for CMB polarization? */
  short has_cl_cmb_lensing_potential; /**< do we need Cl's for CMB lensing potential? */
  short has_cl_lensing_potential;     /**< do we need Cl's for galaxy lensing potential? */
  short has_cl_number_count;          /**< do we need Cl's for density number count? */
  short has_pk_matter;                /**< do we need matter Fourier spectrum? */
  short has_density_transfers;        /**< do we need to output individual matter density transfer functions? */
  short has_velocity_transfers;       /**< do we need to output individual matter velocity transfer functions? */

  short has_nl_corrections_based_on_delta_m;  /**< do we want to compute non-linear corrections with an algorithm relying on delta_m (like halofit)? */

  short has_nc_density;  /**< in dCl, do we want density terms ? */
  short has_nc_rsd1;     /**< in dCl, do we want redshift space distorsion terms (usual Kaiser term) ? */
  short has_nc_rsd2;     /**< in dCl, do we want redshift space distorsion terms (term labelled D2 in CLASSgal's paper ? */
  short has_nc_rsd3;     /**< in dCl, do we want redshift space distorsion terms (term labelled D1 in CLASSgal's paper) ? */
  short has_nc_lens;     /**< in dCl, do we want lensing terms ? */
  short has_nc_gr1;       /**< in dCl, do we want gravity terms (term labelled G2 in CLASSgal's paper)? */
  short has_nc_gr2;       /**< in dCl, do we want gravity terms (term labelled G1 in CLASSgal's paper)? */
  short has_nc_gr3;       /**< in dCl, do we want gravity terms (term labelled G3 in CLASSgal's paper)? */
  short has_nc_gr4;       /**< in dCl, do we want gravity terms (term labelled G4 in CLASSgal's paper)? */
  short has_nc_gr5;       /**< in dCl, do we want gravity terms (term labelled G5 in CLASSgal's paper)? */

  int l_scalar_max; /**< maximum l value for CMB scalars C_ls */
  int l_vector_max; /**< maximum l value for CMB vectors C_ls */
  int l_tensor_max; /**< maximum l value for CMB tensors C_ls */
  int l_lss_max; /**< maximum l value for LSS C_ls (density and lensing potential in  bins) */
  double k_max_for_pk; /**< maximum value of k in 1/Mpc in P(k) (if C_ls also requested, overseeded by value kmax inferred from l_scalar_max if it is bigger) */

  double selection_mean_min;

  int switch_sw;   /**< in temperature calculation, do we want to include the intrinsic temperature + Sachs Wolfe term? */
  int switch_eisw; /**< in temperature calculation, do we want to include the early integrated Sachs Wolfe term? */
  int switch_lisw; /**< in temperature calculation, do we want to include the late integrated Sachs Wolfe term? */
  int switch_dop;  /**< in temperature calculation, do we want to include the Doppler term? */
  int switch_pol;  /**< in temperature calculation, do we want to include the polarisation-related term? */
  double eisw_lisw_split_z; /**< at which redshift do we define the cut between eisw and lisw ?*/

  int k_output_values_num;       /**< Number of perturbation outputs (default=0) */
  double k_output_values[_MAX_NUMBER_OF_K_FILES_];    /**< List of k values where perturbation output is requested. */
  int index_k_output_values[_MAX_NUMBER_OF_K_FILES_]; /**< List of indices corresponding to k-values close to k_output_values */
  FileName root; /**< Same as root in output structure, for writing perturbations.*/


  //@}

  /** @name - useful flags infered from the ones above */

  //@{

  short has_cmb; /**< do we need CMB-related sources (temperature, polarization) ? */
  short has_lss; /**< do we need LSS-related sources (lensing potential, ...) ? */

  //@}

  /** @name - gauge in which to perform the calculation */

  //@{

  enum possible_gauges gauge;

  //@}

  /** @name - indices running on modes (scalar, vector, tensor) */

  //@{

  int index_md_scalars; /**< index value for scalars */
  int index_md_tensors; /**< index value for tensors */
  int index_md_vectors; /**< index value for vectors */

  int md_size; /**< number of modes included in computation */

  //@}

  /** @name - indices running on initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) */

  //@{

  int index_ic_ad; /**< index value for adiabatic */
  int index_ic_cdi; /**< index value for CDM isocurvature */
  int index_ic_bi; /**< index value for baryon isocurvature */
  int index_ic_nid; /**< index value for neutrino density isocurvature */
  int index_ic_niv; /**< index value for neutrino velocity isocurvature */
  int index_ic_ten; /**< index value for unique possibility for tensors */

  int * ic_size;       /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */

  //@}

  /** @name - flags and indices running on types (temperature, polarization, lensing, ...) */

  //@{

  short has_source_t;  /**< do we need source for CMB temperature? */
  short has_source_p;  /**< do we need source for CMB polarisation? */
  short has_source_delta_m;   /**< do we need source for delta of total matter? */
  short has_source_delta_g;    /**< do we need source for delta of gammas? */
  short has_source_delta_b;    /**< do we need source for delta of baryons? */
  short has_source_delta_cdm;  /**< do we need source for delta of cold dark matter? */
  short has_source_delta_dcdm; /**< do we need source for delta of DCDM? */
  short has_source_delta_fld;  /**< do we need source for delta of dark energy? */
  short has_source_delta_dr; /**< do we need source for delta of decay radiation? */
  short has_source_delta_ur; /**< do we need source for delta of ultra-relativistic neutrinos/relics? */
  short has_source_delta_ncdm; /**< do we need source for delta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_theta_m;    /**< do we need source for theta of total matter? */
  short has_source_theta_g;    /**< do we need source for theta of gammas? */
  short has_source_theta_b;    /**< do we need source for theta of baryons? */
  short has_source_theta_cdm;  /**< do we need source for theta of cold dark matter? */
  short has_source_theta_dcdm; /**< do we need source for theta of DCDM? */
  short has_source_theta_fld;  /**< do we need source for theta of dark energy? */
  short has_source_theta_dr; /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_ur; /**< do we need source for theta of ultra-relativistic neutrinos/relics? */
  short has_source_theta_ncdm; /**< do we need source for theta of all non-cold dark matter species (e.g. massive neutrinos)? */
  short has_source_phi;          /**< do we need source for metric fluctuation phi? */
  short has_source_phi_prime;    /**< do we need source for metric fluctuation phi'? */
  short has_source_phi_plus_psi; /**< do we need source for metric fluctuation (phi+psi)? */
  short has_source_psi;          /**< do we need source for metric fluctuation psi? */

  /* remember that the temperature source function includes three
     terms that we call 0,1,2 (since the strategy in class v > 1.7 is
     to avoid the integration by part that would reduce the source to
     a single term) */
  int index_tp_t0; /**< index value for temperature (j=0 term) */
  int index_tp_t1; /**< index value for temperature (j=1 term) */
  int index_tp_t2; /**< index value for temperature (j=2 term) */
  int index_tp_p; /**< index value for polarization */
  int index_tp_delta_m; /**< index value for delta tot */
  int index_tp_delta_g;   /**< index value for delta of gammas */
  int index_tp_delta_b;   /**< index value for delta of baryons */
  int index_tp_delta_cdm; /**< index value for delta of cold dark matter */
  int index_tp_delta_dcdm;/**< index value for delta of DCDM */
  int index_tp_delta_fld;  /**< index value for delta of dark energy */
  int index_tp_delta_dr; /**< index value for delta of decay radiation */
  int index_tp_delta_ur; /**< index value for delta of ultra-relativistic neutrinos/relics */
  int index_tp_delta_ncdm1; /**< index value for delta of first non-cold dark matter species (e.g. massive neutrinos) */
  int index_tp_perturbed_recombination_delta_temp;		/* Gas temperature perturbation */
  int index_tp_perturbed_recombination_delta_chi;			/* Inionization fraction perturbation */

  int index_tp_theta_m;   /**< index value for theta tot */
  int index_tp_theta_g;   /**< index value for theta of gammas */
  int index_tp_theta_b;   /**< index value for theta of baryons */
  int index_tp_theta_cdm; /**< index value for theta of cold dark matter */
  int index_tp_theta_dcdm;/**< index value for theta of DCDM */
  int index_tp_theta_fld;  /**< index value for theta of dark energy */
  int index_tp_theta_ur; /**< index value for theta of ultra-relativistic neutrinos/relics */
  int index_tp_theta_dr; /**< index value for F1 of decay radiation */
  int index_tp_theta_ncdm1; /**< index value for theta of first non-cold dark matter species (e.g. massive neutrinos) */

  int index_tp_phi;          /**< index value for metric fluctuation phi */
  int index_tp_phi_prime;    /**< index value for metric fluctuation phi' */
  int index_tp_phi_plus_psi; /**< index value for metric fluctuation phi+psi */
  int index_tp_psi;          /**< index value for metric fluctuation psi */

  int * tp_size; /**< number of types tp_size[index_md] included in computation for each mode */

  //@}

  /** @name - list of k values for each mode */

  //@{

  int * k_size_cmb;  /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */

  int * k_size_cl;  /**< k_size_cl[index_md] number of k values used
                       for non-CMB Cl calculations, requering a coarse
                       sampling in k-space. */

  int * k_size;     /**< k_size[index_md] = total number of k
                       values, including those needed for P(k) but not
                       for Cl's */

  double ** k;      /**< k[index_md][index_k] = list of values */

  double k_min;     /**< minimum valut (over all modes) */
  double k_max;     /**< maximum valut (over all modes) */

  //@}

  /** @name - list of conformal time values in the source table
      (common to all modes and types) */

  //@{

  int tau_size;          /**< tau_size = number of values */

  double * tau_sampling; /**< tau_sampling[index_tau] = list of tau values */

  double selection_min_of_tau_min; /**< used in presence of selection functions (for matter density, cosmic shear...) */
  double selection_max_of_tau_max; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double selection_delta_tau; /**< used in presence of selection functions (for matter density, cosmic shear...) */

  double * selection_tau_min; /**< value of conformal time below which W(tau) is considered to vanish for each bin */
  double * selection_tau_max; /**< value of conformal time above which W(tau) is considered to vanish for each bin */
  double * selection_tau; /**< value of conformal time at the center of each bin */
  double * selection_function; /** selection function W(tau), normalized to \int W(tau) dtau=1, stored in selection_function[bin*ppt->tau_size+index_tau] */

  //@}

  /** @name - source functions interpolation table */

  //@{

  double *** sources; /**< Pointer towards the source interpolation table
                         sources[index_md]
                         [index_ic * ppt->tp_size[index_md] + index_type]
                         [index_tau * ppt->k_size + index_k] */

  int do_f_nl;
  double *f_nl_kdep;
  double inv_growth_0;

  //@}

  /** @name - technical parameters */

  //@{

  short perturbations_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 */

struct perturb_vector
{
  int index_pt_delta_g;   /**< photon density */
  int index_pt_theta_g;   /**< photon velocity */
  int index_pt_shear_g;   /**< photon shear */
  int index_pt_l3_g;      /**< photon l=3 */
  int l_max_g;            /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_pol0_g;    /**< photon polarization, l=0 */
  int index_pt_pol1_g;    /**< photon polarization, l=1 */
  int index_pt_pol2_g;    /**< photon polarization, l=2 */
  int index_pt_pol3_g;    /**< photon polarization, l=3 */
  int l_max_pol_g;        /**< max momentum in Boltzmann hierarchy (at least 3) */
  int index_pt_delta_b;   /**< baryon density */
  int index_pt_theta_b;   /**< baryon velocity */
  int index_pt_delta_cdm; /**< cdm density */
  int index_pt_theta_cdm; /**< cdm velocity */
  int index_pt_delta_dcdm; /**< dcdm density */
  int index_pt_theta_dcdm; /**< dcdm velocity */
  int index_pt_delta_fld;  /**< dark energy density */
  int index_pt_theta_fld;  /**< dark energy velocity */
  int index_pt_delta_ur; /**< density of ultra-relativistic neutrinos/relics */
  int index_pt_theta_ur; /**< velocity of ultra-relativistic neutrinos/relics */
  int index_pt_shear_ur; /**< shear of ultra-relativistic neutrinos/relics */
  int index_pt_l3_ur;    /**< l=3 of ultra-relativistic neutrinos/relics */
  int l_max_ur;          /**< max momentum in Boltzmann hierarchy (at least 3) */
/* perturbed recombination */
  int index_pt_perturbed_recombination_delta_temp;		/* Gas temperature perturbation */
  int index_pt_perturbed_recombination_delta_chi;			/* Inionization fraction perturbation */

  /** The index to the first Legendre multipole of the DR expansion. Not
      that this is not exactly the usual delta, see Kaplinghat et al.,
      astro-ph/9907388. */
  int index_pt_F0_dr;
  int l_max_dr;          /**< max momentum in Boltzmann hierarchy for dr) */
  int index_pt_psi0_ncdm1;
  int N_ncdm;
  int* l_max_ncdm;
  int* q_size_ncdm;

  int index_pt_eta;       /**< synchronous gauge metric perturbation eta*/
  int index_pt_phi;
  int index_pt_hv_prime;  /**< vector metric perturbation h_v' in synchronous gauge */
  int index_pt_V;         /**< vector metric perturbation V in Newtonian gauge */

  int index_pt_gw;        /**< tensor metric perturbation h (gravitational waves) */
  int index_pt_gwdot;     /**< its time-derivative */
  int pt_size;            /**< size of perturbation vector */

  double * y;             /**< vector of perturbations to be integrated */
  double * dy;            /**< time-derivative of the same vector */

  int * used_in_sources; /**< boolean array specifying which
                            perturbations enter in the calculation of
                            source functions */

};


/**
 * Workspace containing, among other things, the value at a given time
 * of all background/perturbed quantitites, as well as their indices.
 *
 * There will be one such structure created for each mode
 * (scalar/.../tensor) and each thread (in case of parallel computing)
 */

struct perturb_workspace
{

  /** @name - all possible useful indices for those metric
      perturbations which are not integrated over time, but just
      inferred from Einstein equations. "_mt_" stands for "metric".*/

  //@{

  int index_mt_psi;           /**< psi in longitudinal gauge */
  int index_mt_phi_prime;     /**< (d phi/d conf.time) in longitudinal gauge */
  int index_mt_h_prime;       /**< h' (wrt conf. time) in synchronous gauge */
  int index_mt_h_prime_prime; /**< h'' (wrt conf. time) in synchronous gauge */
  int index_mt_eta_prime;     /**< eta' (wrt conf. time) in synchronous gauge */
  int index_mt_alpha;         /**< \alpha = (h' + 6 \eta') / (2 k^2) \f$ in synchronous gauge */
  int index_mt_alpha_prime;   /**< alpha' wrt conf. time) in synchronous gauge */
  int index_mt_gw_prime_prime;/**< second derivative wrt confromal time of gravitational wave field, often called h */
  int index_mt_V_prime;       /**< derivative of Newtonian gauge vector metric perturbation V */
  int index_mt_hv_prime_prime;/**< Second derivative of Synchronous gauge vector metric perturbation h_v */
  int mt_size;                /**< size of metric perturbation vector */

  //@}

  /** @name - value at a given time of all background/perturbed
      quantitites
  */

  //@{

  double * pvecback;          /**< background quantitites */
  double * pvecthermo;        /**< thermodynamics quantitites */
  double * pvecmetric;        /**< metric quantitites */
  struct perturb_vector * pv; /**< pointer to vector of integrated
                                 perturbations and their
                                 time-derivatives */

  double delta_rho;
  double rho_plus_p_theta;
  double rho_plus_p_shear;
  double delta_p;
  double gw_source;
  double vector_source_pi;
  double vector_source_v;

  double tca_shear_g; /**< photon shear in tight-coupling approximation */
  double tca_slip;    /**< photon-baryon slip in tight-coupling approximation */
  double rsa_delta_g; /**< photon density in radiation streaming approximation */
  double rsa_theta_g; /**< photon velocity in radiation streaming approximation */
  double rsa_delta_ur; /**< photon density in radiation streaming approximation */
  double rsa_theta_ur; /**< photon velocity in radiation streaming approximation */

  double * delta_ncdm;
  double * theta_ncdm;
  double * shear_ncdm;

  double delta_m;
  double theta_m;

  FILE * perturb_output_file; /**< filepointer to output file*/

  //@}

  /** @name - indices useful for searching background/termo quantitites in tables */

  //@{

  short inter_mode;

  int last_index_back;   /**< the background interpolation function background_at_tau() keeps memory of the last point called through this index */
  int last_index_thermo; /**< the thermodynamics interpolation function thermodynamics_at_z() keeps memory of the last point called through this index */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int index_ap_tca; /**< index for tight-coupling approximation */
  int index_ap_rsa; /**< index for radiation streaming approximation */
  int index_ap_ufa; /**< index for ur fluid approximation */
  int index_ap_ncdmfa; /**< index for ncdm fluid approximation */
  int ap_size;      /**< number of relevant approximations for a given mode */

  int * approx;     /**< array of approximation flags holding at a given time: approx[index_ap] */

  //@}

  /** @name - approximations used at a given time */

  //@{

  int max_l_max;    /**< maximum l_max for any multipole */
  double * s_l;     /**< array of freestreaming coefficients s_l = sqrt(1-K*(l^2-1)/k^2) */

  //@}

};

/**
 * Structure pointing towards all what the function that perturb_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
 */

struct perturb_parameters_and_workspace {

  struct precision * ppr;         /**< pointer to the precision structure */
  struct background * pba;        /**< pointer to the background structure */
  struct thermo * pth;            /**< pointer to the thermodynamics structure */
  struct perturbs * ppt;          /**< pointer to the precision structure */
  int index_md;                 /**< index of mode (scalar/.../vector/tensor) */
  int index_ic;
  int index_k;
  double k;
  struct perturb_workspace * ppw; /**< worspace defined above */

};

/*************************************************************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int perturb_sources_at_tau(
                             struct perturbs * ppt,
                             int index_md,
                             int index_ic,
                             int index_type,
                             double tau,
                             double * pvecsources
                             );

  int perturb_init(
                   struct precision * ppr,
                   struct background * pba,
                   struct thermo * pth,
                   struct perturbs * ppt
                   );

  int perturb_free(
                   struct perturbs * ppt
                   );

  int perturb_indices_of_perturbs(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct thermo * pth,
                                  struct perturbs * ppt
                                  );

  int perturb_timesampling_for_sources(
                                       struct precision * ppr,
                                       struct background * pba,
                                       struct thermo * pth,
                                       struct perturbs * ppt
                                       );
  int perturb_get_k_list(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt
                         );

  int perturb_workspace_init(
                             struct precision * ppr,
                             struct background * pba,
                             struct thermo * pth,
                             struct perturbs * ppt,
                             int index_md,
                             struct perturb_workspace * ppw
                             );

  int perturb_workspace_free(
                             struct perturbs * ppt,
                             int index_md,
                             struct perturb_workspace * ppw
                             );

  int perturb_solve(
                    struct precision * ppr,
                    struct background * pba,
                    struct thermo * pth,
                    struct perturbs * ppt,
                    int index_md,
                    int index_ic,
                    int index_k,
                    struct perturb_workspace * ppw
                    );

  int perturb_find_approximation_number(
                                        struct precision * ppr,
                                        struct background * pba,
                                        struct thermo * pth,
                                        struct perturbs * ppt,
                                        int index_md,
                                        double k,
                                        struct perturb_workspace * ppw,
                                        double tau_ini,
                                        double tau_end,
                                        int * interval_number,
                                        int * interval_number_of
                                        );

  int perturb_find_approximation_switches(
                                          struct precision * ppr,
                                          struct background * pba,
                                          struct thermo * pth,
                                          struct perturbs * ppt,
                                          int index_md,
                                          double k,
                                          struct perturb_workspace * ppw,
                                          double tau_ini,
                                          double tau_end,
                                          double precision,
                                          int interval_number,
                                          int * interval_number_of,
                                          double * interval_limit,
                                          int ** interval_approx
                                          );

  int perturb_vector_init(
                          struct precision * ppr,
                          struct background * pba,
                          struct thermo * pth,
                          struct perturbs * ppt,
                          int index_md,
                          int index_ic,
                          double k,
                          double tau,
                          struct perturb_workspace * ppw,
                          int * pa_old
                          );

  int perturb_vector_free(
                          struct perturb_vector * pv
                          );

  int perturb_initial_conditions(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct perturbs * ppt,
                                 int index_md,
                                 int index_ic,
                                 double k,
                                 double tau,
                                 struct perturb_workspace * ppw
                                 );

  int perturb_approximations(
                             struct precision * ppr,
                             struct background * pba,
                             struct thermo * pth,
                             struct perturbs * ppt,
                             int index_md,
                             double k,
                             double tau,
                             struct perturb_workspace * ppw
                             );

  int perturb_timescale(
                        double tau,
                        void * parameters_and_workspace,
                        double * timescale,
                        ErrorMsg error_message
                        );

  int perturb_einstein(
                       struct precision * ppr,
                       struct background * pba,
                       struct thermo * pth,
                       struct perturbs * ppt,
                       int index_md,
                       double k,
                       double tau,
                       double * y,
                       struct perturb_workspace * ppw
                       );

  int perturb_total_stress_energy(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct thermo * pth,
                                  struct perturbs * ppt,
                                  int index_md,
                                  double k,
                                  double * y,
                                  struct perturb_workspace * ppw
                                  );

  int perturb_sources(
                      double tau,
                      double * pvecperturbations,
                      double * pvecderivs,
                      int index_tau,
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      );

  int perturb_print_variables(
                              double tau,
                              double * y,
                              double * dy,
                              void * parameters_and_workspace,
                              ErrorMsg error_message
                              );

  int perturb_derivs(
                     double tau,
                     double * y,
                     double * dy,
                     void * parameters_and_workspace,
                     ErrorMsg error_message
                     );

  int perturb_tca_slip_and_shear(
                                 double * y,
                                 void * parameters_and_workspace,
                                 ErrorMsg error_message
                                 );

  int perturb_rsa_delta_and_theta(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct thermo * pth,
                                  struct perturbs * ppt,
                                  double k,
                                  double * y,
                                  double a_prime_over_a,
                                  double * pvecthermo,
                                  struct perturb_workspace * ppw
                                  );

  int perturb_prepare_output_file(struct background * pba,
                                  struct perturbs * ppt,
                                  struct perturb_workspace * ppw,
                                  int index_ikout,
                                  int index_md);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
