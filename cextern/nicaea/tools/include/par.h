#ifndef __PAR_H
#define __PAR_H

#include <stdio.h>

#include "errorlist.h"
#include "io.h"
#include "config.h"

/* Parameters */
typedef enum {
  p_Omegam, p_Omegab, p_Omegade, p_h100, p_Omeganumass,
  p_Omegac, p_OmegaK,
  p_omegam, p_omegab, p_100_omegab, p_omegade, p_omeganumass, p_omegac, p_omegaK,
  p_w0de, p_w1de, p_wpolyde, p_sigma8, p_Delta2R, p_ns, p_Neffnumass,
  p_tau, p_alphas, p_Neffnu0, p_nt, p_r, p_lnr,
  p_a_ymmk, p_b_ymmk, p_c_ymmk, p_c0_ymmk0const, p_a_jonben, p_b_jonben, p_c_jonben,
  p_alpha_ludo, p_beta_ludo, p_z0_ludo, p_M, p_alpha, p_beta, p_beta_z, p_logbeta,
  p_stretch, p_color,
  p_theta10, p_theta11, p_theta12, p_theta13, p_theta14, p_theta15, p_theta16, p_theta17, p_beta_d,
  p_mring_a, p_mring_b, p_Mmin, p_log10Mmin, p_M1, p_log10M1, p_M0, p_log10M0, p_sigma_log_M,
  p_log10Mhalo,
  p_Mstar0, p_log10Mstar0, p_delta, p_gamma, p_B_cut, p_B_sat, p_beta_cut, p_beta_sat,
  p_alpha_halo, p_eta, p_fcen1, p_fcen2, p_Mcut, p_Ng_min, p_Ng_max, p_A_fromz,
  p_z_rescale, p_Sigma, p_qacc, p_Mhalo_av, p_log10Mhalo_av, p_bgal_av, p_Ngal_av, p_fr_sat, p_ngal_den,
  p_log10ngal_den,
  p_A_SZ, p_A_ia, p_A_GGI, p_theta_GGI, p_A_GII, p_theta_GII, p_b_slc, p_gamma_slc,
  p_kb, p_kr, p_wa,
  p_bias_g_amp, p_bias_g_zexp, p_bias_g_Lexp, p_bias_IA_amp, p_bias_IA_zexp, p_bias_IA_Lexp,
  p_dummy,
} par_t;
#define Npar_t 108

#define spar_t(i) ( \
  i==p_Omegam ? "Omega_m" : \
  i==p_Omegab ? "Omega_b" : \
  i==p_Omegade ? "Omega_de" : \
  i==p_h100 ? "h_100" : \
  i==p_Omeganumass ? "Omega_nu_mass" : \
  i==p_Omegac ? "Omega_c" : \
  i==p_OmegaK ? "Omega_K" : \
  i==p_omegam ? "omega_m" : \
  i==p_omegab ? "omega_b" : \
  i==p_100_omegab ? "100_omega_b" : \
  i==p_omegade ? "omega_de" : \
  i==p_omeganumass ? "omega_nu_mass" : \
  i==p_omegac ? "omega_c" : \
  i==p_omegaK ? "omega_K" : \
  i==p_w0de ? "w_0_de" : \
  i==p_w1de ? "w_1_de" : \
  i==p_wpolyde ? "w_poly_de" : \
  i==p_sigma8 ? "sigma_8" : \
  i==p_Delta2R ? "Delta_2_R" : \
  i==p_ns ? "n_s" : \
  i==p_Neffnumass ? "N_eff_nu_mass" : \
  i==p_tau ? "tau" : \
  i==p_alphas ? "alpha_s" : \
  i==p_Neffnu0 ? "N_eff_nu_0" : \
  i==p_nt ? "n_t" : \
  i==p_r ? "r" : \
  i==p_lnr ? "ln_r" : \
  i==p_a_ymmk ? "a_ymmk" : \
  i==p_b_ymmk ? "b_ymmk" : \
  i==p_c_ymmk ? "c_ymmk" : \
  i==p_c0_ymmk0const ? "c0_ymmk0const" : \
  i==p_a_jonben ? "a_jonben" : \
  i==p_b_jonben ? "b_jonben" : \
  i==p_c_jonben ? "c_jonben" : \
  i==p_alpha_ludo ? "alpha_ludo" : \
  i==p_beta_ludo ? "beta_ludo" : \
  i==p_z0_ludo ? "z0_ludo" : \
  i==p_M ? "M" : \
  i==p_alpha ? "alpha" : \
  i==p_beta ? "beta" : \
  i==p_beta_z ? "beta_z" : \
  i==p_logbeta ? "log_beta" : \
  i==p_stretch ? "stretch" : \
  i==p_color ? "color" : \
  i==p_theta10 ? "theta_10" : \
  i==p_theta11 ? "theta_11" : \
  i==p_theta12 ? "theta_12" : \
  i==p_theta13 ? "theta_13" : \
  i==p_theta14 ? "theta_14" : \
  i==p_theta15 ? "theta_15" : \
  i==p_theta16 ? "theta_16" : \
  i==p_theta17 ? "theta_17" : \
  i==p_beta_d  ? "beta_d"  : \
  i==p_mring_a ? "mring_a" : \
  i==p_mring_b ? "mring_b" : \
  i==p_Mmin    ? "Mmin"    : \
  i==p_log10Mmin ? "log10Mmin" : \
  i==p_M1      ? "M1"      : \
  i==p_log10M1 ? "log10M1" : \
  i==p_M0      ? "M0"      : \
  i==p_log10M0 ? "log10M0" : \
  i==p_sigma_log_M ? "sigma_log_M" : \
  i==p_log10Mhalo ? "log10Mhalo" : \
  i==p_alpha_halo  ? "alpha_halo" : \
  i==p_eta  ? "eta" : \
  i==p_fcen1  ? "fcen1" : \
  i==p_fcen2  ? "fcen2" : \
  i==p_Mcut    ? "M_cut"    :  \
  i==p_Mstar0    ? "Mstar0"   :  \
  i==p_log10Mstar0 ? "log10Mstar0"   :  \
  i==p_delta    ? "delta"    :   \
  i==p_gamma    ? "gamma"    :	    \
  i==p_B_cut    ? "B_cut"    :	    \
  i==p_B_sat    ? "B_sat"    :	    \
  i==p_beta_cut ? "beta_cut"    :   \
  i==p_beta_sat ? "beta_sat"  :   \
  i==p_Ng_min  ? "Ng_min"  : \
  i==p_Ng_max  ? "Ng_max"  : \
  i==p_A_fromz ? "A_from_z" : \
  i==p_z_rescale   ? "z_rescale" : \
  i==p_Sigma   ? "Sigma" : \
  i==p_qacc    ? "q_acc" : \
  i==p_Mhalo_av ? "Mhalo_av" : \
  i==p_log10Mhalo_av ? "log10Mhalo_av" : \
  i==p_bgal_av ? "bgal_av" : \
  i==p_Ngal_av ? "Ngal_av" : \
  i==p_fr_sat   ? "fr_sat" : \
  i==p_ngal_den ? "ngal_den" : \
  i==p_log10ngal_den ? "log10ngal_den" : \
  i==p_A_SZ          ? "A_SZ" : \
  i==p_A_ia          ? "A_ia" : \
  i==p_A_GGI         ? "A_GGI" :     \
  i==p_theta_GGI     ? "theta_GGI" : \
  i==p_A_GII         ? "A_GII" :	 \
  i==p_theta_GII     ? "theta_GII" :     \
  i==p_b_slc         ? "b_slc" :	     \
  i==p_gamma_slc     ? "gamma_slc" :	     \
  i==p_kb            ? "kb" : \
  i==p_kr            ? "kr"  : \
  i==p_wa            ? "wa" : \
  i==p_bias_g_amp            ? "bias_g_amp" : \
  i==p_bias_g_zexp           ? "bias_g_zexp" : \
  i==p_bias_g_Lexp           ? "bias_g_Lexp" : \
  i==p_bias_IA_amp            ? "bias_IA_amp" : \
  i==p_bias_IA_zexp           ? "bias_IA_zexp" : \
  i==p_bias_IA_Lexp           ? "bias_IA_Lexp" : \
  i==p_dummy         ? "dummy"   :   \
  "" )

par_t *copy_par_t(const par_t *par, int npar, error **err);
void spar_to_par(par_t **par, int npar, const char *spar[], error **err);


#endif
