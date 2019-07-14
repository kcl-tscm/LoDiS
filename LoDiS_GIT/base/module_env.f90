MODULE ENVIRONMENT

 !===================================
 ! Global variables for environment
 !===================================

 IMPLICIT NONE
  
 REAL(8), ALLOCATABLE :: env_cn(:), env_s(:,:), env_dS(:,:) 
 REAL(8), ALLOCATABLE :: env_dS_dx(:,:), env_dS_dy(:,:), env_dS_dz(:,:)  
 REAL(8), ALLOCATABLE :: env_dcn_dx(:), env_dcn_dy(:), env_dcn_dz(:)  
 REAL(8), ALLOCATABLE :: ener_env_atom(:)
 REAL(8), ALLOCATABLE :: env_fx(:), env_fy(:), env_fz(:)
 REAL(8)      :: ener_env                                              ! potential energy environment
 REAL(8)      :: env_n_pwr =  6
 REAL(8)      :: env_m_pwr = 12
 REAL(8)      :: env_rzero = 0.147d0
 REAL(8)      :: env_ratio_ton
 REAL(8)      :: env_ratio_tom
 REAL(8)      :: env_rnp
 REAL(8)      :: env_rmp
 REAL(8)      :: env_num
 REAL(8)      :: env_comfactor
 REAL(8)      :: env_ratio
 REAL(8)      :: alpha_env = 30.5d0
 REAL(8)      :: rc_env = 3.5d0
 REAL(8)      :: beta_env, beta_env_angstrom

END MODULE ENVIRONMENT
