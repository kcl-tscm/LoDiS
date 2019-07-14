MODULE module_function_var
 !===================================
 ! Global variables for Force
 !===================================
 IMPLICIT NONE
 
 REAL(8) :: funct1, funct2, funct1_dx, funct2_dx, funct1_dy, funct2_dy, exp_term1, exp_term2, exp_term3
 REAL(8), DIMENSION(3)   :: exp_term
 REAL(8), DIMENSION(3)   :: a_, da_dx, da_dy, da_dz
 REAL(8), DIMENSION(3,3) :: b, db_dx, db_dy

END MODULE module_function_var
