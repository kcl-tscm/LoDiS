MODULE module_sticky
 !===================================
 ! Global variables for sticky atoms
 !===================================
 IMPLICIT NONE
  
 REAL(8)         :: sticky_k
 LOGICAL         :: sticky_atoms_wanted
 
 INTEGER, ALLOCATABLE :: sticky_labels(:) 
 REAL(8), ALLOCATABLE :: sticky_fx(:), sticky_fy(:), sticky_fz(:)
 
END MODULE module_sticky
