module distance
   ! LP
   ! Global variables storing pair distances

   real(8),    allocatable :: pair_dist (:, :)
   INTEGER(4), ALLOCATABLE :: pairkindmat (:, :)

   ! xi-xj, yi-yj, zi-zj in array form of dimension (3, natom, natom)
   REAL(8),    ALLOCATABLE :: xyz_dist (:, :, :) 
   
   LOGICAL :: cutoffmat_req = .FALSE. 
   LOGICAL, ALLOCATABLE, DIMENSION (:,:) :: cutoffmat
   REAL(8), DIMENSION(3) :: dcutoff
    
end module distance
