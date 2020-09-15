MODULE META
!==================================
! Global variables for metadynamics
!==================================
  IMPLICIT NONE
  
  ! Input Parameters
  integer :: m_pwr, n_pwr, num_cv
  integer :: m2n_pwr, n2n_pwr, msf_pwr, nsf_pwr
  integer :: metaperiod, NG, metaframe
  real(8) :: gheight, gwidth, g2nwidth, gsfwidth
  real(8) :: rzero               ! parameter of NN distance distribution [bulk latt ref]
  real(8) :: r2n, rsf, d2n, dsf  ! parameters [bulk latt ref]
  logical :: collvar_wanted      ! to calculate CV also during a microcan
  INTEGER :: cn_n_pwr, cn_m_pwr  ! cn bim 
  REAL(8) :: cn_rzero, cn_gwidth ! cn bim
  LOGICAL :: cn_aa, cn_bb, cn_ab ! cn bim
  LOGICAL :: yesno_wanted        ! cn bim
  REAL(8) :: d_com_gwidth        ! d_com
  INTEGER :: natoma, natomb      ! d_com
  LOGICAL :: d_com_wanted        ! d_com
  character(LEN=20), dimension(2) :: collvar_name ! deps on number of cv -- max 2
  integer, dimension(2)           :: collvar_case ! to associate a CV to a number (for select statement)

  ! Common Neighbour Number
  REAL(8) :: cnf_num, gcnnwidth, inflim, suplim
  REAL(8), DIMENSION(14):: cnf_in
  REAL(8), DIMENSION(5) :: win_in
   
  ! arrays used as local, but defined global not to reallocate it
  double precision, allocatable :: deno (:), halfdeno(:)
  double precision, allocatable :: dS(:,:)
  double precision, allocatable :: dVg_ds(:)
  double precision, allocatable :: expoc(:)
  
  double precision, allocatable :: gausswidth(:)
  double precision, ALLOCATABLE :: s_of_t(:,:)
  double precision, ALLOCATABLE :: coll_var(:)   ! vector of current values of collective variables (num_cv)

  ! Meta-force Computation CN
  double precision, ALLOCATABLE :: s(:,:)        ! s(atom-i, atom-j)
  double precision, ALLOCATABLE :: CNatom(:)     ! coord. no. of every atoms
  double precision, ALLOCATABLE :: fgx(:),fgy(:),fgz(:) ! to natom
  DOUBLE PRECISION, ALLOCATABLE :: dS_dx(:), dS_dy(:), dS_dz(:) ! to natom
  
  ! Meta-force Computation C2N
  DOUBLE PRECISION, ALLOCATABLE :: s2n(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: C2Natom(:)
  DOUBLE PRECISION, ALLOCATABLE :: fg2nx(:), fg2ny(:), fg2nz(:)
  DOUBLE PRECISION, ALLOCATABLE :: dS2n_dx(:), dS2n_dy(:), dS2n_dz(:)

  ! Meta-force Computation SFN
  DOUBLE PRECISION, ALLOCATABLE :: ssf(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: SFNatom(:)
  DOUBLE PRECISION, ALLOCATABLE :: fgsfx(:), fgsfy(:), fgsfz(:)
  DOUBLE PRECISION, ALLOCATABLE :: dSsf_dx(:), dSsf_dy(:), dSsf_dz(:)

  ! Meta-force Computation CNF
  DOUBLE PRECISION, ALLOCATABLE :: fgcnnx(:), fgcnny(:), fgcnnz(:)
  DOUBLE PRECISION, ALLOCATABLE :: dScnn_dx(:), dScnn_dy(:), dScnn_dz(:)

  ! Meta-force Computation CN_bim
  double precision, ALLOCATABLE :: cn_s(:,:)        ! s(atom-i, atom-j)
  double precision, ALLOCATABLE :: cn_CNatom(:)     ! coord. no. of every atoms
  double precision, ALLOCATABLE :: cn_fgx(:),cn_fgy(:),cn_fgz(:) ! to natom
  DOUBLE PRECISION, ALLOCATABLE :: cn_dS_dx(:), cn_dS_dy(:), cn_dS_dz(:) ! to natom 
  LOGICAL, ALLOCATABLE, DIMENSION (:,:) :: pair_yesno_mat

  ! Meta-force Computation d_com
  DOUBLE PRECISION              :: d_com_cv 
  DOUBLE PRECISION, ALLOCATABLE :: d_com_dcv_dx(:), d_com_dcv_dy(:), d_com_dcv_dz(:) 
  DOUBLE PRECISION, ALLOCATABLE :: d_com_fgx(:), d_com_fgy(:), d_com_fgz(:) 

  CONTAINS

  !----------------------------------------------------------------------------
  ! Sigmoid function with and without derivatives
  !----------------------------------------------------------------------------
  ELEMENTAL FUNCTION f_sigmoid(dummy, pairkind)
     ! This function is elemental because is called using array arguments
     ! inside a WHERE statement (in comneighfun_no_light)
        
     ! cnf_in(1)  = d0
     ! cnf_in(2)  = r0
     ! cnf_in(3)  = c0
     ! cnf_in(4)  = m
     ! cnf_in(5)  = n
     ! cnf_in(6)  = d0 A-A pairkind = 1    
     ! cnf_in(7)  = d0 B-B pairkind = 2
     ! cnf_in(8)  = d0 A-B pairkind = 3
     ! cnf_in(9)  = r0 A-A pairkind = 1
     ! cnf_in(10) = r0 B-B pairkind = 2
     ! cnf_in(11) = r0 A-B pairkind = 3
     ! cnf_in(12) = c0 A-A pairkind = 1
     ! cnf_in(13) = c0 B-B pairkind = 2
     ! cnf_in(14) = c0 A-B pairkind = 3

     ! cnf_calc(1) = rij 
     ! cnf_calc(2) = rij-d0 
     ! cnf_calc(3) = h0
     ! cnf_calc(4) = (rij-d0)/r0
     ! cnf_calc(4) = (rij-d0-c0)/(r0-c0)

     IMPLICIT NONE
     REAL(8) :: f_sigmoid
     REAL(8), DIMENSION(4) :: cnf_calc 
     REAL(8), INTENT(IN) :: dummy
     INTEGER(4), INTENT(IN) :: pairkind
     

     cnf_calc(1) = dummy
     cnf_calc(2) = cnf_calc(1) - cnf_in(pairkind+5) ! cnf_calc(2) = (rij - d0)
     
     IF (cnf_calc(2).LE.0.d0) THEN
        f_sigmoid = 1.d0
     ELSE
        cnf_calc(3) = 1/(1-(cnf_in(5)*cnf_in(pairkind+8)/ &
                   & ( cnf_in(4)*(cnf_in(pairkind+8)-cnf_in(pairkind+11)) )) )

        IF (    (0.d0.LT.cnf_calc(2)) .AND.(cnf_calc(2).LE.cnf_in(pairkind+8))) THEN
           cnf_calc(4) = cnf_calc(2)/cnf_in(pairkind+8)
           f_sigmoid = (cnf_calc(3)-1.0d0)* cnf_calc(4)**cnf_in(4) + 1.0d0
     
        ELSE  ! (Note: this will work in combination with a cutoff matrix)
           cnf_calc(4) = (cnf_calc(2)-cnf_in(pairkind+11))/(cnf_in(pairkind+8)-cnf_in(pairkind+11))
           f_sigmoid = cnf_calc(3)* cnf_calc(4)**cnf_in(5)    
     
        END IF
     END IF 
  END FUNCTION f_sigmoid


  FUNCTION f_sigmoid_der(dummy, pairkind)
        
     ! cnf_in(1)  = d0
     ! cnf_in(2)  = r0
     ! cnf_in(3)  = c0
     ! cnf_in(4)  = m
     ! cnf_in(5)  = n
     ! cnf_in(6)  = d0 A-A pairkind = 1    
     ! cnf_in(7)  = d0 B-B pairkind = 2
     ! cnf_in(8)  = d0 A-B pairkind = 3
     ! cnf_in(9)  = r0 A-A pairkind = 1
     ! cnf_in(10) = r0 B-B pairkind = 2
     ! cnf_in(11) = r0 A-B pairkind = 3
     ! cnf_in(12) = c0 A-A pairkind = 1
     ! cnf_in(13) = c0 B-B pairkind = 2
     ! cnf_in(14) = c0 A-B pairkind = 3

     ! cnf_calc(1) = rij 
     ! cnf_calc(2) = rij-d0 
     ! cnf_calc(3) = h0
     ! cnf_calc(4) = (rij-d0)/r0
     ! cnf_calc(4) = (rij-d0-c0)/(r0-c0)
     ! cnf_calc(5) = (h0-1) * ((rij-d0)/r0)**(m-1)
     ! cnf_calc(5) = h0 * ((rij-d0-c0)/(r0-c0))**(n-1)

     ! f_sigmoid_der(1) = f_sigmoid
     ! f_sigmoid_der(2) = derivative of f_sigmoid with respect to the pair distance rij

     IMPLICIT NONE
     REAL(8), DIMENSION(2) :: f_sigmoid_der
     REAL(8), DIMENSION(5) :: cnf_calc 
     REAL(8)               :: dummy
     INTEGER(4)            :: pairkind
     

     cnf_calc(1) = dummy
     cnf_calc(2) = cnf_calc(1) - cnf_in(pairkind+5) ! cnf_calc(2) = (rij - d0)
     
     IF (cnf_calc(2).LE.0.d0) THEN
        f_sigmoid_der(1) = 1.d0
        f_sigmoid_der(2) = 0.d0
     ELSE
        cnf_calc(3) = 1/(1-(cnf_in(5)*cnf_in(pairkind+8)/ &
                   & ( cnf_in(4)*(cnf_in(pairkind+8)-cnf_in(pairkind+11)) )) )

        IF (    (0.d0.LT.cnf_calc(2)) .AND.(cnf_calc(2).LE.cnf_in(pairkind+8))) THEN
           cnf_calc(4) = cnf_calc(2)/cnf_in(pairkind+8)
           cnf_calc(5) = (cnf_calc(3)-1.d0) * cnf_calc(4)**(cnf_in(4)-1)

           f_sigmoid_der(1) = cnf_calc(4) * cnf_calc(5) +1.0d0
           f_sigmoid_der(2) = cnf_in(4) * (1/cnf_in(pairkind+8)) * cnf_calc(5) 
     
        ELSE   ! (Note: this will work in combination with a cutoff matrix)
           cnf_calc(4) = (cnf_calc(2)-cnf_in(pairkind+11))/(cnf_in(pairkind+8)-cnf_in(pairkind+11))
           cnf_calc(5) = cnf_calc(3) * cnf_calc(4)**(cnf_in(5)-1)

           f_sigmoid_der(1) = cnf_calc(4) * cnf_calc(5)   
           f_sigmoid_der(2) = cnf_in(5) * (1/(cnf_in(pairkind+8)-cnf_in(pairkind+11))) * cnf_calc(5)   

        END IF
     END IF 
  END FUNCTION f_sigmoid_der

  !----------------------------------------------------------------------------
  ! Window function with and without derivatives
  !----------------------------------------------------------------------------
  FUNCTION f_window(dummy)
     
     ! win_in(1)   = d0
     ! win_in(2)   = r0
     ! win_in(3)   = c0
     ! win_in(4)   = m     
     ! win_in(5)   = n
     ! 
     ! win_calc(1) = rij 
     ! win_calc(2) = rij-d0 
     ! win_calc(3) = h0
     ! win_calc(4) = (rij-d0)/r0
     ! win_calc(4) = (-rij+d0)/r0
     ! win_calc(4) = (rij-d0-c0)/(r0-c0)
     ! win_calc(4) = (-rij+d0-c0)/(r0-c0)

     IMPLICIT NONE
     REAL(8) :: f_window, dummy
     REAL(8), DIMENSION(4) :: win_calc
     

     win_calc(1) = dummy
     win_calc(2) = win_calc(1) - win_in(1) ! win_calc(2) = (rij - d0)
     
     IF ((win_calc(2).LE.-win_in(3)).OR.(win_calc(2).GE.win_in(3))) THEN
        f_window = 0.d0
     ELSE
        win_calc(3) = 1/(1-(win_in(5)*win_in(2)/ &
                      & ( win_in(4)*(win_in(2)-win_in(3)) )) )
        IF (    (0.d0.LE.win_calc(2)) .AND.(win_calc(2).LE.win_in(2))) THEN
           win_calc(4) = win_calc(2)/win_in(2)
           f_window = (win_calc(3)-1.0d0)* win_calc(4)**win_in(4) + 1.0d0
     
        ELSE IF ((-win_in(2).LE.win_calc(2)).AND.(win_calc(2).LT. 0.d0)) THEN
           win_calc(4) = -win_calc(2)/win_in(2)
           f_window = (win_calc(3)-1.0d0)* (win_calc(4))**win_in(4) + 1.0d0
     
        ELSE IF ((win_in(2).LT.win_calc(2)).AND.(win_calc(2).LT. win_in(3))) THEN
           win_calc(4) = (win_calc(2)-win_in(3))/(win_in(2)-win_in(3))
           f_window = win_calc(3)* win_calc(4)**win_in(5)
     
        ELSE
           win_calc(4) = (-win_calc(2)-win_in(3))/(win_in(2)-win_in(3))
           f_window = win_calc(3)* win_calc(4)**win_in(5)        
        END IF
     END IF
  END FUNCTION f_window


  FUNCTION f_window_der(dummy)
     
     ! win_in(1)   = d0
     ! win_in(2)   = r0
     ! win_in(3)   = c0
     ! win_in(4)   = m     
     ! win_in(5)   = n
     ! 
     ! win_calc(1) = rij 
     ! win_calc(2) = rij-d0 
     ! win_calc(3) = h0
     ! win_calc(4) = (rij-d0)/r0
     ! win_calc(4) = (-rij+d0)/r0
     ! win_calc(4) = (rij-d0-c0)/(r0-c0)
     ! win_calc(4) = (-rij+d0-c0)/(r0-c0)
     ! win_calc(5) = (h0-1)* ((rij-d0)/r0)**(m-1)
     ! win_calc(5) = (h0-1)* ((-rij+d0)/r0)**(m-1)
     ! win_calc(5) = h0 * (( rij-d0-c0)/(r0-c0))**(n-1)
     ! win_calc(5) = h0 * ((-rij+d0-c0)/(r0-c0))**(n-1)

     ! f_window_der(1) = f_window
     ! f_window_der(2) = derivative of f_window with respect of l_ij

     IMPLICIT NONE
     REAL(8) :: dummy
     REAL(8), DIMENSION(2) :: f_window_der
     REAL(8), DIMENSION(5) :: win_calc
     

     win_calc(1) = dummy
     win_calc(2) = win_calc(1) - win_in(1) ! win_calc(2) = (rij - d0)
          
     IF ((win_calc(2).LE.-win_in(3)).OR.(win_calc(2).GE.win_in(3))) THEN
        f_window_der(1:2) = 0.d0
     ELSE
        win_calc(3) = 1/(1-(win_in(5)*win_in(2)/ &
                   & ( win_in(4)*(win_in(2)-win_in(3)) )) )

        IF (    (0.d0.LE.win_calc(2)) .AND.(win_calc(2).LE.win_in(2))) THEN
           win_calc(4) = win_calc(2)/win_in(2)
           win_calc(5) = (win_calc(3)-1.0d0)* win_calc(4)**(win_in(4)-1)
           
           f_window_der(1) = win_calc(4) * win_calc(5) + 1.0d0
           f_window_der(2) = (1/win_in(2)) * win_in(4) * win_calc(5) 
     
        ELSE IF ((-win_in(2).LE.win_calc(2)).AND.(win_calc(2).LT. 0.d0)) THEN
           win_calc(4) = -win_calc(2)/win_in(2)
           win_calc(5) = (win_calc(3)-1.0d0)* (win_calc(4))**(win_in(4)-1)
           
           f_window_der(1) = win_calc(4) * win_calc(5) + 1.0d0
           f_window_der(2) = -(1/win_in(2)) * win_in(4) * win_calc(5) 
     
        ELSE IF ((win_in(2).LT.win_calc(2)).AND.(win_calc(2).LT. win_in(3))) THEN
           win_calc(4) = (win_calc(2)-win_in(3))/(win_in(2)-win_in(3))
           win_calc(5) = win_calc(3)* win_calc(4)**(win_in(5)-1)
           
           f_window_der(1) = win_calc(4) * win_calc(5) 
           f_window_der(2) = (1/(win_in(2)-win_in(3))) * win_in(5) * win_calc(5) 
     
        ELSE
           win_calc(4) = (-win_calc(2)-win_in(3))/(win_in(2)-win_in(3))
           win_calc(5) = win_calc(3)* win_calc(4)**(win_in(5)-1)
           
           f_window_der(1) = win_calc(4) * win_calc(5)
           f_window_der(2) = -(1/(win_in(2)-win_in(3))) * win_in(5) * win_calc(5) 
        
        END IF     
     END IF  
          
  END FUNCTION f_window_der


END MODULE META



