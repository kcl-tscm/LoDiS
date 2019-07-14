SUBROUTINE force_choice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Depending on the potential a different subroutine will be called
!
!RGL potential for metallic systems (noble and quasi-noble)
!V. Rosato, M. Guillope, and B. Legrand, 1989, Philos. Mag. A 59, 321.
!A detailed ddescription in C. Mottet PhD thesis
!A very short descritpion in F. Baletto PhD thesis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Lj1 == Lennard-Jones potential, the parameters are for Argon
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Semi-empirical Description of C60 molecules through 
!Gir == Girifalco, L.A. Girifalco, 1992, J. Phys. Chem. 96, 858
!Par == Pacheco, J. Pacheco et al. 1997, Phys. Rev. Lett. 79, 3873
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE CLUSTER, only : type_potential     

Implicit None

if(type_potential=='rgl') then
 call force_rgl
elseif(type_potential=='lj1') then
 call force_lj1
elseif(type_potential=='gir') then
 call force_gir
elseif(type_potential=='par') then
 call force_par
endif

END SUBROUTINE force_choice

