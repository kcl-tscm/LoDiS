SUBROUTINE sub_force_choice

  USE PARACLUSTER
  USE CLUSTER 
  USE POTENTIAL
  USE ENFORCE
  USE DISTANCE
  USE SUBSTRATE

  IMPLICIT NONE
  
  INTEGER :: sub_type

  sub_type = sub_geom

SELECT CASE (sub_type) 
    	CASE (1) 
	if(ipas==1) Write(*,*) "Substrate chosen with a 1-geometry type"
   	   CALL sub_dsq_force
        if(ipas==1) Write(*,*) "After calling dsq_force"
	CASE (2) 
	if(ipas==1) Write(*,*) "Substrate chosen with a 2-geometry type"
   	   CALL sub_hex_force
END SELECT

END SUBROUTINE sub_force_choice
