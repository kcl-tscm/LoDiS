Subroutine subnomi
!
USE PARACLUSTER
USE CLUSTER     !variabili per il file del cluster            
USE ENFORCE
!
Implicit None
Integer :: i,nd

!!if(quenching=='ya') nd=1
!if (metadynamics=='ya') nd=1
nd=nd_proc
!!        ndasaltare=1000*nd    !se nd<9 =>     1000<i<9100 oppure 9300
                              !se nd<99 =>   10000<i<99100 oppure 99300
                              !se nd<999 => 100000<i<999100 oppure 999300

!!to use the cna written by G. Rossi
!!i want out#.xyz
!where #=(nd-1)*npas/scrivo+nfile

        i=(nd-1)*npas/scrivo+nfile
               
         if((i.le.9).and.(i.ge.1))then
        filename(nd*nfile)='out0.xyz'
        write(filename(nd*nfile)(4:4),'(i1)') i 
         elseif((i.le.99).and.(i.ge.10)) then
        filename(nd*nfile)='out00.xyz'
        write(filename(nd*nfile)(4:5),'(i2)') i          
         elseif((i.le.999).and.(i.ge.100)) then                
         filename(nd*nfile)='out000.xyz'
        write(filename(nd*nfile)(4:6),'(i3)') i
         elseif((i.le.9999).and.(i.ge.1000)) then                
         filename(nd*nfile)='out0000.xyz'
        write(filename(nd*nfile)(4:7),'(i4)') i
         elseif((i.le.99999).and.(i.ge.10000)) then                
         filename(nd*nfile)='out00000.xyz'
        write(filename(nd*nfile)(4:8),'(i5)') i
         elseif((i.le.999999).and.(i.ge.100000)) then                
         filename(nd*nfile)='out000000.xyz'
        write(filename(nd*nfile)(4:8),'(i6)') i
        elseif(i>1000000) then
         write(*,*) 'increase number of output file'
         STOP
         endif

End Subroutine Subnomi
