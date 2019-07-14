subroutine writing_coal
    use paracluster
    use cluster
    use potential
    use enforce
    use distance
    use substrate
!    use output

    implicit none

    integer :: i, j, k, unitp, unitc
    integer :: icolor
    integer :: nbond, nbond2, nbond3
    real :: xuar, yvar, zwar, x2cdm, y2cdm, z2cdm, x1cdm, y1cdm, z1cdm, vzsum, vz2sum
    real :: emedia, distcom

    !unit numbers
    unitp = 50
    unitc = 70

    !energy.out contains the energy output
    !coalescing.out contains variables relating coalescence
    !movie.xyz/hdf5 contains positions

        open(unite,file = 'energy.out',status = 'unknown',position = 'append')
        open(unitm, file = 'movie.xyz', status = 'unknown', access = 'append')
        open(unitc, file = 'coalescing.out', status = 'unknown', access = 'append')
        if(ipas == 1) then
           if(mgo_substrate) then
                write(unite, *) '## time[ps]  epot[ev]  etot[ev]  ekin[ev]  <edelta>  <etot>  <t[k]>   esub[ev]'
                write(unitc, *) '## time[ps] comd[a] #bimbonds[a] #monbonds[a] #intclustbond driftv1 driftv2'
            else
                write(unite, *) '## time[ps]  epot[ev]  etot[ev]  ekin[ev]  <edelta>  <etot>  <t[k]>'
                write(unitc, *) '## time[ps] comd[a] #bimbonds[a] #monbonds[a] #intclustbond driftv1 driftv2'
            end if
        end if

    !average quantities
    tpar = tpar + temp
    edelta = edelta + etot
    !emedia = emedia + etot !modified by lp: emedia is not initialized and this value is not used.

    !====================================
    ! writing each scrivo-step
    !====================================
    if(mod(ipas, scrivo) == 0) then
        ! average energy
        emedia = edelta / scrivo !at this point edelta is the sum of etot over a growing number of time steps
        !now emedia is etot averaged over the number of scrivo steps
        !edelta calculation
        edelta = edelta / scrivo ! now edelta becomes etot averaged over the number of scrivo steps
        if(sys == 'mon') then
            !now edelta is really edelta averaged over scrivo
            edelta = (edelta + natom * ecoh(1)) / (natom**(2.d0 / 3.d0))
        else
            edelta = (edelta + ntipo1 * ecoh(1) + ntipo2 * ecoh(2)) / ((natom)**(2.d0 / 3.d0))
        end if

        !average temperature calculation
        tpar = tpar / (scrivo)
        if(mgo_substrate) then
            write(unite, '(1f15.4,1x,3f14.5,1x,4f14.5)') &
            tempo, ener, etot, ecin, edelta, emedia, tpar, ener_sub
            write(*, '(1f15.4,1x,3f14.5,1x,4f14.5)') & !writing on terminal
            tempo, ener, etot, ecin, edelta, emedia, tpar, ener_sub
        else
            write(unite, '(1f15.4,1x,3f14.5,1x,3f14.5)')&
            tempo, ener, etot, ecin, edelta, emedia, tpar
            write(*, '(1f15.4,1x,3f14.5,1x,3f14.5)') & !writing on terminal
            tempo, ener, etot, ecin, edelta, emedia, tpar
        end if

        emedia = 0.d0
        edelta = 0.d0
        tpar = 0.d0

        write(unitm, *) natom
        write(unitm, '(i8, 2(a3,i4), 1x, 1f14.5, 1f7.3, a4)') &
            ipas, elem1, ntipo1, elem2, ntipo2, etot, temp

        do i = 1, natom
            icolor = itype(i)
            xuar = (x(i) + u(i)) * aretebim * fattor
            yvar = (y(i) + v(i)) * aretebim * fattor
            zwar = (z(i) + w(i)) * aretebim * fattor
            !bufferpos(1, i, 1) = xuar
            !bufferpos(2, i, 1) = yvar
            !bufferpos(3, i, 1) = zwar
            !!if (output_xyz) then
                write(unitm,'(a3,1x,3f16.5,i4)') elem(i), xuar, yvar, zwar, icolor
            !!end if
        end do

        !call write_movie

        !=============start of coalescing.out===============!
        !checking hetero bonds
        nbond = 0
        nbond2 = 0
        nbond3 = 0

        !here we check how many nn of target cluster atoms belonged to bullet clusters
        do i = 1, natom - natom2
            do j = 1, nvois(i)
                k = ivois(j, i)
                if ( k > (natom - natom2) ) then
                    nbond = nbond + 1
                end if
            end do
        end do

        do i = 1, natom !monometallic and bimeticallic, counting twice
            do j = 1, nvois(i)
                k = ivois(j, i)
                if (elem(k) /= elem(i)) then
                    nbond2 = nbond2 + 1
                else if (elem(k) == elem(i)) then
                    nbond3 = nbond3 + 1
                end if
            end do
        end do

        nbond2 = nbond2 / 2
        nbond3 = nbond3 / 2

        !COM calculations
        x1cdm = (aretebim * fattor / float(natom - natom2)) * &
            sum(x(1:natom - natom2) + u(1:natom - natom2))
        y1cdm = (aretebim * fattor / float(natom - natom2)) * &
            sum(y(1:natom - natom2) + v(1:natom - natom2))
        z1cdm = (aretebim * fattor / float(natom - natom2)) * &
            sum(z(1:natom - natom2) + w(1:natom - natom2))
        x2cdm = (aretebim * fattor / float(natom2)) * &
            sum(x(natom - natom2 + 1:natom) + u(natom - natom2 + 1:natom))
        y2cdm = (aretebim * fattor / float(natom2)) * &
            sum(y(natom - natom2 + 1:natom) + v(natom - natom2 + 1:natom))
        z2cdm = (aretebim * fattor / float(natom2)) * &
            sum(z(natom - natom2 + 1:natom) + w(natom - natom2 + 1:natom))

        distcom = sqrt((x1cdm - x2cdm)**2 + (y1cdm - y2cdm)**2 + (z1cdm - z2cdm)**2)

        !Drift Velocity
        vzsum = (sum(vz(1:natom - natom2)) / float(natom - natom2)) * (aretebim * fattor * 1.0d-10 / tstep)
        vz2sum = (sum(vz(natom - natom2 + 1:natom)) / float(natom2)) * (aretebim * fattor * 1.0d-10 / tstep)

        write(unitc,*)tempo, distcom, nbond2, nbond3, nbond, vzsum, vz2sum
    end if

    if(ipas == npas) then
        close(unite)
        !!if(output_xyz)
        close(unitm)
    end if

end subroutine
