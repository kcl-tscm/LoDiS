subroutine vel
    use paracluster !using the following modules
    use cluster
    use potential
    use enforce
    use distance

    implicit none

    !local variables
    real(8) :: dlaran
    real(8) :: qx(nsiz), qy(nsiz), qz(nsiz)
    real(8) :: qxx, qyy, qzz
    real(8) :: pcmx, pcmy, pcmz, pcmz2, vxsum, vysum
    real(8) :: versx, versy, versz, versmod, velmod, vel2
    real(8) :: ecinet, poguy, vargau
    real(8) :: vdfa
    real(8) :: vdfae
    integer :: nat3
    integer :: ir(4)
    integer :: i, k, j
    logical :: clustercheck

    !initialization
    vdfa = vnu * tstep
    if(ipas == 1) ir(:) = irand_seed(:)
    ecinet = 0.d0
    temp  = 0.d0
    poguy = 0.d0
    pcmx  = 0.d0
    pcmy  = 0.d0
    pcmz  = 0.d0
    if (coalescence == 'ya') pcmz2 = 0.d0 !for cluster two
    nat3 = natom

    !for deposition type process
    if(deposizione == 'ya') then
        nat3 = natom - 1
        if(nvois(natom) == 0) then
            pcmx = sum(x(1:nat3) + u(1:nat3)) / dble(nat3) !sum using splicing
            pcmy = sum(y(1:nat3) + v(1:nat3)) / dble(nat3)
            pcmz = sum(z(1:nat3) + w(1:nat3)) / dble(nat3)
            versx = pcmx - x(natom) - u(natom)
            versy = pcmy - y(natom) - v(natom)
            versz = pcmz - z(natom) - w(natom)
            versmod = sqrt(versx*versx + versy*versy + versz*versz)
            versx = versx / versmod
            versy = versy / versmod
            versz = versz / versmod
            vel2 = vx(natom)*vx(natom) + vy(natom)*vy(natom) + vz(natom)*vz(natom)
            velmod = sqrt(vel2)
            vx(natom) = velmod * versx
            vy(natom) = velmod * versy
            vz(natom) = velmod * versz
        else if(nvois(natom) /= 0) then
            vx(natom) = (vx(natom) + t2m(itype(natom)) * (fx(natom) + dfx(natom)))
            vy(natom) = (vy(natom) + t2m(itype(natom)) * (fy(natom) + dfy(natom)))
            vz(natom) = (vz(natom) + t2m(itype(natom)) * (fz(natom) + dfz(natom)))
        end if
    end if

    if (coalescence == 'ya') then
        !terminate once near
        clustercheck = .false.
        outer:  do i = 1, natom - natom2  !loop over atoms of cluster two to see if they interact with any atoms in cluster one
            inner:  do j = 1, nvois(i)
                k = ivois(j, i)
                if ( k > (natom - natom2) ) then
                    clustercheck = .true.
                    exit outer
                end if
            end do inner
        end do outer

        if (.not. clustercheck) then !check if no atom in cluster one is interacting with cluster two
            pcmz = sum(z(1:natom - natom2) + w(1:natom - natom2)) / dble(natom - natom2)
            pcmz2 = sum(z(natom - natom2 + 1:natom) + w(natom - natom2 + 1:natom)) / dble(natom2)
            versz = pcmz2 - pcmz
            versmod = sqrt(versz * versz)
            versz = versz / versmod
            vxsum = sum(vx(1:natom - natom2)) / dble(natom - natom2) !supressing drift velocity
            vysum = sum(vy(1:natom - natom2)) / dble(natom - natom2) !in x-y directions for
            vx(1:natom - natom2) = vx(1:natom - natom2) - vxsum      !both cluster one and two
            vy(1:natom - natom2) = vy(1:natom - natom2) - vysum
            vxsum = sum(vx(natom - natom2 + 1:natom)) / dble(natom2)
            vysum = sum(vy(natom - natom2 + 1:natom)) / dble(natom2)
            vx(natom - natom2 + 1:natom) = vx(natom - natom2 + 1:natom) - vxsum
            vy(natom - natom2 + 1:natom) = vy(natom - natom2 + 1:natom) - vysum

            do i = 1, natom - natom2
                velmod = sqrt(vz(i) * vz(i))
                vz(i) = velmod * versz
            end do
            versz = -versz
            do i = natom - natom2 + 1, natom
                velmod = sqrt(vz(i) * vz(i))

                vz(i) = velmod * versz
            end do
            
        end if
    end if

    do i = 1, nat3
        qx(i) = (vx(i) + t2m(itype(i)) * (fx(i) + dfx(i)))
        qy(i) = (vy(i) + t2m(itype(i)) * (fy(i) + dfy(i)))
        qz(i) = (vz(i) + t2m(itype(i)) * (fz(i) + dfz(i)))
        if(quenching =='no') then
            !thermalization of atoms with neighbours with freq vdfa
            vdfae = dlaran(irand_seed)
            if (vdfae > vdfa) then
                vx(i) = qx(i)
                vy(i) = qy(i)
                vz(i) = qz(i)
            else
                !call to the andersen thermostat
                vargau = sqrt(2.d0 * t2m(itype(i)) * cbol * tfin)
                call gauss
                vx(i) = g1 * vargau
                vy(i) = g2 * vargau
                vz(i) = g3 * vargau
            end if  !decision of thermostat call
        else if(quenching == 'ya') then
            vx(i) = qx(i)
            vy(i) = qy(i)
            vz(i) = qz(i)
            if(ipas > itremp) then
                qxx = qx(i) * fx(i)
                qyy = qy(i) * fy(i)
                qzz = qz(i) * fz(i)
                if(qxx < 0.d0) vx(i) = 0.d0
                if(qyy < 0.d0) vy(i) = 0.d0
                if(qzz < 0.d0) vz(i) = 0.d0
            end if
        end if

        !on activation of quenching procedure/thermostat
        !equipartition theorem to calculate the temperature of the system
        if(canonical == 'no') then
            vel2 = vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
            poguy = vel2 / (4.d0 * t2m(itype(i)))
            ecinet = ecinet + poguy !is calculted in ha e^2/bohr(just check) in ev
        end if
    end do   !end i'th atom

    if(canonical == 'ya') then !for vacf
        open(unit = 420, file = "velocity.out", status = "unknown", access = "append")
        call suppress
        do i = 1, natom
            vel2 = vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
            poguy = vel2 / (4.d0 * t2m(itype(i)))
            ecinet = ecinet + poguy
            write(420, *)vx(i), ",", vy(i), ",", vz(i)
        end do
        close(420)
    end if

    ecin = ecinet
    etot = ener + ecin  !ener = potential energy; etot = total energy
    temp = 2.d0 * ecin / (3.d0 * nat3 * cbol) !equipartition theorem

end subroutine
