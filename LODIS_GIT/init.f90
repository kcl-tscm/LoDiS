subroutine init
    use paracluster  !use module where parameters are defined
    use cluster
    use potential
    use enforce
    use distance

    implicit none

    real :: dlaran !function

    !local variables
    integer :: is
    integer :: ir(4)
    real :: facguy
    real :: vxsum, vysum, vzsum
    real :: vxsum2, vysum2, vzsum2
    real :: vxsumcoal, vysumcoal, vzsumcoal !cluster two for coalescense
    real :: vxsum2coal, vysum2coal, vzsum2coal !cluster two for coalescence

    ir(1:4) = irand_seed(1:4)

    !initialization of displacements(old and new) and forces(old and new)
    u(1:natom) = 0.d0
    v(1:natom) = 0.d0
    w(1:natom) = 0.d0
    du(1:natom) = 0.d0
    dv(1:natom) = 0.d0
    dw(1:natom) = 0.d0
    fx(1:natom) = 0.d0
    fy(1:natom) = 0.d0
    fz(1:natom) = 0.d0
    dfx(1:natom) = 0.d0
    dfy(1:natom) = 0.d0
    dfz(1:natom) = 0.d0

    !initial velocities
    do is = 1, natom
        vx(is) = dble(2.d0 * dlaran(ir) - 1.d0)
        vy(is) = dble(2.d0 * dlaran(ir) - 1.d0)
        vz(is) = dble(2.d0 * dlaran(ir) - 1.d0)
    enddo

    !drift velocity of the com
    if (coalescence == 'ya') then
        vxsum = sum(vx(1:natom - natom2)) / float(natom - natom2)
        vysum = sum(vy(1:natom - natom2)) / float(natom - natom2)
        vzsum = sum(vz(1:natom - natom2)) / float(natom - natom2)
        vxsumcoal = sum(vx(natom - natom2 + 1:natom)) / float(natom2)
        vysumcoal = sum(vy(natom - natom2 + 1:natom)) / float(natom2)
        vzsumcoal = sum(vz(natom - natom2 + 1:natom)) / float(natom2)
    else
        vxsum = sum(vx(1:natom)) / float(natom)
        vysum = sum(vy(1:natom)) / float(natom)
        vzsum = sum(vz(1:natom)) / float(natom)
    end if

    write(*,*)vxsum, vysum, vzsum, 'initial drift of cluster (one in the case of coalescence)'
    if (coalescence == 'ya') then
        write(*,*)vxsumcoal, vysumcoal, vzsumcoal, 'drift of cluster two in case of coalescence'
    end if

    !suppressing the drift
    if (coalescence == 'ya') then
        !using splicing opposed to do loops
        vx(1:natom - natom2) = vx(1:natom - natom2) - vxsum
        vy(1:natom - natom2) = vy(1:natom - natom2) - vysum
        vz(1:natom - natom2) = vz(1:natom - natom2) - vzsum
        vx(natom - natom2 + 1:natom) = vx(natom - natom2 + 1:natom) - vxsumcoal
        vy(natom - natom2 + 1:natom) = vy(natom - natom2 + 1:natom) - vysumcoal
        vz(natom - natom2 + 1:natom) = vz(natom - natom2 + 1:natom) - vzsumcoal
    else
        vx(1:natom) = vx(1:natom) - vxsum
        vy(1:natom) = vy(1:natom) - vysum
        vz(1:natom) = vz(1:natom) - vzsum
    end if

    !finding scaling terms
    if (coalescence == 'ya') then
        vxsum2 = sum(vx(1:natom - natom2) * vx(1:natom - natom2)) / float(natom - natom2)
        vysum2 = sum(vy(1:natom - natom2) * vy(1:natom - natom2)) / float(natom - natom2)
        vzsum2 = sum(vz(1:natom - natom2) * vz(1:natom - natom2)) / float(natom - natom2)
        vxsum2coal = sum(vx(natom - natom2 + 1:natom) * vx(natom - natom2 + 1:natom)) / float(natom2)
        vysum2coal = sum(vy(natom - natom2 + 1:natom) * vy(natom - natom2 + 1:natom)) / float(natom2)
        vzsum2coal = sum(vz(natom - natom2 + 1:natom) * vz(natom - natom2 + 1:natom)) / float(natom2)
    else
        vxsum2 = sum(vx(1:natom) * vx(1:natom)) / float(natom)
        vysum2 = sum(vy(1:natom) * vy(1:natom)) / float(natom)
        vzsum2 = sum(vz(1:natom) * vz(1:natom)) / float(natom)
    end if

    write(*,*)vxsum2, vysum2, vzsum2,'scaling for velocity (for cluster one in case of coalescence)'
    if (coalescence == 'ya') then
        write(*,*)vxsum2coal, vysum2coal, vzsum2coal, 'scaling of cluster two in case of coalescence'
    end if

    if (coalescence == 'ya') then
        do is = 1, natom - natom2
            facguy = sqrt(2.d0 * t2m(itype(is)) * cbol * tinit)
            vx(is) = vx(is) * facguy * sqrt(1.d0 / (3.d0 * vxsum2))
            vy(is) = vy(is) * facguy * sqrt(1.d0 / (3.d0 * vysum2))
            vz(is) = vz(is) * facguy * sqrt(1.d0 / (3.d0 * vzsum2))
        enddo
        do is = natom - natom2 + 1, natom
            facguy = sqrt(2.d0 * t2m(itype(is)) * cbol * tinit)
            vx(is) = vx(is) * facguy * sqrt(1.d0 / (3.d0 * vxsum2coal))
            vy(is) = vy(is) * facguy * sqrt(1.d0 / (3.d0 * vysum2coal))
            vz(is) = vz(is) * facguy * sqrt(1.d0 / (3.d0 * vzsum2coal))
        end do
    else
        do is = 1, natom
            facguy = sqrt(2.d0 * t2m(itype(is)) * cbol * tinit)
            vx(is) = vx(is) * facguy * sqrt(1.d0 / (3.d0 * vxsum2))
            vy(is) = vy(is) * facguy * sqrt(1.d0 / (3.d0 * vysum2))
            vz(is) = vz(is) * facguy * sqrt(1.d0 / (3.d0 * vzsum2))
        end do
    end if

    !calclating drift again
    if (coalescence == 'ya') then
        vxsum = sum(vx(1:natom - natom2)) / float(natom - natom2)
        vysum = sum(vy(1:natom - natom2)) / float(natom - natom2)
        vzsum = sum(vz(1:natom - natom2)) / float(natom - natom2)
        vxsumcoal = sum(vx(natom - natom2 + 1:natom)) / float(natom2)
        vysumcoal = sum(vy(natom - natom2 + 1:natom)) / float(natom2)
        vzsumcoal = sum(vz(natom - natom2 + 1:natom)) / float(natom2)
    else
        vxsum = sum(vx(1:natom)) / float(natom)
        vysum = sum(vy(1:natom)) / float(natom)
        vzsum = sum(vz(1:natom)) / float(natom)
    end if

    write(*,*)'init.f90 > facguy = ', facguy
end subroutine
