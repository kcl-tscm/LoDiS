module cluster
    use paracluster

    implicit none

    !system variables
    integer :: irand, sticky_atoms
    integer :: irand_seed(4) !random number seed
    integer :: npas, npassifreezing
    integer :: ipas
    integer :: nome
    integer :: scrivo, error_counter
    integer :: npast !time steps for thermalization
    real(8) :: tempo
    character(len = 59) :: filepos
    character(len = 59) :: filepot
    character(len = 21) :: type_process
    logical :: output_xyz
    logical :: mgo_substrate, metal_on_top, writeheader, cn_cal
    character(len = 59) :: mgo_pot
    logical :: impl_env
    real :: pot_a, pot_b, eta_a, eta_b    

    !for restart
    logical, save :: restart_mode, save_progress
    character(len = 60) restart_from_file
    integer, save :: start_from
    integer, save :: start_from_nd_proc
    character(len = 60) :: progress_file_label
    logical, save :: periodic_save
    integer, save :: periodic_save_every
    integer, save :: periodic_save_last
    integer, save :: periodic_i = 0
    logical, save :: starting_from_midpoint = .false. !this flag indicates if
                                                      !program restarts from a
                                                      !crashed simulation

    !Note: there are some suggested values for some variables
    real(8) :: tstep !=5.0d-15    !time step
    real(8) :: vnu   !=5.0d11     !frequency of the thermostat
    real(8) :: tinit !initial temperature

    !procedures
    character(len = 2) :: quenching
    integer :: itremp
    real(8)    :: tmin !=1.d0-6      !minimum temperature for the quenching

    character(len = 2) :: caloric
    real(8)  :: tcaloric !end time for caloric temperature
    real(8)  :: deltat  !temperature variation

    character(len = 2) :: canonical
    logical :: vel_af
    real(8)  :: tfin    !temperatura running

    character(len = 2) ::coalescence
    real(8) :: somedist !arbitrary distance to ensure atoms don't start clustered
    integer :: natom2
    character(len = 59) :: filepos2

    character(len = 2) :: deposizione
    integer :: natinizio
    integer :: ntipo1, ntipo2
    integer :: ndepmax !number of atoms will be deposited
    integer :: lcs     !biatomic growth = 2; monoatomic growth = 1, c/s = 3
    integer :: at_tipo2 !lcs = 2 with prob deposition of elem2 uo to at_tipo2
    real(8) :: rad     !=5.d0        !distance for deposition
    real(8) :: prob    !only for biatomic growth, to choose the atom type
    real(8) :: tsorg !=1500.d0   !temperature source

    character(len = 2) :: metadyn !see meta_module for more metadynamic variables

    !more system varuables
    character(len = 3) :: type_potential !!stringa per il potenziale di interazione
    character(len = 2) :: elem1          !!stringa per sapere che materiale
    character(len = 2) :: elem2          !!stringa per sapere che materiale
    character(len = 2) :: elem(nsiz)
    integer :: imet1         !!scelta del materiale
    integer :: imet2         !!scelta del materiale
    integer :: itype(nsiz)
    integer :: natom, initntipo2         !!number of atoms, initial number of atoms of species 2 (for growth lcs=2) 

    !periodicity
    logical :: clusters, wires, surface, bulk, substrate
    real(8) :: dilat  !thermal dilatation at finte temperature of the lattice parameters
    real(8) :: pbcx, pbcy, pbcz !only if you have a periodic boundary condition
                                !on z for wire-direction for wires
                                !on x, y for surface
                                !on x, y, z for bulk
    integer :: nd_proc  !it is related to the type of process chosen
                        !number of deposited atoms for deposition
                        !number of times when t is changed
                        !nd = 1

    !conversion variables
    real(8) :: fraz_ecoh, fraz_mass, fraz_radius
    real(8) :: aretebim, fattor
    real(8) :: g1, g2, g3

end module cluster
