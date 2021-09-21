program Lotus
  use preciceMod
  use fieldMod,   only: field
  use flexMod,    only: flexBody
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use imageMod,   only: display
  use geom_global,only: pi
  use gridMod
  use geom_shape
 implicit none
!
! -- Define parameter, declare variables
  real,parameter         :: D = 64,   Re = 1490          ! diameter and Reynolds number
  integer                :: n(3)                         ! numer of points
  integer                :: b(3) = (/2,1,1/)             ! MPI blocks (product must equal n_procs)
  logical                :: root                         ! root processor
  real                   :: dt_precice,dt_lotus          ! time step
  real,allocatable       :: displacements(:),forces(:)   ! arrays for transfer
  integer                :: rank,commsize,ndims=2        ! mpi utils
  integer                :: iter=0,every=4
  type(fluid)            :: flow                         ! viscous flow
  type(flexBody)         :: geom                         ! flexible geometry
  type(precice_coupling) :: precice                      ! coupling
  real                   :: force(3)                     ! global forces on the flag
  !
  eps2 =0.5
  !
  ! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(ndims,set_blocks=b)
  root = mympi_rank()==0
  rank=mympi_rank()
  commsize=product(b)
!
  if(root) print *,'       Setting up the grid, body and fluid      '
  if(root) print *,'------------------------------------------------'

  n = composite(D/2*(/8,3,3/),prnt=root)
  call xg(1)%stretch(n(1), -5.*D, -0.2*D, .8*D, 12.*D, prnt=root)
  call xg(2)%stretch(n(2), 0., 0., 1.2*D, 4*D, prnt=root)
  if(ndims>2)then
    call xg(3)%stretch(n(3), -0.0*D, -0.0*D, 0.6*D, 4*D, prnt=root)
  else
    n(3) = 1
  end if

! build geometry, no scaling can be applied!
  call geom%init(fname='../Solid/geom.inp')
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=D/Re)
  call precice%init(geom,rank,commsize)
  dt_precice = precice%timestep()
  if(root) print*, 'y+: ', sqrt(0.026/Re**(1./7.)/2.)/(D/Re)

  ! allocate array for communication
  allocate(displacements(geom%numberOfVertices*precice%dimensions))
  allocate(forces(geom%numberOfFacets*precice%dimensions))

  if(root) print *,'            Starting time update loop           '
  if(root) print *,'------------------------------------------------'

  ! writing initial data
  if(precice%is_action_required(writeInitialData)) then
    call flow%writeItCheckp
    call geom%writeItCheckp
    call precice%mark_action_fulfilled(writeInitialData)
  end if
  call precice%initialize_data

  ! solve  
  do while(precice%ongoing())

    ! check if we need to read an iteration checkpoint
    if(precice%is_action_required(writeItCheckp)) then
      call flow%writeItCheckp()
      call geom%writeItCheckp()
      call precice%mark_action_fulfilled(writeItCheckp)
    end if
    
    ! read the displacements if they have been written
    if(precice%is_read_data_available()) then
      call precice%read(displacements)
      call geom%updateFlex(displacements,dt_precice)
    end if
    
    ! set timestep, here we limit the CFL of the fluid
    dt_lotus = min(dt_precice,flow%dt)
    flow%dt = dt_lotus

    ! update flow
    call flow%update(geom)
    force = 2*geom%getForces(flow)/(0.1*D**2)
    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt,force
    
    ! get the fores for structure if required
    if(precice%is_write_data_required(dt_lotus)) then
      forces = geom%getFacesForces(flow)
      call precice%write(forces)
    end if

    ! advance the coupling, this is used to communicate
    ! the Lotus time-step length to preCICE
    dt_precice = dt_lotus
    call precice%advance(dt_precice)

    ! check if we need to read an iteration checkpoint
    if(precice%is_action_required(readItCheckp)) then
      call flow%readItCheckp()
      call geom%readItCheckp()
      call precice%mark_action_fulfilled(readItCheckp)
    else
      !if we do not need, we can finalize the simulations
      call flow%endTimeStep()
      iter = iter + 1
    end if

    ! write at the end of every-n time window
    if(precice%is_time_window_complete()) then
      call geom%writeFlex(flow,flow%time)
      call flow%write(geom,lambda=.true.,write_vti=.false.)
      call display(flow%velocity%vorticity_Z(),'out_vort',lim=0.5)
    end if

  end do

  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'------------------------------------------------'

  ! finish
  call precice%finalize
  call mympi_end

end program Lotus