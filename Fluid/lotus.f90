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
  real,parameter         :: D = 64,   Re = 1600         ! diameter and Reynolds number
  integer                :: n(3)                        ! numer of points
  integer                :: b(3) = (/1,1,1/)            ! MPI blocks (product must equal n_procs)
  logical                :: root                        ! root processor
  real                   :: dt_precice,dt_lotus         ! time step
  real,allocatable       :: displacements(:),forces(:)  ! arrays for transfer
  integer                :: rank,commsize,k=0           ! mpi utils
  type(fluid)            :: flow                        ! viscous flow
  type(flexBody)         :: geom                        ! flexible geometry
  type(precice_coupling) :: precice                     ! coupling
  real                   :: force(3),U,a0=0.5,t         ! global forces on the flag
  type(field)            :: f
  real,pointer           :: p(:,:,:)

  eps2 =0.5
  !
  ! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(3,set_blocks=b)
  root = mympi_rank()==0
  rank=mympi_rank()
  commsize=product(b)
!
  if(root) print *,'       Setting up the grid, body and fluid      '
  if(root) print *,'------------------------------------------------'

  n = composite(D*(/3,2,1/),prnt=root)
  call xg(1)%stretch(n(1), -4.5*D, -0.2*D, .8*D, 7.5*D, prnt=root)
  call xg(2)%stretch(n(2), 0., 0., 1.2*D, 5*D, prnt=root)
  n(3)=1
  
! build geometry, no scaling can be appplied!
  call geom%init(fname='../Solid/geom.inp',thick=D/25)
  call flow%init(n/b,geom,V=(/0.,0.,0./),nu=D/Re)
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
   if(root) write(*,*) 'Writing initial data'
  end if

  do while(precice%ongoing())

    ! conversion between dble and real
    dt_lotus =min(0.5, flow%dt)

    ! check if we need to read an iteration checkpoint
    if(precice%is_action_required(writeItCheckp)) then
      if(root) write(*,*) 'Writing iteration checkpoint'
      call flow%writeItCheckp
      call geom%writeItCheckp
      call precice%mark_action_fulfilled(writeItCheckp)
    end if

    ! read the displacements if they have been written
    if(precice%is_read_data_available()) then
      call precice%read(displacements)
      call geom%updateFlex(displacements,dt_precice)
    end if

    ! set time step
    dt_lotus = min(dt_precice, dt_lotus)
    flow%dt = dt_lotus

    ! update flow
    U = velocity((flow%time+flow%dt)/D,0.5)
    call flow%update(geom,V=(/U,0.,0./))
    force = 2*geom%getForces(flow)/(0.2*D**2)
    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt,force
    
    ! get the fores for structure if required
    if(precice%is_write_data_required(dt_lotus)) then
      forces = geom%getFacesForces(flow)
      call precice%write(forces)
      if(mod(k,5)==0) then
        call geom%writeFlex(flow,flow%time)
        call flow%write(geom,lambda=.true.,write_vti=.false.)
      end if
      if(root) print '(f6.1,",",f6.3,",",f10.6,",",f10.6)',flow%time/D,flow%dt,U,force(1)
      f = flow%velocity%vorticity_Z()
      p => f%point()
      call display(p(:,:,:),'out_vort',lim=0.2)
      k=k+1
    end if

    ! advance the coupling
    dt_precice = dt_lotus
    call precice%advance(dt_precice)

    ! check if we need to read an iteration checkpoint
    if(precice%is_action_required(readItCheckp)) then
      if(root) write(*,*) 'Reading iteration checkpoint'
      call flow%readItCheckp
      call geom%readItCheckp
      call precice%mark_action_fulfilled(readItCheckp)
    end if
  end do

  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'------------------------------------------------'

  ! finish
  call geom%writeFlex(flow,flow%time)
  call flow%write(geom,lambda=.true.,write_vti=.false.) 
  call precice%finalize
  call mympi_end
  if(root) write (*,*) 'Lotus: Closing Fortran solver...'

 contains

  real function velocity(time,a0)
    real,intent(in) :: time,a0
    velocity = a0*time+(1+tanh(31.4*(time-1./a0)))/2.*(1-a0*time)
  end function velocity
end program Lotus