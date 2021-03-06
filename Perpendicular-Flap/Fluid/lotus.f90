program Lotus_preCICE
  use preciceMod
  use flexMod,    only: flexBody
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg
  use imageMod,   only: display
  use propMod
 implicit none
!
! -- Define parameter, declare variables
  real,parameter         :: D = 64, Re = 1500   ! diameter and Reynolds number
  integer                :: n(3) = D*(/4,2,2/)  ! numer of points
  integer                :: b(3) = (/1,1,1/)    ! MPI blocks (product must equal n_procs)
  logical                :: root                ! root processor
  type(fluid)            :: flow
  type(flexBody)         :: geom
  real                   :: force(3),time,U,t,a0=0.25
  type(precice_coupling) :: precice 
  integer                :: ndims=2,iter=0,every=32
!
! thin bodies
  eps=0.5
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(ndims,set_blocks=b)
  root = mympi_rank()==0
  rank=mympi_rank()
  commsize=product(b)
!
  if(root) print *,'       Setting up the grid, body and fluid      '
  if(root) print *,'------------------------------------------------'
  call xg(1)%stretch(n(1), -8.*D, -0.2*D, .8*D, 16.*D, prnt=root)
  call xg(2)%stretch(n(2), -0.5*D, -0.5*D, 0.5*D, 8*D, prnt=root)
  if(ndims>2)then
    call xg(3)%stretch(n(3), -0.0*D, -0.0*D, 0.6*D, 4*D, prnt=root)
  else
    n(3) = 1
  end if
!
! build geometry
  call geom%init(fname='../Solid/geom.inp')
  call flow%init(n/b,geom,V=(/0.,0.,0./),nu=D/Re,endStep=.false.) 
  call precice%init(geom,rank,commsize)
  dt_precice = precice%timestep()
!
! - start 
  if(root) print*, 'y+: ', sqrt(0.026/Re**(1./7.)/2.)/(D/Re)
  if(root) print *,'            Setting up the Coupling             '
  if(root) print *,'------------------------------------------------'
!
  flow%dt = 0.1
  call flow%write(geom,write_vti=.false.)
  call geom%writeFlex(flow,flow%time)
!
! -- writing initial data
  if(precice%is_action_required(writeInitialData)) then
    call flow%update(geom)
    forces = geom%getFacesForces(flow)
    call precice%write(forces)
    call precice%mark_action_fulfilled(writeInitialData)
  end if
  call precice%initialize_data
  
  ! read initialized displacements
  if(precice%is_read_data_available()) then
    call precice%read(displacements)
    call geom%updateFlex(displacements,flow%dt)
  end if
  
  if(root) print *,'            Starting time update loop           '
  if(root) print *,'------------------------------------------------'

  do while(precice%ongoing())

    ! set time-step, here we limit the CFL of the fluid
    dt_lotus = min(dt_precice,flow%dt)
    flow%dt = dt_lotus

    ! check if we need to read an iteration checkpoint
    if(precice%is_action_required(writeItCheckp)) then
      call flow%writeItCheckp()
      call geom%writeItCheckp()
      call precice%mark_action_fulfilled(writeItCheckp)
    end if

    ! read the actual displacements
    if(precice%is_read_data_available()) then
      call precice%read(displacements)
      call geom%updateFlex(displacements,dt_precice)
    end if

    U = velocity((flow%time+real(dt_lotus))/D) 
    call flow%update(geom,V=[U,0.,0.])

    ! get the forces on the structure
    if(precice%is_write_data_required(dt_lotus)) then
      forces = geom%getFacesForces(flow)
      ! forces = merge(forces*flow%time/D,forces,flow%time<D)
      call precice%write(forces)
    end if

    ! advance the coupling, this is used to communicate
    ! the Lotus time-step length to preCICE
    dt_precice = flow%dt
    call precice%advance(dt_precice)
    
    if(precice%is_action_required(readItCheckp)) then
      call flow%readItCheckp()
      call geom%readItCheckp()
      call precice%mark_action_fulfilled(readItCheckp)
    else
      ! if we do not need, we can finalize the simulations
      call flow%endTimeStep()
      iter = iter + 1
    end if
    
    if(precice%is_time_window_complete()) then
      force = geom%getForces(flow)/D**2
      write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt,force
      if(mod(iter,every).eq.0) call flow%write(geom,write_vti=.false.)
      if(mod(iter,every).eq.0) call geom%writeFlex(flow,flow%time)
    end if
  end do

  if(root) print *,'                      Done                      '
  if(root) print *,'------------------------------------------------'

  call precice%finalize
  call mympi_end

 contains
   
  function velocity(t) result(omega)
    implicit none
    real,intent(in) :: t
    real            :: omega
    ! omega = 0.5*(1-cos(pi*time))
    omega = 0.25*t+(1+tanh(31.4*(t-1./0.25)))/2.*(1-0.25*t)
  end function velocity

end program Lotus_preCICE
