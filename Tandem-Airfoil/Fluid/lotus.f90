program Lotus_preCICE
  use preciceMod
  use flexMod,    only: flexBody
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg
  use imageMod,   only: display
  ! use geom_global
  use propMod
 implicit none
!
! -- Define parameter, declare variables
  real,parameter         :: D = 64, Re = 1000    ! diameter and Reynolds number
  integer                :: n(3) = D*(/4,3,2/)  ! numer of points
  integer                :: b(3) = (/1,1,1/)    ! MPI blocks (product must equal n_procs)
  logical                :: root                ! root processor
  type(fluid)            :: flow                ! a fluid
  type(flexBody)         :: geom                ! a flexible geometry
  real                   :: force(3),time       ! helpers
  type(precice_coupling) :: precice             ! the coupling
  integer                :: ndims=2,iter=0,every=40
  real,allocatable       :: vertex(:)
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
  call xg(1)%stretch(n(1), -8.*D, -0.5*D, 2.0*D, 16.*D, prnt=root)
  call xg(2)%stretch(n(2), -8*D, -0.8*D, 0.8*D, 8*D, prnt=root)
  if(ndims>2)then
    call xg(3)%stretch(n(3), -0.0*D, -0.0*D, 0.6*D, 4*D, prnt=root)
  else
    n(3) = 1
  end if

  ! build geometry
  call geom%init(fname='../Solid/geom.inp')
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=D/Re,endStep=.false.) 
  call precice%init(geom,rank,commsize)
  call geom%getVertexCoords(vertex)
  dt_precice = precice%timestep()
!
! initalize finishes
  if(root) print*, 'y+: ', sqrt(0.026/Re**(1./7.)/2.)/(D/Re)
  if(root) print *,'            Setting up the Coupling             '
  if(root) print *,'------------------------------------------------'
!
! write initial fields
  flow%dt = 0.25
  call flow%write(geom,write_vti=.false.)
  call geom%writeFlex(flow,flow%time)
!
! writing initial data
  if(precice%is_action_required(writeInitialData)) then
    call flow%update(geom)
    forces = geom%getFacesForces(flow)
    call precice%write(forces)
    call precice%mark_action_fulfilled(writeInitialData)
  end if
  call precice%initialize_data
!
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
      call motion(vertex,displacements,flow%time)
      call geom%updateFlex(displacements,dt_precice)
    end if

    call flow%update(geom)

    ! get the forces on the structure
    if(precice%is_write_data_required(dt_lotus)) then
      forces = geom%getFacesForces(flow)
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

  subroutine motion(vertex,displ,time)
    real,intent(inout) :: vertex(:),displ(:)
    real,intent(in) :: time
    real :: r,x(3)
    integer :: i,k,len
    ! this set what part of the array we move
    len = 128
    do i=1,len
      k = 3*(i-1)
      x = vertex(k+1:k+3)
      displ(k+1:k+3) = rotate(x,[0.,0.,0.],alpha(time),3) - x
    end do
  end subroutine motion

  real function alpha(time)
    real,intent(in)    :: time
    real :: C
    C = 1 - exp(-time/D/8)
    alpha = C*pi/8*sin(0.25*pi*time/D)
  end function

  function rotate(a,cen,alpha,axis) result(b)
    real,intent(in)    :: a(3),cen(3),alpha
    integer,intent(in) :: axis
    integer            :: j,k
    real               :: cs,ss,b(3),x(3)
    x = a-cen; b = a
    j=mod(axis,3)+1; k=mod(j,3)+1
    b(j) = cos(alpha)*x(j)-sin(alpha)*x(k)
    b(k) = cos(alpha)*x(k)+sin(alpha)*x(j)
    b = b + cen
  end function
  
end program Lotus_preCICE
