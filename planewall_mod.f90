!=========================================================================
!
! Title         : planewall_mod.f90
! Application   : Planewall
! Copyright     : Tata Consultancy Services
! Author        : Chetan Malhotra
! Creation Date : 06/05/99
! Requirements  : None
! Description   : This file contains the module 'planewall_mod' which 
!				  calculates the transient/steady-state temperature 
!				  distribution in a composite plane wall with variable 
!				  properties.
! Limitations   : No heat generation within the wall
! Dependencies  : Fortran 90
! Modifications :       WHO           WHEN         WHY
!				  Chetan Malhotra   23-2-2000   Added and tested variable
!												properties and formatted
!												as per TCS-QMS
!=========================================================================

module planewall_mod

contains

!=========================================================================
!
! Name              : planewall
! Description       : This subroutine calculates the transient/
!					  steady-state temperature distribution in a composite 
!					  plane wall with variable properties but without heat 
!					  generation. Constant temperature or (heat input and/
!					  or convection) boundary conditions may be used at 
!					  inner or outer surfaces 
! Parameters        : [I] N_layers
!					  [I] N_nodes_layer
!					  [I] N_nodes
!					  [I] Therm_cond
!					  [I] thick
!					  [I] Area
!					  [I] ibc_flag
!					  [I] h_in
!					  [I] T_in
!					  [I] heat_inner
!					  [I] obc_flag
!					  [I] h_out
!					  [I] T_out
!					  [I] heat_outer
!					  [I] Transient_flag
!					  [I] delta_t
!					  [I] Temp_old
!					  [I] rho
!					  [I] Cp
!					  [O] X
!                     [O] Temp
! Limitations       : No heat generation within the wall
! Globals updated   : None                     
! Files Read		: None
! Files updated     : None
! Calls             : subroutine tridiag
!
!=========================================================================

subroutine planewall(N_layers,N_nodes_layer,N_nodes,Therm_cond,thick, &
                     Area,ibc_flag,h_in,T_in,heat_inner,obc_flag, &
					 h_out,T_out,heat_outer,Transient_flag,delta_t, &
					 Temp_old,rho,Cp,X,Temp)

! =========
! Externals
! =========

!use mo_interpolate_mod
!use InterpolateLinear_mod

! =====================
! Internal Declarations
! =====================

implicit none

! Parameters

integer, parameter :: max_iter=50	! Max number of iterations
real*8, parameter :: eps_Temp=1.e-4 ! Convergence criterion for
									! temperature

! Passed variables

integer N_layers, &                 ! Number of layers in wall
		N_nodes, &					! Number of temperature nodes 
									! = Total no. of internal nodes + 
									!   no. of layers + 1 
		ibc_flag, &                 ! Flag=1/2 indicating constant 
									! temperature/convection boundary 
									! condition at inner surface of wall
		obc_flag, &                 ! Flag=1/2 indicating constant 
									! temperature/convection boundary 
									! condition at outer surface of wall
		Transient_flag	, &			! Flag=1/2 indicating Steady/Transient
									! run
	    preq,order_flag
integer, dimension(N_layers)  :: &
    N_nodes_layer                   ! Number of internal nodes in each 
									! layer (minimum=1, for each layer)
real*8 h_in, &						! Heat transfer coefficient at inner 
									! surface
									! Not required if ibc_flag=1 
									! Produces cumulative effect with heat 
									! addition
									! Enter=0 if only pure heat addition 
									! is required
       T_in, &                      ! Temperature at inner surface
									! = temperature of gas if ibc_flag=2
									! Not used for pure heat addition 
									! problems
	   heat_inner, &				! Heat added in addition to convection
									! at inner surface 
									! Not required if ibc_flag=1
	   h_out, &						! Heat transfer coefficient at outer
									! surface
									! Not required if obc_flag=1 
									! Produces cumulative effect with heat 
									! addition
									! Enter=0 if only pure heat addition 
									! is required
       T_out, &                     ! Temperature at outer surface
									! = temperature of gas if obc_flag=2
									! Not used for pure heat addition 
									! problems
	   heat_outer, &				! Heat added in addition to convection
									! at outer surface 
									! Not required if obc_flag=1
	   delta_t, &						! Time stepsize
									! Not required for steady-state run
	
	Therm_cond, &
	rho, Cp
real*8, dimension(N_layers) :: &
    thick, &                        ! Thickness of each layer
    Area							! Area of each layer
real*8, dimension(N_Nodes) :: &
	X, &							! Array containing location of nodes
	Temp_old, &						! Temperatures at pervious time step
	Temp							! Temperatures at steady-state/
									! current time step

real*8, dimension(N_Nodes-1) :: &
    X_Temp, &                       ! Array containing location of nodes at interfaces
    HF_Temp                         ! Temperatures at steady-state/
                                    ! current time step

real*8 xreqd(2),qreqd(2)


! Subroutine variables

integer i, &						! Node counter
		j, &						! Layer counter
		k, &						! Node counter within layer
		niter, &					! Iteration counter
		N_nodes_test, &				! Variable to test whether number of 
									! nodes passed is consistent with the
									! number of nodes calculated by adding
									! internal and interfacial nodes
		layer_no, &					! Dynamically calculated layer number
		layer_no_cont_vol_interface, & ! Layer number at control volume
									   ! interface
        layer_left, &				! Layer to the left of a node
		layer_right					! Layer to the right of a node
integer, dimension(N_Nodes) :: &
	point_map						! Array which maps nodes to layers
									! A 0 value indicates that the node
									! lies at the interface of two layers
real*8 conductivity, &				! Dynamically calculated conductivity
	   density, &					! Dynamically calculated density
	   specificheat, &				! Dynamically calculated specific heat
									! capacity
	   conductivity_left, &			! Conductivity to the left of a node
	   conductivity_right, &		! Conductivity to the right of a node
	   density_left, &				! Density to the left of a node
	   density_right, &				! Density to the right of a node
	   specificheat_left, &			! Specific heat capacity to the left
									! of the node
	   specificheat_right, &		! Specific heat capacity to the right
									! of the node
	   alpha_left, &				! Value of alpha to the left of the 
									! node (for definition of alpha, refer
									! to code manual)
	   alpha_right, &				! Value of alpha to the right of the 
									! node (for definition of alpha, refer
									! to code manual)
	   beta_left, &					! Value of beta to the left of the 
									! node (for definition of beta, refer
									! to code manual)
	   beta_right, &				! Value of beta to the right of the 
									! node (for definition of beta, refer
									! to code manual)
	   alpha_1, &					! Value of alpha for the innermost
									! node (for definition of alpha, refer
									! to code manual)
	   beta_1, &					! Value of beta for the innermost
									! node (for definition of beta, refer
									! to code manual)
	   alpha_N_nodes, &				! Value of alpha for the outermost
									! node (for definition of alpha, refer
									! to code manual)
	   beta_N_nodes, &				! Value of beta for the outermost
									! node (for definition of beta, refer
									! to code manual)
	   eps, &						! Difference in temperature values
									! after each iteration
	   epsmax						! Max value of eps among all nodes for
									! the iteration
real*8, dimension(N_layers) :: &
	dx								! Array containing the thickness of 
									! control volumes in every layer
real*8, dimension(N_Nodes-1) :: &
	Tbar, &							! Array containing the temperatures at
									! control volume interfaces
	alpha, &						! For definition of alpha, refer to
									! code manual
	beta							! For definition of beta, refer to
									! code manual
real*8, dimension(N_Nodes) :: &
	TDA, &							! First diagonal of the tridiagonal 
									! matrix
	TDB, &							! Second diagonal of the tridiagonal 
									! matrix
	TDC, &							! Third diagonal of the tridiagonal 
									! matrix
	TDD, &							! R.H.S. of the matrix equation
    TDX, &							! Solution of the matrix equation
    Temp_prev_iter					! Temperature values at the previous
									! iteration

! =======
! Intents
! =======

intent(in)   :: N_layers, N_nodes_layer, N_nodes, Therm_cond, thick, &
				Area, ibc_flag, h_in, T_in, heat_inner, obc_flag, h_out, &
				T_out, heat_outer, Transient_flag, delta_t, Temp_old, &
				rho, Cp
intent(out)  :: X, Temp

! =====================
! Executable Statements
! =====================

! Input Checks 

if(ibc_flag.ne.1.and.ibc_flag.ne.2)then
  write(*,*) 'Subroutine planewall called with bad ibc_flag'
  write(*,*) 'Exiting'
  stop
endif
if(obc_flag.ne.1.and.obc_flag.ne.2)then
  write(*,*) 'Subroutine planewall called with bad obc_flag'
  write(*,*) 'Exiting'
  stop
endif
if(transient_flag.ne.1.and.transient_flag.ne.2)then
  write(*,*) 'Subroutine planewall called with bad transient_flag'
  write(*,*) 'Exiting'
  stop
endif
N_nodes_test=0
do j=1,N_layers
  if(N_nodes_layer(j).lt.1) then
    write(*,'(" Number of nodes in layer ",I2," is ",I3, &
		  " which is less than minimum required (1)")') &
		  j,N_nodes_layer(j)
	write(*,*)'Exiting'
	stop
  endif
  N_nodes_test=N_nodes_test+N_nodes_layer(j)
enddo
N_nodes_test=N_nodes_test+N_layers+1
if(N_nodes.ne.N_nodes_test)then
  write(*,*) 'Value of N_nodes passed to subroutine planewall invalid'
  write(*,*) "N_nodes_test,N_nodes",N_nodes_test,N_nodes
  write(*,*) 'Should be = Total no. of internal nodes + no. of layers + 1' 
  write(*,*) 'Exiting'
  stop
endif

! Fill location array : point_map

point_map=0
i=1
point_map(1)=1
do j=1,N_layers
  do k=1,N_nodes_layer(j)
    i=i+1
    point_map(i)=j
  enddo
  i=i+1
enddo
point_map(N_nodes)=N_layers

! Fill grid arrays : dx,X 

do j=1,N_layers
  dx(j)=thick(j)/(N_nodes_layer(j)+1)
enddo
X(1)=0.
do i=2,N_nodes
  layer_no=point_map(i)
  if(layer_no==0) layer_no=point_map(i-1)
  X(i)=X(i-1)+dx(layer_no)
enddo

! Guess temperatures

  if(transient_flag==1) then
      Temp=25.0
  elseif(transient_flag==2)then
    Temp=Temp_old
  endif

! Initialize beta's to 0.

  beta_1=0.
  do i=1,N_nodes-1
    beta(i)=0.
  enddo
  beta_left=0.
  beta_right=0.
  beta_N_nodes=0.


! Iteration loop

do niter=1,max_iter
  
  ! Find control volume interface temperatures (Tbar's)

    do i=1,N_nodes-1
	  Tbar(i)=(Temp(i)+Temp(i+1))/2.
	enddo

  ! Find alpha's at control volume interfaces (and beta's at interior 
  ! nodes in transient run) 

    do i=1,N_nodes-1
	  layer_no=point_map(i)
	  if(layer_no.ne.0) then
	    layer_no_cont_vol_interface=layer_no
	  else
	    layer_no_cont_vol_interface=point_map(i+1)
	  endif


	  conductivity=Therm_cond


	  alpha(i)=conductivity/dx(layer_no_cont_vol_interface)**2
	  if(transient_flag==2.and.layer_no.ne.0) then
	    density=rho
		specificheat=Cp

	    beta(i)=density*specificheat/delta_t
	  endif
	enddo
	 
  ! Fill tridiagonal vectors

    ! Fill 1st row

	if(ibc_flag==1) then
      TDB(1)=1.
      TDC(1)=0.
      TDD(1)=T_in
	else

	  conductivity=Therm_cond
	  
	  alpha_1=conductivity/dx(1)
	  if(transient_flag==2)then

	    density=rho
		specificheat=Cp
		
	    beta_1=density*dx(1)*specificheat/2./delta_t
	  endif
	  TDB(1)=h_in + beta_1 + alpha_1
	  TDC(1)=-alpha_1
	  if (Area(1).gt.1d-28) then
	    TDD(1)=beta_1*Temp_old(1) + h_in*T_in + heat_inner/Area(1)
	  else
        TDD(1)=beta_1*Temp_old(1) + h_in*T_in
      endif
	endif
	  
	! Fill i'th row
	 
    do i=2,N_nodes-1
	  if(point_map(i).ne.0) then

	  ! Fill coefficients of interior nodes

	    TDA(i)=-alpha(i-1)
		TDB(i)=alpha(i-1)+beta(i)+alpha(i)
		TDC(i)=-alpha(i)
		TDD(i)=beta(i)*Temp_old(i)
	  else

	  ! Fill coefficients of interface nodes

		layer_left=point_map(i-1)
		layer_right=point_map(i+1)


		conductivity_left=Therm_Cond
		conductivity_right=Therm_Cond

		alpha_left=conductivity_left*Area(layer_left)/dx(layer_left)
		alpha_right=conductivity_right*Area(layer_right)/dx(layer_right)
		if(transient_flag==2)then

		     density_left=rho
			 density_right=rho
		     specificheat_left=Cp
			 specificheat_right=Cp

		  beta_left=density_left*Area(layer_left)*dx(layer_left)* &
		    specificheat_left/2./delta_t
		  beta_right=density_right*Area(layer_right)*dx(layer_right)* &
		    specificheat_right/2./delta_t
		endif
		TDA(i)=-alpha_left
		TDB(i)=alpha_left + beta_left + beta_right + alpha_right
		TDC(i)=-alpha_right
		TDD(i)=(beta_left + beta_right)*Temp_old(i)
	  endif

	enddo
    
	! Fill last row

	if(obc_flag==1) then
	  TDA(N_nodes)=0.
	  TDB(N_nodes)=1.
	  TDD(N_Nodes)=T_out
	else

		 conductivity=Therm_Cond

	  alpha_N_nodes=conductivity/dx(N_layers)
	  if(transient_flag==2)then

		   density=rho
		   specificheat=Cp

	    beta_N_nodes=density*dx(N_layers)*specificheat/2./delta_t
	  endif
	  TDA(N_nodes)=-alpha_N_nodes
	  TDB(N_nodes)=alpha_N_nodes + beta_N_nodes + h_out
	  if (Area(N_layers).gt.1d-28) then
	    TDD(N_nodes)=beta_N_nodes*Temp_old(N_nodes) + h_out*T_out + &
	        heat_outer/Area(N_layers)
	  else
        TDD(N_nodes)=beta_N_nodes*Temp_old(N_nodes) + h_out*T_out
      endif
	endif

  ! Solve for Temperature values

  call tridiag(TDA,TDB,TDC,TDD,TDX,N_nodes)
  Temp_prev_iter=Temp
  Temp=TDX

  ! Check for convergence

  epsmax=0.
  do i=2,N_nodes
    eps=Temp(i)-Temp_prev_iter(i)
	if(abs(eps).gt.epsmax) epsmax=abs(eps)
  enddo
  if(epsmax.lt.eps_Temp) then

   return
  endif
enddo
write(*,*) 'Maximum iterations (50) in subroutine planewall exceeded'
write(*,*) 'Exiting'
stop

end subroutine planewall

! -------------------------
! Tridiagonal matrix solver
! -------------------------

subroutine tridiag(A,B,C,R,U,N)

implicit none

integer,parameter :: NMAX=1000
integer N,j
double precision GAM(NMAX),A(N),B(N),C(N),R(N),U(N),BET

if(B(1).eq.0.)pause
BET=B(1)
U(1)=R(1)/BET
do j=2,N
  GAM(j)=C(j-1)/BET
  BET=B(j)-A(j)*GAM(j)
  if(BET.eq.0.)pause
  U(j)=(R(j)-A(j)*U(j-1))/BET
enddo
do j=N-1,1,-1
  U(j)=U(j)-GAM(j+1)*U(j+1)
enddo

end subroutine tridiag

end module planewall_mod
