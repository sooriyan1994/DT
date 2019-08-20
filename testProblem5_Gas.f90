program Test_Problem

use MatrixSolve_mod
use planewall_mod

implicit none

integer, parameter:: &
	SC_NoOfNodes_i=10, &  !SC_NoOfNodes_i denotes the no. of nodes in the slab conduction
	nWalls = 9
	
integer :: i,j,k,z, &
	SC_NoLayers_i, &
	SC_IbcFlag_i, SC_ObcFlag_i, &
	SC_TransientFlag_i, &
	no_walls

real*8 :: sigma_d, &		  
		  SC_Hin_d, SC_Tin_d, SC_HeatInner_d, &
		  SC_Hout_d, SC_Tout_d, SC_HeatOuter_d, &
		  SC_DeltaT_d, simTime_d, endTime_d, &
		  SC_ConvCriterion_d, &
		  StoredEnergy_d, netOutHeat_d, &
		  Qin_plate, &
		  num,den
		  
real *8, dimension(SC_NoOfNodes_i) :: &
	SC_TempOld_ad, &
	SC_Temp_ad, &
	SC_X_ad
	
real *8, dimension(nWalls,SC_NoOfNodes_i) :: &
	Wall_TempOld_ad, &
	Wall_Temp_ad, &
	Wall_X_ad

real*8, dimension(nWalls) :: &
	wallTemp_ad, &
	wallEmiss_ad, &
	matrixC_ad, &
	wallHeatFlux_ad,&
	wallHeatFluxOut_ad,&
	SC_ThermCond_ad,SC_Rho_ad,SC_Cp_ad, &
	heat_inner_walls, &
	Walls_StoredEnergy_ad,&
	Walls_outHeat_ad, &
	GasHTC_ad, &
	Wall_Area_ad,&
	Wall_thk_ad
	
real*8, dimension(nWalls,nWalls) :: &
	shapeFactors_ad, &
	matrixA_ad

integer, dimension(1) :: SC_NoNodesLayer_ai

real*8, dimension(1) :: &
	SC_Thick_ad, &
	SC_Area_ad, &
	GasTemp_ad, &
	GasTempPrev_ad
	
integer, dimension(nWalls,nWalls) :: &	
	deltaKronecker_ai

!Simulation Time variable
simTime_d=0.0
sigma_d=5.67E-08

!Setting up Kronecker's delta
do i=1,nWalls
	do j=1,nWalls
		if (j==i) then
			deltaKronecker_ai(i,j)=1
		else 
			deltaKronecker_ai(i,j)=0
		end if
	end do
end do

!Files to write Gas temperatures and results
open(200, file='Results.dat', ACTION='WRITE')
open(201, file='WallTemp.dat', ACTION='WRITE')
open(202, file='GasTemp.dat', ACTION='WRITE')

write(200,*) "Time, Temperature, HeatIn, HeatOut, HeatLossConvection, EnergyStored"
write(201,*) "Simulation_Time, Wall_2, Wall_3, Wall_4, Wall_5, Wall_6, Wall_7, Wall_8, Wall_9"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading data from input file
open(2, file = 'inputs3.dat')

read(2,*) ! Wall Temperatures
read(2,*)(wallTemp_ad(i),i=1,nWalls)

read(2,*) ! Wall Emissivities
read(2,*)(wallEmiss_ad(i),i=1,nWalls)

read(2,*) ! Wall HeatFlux
read(2,*)(wallHeatFlux_ad(i),i=2,nWalls)

read(2,*) ! Shape Factors
do j=1,nWalls
	read(2,*)(shapeFactors_ad(j,i),i=1,nWalls)
end do

! read simulation end Time
read(2,*)
read(2,*) endTime_d

! read Time step
read(2,*)
read(2,*) SC_DeltaT_d

! read material properties: k, rho, Cp
read(2,*)
do i=1,nWalls
	read(2,*) SC_ThermCond_ad(i), SC_Rho_ad(i), SC_Cp_ad(i)
end do 

! read area
read(2,*)
read(2,*)(Wall_Area_ad(i),i=1,nWalls)

! read thickness
read(2,*)
read(2,*)(Wall_thk_ad(i),i=1,nWalls)

! read outer surface BC:: h(W/m2k), T_ambient (K)
read(2,*)
read(2,*) SC_Hout_d, SC_Tout_d 

! read Convergence criterion
read(2,*)
read(2,*) SC_ConvCriterion_d

! read Qin plate :
read(2,*)
read(2,*) Qin_plate

! read Gas temperature :
read(2,*)
read(2,*) GasTemp_ad

! read HTC for gas convection on inner sides
read(2,*)
read(2,*)(GasHTC_ad(i),i=1,nWalls)

close(2)

! Finished reading inputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Heat flux due to the radiant plate
wallHeatFlux_ad(1)=Qin_plate/Wall_Area_ad(1)

do i=1, nWalls
	do j=1, SC_NoOfNodes_i
		Wall_TempOld_ad(i,j)=wallTemp_ad(i)
		Wall_Temp_ad(i,j)=wallTemp_ad(i)
	end do	
end do

GasTempPrev_ad = GasTemp_ad

simTime_d = simTime_d + SC_DeltaT_d
	
! creating matrix A
20 do while(simTime_d<=endTime_d)
	z=0

	10  z=z+1
	
	! Kronecker's Delta
	do i=1,nWalls
		do j=1,nWalls
			if (i.eq.1) then
				matrixA_ad(i,j)=deltaKronecker_ai(i,j)-shapeFactors_ad(i,j)
			else
				matrixA_ad(i,j)=deltaKronecker_ai(i,j) - &
							(1.0d0-wallEmiss_ad(i))*shapeFactors_ad(i,j)
			end if
		end do
	end do

	! creating matrix C

	do i=1,nWalls
		if (i.eq.1) then
			matrixC_ad(i)=wallHeatFlux_ad(i)
		else
			matrixC_ad(i)=wallEmiss_ad(i)*sigma_d*wallTemp_ad(i)**4
		end if
	end do
	
	! Solving AX = C
	call matrixsolver(nWalls,matrixA_ad,matrixC_ad,wallHeatFluxOut_ad)


	!"-----------------------   RESULTS --------------------------"

	do i=2,nWalls
		heat_inner_walls(i) = wallEmiss_ad(i)/(1.0d0-wallEmiss_ad(i))* &
				(sigma_d * wallTemp_ad(i)**4 - wallHeatFluxOut_ad(i)) + &
				(GasHTC_ad(i)*(wallTemp_ad(i) - GasTemp_ad(1)))
		!To compare the heat input and heat output, divide by area
		heat_inner_walls(i)=heat_inner_walls(i)*Wall_Area_ad(i) 
	end do	


	!"-----------------------   RESULTS --------------------------"

	do i=1,5
		if (i.eq.1) then
		wallTemp_ad(i)=(1/sigma_d * (wallHeatFlux_ad(i)* (1-wallEmiss_ad(i))/wallEmiss_ad(i) &
							+ wallHeatFluxOut_ad(i)))**0.25
		end if
	end do

	do i=2,nWalls
		SC_NoLayers_i=1
		SC_NoNodesLayer_ai(1) = SC_NoOfNodes_i - 2
		SC_Thick_ad(1)=Wall_thk_ad(i)
		SC_Area_ad(1)=Wall_Area_ad(i)
		SC_IbcFlag_i=2
		SC_Hin_d=0.0
		SC_Tin_d=0.0
		SC_HeatInner_d=-heat_inner_walls(i)
		SC_ObcFlag_i=2
		SC_HeatOuter_d=0.0
		SC_TransientFlag_i=2
		
		do k = 1, SC_NoOfNodes_i
			SC_TempOld_ad(k)=Wall_TempOld_ad(i,k)
		end do	
		
		! calling slab conduction 

		call planewall(SC_NoLayers_i, &
						 SC_NoNodesLayer_ai, &
						 SC_NoOfNodes_i, &
						 SC_ThermCond_ad(i), &
						 SC_Thick_ad, &
						 SC_Area_ad, &
						 SC_IbcFlag_i, &
						 SC_Hin_d, &
						 SC_Tin_d, &
						 SC_HeatInner_d, &
						 SC_ObcFlag_i, &
						 SC_Hout_d, &
						 SC_Tout_d, &
						 SC_HeatOuter_d, &
						 SC_TransientFlag_i, &
						 SC_DeltaT_d, &
						 SC_TempOld_ad, &
						 SC_Rho_ad(i), &
						 SC_Cp_ad(i), &
						 SC_X_ad, &
						 SC_Temp_ad)
						 
		
		do j = 1, SC_NoOfNodes_i
			Wall_Temp_ad(i,j) = SC_Temp_ad(j)
			Wall_X_ad(i,j) = SC_X_ad(j)
		end do	
			
	end do	
	
	! Calculating Gas Temperature
	num=0.0
	den=0.0
	do i=1,nWalls
		num=num + GasHTC_ad(i)*Wall_Area_ad(i)*wallTemp_ad(i)
		den=den + GasHTC_ad(i)*Wall_Area_ad(i)
	end do
	GasTemp_ad(1)=num/den

	if ((dabs(Wall_Temp_ad(2,1) - wallTemp_ad(2)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(3,1) - wallTemp_ad(3)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(4,1) - wallTemp_ad(4)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(5,1) - wallTemp_ad(5)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(6,1) - wallTemp_ad(6)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(7,1) - wallTemp_ad(7)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(8,1) - wallTemp_ad(8)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(9,1) - wallTemp_ad(9)).gt.SC_ConvCriterion_d).and. &
		(dabs(GasTempPrev_ad(1) - GasTemp_ad(1)).gt.SC_ConvCriterion_d))then 
		do i=2,nWalls
			wallTemp_ad(i)=Wall_Temp_ad(i,1)
		end do
		GasTempPrev_ad(1) = GasTemp_ad(1)
		goto 10	
	else
		do i=2,nWalls
			do k = 1, SC_NoOfNodes_i
				SC_TempOld_ad(k)=Wall_TempOld_ad(i,k)
				SC_Temp_ad(k)=Wall_Temp_ad(i,k)
				SC_X_ad(k)=Wall_X_ad(i,k)
			end do
			SC_Area_ad(1)=Wall_Area_ad(i)	
			call energyStored(SC_NoOfNodes_i,&
							  SC_DeltaT_d,&
							  SC_X_ad,&
							  SC_TempOld_ad,&
							  SC_Temp_ad,&
							  SC_Rho_ad(i),&
							  SC_Cp_ad(i),&
							  SC_Area_ad(1),&
							  StoredEnergy_d)
			Walls_StoredEnergy_ad(i)=StoredEnergy_d				  
		end do
						  
		Wall_TempOld_ad=Wall_Temp_ad
		do i=2,nWalls
			Walls_outHeat_ad(i)=SC_Hout_d*Wall_Area_ad(i)*(Wall_Temp_ad(i,SC_NoOfNodes_i)-SC_Tout_d)
		end do	
		
		netOutHeat_d = 0.0
		do i=2,nWalls
			netOutHeat_d = netOutHeat_d + Walls_outHeat_ad(i) + Walls_StoredEnergy_ad(i)
		end do
					
		write(200,*) simTime_d, Qin_plate,&
					netOutHeat_d, Walls_outHeat_ad(2), Walls_outHeat_ad(3),&
					Walls_outHeat_ad(4),Walls_outHeat_ad(5), Walls_outHeat_ad(6), &
					Walls_outHeat_ad(7), Walls_outHeat_ad(8),Walls_outHeat_ad(9), &
					Walls_StoredEnergy_ad(2), Walls_StoredEnergy_ad(3), &
					Walls_StoredEnergy_ad(4), Walls_StoredEnergy_ad(5), &
					Walls_StoredEnergy_ad(6), Walls_StoredEnergy_ad(7), &
					Walls_StoredEnergy_ad(8), Walls_StoredEnergy_ad(9)
					

		write(201,*) simTime_d, (Wall_Temp_ad(i,1), i=2,nWalls), (Wall_Temp_ad(i,2), i=2,nWalls), &
			(Wall_Temp_ad(i,3), i=2,nWalls), (Wall_Temp_ad(i,4), i=2,nWalls), &
			(Wall_Temp_ad(i,5), i=2,nWalls), (Wall_Temp_ad(i,6), i=2,nWalls), &
			(Wall_Temp_ad(i,7), i=2,nWalls), (Wall_Temp_ad(i,8), i=2,nWalls), &
			(Wall_Temp_ad(i,9), i=2,nWalls), (Wall_Temp_ad(i,10), i=2,nWalls)
		
		write(202,*)simTime_d, (Wall_Temp_ad(2,1)+Wall_Temp_ad(3,1)&
			+Wall_Temp_ad(4,1)+Wall_Temp_ad(5,1)+Wall_Temp_ad(6,1)+Wall_Temp_ad(7,1)&
			+Wall_Temp_ad(8,1)+Wall_Temp_ad(9,1))/4, GasTemp_ad
		
		simTime_d = simTime_d + SC_DeltaT_d
		goto 20	
	end if	
end do
print *, "FINAL results"

close(200)
close(201)
close(202)
	
end program Test_Problem

! --------------------------------------------------!
!   Energy Stored Subroutine (thermal inertia)      !
! --------------------------------------------------!

subroutine energyStored(noOfNodes_i,timeStep_d,loc_ad,prevTemp_ad,currTemp_ad,&
						rho,Cp,area,enStored_d)

implicit none

integer noOfNodes_i, ES_i

real*8, dimension(noOfNodes_i) :: &
	loc_ad,prevTemp_ad,currTemp_ad
	
real *8 enStored_d, ES_X, ES_TempOld, ES_TempNew, &
		rho,Cp,area, timeStep_d


enStored_d=0.0

do ES_i=1,(noOfNodes_i-1)
	ES_X=loc_ad(ES_i+1)-loc_ad(ES_i)
	ES_TempOld=(prevTemp_ad(ES_i+1)+prevTemp_ad(ES_i))/2.0d0
	ES_TempNew=(currTemp_ad(ES_i+1)+currTemp_ad(ES_i))/2.0d0
	
	enStored_d=enStored_d+(rho*ES_X*area)*Cp*(ES_TempNew-ES_TempOld)/timeStep_d
	!write(202,*)ES_i,ES_X,ES_TempOld,ES_TempNew,enStored_d
end do
!print *,"in subroutine------------------------",loc_ad
!print *, prevTemp_ad
!print *, currTemp_ad


end subroutine energyStored

