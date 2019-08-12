program Test_Problem

use MatrixSolve_mod
use planewall_mod

implicit none

integer, parameter:: &
	nWalls=5, &
	SC_NoOfNodes_i=10
	
integer :: i,j,k,z, &
	SC_NoLayers_i, &
	SC_IbcFlag_i, SC_ObcFlag_i, &
	SC_TransientFlag_i

real*8 :: sigma_d, temp_c, heat_1, heat_2, &		  
		  SC_Hin_d, SC_Tin_d, SC_HeatInner_d, &
		  SC_Hout_d, SC_Tout_d, SC_HeatOuter_d, &
		  SC_DeltaT_d, simTime_d, endTime_d, &
		  outHeat_d, SC_ConvCriterion_d, &
		  StoredEnergy_d, netOutHeat_d, &
		  Cy_ThermCond_d,Cy_Rho_d,Cy_Cp_d, &
		  Qin_Cy_d, Cy_TempOld_ad, Cy_Temp_ad, &
		  	Cy_Thick_ad, &
			Cy_Dia_ad,Cy_StoredEnergy_d, &
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
	GasHTC_ad
	

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

simTime_d=0.0
sigma_d=5.67E-08

do i=1,nWalls
	do j=1,nWalls
		if (j==i) then
			deltaKronecker_ai(i,j)=1
		else 
			deltaKronecker_ai(i,j)=0
		end if
	end do
end do


open(200, file='Results.dat', ACTION='WRITE')
open(201, file='WallTemp2.dat', ACTION='WRITE')
open(202, file='WallTemp3.dat', ACTION='WRITE')
open(203, file='WallTemp4.dat', ACTION='WRITE')
open(204, file='WallTemp5.dat', ACTION='WRITE')
open(205, file='GasTemp.dat', ACTION='WRITE')
write(200,*) "Time, Temperature, HeatIn, HeatOut, HeatLossConvection, EnergyStored"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading data from input file
open(2, file='inputs3.dat')
read(2,*) ! Wall Temperatures
read(2,*)(wallTemp_ad(i),i=1,nWalls)
print *,wallTemp_ad

read(2,*) ! Wall Emissivities
read(2,*)(wallEmiss_ad(i),i=1,nWalls)

read(2,*) ! Wall HeatFlux
read(2,*)(wallHeatFlux_ad(i),i=1,nWalls)
print *, "WallHeatFlux ", wallHeatFlux_ad
!print *,wallEmiss_ad
read(2,*) ! Shape Factors
do j=1,nWalls
	read(2,*)(shapeFactors_ad(j,i),i=1,nWalls)
end do

do i=1,nWalls
	print *, shapeFactors_ad(i,:)
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
! read outer surface BC:: h(W/m2k), T_ambient (K)
read(2,*)
read(2,*) 	SC_Hout_d, SC_Tout_d 

! read Convergence criterion
read(2,*)
read(2,*) 	SC_ConvCriterion_d

! read cylinder material properties: k, rho, Cp
read(2,*)
read(2,*) Cy_ThermCond_d, Cy_Rho_d, Cy_Cp_d

! read cylinder : D, T
read(2,*)
read(2,*) Cy_Thick_ad, Cy_Dia_ad	

! read Qin cylinder :
read(2,*)
read(2,*) Qin_Cy_d

! read HTC for gas convection on inner sides
read(2,*)
read(2,*) (GasHTC_ad(i), i=1,nWalls)

 close(2)
! Finished reading inputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


wallHeatFlux_ad(1)=Qin_Cy_d/(3.14*Cy_Dia_ad)

do i=1, nWalls
	do j=1, SC_NoOfNodes_i
		Wall_TempOld_ad(i,j)=wallTemp_ad(i)
		Wall_Temp_ad(i,j)=wallTemp_ad(i)
	end do	
end do

GasTemp_ad(1)=750.0
GasTempPrev_ad = GasTemp_ad
!do i=1, nWalls
!	write(200,*) i, simTime_d, SC_Temp_ad(i,1)
!end do

simTime_d=simTime_d+SC_DeltaT_d
	
!	Cy_TempOld_ad = wallTemp_ad(1) 


! creating matrix A
20 do while(simTime_d<=endTime_d)
!	simTime_d=simTime_d+SC_DeltaT_d
	z=0
	print *, "simTime_d", simTime_d

	10 print *,"		Iteration: ", z
	

	z=z+1
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

	!print *, "Matrix A"
	!do i=1,nWalls
	!	print *, matrixA_ad(i,:)
	!end do

	!print *
	! creating matrix C

	do i=1,nWalls
		if (i.eq.1) then
			matrixC_ad(i)=wallHeatFlux_ad(i)
		else
			matrixC_ad(i)=wallEmiss_ad(i)*sigma_d*wallTemp_ad(i)**4
		end if
	end do

	!print *, "Matrix C"
	!print *, matrixC_ad

		call matrixsolver(nWalls,matrixA_ad,matrixC_ad,wallHeatFluxOut_ad)


	!print *
	!print *, "qo"
	!print *, wallHeatFluxOut_ad
	!print *


	!print *,"-----------------------   RESULTS --------------------------"


!	heat_1=wallEmiss_ad(1)/(1-wallEmiss_ad(1))* &
!				(sigma_d*wallTemp_ad(1)**4-wallHeatFluxOut_ad(1))
!	heat_1=heat_1*2.0*22.0/7.0*0.5
	!print *, "Heat in 1st surface = ", heat_1
	do i=2,nWalls
		heat_inner_walls(i)=wallEmiss_ad(i)/(1.0d0-wallEmiss_ad(i))* &
				(sigma_d*wallTemp_ad(i)**4-wallHeatFluxOut_ad(i)) + &
				(GasHTC_ad(i)*(wallTemp_ad(i)-GasTemp_ad(1)))
		heat_inner_walls(i)=heat_inner_walls(i)*2.0
	enddo	

	!print *, "Heat in 2nd surface = ", heat_2

	!print *,"-----------------------   RESULTS --------------------------"

	do i=1,5
		if (i.eq.1) then
		temp_c=(1/sigma_d * (wallHeatFlux_ad(i)* wallEmiss_ad(i)/(1-wallEmiss_ad(i)) &
							+ wallHeatFluxOut_ad(i)))**0.25
		wallTemp_ad(i)=temp_c
	!	write(304,*)simTime_d, i, temp_c
		end if
	!	print *, "Temp of wall" ,i, temp_c 
	end do

	do i=2,nWalls
		SC_NoLayers_i=1
		SC_NoNodesLayer_ai(1) = SC_NoOfNodes_i - 2
		SC_Thick_ad(1)=0.33
		SC_Area_ad(1)=2.0
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
						 
				!write(303,*) "old",Cy_TempOld_ad
!				Cy_Temp_ad = ((Qin_Cy_d*SC_DeltaT_d)/ &
!				(Cy_Rho_d*Cy_Cp_d*3.14*Cy_Thick_ad*((Cy_Dia_ad/2)**2)))+Cy_TempOld_ad	
			
			!write(303,*) Cy_Temp_ad
		do j = 1, SC_NoOfNodes_i
			Wall_Temp_ad(i,j) = SC_Temp_ad(j)
			Wall_X_ad(i,j) = SC_X_ad(j)
		end do	
			
	end do	
	
		! Calculating Gas Temperature
	num=0.0
	den=0.0
	do i=1,nWalls
		if (i.eq.1) then
			num=num + GasHTC_ad(i)*(3.14*Cy_Dia_ad)*wallTemp_ad(i)
			den=den + GasHTC_ad(i)*(3.14*Cy_Dia_ad)
		else
			num=num + GasHTC_ad(i)*2.0*wallTemp_ad(i)
			den=den + GasHTC_ad(i)*2.0	
		endif
	end do
	GasTemp_ad(1)=num/den
!	print *, z, GasTemp_ad	
	
	 !print *,"X array ", SC_X_ad
	 !print *, "Temperatures ", SC_Temp_ad
	 !outHeat_d=SC_Hout_d*SC_Area_ad(1)*(SC_Temp_ad(SC_NoOfNodes_i)-SC_Tout_d)
	 
	!print *, "Error ", dabs(SC_Temp_ad(1) - wallTemp_ad(2))
	if ((dabs(Wall_Temp_ad(2,1) - wallTemp_ad(2)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(3,1) - wallTemp_ad(3)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(4,1) - wallTemp_ad(4)).gt.SC_ConvCriterion_d).and. &
		(dabs(Wall_Temp_ad(5,1) - wallTemp_ad(5)).gt.SC_ConvCriterion_d).and. &
		(dabs(GasTempPrev_ad(1) - GasTemp_ad(1)).gt.SC_ConvCriterion_d))then 
		wallTemp_ad(2)=Wall_Temp_ad(2,1)
		wallTemp_ad(3)=Wall_Temp_ad(3,1)
		wallTemp_ad(4)=Wall_Temp_ad(4,1)
		wallTemp_ad(5)=Wall_Temp_ad(5,1)
		GasTempPrev_ad(1) = GasTemp_ad(1)
		goto 10	
	else
		do i=2,nWalls
			do k = 1, SC_NoOfNodes_i
				SC_TempOld_ad(k)=Wall_TempOld_ad(i,k)
				SC_Temp_ad(k)=Wall_Temp_ad(i,k)
				SC_X_ad(k)=Wall_X_ad(i,k)
			end do
		SC_Area_ad(1)=2.0	
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
						  
	!Cy_StoredEnergy_d = Cy_Rho_d*Cy_Cp_d*3.14*Cy_Thick_ad*((Cy_Dia_ad/2)**2)&
	!					*(Cy_Temp_ad-Cy_TempOld_ad)/SC_DeltaT_d
						
		Wall_TempOld_ad=Wall_Temp_ad
		!Cy_TempOld_ad = Cy_Temp_ad
		do i=2,nWalls
			Walls_outHeat_ad(i)=SC_Hout_d*SC_Area_ad(1)*(Wall_Temp_ad(i,SC_NoOfNodes_i)-SC_Tout_d)
		end do	
		netOutHeat_d=Walls_outHeat_ad(2)+Walls_StoredEnergy_ad(2)+&
					Walls_outHeat_ad(3)+Walls_StoredEnergy_ad(3)+&
					Walls_outHeat_ad(4)+Walls_StoredEnergy_ad(4)+&
					Walls_outHeat_ad(5)+Walls_StoredEnergy_ad(5)
					
		write(200,*) simTime_d, Qin_Cy_d,&
					netOutHeat_d, Walls_outHeat_ad(2), Walls_outHeat_ad(3),&
					Walls_outHeat_ad(4),Walls_outHeat_ad(5),Walls_StoredEnergy_ad(2),&
					Walls_StoredEnergy_ad(3),&
					Walls_StoredEnergy_ad(4),&
					Walls_StoredEnergy_ad(5)
					

			write(201,*) simTime_d, (Wall_Temp_ad(2,i), &
				i=1,SC_NoOfNodes_i)
							write(202,*) simTime_d, (Wall_Temp_ad(3,i), &
				i=1,SC_NoOfNodes_i)
							write(203,*) simTime_d, (Wall_Temp_ad(4,i), &
				i=1,SC_NoOfNodes_i)
							write(204,*) simTime_d, (Wall_Temp_ad(5,i), &
				i=1,SC_NoOfNodes_i)
			
			!write(205,*) simTime_d, GasTemp_ad
			write(205,*)simTime_d, (Wall_Temp_ad(2,1)+Wall_Temp_ad(3,1)&
			+Wall_Temp_ad(4,1)+Wall_Temp_ad(5,1))/4, GasTemp_ad

		simTime_d=simTime_d+SC_DeltaT_d
		goto 20	
	end if	
end do
 print *, "FINAL results"
! print *,"X array ", SC_X_ad
 !print *, "Temperatures ", SC_Temp_ad
close(200)
close(201)
close(202)
close(203)
close(204)
close(205)
	
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

!print *,"Out subroutine------------------------"


end subroutine energyStored

