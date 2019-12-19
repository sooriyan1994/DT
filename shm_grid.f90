!=========================================================================
!
! Title         : shm_grid.f90
! Application   : Shape
! Copyright     : Tata Consultancy Services
! Author        : Niranjan Pedanekar, Chetan Malhotra
! Creation Date : 19/07/99
! Requirements  : Module "shm_glob_vars_mod"
! Description   : This file contains the subroutine which divides the 
!                 various inner surfaces of the reheating furnace into
!                 sections and stores their coordinates in structures
! Limitations   : 
! Dependencies  : Fortran 90
! Modifications :		WHO			WHEN		WHY
!				  Chetan Malhotra  20-8-99    Added gas
!				  Chetan Malhotra  9-3-2000   Made conducive to RHF model
!				  Chetan Malhotra  24-3-2000  Split wall into above and 
!											  below gas for gas 
!											  calculations
!=========================================================================
!
! Name              : shm_grid
! Description       : Divides the furnace geometry into arrays of surfaces
!                     and divides the surfaces into sections and stores 
!                     their coordinates into structures
! Parameters        : None
! Limitations		:
! Globals updated   : gv_roof_as,
!					  gv_wtf_as,gv_wtb_as,
!					  gv_wtfag_as,gv_wtbag_as,gv_wtfbg_as,gv_wtbbg_as,
!					  gv_st_as,
!					  gv_gtu_as,gv_gtd_as
!                     gv_floor_as,
!					  gv_wbf_as,gv_wbb_as,
!					  gv_wbfag_as,gv_wbbag_as,gv_wbfbg_as,gv_wbbbg_as,
!					  gv_sb_as,
!					  gv_gbu_as,gv_gbd_as
! Files updated     : None
! Calls             : None
!
!=========================================================================

subroutine shm_grid

use shm_glob_vars_mod
implicit none

integer gr_i_i, &                   ! Surface counter
        gr_j_i, &                   ! Division counter
        gr_k_i, &                   ! Flag/Local counter
        gr_m_i, &                   ! Section counter
		gr_l_i, &					! Division counter
		gr_n_i, &					! Local Counter
		gr_temp_i					! Temporary counter

real*8  gr_sect_w_d, &              ! Section width
		gr_sect_height_d, &         ! Section height/length
        gr_slope_d, &               ! Slope of surface in xy plane
        gr_theta_d                  ! Angle made by surface with x axis

real*8	gr_tempbp_ad(gv_MAX_SECTIONS)

!---------------------------------------------
! Segment Coordinates for Roof Active faces     
!---------------------------------------------
do gr_i_i = 1,gv_NoOfTopSurf_i
	gv_roof_as(gr_i_i)%st_x_d = gv_rbpx_ad(gr_i_i)
	gv_roof_as(gr_i_i)%st_y_d = gv_rbpy_ad(gr_i_i)
	gv_roof_as(gr_i_i)%end_x_d = gv_rbpx_ad(gr_i_i+1)
	gv_roof_as(gr_i_i)%end_y_d = gv_rbpy_ad(gr_i_i+1)
	
	! Dividing the surface along the front view and calculating the 
	! width of sections within the surface
	gr_j_i = 0
	gr_k_i = 1
	do while(gr_k_i.eq.1)
		gr_j_i = gr_j_i+1
		gr_sect_w_d = (((gv_roof_as(gr_i_i)%st_x_d- &
                   gv_roof_as(gr_i_i)%end_x_d)**2.0 + &
                  (gv_roof_as(gr_i_i)%st_y_d- &
                   gv_roof_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i   
		if(gr_sect_w_d.le.gv_sect_max_w_d) then
			gv_roof_as(gr_i_i)%n_sections_i = gr_j_i
			gv_roof_as(gr_i_i)%sect_width_d = gr_sect_w_d
			gr_k_i = 0
		endif
	end do
	
	! Dividing the surface along the plan view and calculating the 
	! height/length of sections within the surface
	do gr_m_i = 1, gv_roof_as(gr_i_i)%n_sections_i
		gr_k_i = 1
		gr_l_i = 0
		do while(gr_k_i.eq.1)
			gr_l_i = gr_l_i + 1
			gr_sect_height_d = gv_fwidth_d/gr_l_i
			if(gr_sect_height_d.le.gv_sect_max_len_d) then
				gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = gr_l_i
				gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
														gr_sect_height_d
				gr_k_i = 0	  	 
			endif
		enddo
	enddo
	
	! Calculating the slope of the surface
	if(dabs(gv_roof_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		gr_slope_d = ((gv_roof_as(gr_i_i)%end_y_d)- &
                (gv_roof_as(gr_i_i)%st_y_d))/ &
				((gv_roof_as(gr_i_i)%end_x_d)- &
                (gv_roof_as(gr_i_i)%st_x_d))
	else
		gr_slope_d = 0.0
	endif

	! Coordinate points
	gr_theta_d = datan(gr_slope_d)
	do gr_m_i = 1, gv_roof_as(gr_i_i)%n_sections_i
		do gr_n_i = 1, gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_roof_as(gr_i_i)%st_x_d + (gr_m_i-1) * &
					gv_roof_as(gr_i_i)%sect_width_d * dcos(gr_theta_d)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_roof_as(gr_i_i)%st_y_d + (gr_m_i-1) * &
					gv_roof_as(gr_i_i)%sect_width_d * dsin(gr_theta_d)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1) = - gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * (gr_n_i - 1)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(1)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2) = - gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * gr_n_i
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3) = gv_roof_as(gr_i_i)%st_x_d+ &
					gr_m_i * gv_roof_as(gr_i_i)%sect_width_d * dcos(gr_theta_d)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3) = gv_roof_as(gr_i_i)%st_y_d+ &
					gr_m_i * gv_roof_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(2)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(3)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(3)
			gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4) = gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(1)
		enddo
	end do
enddo

! write(10,*)"printing Roof coords"
! do gr_i_i = 1,gv_NoOfTopSurf_i
	! write(10,*)"Roof",gr_i_i 
	! do gr_m_i =1,gv_roof_as(gr_i_i)%n_sections_i
		! write(10,*)"section_x",gr_m_i
		! do gr_n_i =1,gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			! write(10,*)"section_y",gr_n_i
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(1)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(1)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(1)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(2)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(2)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(2)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(3)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(3)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(3)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(4)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(4)
		! write(10,*)gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(4)
		! enddo
	! enddo
! enddo

!------------------------------------------------------------
! Segment Coordinates for Wall Top Front / Back Active sfaces
!------------------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
	gv_wtf_as(gr_i_i)%st_x_d = gv_rbpx_ad(gr_i_i)
	gv_wtf_as(gr_i_i)%st_y_d = gv_rbpy_ad(gr_i_i)
	gv_wtf_as(gr_i_i)%end_x_d = gv_rbpx_ad(gr_i_i+1)
	gv_wtf_as(gr_i_i)%end_y_d = gv_rbpy_ad(gr_i_i+1)
	gv_wtf_as(gr_i_i)%base_y_d = gv_HearthYCoord_d + gv_SlabThick_d
	
	! Dividing the surface along the front view and calculating the 
	! width of sections within the surface
	gr_j_i=0
	gr_k_i=1
	do while(gr_k_i.eq.1)
		gr_j_i = gr_j_i+1
		gr_sect_w_d = (((gv_wtf_as(gr_i_i)%st_x_d- &
                   gv_wtf_as(gr_i_i)%end_x_d)**2.0 + &
                  (gv_wtf_as(gr_i_i)%st_y_d- &
                   gv_wtf_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i   
		if(gr_sect_w_d.le.gv_sect_max_w_d) then
			gv_wtf_as(gr_i_i)%n_sections_i = gr_j_i
			gv_wtf_as(gr_i_i)%sect_width_d = gr_sect_w_d
			gr_k_i = 0
		endif
	end do
	
	! Calculating the slope of the surface
	if(dabs(gv_wtf_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		gr_slope_d = ((gv_wtf_as(gr_i_i)%end_y_d)- &
                (gv_wtf_as(gr_i_i)%st_y_d))/ &
				((gv_wtf_as(gr_i_i)%end_x_d)- &
                (gv_wtf_as(gr_i_i)%st_x_d))
	else
		gr_slope_d = 0.0
	endif
	gr_theta_d = datan(gr_slope_d)
	
	do gr_m_i = 1, gv_wtf_as(gr_i_i)%n_sections_i
		gr_k_i = 1
		gr_l_i = 0
		do while(gr_k_i.eq.1)
			gr_l_i = gr_l_i + 1
			gr_sect_height_d = gv_fwidth_d/gr_l_i
			if(gr_sect_height_d.le.gv_sect_max_len_d) then
				gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = gr_l_i
				gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
														gr_sect_height_d
				gr_k_i = 0	  	 
			endif
		enddo
	enddo
	
	! Wall Top Back Surfaces
	gv_wtb_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	gv_wtb_as(gr_i_i)%st_y_d = gv_wtf_as(gr_i_i)%st_y_d
	gv_wtb_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	gv_wtb_as(gr_i_i)%end_y_d = gv_wtf_as(gr_i_i)%end_y_d
	gv_wtb_as(gr_i_i)%n_sections_i = gv_wtf_as(gr_i_i)%n_sections_i
	gv_wtb_as(gr_i_i)%sect_width_d = gv_wtf_as(gr_i_i)%sect_width_d
										
	! Coordinate points
	do gr_m_i = 1, gv_wtf_as(gr_i_i)%n_sections_i
		gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = &
				gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = &
				gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		
		do gr_n_i = 1, gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_wtf_as(gr_i_i)%st_x_d + (gr_m_i-1)*&
					gv_wtf_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_wtf_as(gr_i_i)%st_y_d + (gr_m_i-1)*&
					gv_wtf_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1) = 0.0
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2) = gv_wtf_as(gr_i_i)%st_x_d + gr_m_i * &
					gv_wtf_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2) = gv_wtf_as(gr_i_i)%st_y_d + gr_m_i * &
					gv_wtf_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2) = 0.0
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(2)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3) = gv_wtf_as(gr_i_i)%base_y_d 
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3) = 0.0
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4) = gv_wtf_as(gr_i_i)%base_y_d
			gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4) = 0.0
		
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(1)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1) = -gv_fwidth_d
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(4)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(4)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2) = -gv_fwidth_d
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(3)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(3)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3) = -gv_fwidth_d
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(2)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4) = gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(2)
			gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4) = -gv_fwidth_d
		end do
	enddo
enddo
    
!write(10,*)"printing Wall Top Front coords"
!do gr_i_i = 1,gv_NoOfTopSurf_i
	!write(10,*)"Wall",gr_i_i 
	!do gr_m_i =1,gv_wtf_as(gr_i_i)%n_sections_i
		!write(10,*)"section_x",gr_m_i
		!do gr_n_i =1,gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		!	write(10,*)"section_y",gr_n_i
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(1)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(1)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(1)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(2)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(2)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(2)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(3)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(3)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(3)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(4)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(4)
		!write(10,*)gv_wtf_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(4)
		!enddo
	!enddo
!enddo

!write(10,*)"printing Wall Top Back coords"
!do gr_i_i = 1,gv_NoOfTopSurf_i
	!write(10,*)"Wall",gr_i_i 
	!do gr_m_i =1,gv_wtb_as(gr_i_i)%n_sections_i
		!write(10,*)"section_x",gr_m_i
		!do gr_n_i =1,gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		!	write(10,*)"section_y",gr_n_i
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(1)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(1)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(1)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(2)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(2)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(2)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(3)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(3)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(3)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		x_coord_ad(4)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		y_coord_ad(4)
		!write(10,*)gv_wtb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
		!		z_coord_ad(4)
		!enddo
	!enddo
!enddo
!stop

!---------------------------------------------------
! Segment Coordinates for Slab Top Active surfaces
!---------------------------------------------------

do gr_i_i=1,gv_NoOfSlabSurf_i
	gv_st_as(gr_i_i)%st_x_d = (gr_i_i-1)*gv_SlabWidth_d
	gv_st_as(gr_i_i)%st_y_d = gv_HearthYCoord_d + gv_SlabThick_d
	gv_st_as(gr_i_i)%end_x_d = gr_i_i*gv_SlabWidth_d
	gv_st_as(gr_i_i)%end_y_d = gv_HearthYCoord_d + gv_SlabThick_d
	gv_st_as(gr_i_i)%sect_width_d = gv_SlabWidth_d
	gv_st_as(gr_i_i)%n_sections_i = 1
	
	! Dividing the surface along the plan view and calculating the 
	! height/length of sections within the surface
	do gr_m_i = 1, gv_st_as(gr_i_i)%n_sections_i
		gr_k_i = 1
		gr_j_i = 0
		do while(gr_k_i.eq.1)
			gr_j_i = gr_j_i + 1
			gr_sect_height_d = gv_fwidth_d/gr_j_i
			if(gr_sect_height_d.le.gv_sect_max_len_d) then
				gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = gr_j_i
				gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
														gr_sect_height_d
				gr_k_i = 0	  	 
			endif
		enddo
	enddo
	
	! Coordinate points	
	do gr_m_i = 1,gv_st_as(gr_i_i)%n_sections_i
		do gr_n_i = 1, gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_st_as(gr_i_i)%st_x_d
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_st_as(gr_i_i)%st_y_d
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1)= - (gr_n_i - 1) * gv_st_as(gr_i_i)%&
					surf_sect_as(gr_m_i)%height_d 
					!-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2)= gv_st_as(gr_i_i)%end_x_d
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2)= gv_st_as(gr_i_i)%end_y_d
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2)= - (gr_n_i - 1) * gv_st_as(gr_i_i)%&
					surf_sect_as(gr_m_i)%height_d 
					!-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3)= gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(2)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3)= gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(2)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3)= - gr_n_i * gv_st_as(gr_i_i)%&
					surf_sect_as(gr_m_i)%height_d
					!-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4)= gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4)= gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(1)
			gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4)= - gr_n_i * gv_st_as(gr_i_i)%&
					surf_sect_as(gr_m_i)%height_d
					!-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
		enddo
	enddo
enddo

! write(10,*)"printing slab  coords"
! do gr_i_i = 1,gv_NoOfSlabSurf_i
	! write(10,*)"Slab",gr_i_i 
	! do gr_m_i =1,gv_st_as(gr_i_i)%n_sections_i
		! write(10,*)"section_x",gr_m_i
		! do gr_n_i =1,gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			! write(10,*)"section_y",gr_n_i
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(1)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(1)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(1)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(2)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(2)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(2)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(3)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(3)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(3)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! x_coord_ad(4)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! y_coord_ad(4)
		! write(10,*)gv_st_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
				! z_coord_ad(4)
		! enddo
	! enddo
! enddo
! stop

!-------------------------------------------------
! Segment Coordinates for Gas top-up Active sfaces    
!-------------------------------------------------

! do gr_i_i=1,gv_NoOfTopSurf_i
	! gv_gtu_as(gr_i_i)%st_x_d=gv_rbpx_ad(gr_i_i)
	! gv_gtu_as(gr_i_i)%st_y_d=gv_gtbpy_ad(gr_i_i)
	! gv_gtu_as(gr_i_i)%end_x_d=gv_rbpx_ad(gr_i_i+1)
	! gv_gtu_as(gr_i_i)%end_y_d=gv_gtbpy_ad(gr_i_i+1)

	! gr_j_i=0
	! gr_k_i=1
	! do while(gr_k_i.eq.1)
		! gr_j_i = gr_j_i+1
		! gr_sect_w_d = (((gv_gtu_as(gr_i_i)%st_x_d- &
                   ! gv_gtu_as(gr_i_i)%end_x_d)**2.0 + &
                  ! (gv_gtu_as(gr_i_i)%st_y_d- &
                   ! gv_gtu_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i   
		! if(gr_sect_w_d.le.gv_sect_max_w_d) then
			! gv_gtu_as(gr_i_i)%n_sections_i = gr_j_i
			! gv_gtu_as(gr_i_i)%sect_width_d = gr_sect_w_d
			! gr_k_i = 0
		! endif
	! end do
	
	! if(dabs(gv_gtu_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		! gr_slope_d=((gv_gtu_as(gr_i_i)%end_y_d)- &
					! (gv_gtu_as(gr_i_i)%st_y_d))/ &
					! ((gv_gtu_as(gr_i_i)%end_x_d)- &
					! (gv_gtu_as(gr_i_i)%st_x_d))
	! else
		! gr_slope_d=0.0
	! endif
  
	! gr_theta_d=datan(gr_slope_d)
  
	! do gr_m_i=1,gv_gtu_as(gr_i_i)%n_sections_i
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_gtu_as(gr_i_i)%st_y_d+ &
				! (gr_m_i-1)*gv_gtu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_gtu_as(gr_i_i)%st_y_d+ &
				! gr_m_i*gv_gtu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
	! end do  !do for gr_m_i no. of sections
! enddo   


!---------------------------------------------------
! Segment Coordinates for Gas top-down Active sfaces    
!---------------------------------------------------

! do gr_i_i=1,gv_NoOfTopSurf_i
	! gv_gtd_as(gr_i_i)%st_x_d=gv_gtu_as(gr_i_i)%st_x_d
	! gv_gtd_as(gr_i_i)%st_y_d=gv_gtu_as(gr_i_i)%st_y_d
	! gv_gtd_as(gr_i_i)%end_x_d=gv_gtu_as(gr_i_i)%end_x_d
	! gv_gtd_as(gr_i_i)%end_y_d=gv_gtu_as(gr_i_i)%end_y_d
	! gv_gtd_as(gr_i_i)%n_sections_i=gv_gtu_as(gr_i_i)%n_sections_i
	! gv_gtd_as(gr_i_i)%sect_width_d=gv_gtu_as(gr_i_i)%sect_width_d
	
	! do gr_m_i=1,gv_gtd_as(gr_i_i)%n_sections_i
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
			! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
	! end do  !do for gr_m_i no. of sections
! enddo   


!------------------------------------------------------------
! Segment Coordinates for Wall Top Front/Back Above/Below Gas 
! Active sfaces
!------------------------------------------------------------

! do gr_i_i=1,gv_NoOfTopSurf_i

	! gv_wtfag_as(gr_i_i)%st_x_d = gv_rbpx_ad(gr_i_i)
	! gv_wtfag_as(gr_i_i)%st_y_d = gv_rbpy_ad(gr_i_i)
	! gv_wtfag_as(gr_i_i)%end_x_d = gv_rbpx_ad(gr_i_i+1)
	! gv_wtfag_as(gr_i_i)%end_y_d = gv_rbpy_ad(gr_i_i+1)
	! gv_wtfag_as(gr_i_i)%n_sections_i = gv_roof_as(gr_i_i)%n_sections_i
	! gv_wtfag_as(gr_i_i)%sect_width_d = 0.0d0

	! gv_wtbag_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	! gv_wtbag_as(gr_i_i)%st_y_d = gv_wtf_as(gr_i_i)%st_y_d
	! gv_wtbag_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	! gv_wtbag_as(gr_i_i)%end_y_d = gv_wtf_as(gr_i_i)%end_y_d
	! gv_wtbag_as(gr_i_i)%n_sections_i = gv_wtf_as(gr_i_i)%n_sections_i
	! gv_wtbag_as(gr_i_i)%sect_width_d = gv_wtf_as(gr_i_i)%sect_width_d

	! gv_wtfbg_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	! gv_wtfbg_as(gr_i_i)%st_y_d = gv_wtf_as(gr_i_i)%st_y_d
	! gv_wtfbg_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	! gv_wtfbg_as(gr_i_i)%end_y_d = gv_wtf_as(gr_i_i)%end_y_d
	! gv_wtfbg_as(gr_i_i)%n_sections_i = gv_wtf_as(gr_i_i)%n_sections_i
	! gv_wtfbg_as(gr_i_i)%sect_width_d = gv_wtf_as(gr_i_i)%sect_width_d

	! gv_wtbbg_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	! gv_wtbbg_as(gr_i_i)%st_y_d = gv_wtf_as(gr_i_i)%st_y_d
	! gv_wtbbg_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	! gv_wtbbg_as(gr_i_i)%end_y_d = gv_wtf_as(gr_i_i)%end_y_d
	! gv_wtbbg_as(gr_i_i)%n_sections_i = gv_wtf_as(gr_i_i)%n_sections_i
	! gv_wtbbg_as(gr_i_i)%sect_width_d = gv_wtf_as(gr_i_i)%sect_width_d

	! do gr_m_i=1,gv_wtf_as(gr_i_i)%n_sections_i
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0

		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_HearthYCoord_d + gv_SlabThick_d
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_HearthYCoord_d + gv_SlabThick_d
		! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
	

		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d

		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
	
	! end do
! enddo     

!----------------------------------------------
! Segment Coordinates for Floor Active sfaces    
!----------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i
	! gv_floor_as(gr_i_i)%st_x_d=gv_fbpx_ad(gr_i_i)
	! gv_floor_as(gr_i_i)%st_y_d=gv_fbpy_ad(gr_i_i)
	! gv_floor_as(gr_i_i)%end_x_d=gv_fbpx_ad(gr_i_i+1)
	! gv_floor_as(gr_i_i)%end_y_d=gv_fbpy_ad(gr_i_i+1)  
	
	! gr_j_i=0
	! gr_k_i=1
	! do while(gr_k_i.eq.1)
		! gr_j_i=gr_j_i+1
		! gr_sect_w_d=(((gv_floor_as(gr_i_i)%st_x_d-gv_floor_as(gr_i_i)% &
                ! end_x_d)**2.0 + &
				! (gv_floor_as(gr_i_i)%st_y_d-gv_floor_as(gr_i_i)% &
                ! end_y_d)**2.0)**0.5)/gr_j_i  
    
		! if (gr_sect_w_d.le.gv_sect_max_w_d) then
			! gv_floor_as(gr_i_i)%n_sections_i = gr_j_i
			! gv_floor_as(gr_i_i)%sect_width_d = gr_sect_w_d
			! gr_k_i = 0
		! endif
	! end do
  
	! if(dabs(gv_floor_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		! gr_slope_d=((gv_floor_as(gr_i_i)%end_y_d)- &
            ! (gv_floor_as(gr_i_i)%st_y_d))/ &
            ! ((gv_floor_as(gr_i_i)%end_x_d)- &
            ! (gv_floor_as(gr_i_i)%st_x_d))
	! else
		! gr_slope_d=0.0
	! endif
  
	! gr_theta_d=datan(gr_slope_d)
	! do gr_m_i=1,gv_floor_as(gr_i_i)%n_sections_i
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_floor_as(gr_i_i)%st_x_d+ &
				! (gr_m_i-1)*gv_floor_as(gr_i_i)%sect_width_d*dcos(gr_theta_d) 
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_floor_as(gr_i_i)%st_y_d+ &
				! (gr_m_i-1)*gv_floor_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_floor_as(gr_i_i)%st_x_d+ &
				! gr_m_i*gv_floor_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_floor_as(gr_i_i)%st_y_d+ &
				! gr_m_i*gv_floor_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
	! end do  
! enddo   

!---------------------------------------------------------------
! Segment Coordinates for Wall Bottom Front / Back Active sfaces
!---------------------------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i
	! gv_wbf_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
	! gv_wbf_as(gr_i_i)%st_y_d=gv_floor_as(gr_i_i)%st_y_d
	! gv_wbf_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
	! gv_wbf_as(gr_i_i)%end_y_d=gv_floor_as(gr_i_i)%end_y_d
	! gv_wbf_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
	! gv_wbf_as(gr_i_i)%sect_width_d=0.0d0
	
	! gv_wbb_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
	! gv_wbb_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
	! gv_wbb_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
	! gv_wbb_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
	! gv_wbb_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
	! gv_wbb_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d
	
	! do gr_m_i=1,gv_wbf_as(gr_i_i)%n_sections_i
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)=gv_HearthYCoord_d
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)=gv_HearthYCoord_d
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
		
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
	! enddo
! enddo     
 
!----------------------------------------------------
! Segment Coordinates for Slab Bottom Active sfaces
!----------------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i
	! gv_sb_as(gr_i_i)%st_x_d = gv_floor_as(gr_i_i)%st_x_d
	! gv_sb_as(gr_i_i)%st_y_d = gv_HearthYCoord_d
	! gv_sb_as(gr_i_i)%end_x_d = gv_floor_as(gr_i_i)%end_x_d
	! gv_sb_as(gr_i_i)%end_y_d = gv_HearthYCoord_d
	! gv_sb_as(gr_i_i)%n_sections_i = gv_floor_as(gr_i_i)%n_sections_i
	! gv_sb_as(gr_i_i)%sect_width_d = 0.0d0
	
	! do gr_m_i=1,gv_sb_as(gr_i_i)%n_sections_i
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1) = &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1) = gv_HearthYCoord_d
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1) = &
				! -(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2) = &
				! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2) = gv_HearthYCoord_d
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2) = &
				! -(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3) = &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3) = gv_HearthYCoord_d
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3) = &
				! -(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4) = &
				! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4) = gv_HearthYCoord_d
		! gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4) = &
				! -(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
	! enddo
! enddo

!----------------------------------------------------
! Segment Coordinates for Gas bottom-up Active sfaces    
!----------------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i
	! gv_gbu_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
	! gv_gbu_as(gr_i_i)%st_y_d=gv_gbbpy_ad(gr_i_i)
	! gv_gbu_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
	! gv_gbu_as(gr_i_i)%end_y_d=gv_gbbpy_ad(gr_i_i+1)  
	! gv_gbu_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
	! gv_gbu_as(gr_i_i)%sect_width_d=&
			! (((gv_gbu_as(gr_i_i)%st_x_d- gv_gbu_as(gr_i_i)%end_x_d)**2.0 + &
            ! (gv_gbu_as(gr_i_i)%st_y_d- gv_gbu_as(gr_i_i)%end_y_d)**2.0)**0.5)/ &
			! gv_gbu_as(gr_i_i)%n_sections_i
			
	! if(dabs(gv_gbu_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		! gr_slope_d=((gv_gbu_as(gr_i_i)%end_y_d)- &
				! (gv_gbu_as(gr_i_i)%st_y_d))/ &
				! ((gv_gbu_as(gr_i_i)%end_x_d)- &
                ! (gv_gbu_as(gr_i_i)%st_x_d))
	! else
		! gr_slope_d=0.0
	! endif
	
	! gr_theta_d=datan(gr_slope_d)
	! do gr_m_i=1,gv_gbu_as(gr_i_i)%n_sections_i
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
			! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
			! gv_gbu_as(gr_i_i)%st_y_d+ &
			! (gr_m_i-1)*gv_gbu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
			! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
			! gv_gbu_as(gr_i_i)%st_y_d+ &
			! gr_m_i*gv_gbu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
	! end do 
! enddo   

!------------------------------------------------------
! Segment Coordinates for Gas bottom-down Active sfaces    
!------------------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i
	! gv_gbd_as(gr_i_i)%st_x_d=gv_gbu_as(gr_i_i)%st_x_d
	! gv_gbd_as(gr_i_i)%st_y_d=gv_gbu_as(gr_i_i)%st_y_d
	! gv_gbd_as(gr_i_i)%end_x_d=gv_gbu_as(gr_i_i)%end_x_d
	! gv_gbd_as(gr_i_i)%end_y_d=gv_gbu_as(gr_i_i)%end_y_d
	! gv_gbd_as(gr_i_i)%n_sections_i=gv_gbu_as(gr_i_i)%n_sections_i
	! gv_gbd_as(gr_i_i)%sect_width_d=gv_gbu_as(gr_i_i)%sect_width_d
	! do gr_m_i=1,gv_gbd_as(gr_i_i)%n_sections_i
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
	! end do  
! enddo   

!---------------------------------------------------------------
! Segment Coordinates for Wall Bottom Front/Back Above/Below Gas 
! Active sfaces
!---------------------------------------------------------------

! do gr_i_i=1,gv_NoOfBotSurf_i

	! gv_wbfag_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
	! gv_wbfag_as(gr_i_i)%st_y_d=gv_floor_as(gr_i_i)%st_y_d
	! gv_wbfag_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
	! gv_wbfag_as(gr_i_i)%end_y_d=gv_floor_as(gr_i_i)%end_y_d
	! gv_wbfag_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
	! gv_wbfag_as(gr_i_i)%sect_width_d=0.0d0
	
	! gv_wbbag_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
	! gv_wbbag_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
	! gv_wbbag_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
	! gv_wbbag_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
	! gv_wbbag_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
	! gv_wbbag_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

	! gv_wbfbg_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
	! gv_wbfbg_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
	! gv_wbfbg_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
	! gv_wbfbg_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
	! gv_wbfbg_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
	! gv_wbfbg_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

	! gv_wbbbg_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
	! gv_wbbbg_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
	! gv_wbbbg_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
	! gv_wbbbg_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
	! gv_wbbbg_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
	! gv_wbbbg_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

	! do gr_m_i=1,gv_wbf_as(gr_i_i)%n_sections_i
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= gv_HearthYCoord_d
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= gv_HearthYCoord_d
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0

		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0

		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d

		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				! gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d

	! end do
! enddo     

!----------------------------------------------------------------------
! Segment Coordinates for Roof Radiant plate top active faces     
!----------------------------------------------------------------------

do gr_i_i=1,gv_NoOfTopRadPlates_i
	gv_roof_RPt_as(gr_i_i)%st_x_d = gv_r_rpx_st_ad(gr_i_i)
	gv_roof_RPt_as(gr_i_i)%st_y_d = gv_r_rpy_st_ad(gr_i_i)
	gv_roof_RPt_as(gr_i_i)%end_x_d = gv_r_rpx_end_ad(gr_i_i)
	gv_roof_RPt_as(gr_i_i)%end_y_d = gv_r_rpy_end_ad(gr_i_i)
	
	! Dividing the surface along the plan view and calculating the
	! width of sections within the surface
	gr_j_i = 0
	gr_k_i = 1
	do while(gr_k_i.eq.1)
		gr_j_i = gr_j_i+1
		gr_sect_w_d = (((gv_roof_RPt_as(gr_i_i)%st_x_d- &
                gv_roof_RPt_as(gr_i_i)%end_x_d)**2.0 + &
                (gv_roof_RPt_as(gr_i_i)%st_y_d- &
                gv_roof_RPt_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i
    
		if(gr_sect_w_d.le.gv_sect_max_w_d) then
			gv_roof_RPt_as(gr_i_i)%n_sections_i = gr_j_i
			gv_roof_RPt_as(gr_i_i)%sect_width_d = gr_sect_w_d
			gr_k_i = 0	  	 
		endif
	end do
	
	! Dividing the surface along the plan view and calculating the 
	! height/length of sections within the surface
	do gr_m_i = 1, gv_roof_RPt_as(gr_i_i)%n_sections_i
		gr_k_i = 1
		gr_l_i = 0
		do while(gr_k_i.eq.1)
			gr_l_i = gr_l_i + 1
			gr_sect_height_d = gv_fwidth_d/gr_l_i
			if(gr_sect_height_d.le.gv_sect_max_len_d) then
				gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = gr_l_i
				gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
														gr_sect_height_d
				gr_k_i = 0	  	 
			endif
		enddo
	enddo
	
	! Calculating the slope of the surface
	if(dabs(gv_roof_RPt_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		gr_slope_d = ((gv_roof_RPt_as(gr_i_i)%end_y_d)- &
				(gv_roof_RPt_as(gr_i_i)%st_y_d))/ &
				((gv_roof_RPt_as(gr_i_i)%end_x_d)- &
                (gv_roof_RPt_as(gr_i_i)%st_x_d))
	else
		gr_slope_d = 0.0
	endif
	gr_theta_d = datan(gr_slope_d)
	
	! Segment Coordinates for Roof Radiant plate top down active faces  
	gv_roof_RPb_as(gr_i_i)%st_x_d = gv_roof_RPt_as(gr_i_i)%st_x_d
	gv_roof_RPb_as(gr_i_i)%st_y_d = gv_roof_RPt_as(gr_i_i)%st_y_d
	gv_roof_RPb_as(gr_i_i)%end_x_d = gv_roof_RPt_as(gr_i_i)%end_x_d
	gv_roof_RPb_as(gr_i_i)%end_y_d = gv_roof_RPt_as(gr_i_i)%end_y_d
	gv_roof_RPb_as(gr_i_i)%n_sections_i = gv_roof_RPt_as(gr_i_i)%n_sections_i
	gv_roof_RPb_as(gr_i_i)%sect_width_d = gv_roof_RPt_as(gr_i_i)%sect_width_d
	
	! Coordinate points 
	do gr_m_i = 1, gv_roof_as(gr_i_i)%n_sections_i
		gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = &
				gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
				gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d
		
		do gr_n_i = 1, gv_roof_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_roof_RPt_as(gr_i_i)%st_x_d + (gr_m_i-1) * &
					gv_roof_RPt_as(gr_i_i)%sect_width_d * dcos(gr_theta_d)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_roof_RPt_as(gr_i_i)%st_y_d + (gr_m_i-1) * &
					gv_roof_RPt_as(gr_i_i)%sect_width_d * dsin(gr_theta_d)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1) = - gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * (gr_n_i - 1)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(1)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4) = - gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * gr_n_i
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3) = gv_roof_RPt_as(gr_i_i)%st_x_d+ &
					gr_m_i*gv_roof_RPt_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3) = gv_roof_RPt_as(gr_i_i)%st_y_d+ &
					gr_m_i*gv_roof_RPt_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3) = - gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * gr_n_i
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(3)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(3)
			gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2) = - gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					height_d * (gr_n_i - 1)

			! Radiant plate top-down
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(1)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(1)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(1)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(4)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(4)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(4)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(3)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(3)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(3)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%x_coord_ad(2)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%y_coord_ad(2)
			gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4) = gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%&
					sect_as(gr_n_i)%z_coord_ad(2)
		enddo
	end do
enddo

write(10,*)"printing RP coords"
do gr_i_i = 1,gv_NoOfTopRadPlates_i
	write(10,*)"RD",gr_i_i
	do gr_m_i =1,gv_roof_RPt_as(gr_i_i)%n_sections_i
		write(10,*)"section_x",gr_m_i
		do gr_n_i =1,gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			write(10,*)"section_y",gr_n_i
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4)
			write(10,*)gv_roof_RPt_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4)
		enddo
	enddo
enddo

write(10,*)"printing RP down coords"
do gr_i_i = 1,gv_NoOfTopRadPlates_i
	write(10,*)"RD",gr_i_i
	do gr_m_i =1,gv_roof_RPb_as(gr_i_i)%n_sections_i
		write(10,*)"section_x",gr_m_i
		do gr_n_i =1,gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
			write(10,*)"section_y",gr_n_i
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(1)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(1)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(1)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(2)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(2)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(2)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(3)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(3)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(3)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					x_coord_ad(4)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					y_coord_ad(4)
			write(10,*)gv_roof_RPb_as(gr_i_i)%surf_sect_as(gr_m_i)%sect_as(gr_n_i)%&
					z_coord_ad(4)
		enddo
	enddo
enddo

!------------------------------------------------------------
! Segment Coordinates for Wall Top Front/Back Above/Below Radiant plate 
! Active sfaces
!------------------------------------------------------------
do gr_i_i=1,gv_NoOfTopSurf_i

	gv_wtfarp_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	gv_wtfarp_as(gr_i_i)%st_y_d = gv_wtf_as(gr_i_i)%st_y_d
	gv_wtfarp_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	gv_wtfarp_as(gr_i_i)%end_y_d = gv_wtf_as(gr_i_i)%end_y_d
	gv_wtfarp_as(gr_i_i)%sect_width_d = gv_wtf_as(gr_i_i)%sect_width_d
	gv_wtfarp_as(gr_i_i)%n_sections_i = gv_wtf_as(gr_i_i)%n_sections_i
	gv_wtfarp_as(gr_i_i)%base_y_d = gv_wtf_as(gr_i_i)%base_y_d
	
	gr_temp_i = 0
	gr_tempbp_ad = 0.0d0
	do gr_j_i=1, gv_NoOfTopRadPlates_i
		if(MAX(gv_wtfarp_as(gr_i_i)%st_y_d,gv_wtfarp_as(gr_i_i)%end_y_d).lt.&
		MIN(gv_r_rpy_st_ad(gr_j_i),gv_r_rpy_end_ad(gr_j_i))) then
			cycle
		else
			gr_temp_i = gr_temp_i + 1
			gr_tempbp_ad(gr_temp_i) = MIN(gv_r_rpy_st_ad(gr_j_i),&
					gv_r_rpy_end_ad(gr_j_i))
			if(gv_r_rpy_st_ad(gr_j_i).ne.gv_r_rpy_end_ad(gr_j_i).and.&
			MAX(gv_wtfarp_as(gr_i_i)%st_y_d,gv_wtfarp_as(gr_i_i)%end_y_d).gt.&
			MAX(gv_r_rpy_st_ad(gr_j_i),gv_r_rpy_end_ad(gr_j_i))) then
				gr_temp_i = gr_temp_i + 1
				gr_tempbp_ad(gr_temp_i) = MAX(gv_r_rpy_st_ad(gr_j_i),&
					gv_r_rpy_end_ad(gr_j_i))
			endif
		endif
	enddo
	
	if(gr_temp_i.ne.0) then
		gv_wtfarp_as(gr_i_i)%base_y_d = MIN(gr_tempbp_ad)
	
	! Dividing the surface along the plan view and calculating the
	! height/length of sections
	gr_j_i=0
	gr_k_i=1
	do while(gr_k_i.eq.1)
		gr_j_i = gr_j_i+1
		gr_sect_height_d = (MAX(gv_wtfarp_as(gr_i_i)%st_y_d, &
					gv_wtfarp_as(gr_i_i)%end_y_d) - &
					gv_wtfarp_as(gr_i_i)%base_y_d)/gr_j_i   
		if(gr_sect_height_d.le.gv_sect_max_w_d) then
			gv_wtfarp_as(gr_i_i)%n_sections_i = &
					gv_wtf_as(gr_i_i)%n_sections_i * gr_j_i
			gv_wtfarp_as(gr_i_i)%sect_height_d = gr_sect_height_d
			gr_k_i = 0
		endif
	end do
	
	! Dividing the surface along the plan view and calculating the 
	! height/length of sections within the surface
	do gr_m_i = 1, gv_wtfarp_as(gr_i_i)%n_sections_i
		gr_k_i = 1
		gr_l_i = 0
		do while(gr_k_i.eq.1)
			gr_l_i = gr_l_i + 1
			gr_sect_height_d = (MAX(gv_wtfarp_as(gr_i_i)%st_y_d, &
					gv_wtfarp_as(gr_i_i)%end_y_d) - &
					gv_wtfarp_as(gr_i_i)%base_y_d)/gr_j_i
			if(gr_sect_height_d.le.gv_sect_max_len_d) then
				gv_wtfarp_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = gr_l_i
				gv_wtfarp_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
														gr_sect_height_d
				gr_k_i = 0	  	 
			endif
		enddo
	enddo
	
	! Wall Top Back above Radiant plate 
	gv_wtbarp_as(gr_i_i)%st_x_d = gv_wtfarp_as(gr_i_i)%st_x_d
	gv_wtbarp_as(gr_i_i)%st_y_d = gv_wtfarp_as(gr_i_i)%st_y_d
	gv_wtbarp_as(gr_i_i)%end_x_d = gv_wtfarp_as(gr_i_i)%end_x_d
	gv_wtbarp_as(gr_i_i)%end_y_d = gv_wtfarp_as(gr_i_i)%end_y_d
	gv_wtbarp_as(gr_i_i)%n_sections_i = gv_wtfarp_as(gr_i_i)%n_sections_i
	gv_wtbarp_as(gr_i_i)%sect_width_d = gv_wtfarp_as(gr_i_i)%sect_width_d
	gv_wtbarp_as(gr_i_i)%base_y_d = gv_wtfarp_as(gr_i_i)%base_y_d
	
	do gr_m_i = 1, gv_wtf_as(gr_i_i)%n_sections_i
		gv_wtbarp_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i = &
						gv_wtfarp_as(gr_i_i)%surf_sect_as(gr_m_i)%n_sections_i
		gv_wtbarp_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d = &
						gv_wtfarp_as(gr_i_i)%surf_sect_as(gr_m_i)%height_d
		
		do gr_n_i = 1, gr_j_i
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(1)&
					= gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(1)&
					= MIN(gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1), &
					gv_wtfarp_as(gr_i_i)%base_y_d + gr_l_i * &
					gv_wtfarp_as(gr_i_i)%sect_height_d)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(1)&
					= 0.0
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(2)&
					= gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(2)&
					= MIN(gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1), & 
					gv_wtfarp_as(gr_i_i)%base_y_d + gr_l_i * &
					gv_wtfarp_as(gr_i_i)%sect_height_d)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(2)&
					= 0.0
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(3)&
					= gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(3)&
					= gv_wtfarp_as(gr_i_i)%base_y_d + (gr_l_i-1) * &
					gv_wtfarp_as(gr_i_i)%sect_height_d 
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(3)&
					= 0.0
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(4)&
					= gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(4)&
					= gv_wtfarp_as(gr_i_i)%base_y_d + (gr_l_i-1) * &
					gv_wtfarp_as(gr_i_i)%sect_height_d
			gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(4)&
					= 0.0
		
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(1)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(1)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(1)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(1)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(1)&
					= -gv_fwidth_d
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(2)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(4)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(2)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(4)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(2)&
					= -gv_fwidth_d
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(3)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(3)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(3)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(3)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(3)&
					= -gv_fwidth_d
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(4)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(2)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(4)&
					= gv_wtfarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(2)
			gv_wtbarp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(4)&
					= -gv_fwidth_d
		enddo
	enddo

	! WALL TOP FRONT/BOTTOM BELOW RP
	gv_wtfbrp_as(gr_i_i)%st_x_d = gv_wtf_as(gr_i_i)%st_x_d
	gv_wtfbrp_as(gr_i_i)%st_y_d = gv_r_rpy_st_ad(1)
	gv_wtfbrp_as(gr_i_i)%end_x_d = gv_wtf_as(gr_i_i)%end_x_d
	gv_wtfbrp_as(gr_i_i)%end_y_d = gv_r_rpy_end_ad(1)
	gv_wtfbrp_as(gr_i_i)%sect_width_d = gv_wtfbrp_as(gr_i_i)%end_x_d - &
							gv_wtfbrp_as(gr_i_i)%st_x_d
	gv_wtfbrp_as(gr_i_i)%base_y_d = gv_HearthYCoord_d + &
							gv_SlabThick_d
	
	! Dividing the surface along the plan view and calculating the
	! height/length of sections
	gr_j_i=0
	gr_k_i=1
	do while(gr_k_i.eq.1)
		gr_j_i = gr_j_i+1
		gr_sect_height_d = (MAX(gv_wtfbrp_as(gr_i_i)%st_y_d, &
					gv_wtfbrp_as(gr_i_i)%end_y_d) - &
					gv_wtfbrp_as(gr_i_i)%base_y_d)/gr_j_i   
		if(gr_sect_height_d.le.gv_sect_max_len_d) then
			gv_wtfbrp_as(gr_i_i)%n_sections_i = &
					gv_wtf_as(gr_i_i)%n_sections_i * gr_j_i
			gv_wtfbrp_as(gr_i_i)%sect_height_d = gr_sect_height_d
			gr_k_i = 0
		endif
	end do
	
	! Wall Top Back below Radiant plate 
	gv_wtbbrp_as(gr_i_i)%st_x_d = gv_wtfbrp_as(gr_i_i)%st_x_d
	gv_wtbbrp_as(gr_i_i)%st_y_d = gv_wtfbrp_as(gr_i_i)%st_y_d
	gv_wtbbrp_as(gr_i_i)%end_x_d = gv_wtfbrp_as(gr_i_i)%end_x_d
	gv_wtbbrp_as(gr_i_i)%end_y_d = gv_wtfbrp_as(gr_i_i)%end_y_d
	gv_wtbbrp_as(gr_i_i)%n_sections_i = gv_wtfbrp_as(gr_i_i)%n_sections_i
	gv_wtbbrp_as(gr_i_i)%sect_width_d = gv_wtfbrp_as(gr_i_i)%sect_width_d
	gv_wtbbrp_as(gr_i_i)%base_y_d = gv_wtfbrp_as(gr_i_i)%base_y_d
	gv_wtbbrp_as(gr_i_i)%sect_height_d = gv_wtfbrp_as(gr_i_i)%sect_height_d
	
	do gr_m_i = 1, gv_wtf_as(gr_i_i)%n_sections_i
		do gr_l_i = 1, gr_j_i
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(1)&
					= gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(1)&
					= MIN(gv_r_rpy_st_ad(1), &
					gv_wtfbrp_as(gr_i_i)%base_y_d + gr_l_i * &
					gv_wtfbrp_as(gr_i_i)%sect_height_d)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(1)&
					= 0.0
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(2)&
					= gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(2)&
					= MIN(gv_r_rpy_end_ad(1), & 
					gv_wtfbrp_as(gr_i_i)%base_y_d + gr_l_i * &
					gv_wtfbrp_as(gr_i_i)%sect_height_d)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(2)&
					= 0.0
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(3)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(2)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(3)&
					= gv_wtfbrp_as(gr_i_i)%base_y_d + (gr_l_i-1) * &
					gv_wtfbrp_as(gr_i_i)%sect_height_d 
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(3)&
					= 0.0
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(4)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(1)
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(4)&
					= gv_wtfbrp_as(gr_i_i)%base_y_d + (gr_l_i-1) * &
					gv_wtfbrp_as(gr_i_i)%sect_height_d
			gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(4)&
					= 0.0
		
		
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(1)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(1)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(1)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(1)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(1)&
					= -gv_fwidth_d
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(2)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(4)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(2)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(4)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(2)&
					= -gv_fwidth_d
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(3)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(3)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(3)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(3)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(3)&
					= -gv_fwidth_d
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%x_coord_ad(4)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					x_coord_ad(2)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%y_coord_ad(4)&
					= gv_wtfbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%&
					y_coord_ad(2)
			gv_wtbbrp_as(gr_i_i)%sect_as((gr_m_i-1)*gr_j_i+gr_l_i)%z_coord_ad(4)&
					= -gv_fwidth_d
		enddo
	enddo
enddo   
	
! write(*,*)"printing wall front above RP coords"
! do gr_i_i = 1,gv_NoOfTopSurf_i
	! write(10,*)"Surface",gr_i_i 
	! do gr_m_i =1,gv_wtfarp_as(gr_i_i)%n_sections_i
		! write(10,*)"WTFARP",gr_m_i 
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! write(10,*)gv_wtfarp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
	! enddo
! enddo
! stop
! write(10,*)"printing wall front below RP coords"
! do gr_i_i = 1,gv_NoOfTopSurf_i
	! write(10,*)"Roof",gr_i_i 
	! do gr_m_i =1,gv_wtbbrp_as(gr_i_i)%n_sections_i
		! write(10,*)"WTBBRP",gr_m_i 
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
		! write(10,*)gv_wtbbrp_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
	! enddo
! enddo
! stop

!---------------------------------------------
! Segment Coordinates for Floor Radiant plate faces     
!---------------------------------------------
do gr_i_i=1,gv_NoOfBotRadPlates_i
	gv_floor_RP_as(gr_i_i)%st_x_d=gv_f_rpx_st_ad(gr_i_i)
	gv_floor_RP_as(gr_i_i)%st_y_d=gv_f_rpy_st_ad(gr_i_i)
	gv_floor_RP_as(gr_i_i)%end_x_d=gv_f_rpx_end_ad(gr_i_i)
	gv_floor_RP_as(gr_i_i)%end_y_d=gv_f_rpy_end_ad(gr_i_i)
	
	! Calculating the width of sections within the surface
	gr_j_i=0
	gr_k_i=1
	do while(gr_k_i.eq.1)
		gr_j_i=gr_j_i+1
		gr_sect_w_d=(((gv_floor_RP_as(gr_i_i)%st_x_d- &
                gv_floor_RP_as(gr_i_i)%end_x_d)**2.0 + &
				(gv_floor_RP_as(gr_i_i)%st_y_d- &
                gv_floor_RP_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i   
    
		if(gr_sect_w_d.le.gv_sect_max_w_d) then
			gv_floor_RP_as(gr_i_i)%n_sections_i=gr_j_i
			gv_floor_RP_as(gr_i_i)%sect_width_d=gr_sect_w_d
			gr_k_i=0	  	 
		endif
	end do
	
	! Calculating the slope of the surface
	if(dabs(gv_floor_RP_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
		gr_slope_d=((gv_floor_RP_as(gr_i_i)%end_y_d)- &
				(gv_floor_RP_as(gr_i_i)%st_y_d))/ &
				((gv_floor_RP_as(gr_i_i)%end_x_d)- &
                (gv_floor_RP_as(gr_i_i)%st_x_d))
	else
		gr_slope_d=0.0
	endif
	
	! Coordinate points
	gr_theta_d=datan(gr_slope_d)
	do gr_m_i=1,gv_floor_RP_as(gr_i_i)%n_sections_i
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
				gv_floor_RP_as(gr_i_i)%st_x_d+ &
				(gr_m_i-1)*gv_floor_RP_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
				gv_floor_RP_as(gr_i_i)%st_y_d+ &
				(gr_m_i-1)*gv_floor_RP_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
				gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
				gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
				gv_floor_RP_as(gr_i_i)%st_x_d+ &
				gr_m_i*gv_floor_RP_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
				gv_floor_RP_as(gr_i_i)%st_y_d+ &
				gr_m_i*gv_floor_RP_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
				gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
				gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
		gv_floor_RP_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
	end do
enddo

write(*,*) "Inputs completed" 
return

end subroutine shm_grid