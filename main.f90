!=============================================================================
! Fortran 90 code for solving income fluctuation problem
! Created by Alessandro Di Nola (University of Konstanz)
! Many thanks to Haomin Wang, Arnau V.E. and to the excellent website
! www.quantecon.org
! For an explanation of the income fluctuation problem and 
! a simple calibration, see:
! https://python.quantecon.org/ifp.html
!=============================================================================
! WHAT'S NEW: I created a new module, mod_numerical.f90, where I store
! all numerical routines. In the module mode_baselib.f90 I store instead
! procedures related to this specific project (utility function, etc.)
!=============================================================================
    
    
    program main
   
    use mod_param   !parameters and arrays as globals
    use mod_baselib
    use mod_numerical, only: grid
    
    
	implicit none
    
    
	
	!Time
	!call cpu_time(t1)
			
    !=================================================
    ! ECONOMIC PARAMETERS
    !=================================================
    
    !Assign values to parameters and grids
    beta     = 0.96d0
    r        = 0.03d0
    PZ(1,:)  = (/0.6d0, 0.40d0/)
    PZ(2,:)  = (/0.05d0, 0.95d0/)
    z_vals   = (/0.5d0, 1.0d0/)
    b        = 0d0 !borrowing limit
    grid_max = 4.0d0 !upper bound for assets
    s        = 1.0d0  !linear spacing
    
    
    !=================================================
    ! COMPUTATIONAL PARAMETERS
    !=================================================
    
	na             = 1000
    nz             = size(z_vals)
    valtol         = 1.0d-7 !tolerance for VFI
    maxiter        = 500 !max numb iters for VFI
    verbose        = 1 !0-1 display or not VFI errors
	disp_on_screen = 1 !!0-1 display on black screen some output

	
    !=================================================
    ! Initialize
    !=================================================
    
    allocate(v0(na,nz))
    allocate(c0(na,nz))
    allocate(asset_grid(na))
    
    !Inizialize grids calling an internal subroutine
    call initialize()
    
    !=================================================
    ! Print some output on the screen
    !=================================================
    if (disp_on_screen==1) then
        
!    write(*,'(a)') 'asset_grid = ' 
!    write(*,'(F12.4)') (asset_grid(i), i=1,na)
!    write(*,'(a)') '....'
!    write(*,'(a)') 'v0 = ' 
!    write(*,100) ((v0(i,j), j=1,nz), i=1,na)
!100 format (1X,F12.4,1X,F12.4)
!    write(*,'(a)') '....'
!    write(*,'(a)') 'c0 = '
!    write(*,100) ((c0(i,j), j=1,nz), i=1,na)
    
    endif 
    
    
    !=================================================
    ! Value function iteration
    !=================================================
    
    !Allocate
    allocate(V(na,nz))
    allocate(TV(na,nz))
    allocate(EVz(na))
    allocate(EV(na,nz))
    allocate(c_pol(na,nz))
    allocate(a_pol(na,nz))
    
    !outputs: TV (final value function), Policy functions (real variabl): c_pol, a_pol
    call sub_vfi(v0,TV,c_pol,a_pol)
    
    !check
    if (disp_on_screen==1) then
        
    !write(*,'(a)') 'TV = ' 
    !write(*,100) ((TV(i,j), j=1,nz), i=1,na)
    !write(*,'(a)') '....'
    !write(*,'(a)') 'apol = ' 
    !write(*,100) ((a_pol(i,j), j=1,nz), i=1,na)
        
        write(*,'(F12.6)') (a_pol(i,2),i=1,10)
    
    endif
    
    !Deallocate
    deallocate(EVz)
    deallocate(EV)
    deallocate(v0)    
    deallocate(c0)
    deallocate(V)
   
    
    !=================================================
    ! Writing out objects
    !=================================================
    
    !Allocate arrays for output
    !allocate(mu_vec(nagrid_dist*negrid*ntgrid) )
    
    !Writing routine
    call sub_write
    
    ! Deallocate
    deallocate(TV)
    deallocate(c_pol)
    deallocate(a_pol)
    

    !Time
    !call cpu_time(t2)
    
    !write(*,*) "============================================================="
    !write(*,*) 'FORTRAN executable runs for',real(t2-t1),'seconds.'
    !write(*,*) "============================================================="
    
	pause
    
    CONTAINS
    
    !=================================================
    ! INTERNAL PROCEDURES
    !=================================================
    
    subroutine initialize()
    
    
    
    implicit none
    
    !Subroutine "grid" is in the module "mod_baselib"
    call grid(asset_grid,-b,grid_max,s) !generate asset_grid
    
    
    do iz = 1,nz
        do ia = 1,na
            c_max = (1+r)*asset_grid(ia)+z_vals(iz)+b-small
            c0(ia,iz) = c_max
            v0(ia,iz) = (utilfun(c_max))/(1-beta)
        enddo
    enddo
    
    end subroutine initialize
    
    
        	
end program main
