module mod_param
	implicit none
    !Here define type of parameters and arrays (as allocatable)
    !Parameters and arrays are then assigned a numerical value (i.e. they are initialized)
    !in the main file
	
	!Declaring parameters:	
    !integer, parameter :: dp = kind(0.0d0)
    real(8), parameter :: large_negative = -1d100 !-10.0d0**25
    real(8), parameter :: small   = 1.0d-8
	integer, parameter :: par_fortran = 1
    character(len=*), parameter :: savedir='output\'
    
    
    
        
    !Declaring grid dimensions
    integer :: na, nz
    
    !Integers:
	integer :: maxiter, iz, ia, i, j, exitflag_vfi, verbose, disp_on_screen

 
	!Real variables:
    real(8) :: beta, r, b, s, valtol, grid_max, t1, t2, c_max, optim, f_max
    real(8) :: check1, check2, check3, start, finish
    real(8) :: PZ(2,2), z_vals(2), PZ_iz(2)

    
    !Decalring allocatable integer variables:
    !integer, allocatable :: int_params(:)
   
    
                
    !Declaring allocatable real variables:
    real(8), allocatable :: asset_grid(:), c0(:,:), v0(:,:), TV(:,:), EVz(:), c_pol(:,:)
    real(8), allocatable :: a_pol(:,:), V(:,:), EV(:,:)
    
    !Communication variables
    real(8) :: z_today, a_today
    
    !=================================================
    ! OMP DIRECTIVES
    !=================================================
    !stuff to be available to (sub)procedures and which may need to be varied across OpenMP threads
    !$omp threadprivate(z_today,a_today,EVz)
    !!!$omp threadprivate(EVz)
    
    
   
    
       
end module mod_param
