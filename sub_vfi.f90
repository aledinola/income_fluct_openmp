subroutine sub_vfi()
	use mod_param
	use mod_baselib
    use mod_numerical, only: golden_method, linint
    use omp_lib
   
	implicit none
    
    ! Declare local variables for this subroutine
    real(8) :: diff, c_ub
    integer :: itervf
    !IMPORTANT: For openMP, always declare communication variables
    !in a module (mod_param.f90) and not in the subroutine that hosts
    !the internal procedure (to which these comm have to be passed to)
    !real(8) :: z_today, a_today
    
    !=================================================
    ! Value function iteration
    !=================================================
    
   
    write(*,*) "============================================================="
    write(*,*) "VALUE FUNCTION ITERATION..."
    write(*,*) "============================================================="
    
    write(*,*) "available threads = ", omp_get_max_threads()
    
   !$omp parallel
        write(*,*) 'Parallel hello to you!'
    !$omp end parallel
    
    pause

    diff = 10.0d0
    itervf = 0
    ! Initialize value function 
	V = v0
    
    !c_pol = 0.0d0 !Policy function for cons, dim: (na,nz)
    !TV   = 0.0d0 !New value function, dim: (na,nz)

    
    !call the clock to get the start time 
    start = omp_get_wtime()
    
    do while (diff > valtol .and. itervf <= maxiter)
        
        itervf = itervf+1
        
                
        !Expected future value function has dim (na,nz)
        !EV(a',z) = sum_{z'} V(a',z')*PZ(z,z')
        EV = matmul(V,transpose(PZ))
        
        !Loop over exogenous state "z"
        !$omp parallel private(ia,PZ_iz,c_ub,optim,f_max)
        !$omp do collapse(2)
        !EVz,a_today and z_today are global variables defined as "threadprivate"
        !in module mod_param.f90 (see also comment at the start of this file)
        do iz = 1,nz
            !Loop over endogenous state "a"
            do ia = 1,na
                !z_today and EVz could be pulled off the ia loop
                !but when using collapse they cannot
                z_today = z_vals(iz)
                EVz = EV(:,iz) ! this is EV(a') _given_ z
                
                a_today = asset_grid(ia)
                c_ub    = (1+r)*a_today+z_today+b
                !Call golden to max RHS(c), for c in [0,c_max]
                !optim is argmax, f_max is the max
                call golden_method(rhs, small, c_ub-small, optim, f_max,1.0d-5)
                
                c_pol(ia,iz) = optim
                TV(ia,iz)    = f_max
            
            enddo ! ia
            
        enddo ! iz
        !$omp enddo nowait
        !$omp end parallel
        
        !Compute sup-norm b/w two successive iterations
        diff = maxval(abs(TV-V)/(1+abs(V)))
        
        if (verbose==1) then
            write(*,*) "-----------------------------------"
            write(*,'(a,I4)') " iteration = ", itervf
            write(*,'(a,f0.7)') " norm_total = ", diff
            write(*,*) "-----------------------------------"
        endif
        
        !Update the guess for the value function             
        V = TV 
        
    enddo !while loop ends here
    
    !call the clock to get the end time and closes log file
    finish = omp_get_wtime()
    write(*,"(A,F15.1,A)") "Finished completely in ",finish-start," seconds."
    
    write(*,*) "============================================================="
    write(*,*) "VALUE FUNCTION FOUND"
    write(*,*) "============================================================="
   
    
    if (diff <= valtol) then
        write(*,*) 'VFI converged succesfully'
        exitflag_vfi = 0
    else
        write(*,*) 'VFI did not converge succesfully'
        exitflag_vfi = 1
    endif 
    
    !=================================================
    ! Compute policy functions
    !=================================================
    
    ! Policy function for next-period assets
    !Compute a_pol using c_pol
    a_pol = 0.0d0
    do iz= 1,nz
        a_pol(:,iz) = (1+r)*asset_grid+z_vals(iz)-c_pol(:,iz)
    enddo
    
    
    contains
    
    pure function rhs(c) 
        !Evaluates the RHS of the Bellman equation
        !This function could be put instead as a module function
        !in the module mod_baselib.f90
        implicit none
        real(8), intent(in) :: c
        real(8) :: rhs
        !Local variables
        real(8) :: ap_next, v_next
        !a_today, z_today, EVz are communication variables
        !For parallel computing, it is important that they 
        !are defined in a module (mod_param.f90) and NOT
        !in the VFI subroutine that hosts this internal function
        
        
        !Next-period assets
        ap_next = (1+r)*a_today+z_today-c
        !linint(x,y,xi)
        v_next  = linint(asset_grid,EVz,ap_next)
        rhs = utilfun(c) + beta*v_next
    
        !Reverse the sign for the minimizer
        !rhs = -rhs
    
    end function
    

    
end subroutine sub_vfi