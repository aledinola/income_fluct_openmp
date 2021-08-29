module mod_baselib

! This module mod_baselib contains procedures specific to this project
! Numerical routines are contained in another module
    
use mod_param  !Make calibrated variables known to these procedures

    
implicit none
    
!private !Variables and internal proc are not visible outside of this module
    
!public :: grid, interp1, utilfun
    
!Declaring parameters:	
   
    
!Module procedures
CONTAINS
    
  !******************************************************************************!
  
    pure function utilfun(c) result(F)
        implicit none

        !Declare inputs:
        real(8), intent(in) :: c ! consumption
        !Declare outputs:
        real(8) :: F 
        
        F = log(c)
        
    end function utilfun
    !******************************************************************************!
    
    !function rhs(c) 
    !    use mod_numerical, only: linint
    !    !Evaluates the RHS of the Bellman equation
    !    implicit none
    !    real(8), intent(in) :: c
    !    real(8) :: rhs, ap_next, v_next
    !
    !    !Next-period assets
    !    ap_next = (1+r)*a_today+z_today-c
    !    !linint(x,y,xi)
    !    v_next  = linint(asset_grid,EVz,ap_next)
    !    rhs = utilfun(c) + beta*v_next
    !
    !    !Reverse the sign for the minimizer
    !    !rhs = -rhs
    !
    !end function
    
   
    
 
    
    !******************************************************************************!
    !pure function prob_audit(phi,profit,theta,k,n,le,flag) result(pk)
    !     inputs
    !     phi: fraction of misreported business income
    !     profit: business income
    !    implicit none
    !
    !    declare inputs:
    !    real(dp1), intent(in) :: phi, profit, theta, k, n, le
    !    integer, intent(in) :: flag
    !    declare outputs:
    !    real(dp1) :: pk 
    !     declare locals:
    !    real(dp1) :: temp ! 
    !    
    !    select case (flag)
    !        
    !    case (1) !misreported income
    !        temp = phi*profit
    !        
    !    case (2) !f(theta,k,n)
    !        temp = theta*prodfun(k,le,n)
    !        
    !    case (3) !f(k,n)
    !        temp = prodfun(k,le,n)
    !    
    !    case (4) !k
    !        temp = k
    !        
    !    case (5) !n hired labor
    !        temp = n
    !        
    !    end select
    !    
    !    
    !    temp = (1.0d0-phi)*profit 
    !    pn_3 + (1-pn_3)./(1 + pn_1*exp(-pn_2*kgrid));
    !    pk   = pn_3 + (1.0d0-pn_3)/(1.0d0 + pn_1*exp(-pn_2*temp ))
    !    
    !   
    !end function prob_audit
    
   

end module mod_baselib
