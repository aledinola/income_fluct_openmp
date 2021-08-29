module mod_numerical
    
implicit none
    
!private !Variables and internal proc are not visible outside of this module
    
!public :: grid, interp1, utilfun
    
!Declaring parameters:	
   
    
!Module procedures
    CONTAINS
    
    SUBROUTINE grid(x,xmin,xmax,s)
    ! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
    ! s=1		linear spacing
    ! s>1		left skewed grid spacing with power s
    ! 0<s<1		right skewed grid spacing with power s
    ! s<0		geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
    ! s=-1		logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
    ! s=0		logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(OUT) :: x
	REAL(8), INTENT(IN) :: xmin,xmax,s
	REAL(8) :: c ! growth rate of grid subintervals for logarithmic spacing
	INTEGER :: n,i
	n=size(x)
	FORALL(i=1:n) x(i)=(i-1)/real(n-1,8)
	IF (s>0.0d0) THEN
		x=x**s*(xmax-xmin)+xmin
		IF (s==1.0d0) THEN
			PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
		ELSE
			PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
		END IF
	ELSE
		IF (s==-1.0d0) THEN
			c=xmax-xmin+1
!		ELSEIF (s==0.0_WP) THEN
!			IF (xmin>0.0_WP) THEN
!				c=xmax/xmin
!			ELSE
!				STOP 'grid: can not use logarithmic spacing for nonpositive values'
!			END IF
		ELSE
			c=-s
		END IF
		PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
		x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
	END IF
END SUBROUTINE grid


!******************************************************************************!
   PURE FUNCTION linint(x,y,xi)
! linear interpolation of function y on grid x at interpolation point xi
   !To make it pure, cannot use PRINT or STOP statements
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y
	REAL(8), INTENT(IN) :: xi
	REAL(8) :: linint
	REAL(8) :: a,b,d
	INTEGER :: n,i
	n=size(x)
	!IF (size(y)/=n) THEN
	!	PRINT *, 'linint: x and y must be of the same size'
	!	STOP 'program terminated by linint'
	!END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	!IF (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xi)/d
	b=(xi-x(i))/d
	linint=a*y(i)+b*y(i+1)
END FUNCTION linint
    !******************************************************************************!
   PURE FUNCTION locate(xx,x)
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: xx
	REAL(8), INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,il,im,iu
	n=size(xx)
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
		im=(iu+il)/2
		if (x >= xx(im)) then
			il=im
		else
			iu=im
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=il
	end if
END FUNCTION locate
    
    !!******************************************************************************!
    !subroutine bracket(x, xval, l, r)
    !    ! original code by  John Burkardt
    !    real(8), intent(in), dimension(:) :: x
    !    real(8), intent(in) :: xval
    !    integer, intent(out) :: l, r
    !    integer :: i, n
    !
    !    n=size(x)
    !    do i = 2, n - 1
    !       if ( xval < x(i) ) then
    !          l = i - 1
    !          r = i
    !          return
    !       end if
    !    end do
    !    l = n - 1
    !    r = n
    !end subroutine bracket
    !
    !!******************************************************************************!
    
  !===== GOLDEN SEARCH SECTION ====================================================================
  subroutine golden_method(f, a, b, x1, f1, mytol, mymaxit)
  ! Applies Golden-section search to search for the *maximum* of a function in the interval (a, b)
  !
  ! https://en.wikipedia.org/wiki/Golden-section_search
  ! Adapted to Fortran90 from: https://github.com/QuantEcon
    
    !---------------------------------------------------!
    !INPUTS
    interface
        function f(x)
        implicit none
        real(8), intent(in) :: x
        real(8) :: f
        end function f
    end interface
    real(8), intent(in) :: a, b
    !Some optional inputs
    integer,  optional :: mymaxit
    real(8), optional :: mytol
    !OUTPUTS
    real(8), intent(out) :: x1, f1
    !---------------------------------------------------!
    
    !Locals
    integer :: maxit, it
    real(8) :: tol, alpha1, alpha2, d, f2, x2, s
  
    !! Assign default value to maxit if not defined by user
    if (present(mymaxit)) then
        maxit = mymaxit
    else
        maxit = 1000
    end if
    
    ! Assign default value to tol if not defined by user
    if (present(mytol)) then
        tol = mytol
    else
        tol = 1.0d-6
    end if
  
    alpha1 = (3.d0 - sqrt(5.d0)) / 2.d0
    alpha2 = 1.d0 - alpha1
    d = b - a
    x1 = a + alpha1*d
    x2 = a + alpha2*d
    s = 1.d0
    f1 = f(x1)
    f2 = f(x2)
    d = alpha1*alpha2*d
  
    it = 0
  
    do while ((d.gt.tol).and.(it.lt.maxit))
        it = it + 1
        d = d*alpha2
  
        if (f2.gt.f1) then
            x1 = x2
            f1 = f2
            x2 = x1 + s*d
        else
            x2 = x1 - s*d
        end if
  
        s = sign(s, x2-x1)
        f2 = f(x2)
    end do
  
    if (it.ge.maxit) then
        print *, "Golden method: Maximum iterations exceeded"
    end if
  
    if (f2.gt.f1) then
        x1 = x2
        f1 = f2
    end if
  end subroutine golden_method

    

   

end module mod_numerical
