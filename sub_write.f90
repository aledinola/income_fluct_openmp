subroutine sub_write
	use mod_param
	implicit none
	
    ! Write the outputs as txt files (to be imported into matlab) 
    ! METHOD 1: You can first vectorize all arrays and then write into txt using a single loop
    ! METHOD 2: You can directly write into txt the vectorized array using multiple loops (starting from the innermost)
    ! Vectorize all arrays
    !--------- EXAMPLE METHOD 1 --------------------------------!
 !   mu_vec           = reshape(mu, [nagrid_dist*negrid*ntgrid] )
 !   
 !   open(15, file='mu.txt', status='replace')
	!do i=1,nagrid*negrid*ntgrid
	!write(15,*) mu_vec(i)
 !   enddo
	!close(15)
 !   !Note: Method 1 requires that you define array mu_vec in advance !!!!!!!
    
    !--------- END EXAMPLE METHOD 1 --------------------------------!
	
	! 
	call create_directory( savedir )
    
    open(15,file=savedir//'constants.txt', status='replace')
    write(15,*) beta
    write(15,*) r
    write(15,*) b
    write(15,*) valtol
    write(15,*) maxiter
    close(15)
    
    open(15, file=savedir//'z_vals.txt', status='replace') 
    do i=1,nz
        write(15,*) z_vals(i)
    enddo
    close(15)
    
    open(15, file=savedir//'asset_grid.txt', status='replace') 
    do i=1,na
        write(15,*) asset_grid(i)
    enddo
    close(15)
    
    open(15, file=savedir//'TV.txt', status='replace')
    !Vectorize TV in column order (as using (:) in Matlab) 
    do j=1,nz
        do i=1,na
            write(15,*) TV(i,j)
        enddo
    enddo
    close(15)
    
    open(15, file=savedir//'c_pol.txt', status='replace')
	do j=1,nz
        do i=1,na
            write(15,*) c_pol(i,j)
        enddo
    enddo
    close(15)
    
    open(15, file=savedir//'a_pol.txt', status='replace')
	do j=1,nz
        do i=1,na
            write(15,*) a_pol(i,j)
        enddo
    enddo
    close(15)
    
   
              
end subroutine sub_write

subroutine create_directory( newDirPath )
    ! Author:  Jess Vriesema
    ! Date:    Spring 2011
    ! Purpose: Creates a directory at ./newDirPath

    implicit none

    character(len=*), intent(in) :: newDirPath
    character(len=256)           :: mkdirCmd
    logical                      :: dirExists

    ! Check if the directory exists first
!   inquire( file=trim(newDirPath)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
    inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran


    if (dirExists) then
!      write (*,*) "Directory already exists: '"//trim(newDirPath)//"'"
    else
        mkdirCmd = 'mkdir -p '//trim(newDirPath)
        write(*,'(a)') "Creating new directory: '"//trim(mkdirCmd)//"'"
        call system( mkdirCmd )
    endif
end subroutine create_directory
