subroutine output(nmax1,nmax2,nmax3,Mj,Kt,selectquark,loopnumber_loc,outputparameter)
    use numbers
    use basis_info
    implicit none

    integer :: nmax1,nmax2,nmax3,Mj,Kt,selectquark,loopnumber_loc,loopnumbers
    double precision, dimension(7,*) :: outputparameter
    !double precision, dimension(:), allocatable :: x,q2,functions1,functions2
    !integer :: error

    integer :: outputunit,i
    character(len=128) :: outputfile
    character(len=4) :: nmaxchar1,nmaxchar2,nmaxchar3,Mjchar,ktchar
    character(len=16) :: kappa4char,kappatchar,kappalchar,cpchar

    !if (select_distribution_function.eq.1) then

    !    call MPI_Reduce(loopnumber_loc,loopnumbers,1,MPI_INTEGER,MPI_SUM,root_process,comm,ierr)
    !    call MPI_Bcast(loopnumbers, 1, MPI_INTEGER, 0, comm, ierr)

    !else if (select_distribution_function.eq.2) then

        loopnumbers=loopnumber_loc
        
    !endif
    !Print*, loopnumber_loc,loopnumbers,myid

    !allocate(x(loopnumbers),q2(loopnumbers),functions1(loopnumbers),functions2(loopnumbers),stat=error)
    !    if(error.ne.0) then
    !        print *, 'cannot allocate basis array in output.f90'
    !    end if

    !call MPI_gather(outputparameter1,loopnumber_loc,MPI_DOUBLE_PRECISION,x,loopnumber_loc,MPI_DOUBLE_PRECISION,&
    !    & root_process,comm,ierr)
    !call MPI_gather(outputparameter2,loopnumber_loc,MPI_DOUBLE_PRECISION,q2,loopnumber_loc,MPI_DOUBLE_PRECISION,&
    !    & root_process,comm,ierr)
    !call MPI_gather(outputparameter3,loopnumber_loc,MPI_DOUBLE_PRECISION,functions1,loopnumber_loc,MPI_DOUBLE_PRECISION,&
    !    & root_process,comm,ierr)
    !call MPI_gather(outputparameter4,loopnumber_loc,MPI_DOUBLE_PRECISION,functions2,loopnumber_loc,MPI_DOUBLE_PRECISION,&
    !    & root_process,comm,ierr)


    write(nmaxchar1,'(I4)') nmax1
    write(nmaxchar2,'(I4)') nmax2
    write(nmaxchar3,'(I4)') nmax3
    write(Mjchar,'(I4)')   Mj
    write(ktchar,'(I4)')   Kt
    write(kappa4char,'(F16.2)')   kappa4
    write(kappatchar,'(F16.2)')   kappat1
    write(kappalchar,'(F16.2)')   kappal1
    write(cpchar,'(F16.4)') coupling_default3

    if (select_distribution_function.eq.2) then

        if(selectquark.eq.1) then

            outputfile="outputfile"//"_FF"//"_d_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.2) then

            outputfile="outputfile"//"_FF"//"_u_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"
        else if(selectquark.eq.3) then

            outputfile="outputfile"//"_FF"//"_s_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.4) then

            outputfile="outputfile"//"_FF"//"_sea1_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.5) then

            outputfile="outputfile"//"_FF"//"_sea2_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        endif
    else if (select_distribution_function.eq.3) then
        if(selectquark.eq.1) then

            outputfile="outputfile"//"_GPD"//"_d_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.2) then

            outputfile="outputfile"//"_GPD"//"_u_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"
        else if(selectquark.eq.3) then

            outputfile="outputfile"//"_FF"//"_s_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.4) then

            outputfile="outputfile"//"_FF"//"_sea1_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.5) then

            outputfile="outputfile"//"_FF"//"_sea2_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"
        endif

    else if (select_distribution_function.eq.4) then
        if(selectquark.eq.1) then

            outputfile="outputfile"//"_AFF"//"_d_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.2) then

            outputfile="outputfile"//"_AFF"//"_u_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        endif
    else if (select_distribution_function.eq.5) then
        if(selectquark.eq.1) then

            outputfile="outputfile"//"_tGPD"//"_d_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.2) then

            outputfile="outputfile"//"_tGPD"//"_u_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        endif
    else if (select_distribution_function.eq.6) then
        if(selectquark.eq.1) then

            outputfile="outputfile"//"_h"//"_d_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        else if(selectquark.eq.2) then

            outputfile="outputfile"//"_h"//"_u_quark_"//"nmax_"//trim(adjustl(nmaxchar1))//"_"//trim(adjustl(nmaxchar1))&
                & //"_"//trim(adjustl(nmaxchar1))//"_Mj_"//trim(adjustl(Mjchar))&
                & //"_kmax_"//trim(adjustl(ktchar))//"_kappa_"//trim(adjustl(kappa4char))//"_"//trim(adjustl(kappatchar))&
                & //"_"//trim(adjustl(kappalchar))//"_couplingconstant_"//trim(adjustl(cpchar))//".dat"

        endif
    endif


    !Print*, outputfile

        open (18,file=outputfile,status="replace")
        outputunit=18

        write(outputunit,*) loopnumbers

        do i=1,loopnumbers

            !write(outputunit,*) x(i),q2(i),functions1(i),functions2(i)
            write(outputunit,*) outputparameter(1,i),outputparameter(2,i),outputparameter(3,i),outputparameter(4,i)!,&
                !& outputparameter(5,i),outputparameter(6,i),outputparameter(7,i)

        enddo
    
    return
end subroutine output
