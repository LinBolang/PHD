subroutine wfproduction_single(nmax1,nmax2,nmax3,Mj,Kt,method,massoffset,nmcutoff,xcutoff,nestate,ev,evector)

    use numbers
    use basis_info
    use hamiMatrix

        implicit none

  !input parameter
    integer :: nmax1,nmax2,nmax3,Mj,Kt,nestate,method
    double precision :: massoffset,nmcutoff,xcutoff

  !output parameter
    double precision :: wavefunction

  !private parameter
    integer ::  nmpstate1,nmpstate2,nmpstate3
    !double precision :: mu1,md1,mg1,ms1,mu2,md2,mg2,ms2,mu3,md3,mg3,ms3
  !!!!!!!!!! For the diagonalization !!!!!!!!!!!!!!!!!!!!!!!
    double precision,dimension(nestate) :: ev
    double precision,dimension(nestate,dimtot) :: evector
    double precision,dimension(:),allocatable :: initial_vector

    integer :: verbosesw,renorm_flag
    integer :: error
    integer :: test,test2
    double precision,dimension(nestate) :: normalization

  ! set parameter
    ! mu1=Mu1_DEFAULT
    ! md1=Md1_DEFAULT
    ! ms1=MS1_DEFAULT
    ! mg1=mg1_default
    ! mu2=Mu2_DEFAULT
    ! md2=Md2_DEFAULT
    ! ms2=MS2_DEFAULT
    ! mg2=mg2_default
    ! mu3=Mu3_DEFAULT
    ! md3=Md3_DEFAULT
    ! ms3=MS3_DEFAULT
    ! mg3=mg3_default

  !basis enumeration

      !!!!!!!############### generate the mass counter term ##############!!!!!!!!!!!!!!!!


      call renormalization(nmax1,nmax2,nmax3,Mj,Kt,method,1,nestate,deltam_in,deltam_out)

      allocate(deltamq(1),deltamg(1),stat=error)
        if(error.ne.0) then
          print *, 'cannot allocate delta m'
        end if

      deltamq(1)=deltam_out(1)
      deltamg(1)=deltam_out(2)

      If (method.eq.1) then

        call dimtotal1p(nmax2,Mj,Kt,dimtot1)
        call dimtotal3p(nmax2,Mj,Kt,dimtot3)

        dimtot=dimtot1+3*dimtot3

      else if (method.eq.2) then

        call dimtotal1p(nmax2,Mj,Kt,dimtot1)
        call dimtotal2p(nmax2,Mj,Kt,dimtot2)
        call dimtotal3p(nmax2,Mj,Kt,dimtot3)

        dimtot=dimtot1+2*dimtot2+3*dimtot3

      else if (method.eq.3) then

        call dimtotal1p(nmax2,Mj,Kt,dimtot1)
        call dimtotal2q(nmax2,Mj,Kt,dimtot2)

        dimtot=dimtot1+3*dimtot2

      endif

      allocate(s1dat(dimtot),s2dat(dimtot),s3dat(dimtot),n1dat(dimtot),m1dat(dimtot)&
        & ,n2dat(dimtot),m2dat(dimtot),n3dat(dimtot),m3dat(dimtot),k1dat(dimtot)&
        & ,k2dat(dimtot),k3dat(dimtot),stat=error)
      if(error.ne.0) then
          print *, 'cannot allocate basis array'
      end if

      If (method.eq.1.or.method.eq.2) then

        call mpstatesfillquark(nmax2,mj,kt,method,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,&
          & n3dat,m3dat,k1dat,k2dat,k3dat,nmpstate1)

          if(nmpstate1.ne.dimtot)then              !check the basis enumeration and dimension of states
            Print *, 'get the wrong dimension', dimtot, nmpstate1
          end if

      else If (method.eq.3)

        call mpstatesfillPhoton(nmax1,mj,kt,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,&
            & k2dat,k3dat,nmpstate1)

          if(nmpstate1.ne.dimtot)then              !check the basis enumeration and dimension of states
            Print *, 'get the wrong dimension', dimtot, nmpstate1
          end if

      endif

  !!!!!!!!!!!!!!################# the hamiltonian matrix element for single particle ###############!!!!!!!!!!!!!!!!!!!!

    If (method.eq.1) then

      allocate(i_nzk(dimtot+1),j_nzk(dimtot*100),hamiltoniankineticvalue(dimtot*100),stat=error)

        if(error.ne.0) then
          print *, 'cannot allocate kinetic'
        end if

      allocate(i_nz_inst3t5(dimtot1+1),j_nz_inst3t5(dimtot1*1000),hamiltonianinteractinst3t5(dimtot1*1000),stat=error)

            if(error.ne.0) then
              print *, 'cannot allocate interaction OGE 3t5'
            end if

      allocate(i_nz_oge3t5(dimtot1+1),j_nz_oge3t5(dimtot1*1000),hamiltonianinteractogevalue3t5(dimtot1*1000),stat=error)

            if(error.ne.0) then
              print *, 'cannot allocate interaction OGE 3t5'
            end if

      call hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,deltamq(1),deltamg(1),i_nzk,j_nzk,hamiltoniankineticvalue)
      call hamiltonianinteractOGE1to3(nmax2,Mj,Kt,1,deltamq(1),i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
      call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,1,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)

    else If (method.eq.2) then

      allocate(i_nzk(dimtot+1),j_nzk(dimtot*100),hamiltoniankineticvalue(dimtot*100),stat=error)

        if(error.ne.0) then
          print *, 'cannot allocate kinetic'
        end if

      allocate(i_nz_vcvgtqq(dimtot+1),j_nz_vcvgtqq(dimtot*1000),hamiltonianinteractvcvgtqq(dimtot*100),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate Vertex g to qq'
          end if

      allocate(i_nz_vcvqtqg(dimtot+1),j_nz_vcvqtqg(dimtot*1000),hamiltonianinteractvcvqtqg(dimtot*1000),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate Vertex q to qg'
          end if

      allocate(i_nz_oge3t5(dimtot1+1),j_nz_oge3t5(dimtot1*1000),hamiltonianinteractogevalue3t5(dimtot1*1000),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate interaction OGE 3t5'
          end if
      
      call hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,deltamq(1),deltamg(1),i_nzk,j_nzk,hamiltoniankineticvalue)
      call hamiltonianinteractqtoqgSingle(nmax2,Mj,Kt,1,deltamq(1),i_nz_vcvqtqg,j_nz_vcvqtqg,hamiltonianinteractvcvqtqg)
      call hamiltonianinteractgtoqqSingle(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
      call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,1,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)

    else If (method.eq.2) then

      allocate(i_nzk(dimtot+1),j_nzk(dimtot*100),hamiltoniankineticvalue(dimtot*100),stat=error)

        if(error.ne.0) then
          print *, 'cannot allocate kinetic'
        end if

      allocate(i_nz_vcvgtqq(dimtot+1),j_nz_vcvgtqq(dimtot*100),hamiltonianinteractvcvgtqq(dimtot*100),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate Vertex g to qq'
          end if

      call hamiltoniankineticSinglephoton(nmax2,Mj,Kt,deltamg(1),i_nzk,j_nzk,hamiltoniankineticvalue)
      call hamiltonianinteractgtoqqphoton(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)

    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!############ diagonalize the hamiltonian matrix ##################!!!!!!!!!!!!!!!!!!!!!!!!!!!

    verbosesw=0
    renorm_flag=0
    If(method.eq.1) then
      dimtot=dimtot1+3*dimtot3
    else if(method.eq.2) then
      dimtot=dimtot1+2*dimtot2+3*dimtot3
      !dimtot=3*dimtot3
    else if(method.eq.3) then
      ! dimtot=dimtot1+3*dimtot3
      dimtot=dimtot1+2*dimtot2
    endif

    allocate(initial_vector(dimtot),stat=error)
      if(error.ne.0) then
        print *, 'cannot allocate basis array'
      end if

      Print*, dimtot,dimtot1,dimtot2,dimtot3

    If (method.eq.1.or.method.eq.2) then

      call single_diag_arpack(nestate,dimtot,method,2,initial_vector,verbosesw,ev,evector)

    else if (method.eq.3) then

      call single_diag_arpack(nestate,dimtot,method,1,initial_vector,verbosesw,ev,evector)

    endif

    normalization=0.D0

    do test2=1,nestate
    !test2=4
    Print*,"eigenvalue","  ", "eigenvector"
    Print*,Sqrt(ev(test2)*(dble(Kt)+0.5D0))
      ! do test=1,dimtot1
      !   If(m1dat1(test)+m2dat1(test)+m3dat1(test).eq.0) then
      !     Print*,evector(test2,test),test
      !   endif
      ! enddo
      do test=1,dimtot1
        normalization(test2)=normalization(test2)+evector(test2,test)**2
      enddo
      print*,"normal1=",normalization(test2)
      normalization=0.D0
         do test=dimtot1,dimtot1+2*dimtot2
           normalization(test2)=normalization(test2)+evector(test2,test)**2
         enddo
         print*,"normal2=",normalization(test2)
      normalization=0.D0
         do test=dimtot1+2*dimtot2,dimtot1+2*dimtot2+3*dimtot3
           normalization(test2)=normalization(test2)+evector(test2,test)**2
         enddo
         print*,"normal3=",normalization(test2)
    enddo

    deallocate(i_nzk,j_nzk,hamiltoniankineticvalue)
    deallocate(i_nz_inst3t5,j_nz_inst3t5,hamiltonianinteractinst3t5)
    If(method.eq.1) then
      deallocate(i_nz_oge3,j_nz_oge3,hamiltonianinteractogevalue3)
      deallocate(i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    else if(method.eq.2) then
      deallocate(i_nz_vcvqtqg,j_nz_vcvqtqg,hamiltonianinteractvcvqtqg)
      deallocate(i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
      deallocate(i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    else if(method.eq.3) then
      deallocate(i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
    endif

    return
  end subroutine wfproduction_single