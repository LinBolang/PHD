subroutine generate_renorm_table(nmax1,nmax2,nmax3,Mj,Kt,method)
  use basis_info
  implicit none
  integer :: nmax1,nmax2,nmax3,Mj,Kt,method
  integer :: k,n,error,dim,particlenumber,i
  double precision :: deltam_in(2),deltam_out(2)

  dim=nmax2*Kt*3

  allocate(deltamq(dim),deltamg(dim),stat=error)
    if(error.ne.0) then
      print *, 'cannot allocate delta m'
    end if

  i=0
  do k=1,Kt
    do n=1,nmax2
      do particlenumber=1,3
        i=i+1
        call renormalization(n,n,n,Mj,K,method,particlenumber,1,deltam_in,deltam_out)
        deltamq(i)=deltam_out(1)
        deltamg(i)=deltam_out(2)
      enddo
    enddo
  enddo


  return
end subroutine generate_renorm_table

subroutine renormalization(nmax1,nmax2,nmax3,Mj,Kt,method,particlenumber,nestate,deltam_in,deltam_out)

    use numbers
    use basis_info
    use hamiMatrix

    implicit none

  !input parameter
    integer :: nmax1,nmax2,nmax3,Mj,Kt,nestate,method,particlenumber
    double precision :: massoffset,nmcutoff,xcutoff,deltam_in(2)

  !output parameter
    double precision :: wavefunction,deltam_out(2),mg

  !private parameter
    integer ::  nmpstate1,nmpstate2,nmpstate3
    !double precision :: mu1,md1,mg1,ms1,mu2,md2,mg2,ms2,mu3,md3,mg3,ms3
  !!!!!!!!!! For the diagonalization !!!!!!!!!!!!!!!!!!!!!!!
    double precision,dimension(nestate) :: ev
    double precision,dimension(:,:), allocatable :: evector
    double precision,dimension(:),allocatable :: initial_vector

    integer :: verbosesw,renorm_flag
    integer :: error
    integer :: test,test2
    double precision,dimension(nestate) :: normalization
    double precision :: energy,eigenv(2),renormass(3),mass(3)
    integer :: loopnumber,loopnumbermax

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

      !!!!!!!############### renormalization of the photon mass ##############!!!!!!!!!!!!!!!!

      If (method.eq.2.or.method.eq.3) then

        call dimtotal2q(Nmax2,Mj,Kt,dimtot2)
        call dimtotal1p(Nmax2,Mj,Kt,dimtot1)

        dimtot=1+dimtot2
        loopnumbermax=30

        allocate(s1dat(dimtot),s2dat(dimtot),s3dat(dimtot),n1dat(dimtot),m1dat(dimtot)&
          & ,n2dat(dimtot),m2dat(dimtot),n3dat(dimtot),m3dat(dimtot),k1dat(dimtot)&
          & ,k2dat(dimtot),k3dat(dimtot),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

        call mpstatesfillPhoton(nmax2,mj,kt,method,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,&
          & n3dat,m3dat,k1dat,k2dat,k3dat,nmpstate1)

        if(nmpstate1.ne.dimtot)then              !check the basis enumeration and dimension of states
          Print *, 'get the wrong dimension', dimtot, nmpstate1
        end if


        allocate(i_nzk(dimtot+1),j_nzk(dimtot*100),hamiltoniankineticvalue(dimtot*100),stat=error)

        if(error.ne.0) then
          print *, 'cannot allocate kinetic'
        end if

        allocate(i_nz_vcvgtqq(dimtot+1),j_nz_vcvgtqq(dimtot*100),hamiltonianinteractvcvgtqq(dimtot*100),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate Vertex g to qq'
          end if

        allocate(initial_vector(dimtot),evector(nestate,dimtot),stat=error)
          if(error.ne.0) then
            print *, 'cannot allocate basis array'
          end if

        !!!!!!!!!!!!!!!!!!!!!!! initial point for renormalization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        energy=0.D0
        renormass(1)=0.0D0

        call hamiltoniankineticSinglephoton(nmax2,Mj,Kt,renormass(1),i_nzk,j_nzk,hamiltoniankineticvalue)
        !call hamiltonianinteractgtoqqphoton(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)

        verbosesw=0
        initial_vector=0.D0
        initial_vector(1)=1.0D0
        renorm_flag=1

        call single_diag_arpack(nestate,dimtot,method,1,renorm_flag,initial_vector,verbosesw,ev,evector)

        eigenv(1)=ev(1)*Kt

        Print*,"ev=",ev(1)
        do test=1,dimtot
          Print*,evector(1,test)
        enddo

        ! renormass(2)=1.0D0
        ! call hamiltoniankineticSinglephoton(nmax2,Mj,Kt,renormass(2),i_nzk,j_nzk,hamiltoniankineticvalue)
        ! call hamiltonianinteractgtoqqphoton(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)

        ! verbosesw=0
        ! initial_vector=0.D0
        ! renorm_flag=0
        ! call single_diag_arpack(nestate,dimtot,method,1,renorm_flag,initial_vector,verbosesw,ev,evector)
        ! eigenv(2)=ev(1)*Kt

        ! Print*,eigenv(1),eigenv(2)
        !!!!!!!!!!!!!! iteration for the renormalization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! loopnumber=0
        ! Do while (abs(energy-eigenv(2)).gt.epsilon(abs(energy-eigenv(2))).and.abs(eigenv(2)-eigenv(1))&
        !   & .gt.epsilon(abs(eigenv(2)-eigenv(1))).and.loopnumber.le.loopnumbermax)

        !   renormass(3)= (renormass(2)+((energy-eigenv(2))/(eigenv(2)-eigenv(1)))*(renormass(2)-renormass(1)))

        !   call hamiltoniankineticSinglephoton(nmax2,Mj,Kt,sqrt(abs(renormass(3))),i_nzk,j_nzk,hamiltoniankineticvalue)
        !   call hamiltonianinteractgtoqqphoton(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)

        !   verbosesw=0
        !   initial_vector=0.D0
        !   renorm_flag=0
        !   call single_diag_arpack(nestate,dimtot,method,1,renorm_flag,initial_vector,verbosesw,ev,evector)

        !   eigenv(1)=eigenv(2)
        !   eigenv(2)=ev(1)*Kt
        !   renormass(1)=renormass(2)
        !   renormass(2)=renormass(3)
        !   loopnumber=loopnumber+1

        !   Print*,eigenv(2),renormass(2)

        ! enddo

        ! deltam_out(2)=renormass(2)

        deallocate(i_nzk,j_nzk,hamiltoniankineticvalue,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
        deallocate(s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat)
        deallocate(initial_vector)

      endif

      !!!!!!!################# renormalization of the quark mass ###############!!!!!!!!!!!!!!!!!!

      ! If (method.eq.1.or.method.eq.2) then
      If (method.eq.4) then

        call dimtotal1p(nmax2,Mj,Kt,dimtot1)
        call dimtotal2p(nmax2,Mj,Kt,dimtot2)
        call dimtotal3p(nmax2,Mj,Kt,dimtot3)
        mass(1)=Mu1_DEFAULT
        mass(2)=Md1_DEFAULT
        mass(3)=Ms1_DEFAULT

          If(method.eq.1) then
            dimtot=dimtot1+3*dimtot3
          else if(method.eq.2) then
            dimtot=dimtot1+2*dimtot2+3*dimtot3
            !dimtot=3*dimtot3
          else
            ! dimtot=dimtot1+3*dimtot3
            dimtot=dimtot1
          endif

        allocate(s1dat(dimtot),s2dat(dimtot),s3dat(dimtot),n1dat(dimtot),m1dat(dimtot)&
          & ,n2dat(dimtot),m2dat(dimtot),n3dat(dimtot),m3dat(dimtot),k1dat(dimtot)&
          & ,k2dat(dimtot),k3dat(dimtot),stat=error)
        if(nmpstate1.ne.dimtot)then              !check the basis enumeration and dimension of states
          Print *, 'get the wrong dimension', dimtot, nmpstate1
        end if

        call mpstatesfillquark(nmax2,mj,kt,method,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat,nmpstate1)

        allocate(i_nzk(dimtot+1),j_nzk(dimtot*100),hamiltoniankineticvalue(dimtot*100),stat=error)

          if(error.ne.0) then
            print *, 'cannot allocate kinetic'
          end if

        allocate(initial_vector(dimtot),stat=error)
            if(error.ne.0) then
              print *, 'cannot allocate basis array'
            end if

        If (method.eq.2) then

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

        else if (method.eq.1) then

          allocate(i_nz_inst3t5(dimtot1+1),j_nz_inst3t5(dimtot1*1000),hamiltonianinteractinst3t5(dimtot1*1000),stat=error)

            if(error.ne.0) then
              print *, 'cannot allocate interaction OGE 3t5'
            end if

          allocate(i_nz_oge3t5(dimtot1+1),j_nz_oge3t5(dimtot1*1000),hamiltonianinteractogevalue3t5(dimtot1*1000),stat=error)

            if(error.ne.0) then
              print *, 'cannot allocate interaction OGE 3t5'
            end if

        endif

      !!!!!!!!!!!!!!!!!!!!!!! initialized the normalization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        renormass(1)=deltam_in(1)
        mg=deltam_out(2)

        call hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,renormass(1),mg,i_nzk,j_nzk,hamiltoniankineticvalue)

        If(method.eq.2) then

          call hamiltonianinteractqtoqgSingle(nmax2,Mj,Kt,particlenumber,renormass(1),i_nz_vcvqtqg,j_nz_vcvqtqg,&
            & hamiltonianinteractvcvqtqg)
          call hamiltonianinteractgtoqqSingle(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
          call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
            & hamiltonianinteractogevalue3t5)

        else If(method.eq.1) then

          call hamiltonianinteractOGE1to3(nmax2,Mj,Kt,particlenumber,renormass(1),i_nz_oge3t5,j_nz_oge3t5&
            & ,hamiltonianinteractogevalue3t5)
          call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
            & hamiltonianinteractogevalue3t5)

        endif

        verbosesw=0
        initial_vector=0.D0
        renorm_flag=0
        call single_diag_arpack(nestate,dimtot,method,2,renorm_flag,initial_vector,verbosesw,ev,evector)
        eigenv(1)=ev(1)*(dble(Kt)+0.5d0)


        renormass(2)=sqrt((mass(particlenumber)+renormass(1))**2+mass(particlenumber)**2-eigenv(1))

        call hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,renormass(2),mg,i_nzk,j_nzk,hamiltoniankineticvalue)

        If(method.eq.2) then

          call hamiltonianinteractqtoqgSingle(nmax2,Mj,Kt,particlenumber,renormass(2),i_nz_vcvqtqg,j_nz_vcvqtqg,&
            & hamiltonianinteractvcvqtqg)
          call hamiltonianinteractgtoqqSingle(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
          call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
            & hamiltonianinteractogevalue3t5)

        else If(method.eq.1) then

          call hamiltonianinteractOGE1to3(nmax2,Mj,Kt,particlenumber,renormass(2),i_nz_oge3t5,j_nz_oge3t5,&
            & hamiltonianinteractogevalue3t5)
          call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
            & hamiltonianinteractogevalue3t5)

        endif

        verbosesw=0
        initial_vector=0.D0
        renorm_flag=0
        call single_diag_arpack(nestate,dimtot,method,2,renorm_flag,initial_vector,verbosesw,ev,evector)
        eigenv(2)=ev(1)*(dble(Kt)+0.5d0)

        do loopnumber=1,loopnumbermax

          If(abs(eigenv(2)-mass(particlenumber)**2).lt.epsilon(abs(eigenv(2)-mass(particlenumber)**2))) then
            exit
          endif

          renormass(3)=Sqrt((mass(particlenumber)+renormass(1))**2+((mass(particlenumber)+renormass(2))**2-&
            & (mass(particlenumber)+renormass(1))**2)/(eigenv(2)-eigenv(1))*(mass(particlenumber)**2-eigenv(1)))&
            & -mass(particlenumber)

          call hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,renormass(3),mg,i_nzk,j_nzk,hamiltoniankineticvalue)

          If(method.eq.2) then

            call hamiltonianinteractqtoqgSingle(nmax2,Mj,Kt,particlenumber,renormass(3),i_nz_vcvqtqg,j_nz_vcvqtqg,&
              & hamiltonianinteractvcvqtqg)
            call hamiltonianinteractgtoqqSingle(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
            call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
              & hamiltonianinteractogevalue3t5)

          else If(method.eq.1) then

            call hamiltonianinteractOGE1to3(nmax2,Mj,Kt,particlenumber,renormass(3),i_nz_oge3t5,j_nz_oge3t5,&
              & hamiltonianinteractogevalue3t5)
            call hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,&
              & hamiltonianinteractogevalue3t5)

          endif

          verbosesw=0
          initial_vector=0.D0
          renorm_flag=0
          call single_diag_arpack(nestate,dimtot,method,2,renorm_flag,initial_vector,verbosesw,ev,evector)

          renormass(1)=renormass(2)
          renormass(2)=renormass(3)
          eigenv(1)=eigenv(2)
          eigenv(2)=ev(1)*(dble(Kt)+0.5d0)


        enddo

          deltam_out(1)=renormass(2)

          deallocate(i_nzk,j_nzk,hamiltoniankineticvalue)
          deallocate(s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat)
          deallocate(initial_vector)
          If(method.eq.1) then
            deallocate(i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
            deallocate(i_nz_inst3t5,j_nz_inst3t5,hamiltonianinteractinst3t5)
          else if(method.eq.2) then
            deallocate(i_nz_inst3t5,j_nz_inst3t5,hamiltonianinteractinst3t5)
            deallocate(i_nz_vcvqtqg,j_nz_vcvqtqg,hamiltonianinteractvcvqtqg)
            deallocate(i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
          endif
      endif

    return
  end subroutine renormalization