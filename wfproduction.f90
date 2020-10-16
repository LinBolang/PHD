!this subroutine design for producing the wave function of proton
!In this part we include the basis enumeration, total dimension calculation and diagonalization the matrix of Hamiltonian

subroutine wfproduction(nmax1, nmax2, nmax3, Mj, Kt, method, massoffset, nmcutoff, xcutoff, nestate, ev, evector)

  use numbers
  use basis_info
  use hamiMatrix

  implicit none

  !input parameter
  integer:: nmax1, nmax2, nmax3, Mj, Kt, nestate, method
  double precision:: massoffset, nmcutoff, xcutoff

  !output parameter
  double precision:: wavefunction

  !private parameter
  integer ::  nmpstate1, nmpstate2, nmpstate3
  !double precision:: mu1, md1, mg1, ms1, mu2, md2, mg2, ms2, mu3, md3, mg3, ms3
  !!!!!!!!!! For the diagonalization !!!!!!!!!!!!!!!!!!!!!!!
  double precision, dimension(nestate):: ev
  double precision, dimension(nestate, dimtot):: evector
  double precision, dimension(:), allocatable:: initial_vector

  integer:: verbosesw, renorm_flag
  integer:: error
  integer:: test, test2
  double precision, dimension(nestate):: normalization

  ! set parameter
  ! mu1 = Mu1_DEFAULT
  ! md1 = Md1_DEFAULT
  ! ms1 = MS1_DEFAULT
  ! mg1 = mg1_default
  ! mu2 = Mu2_DEFAULT
  ! md2 = Md2_DEFAULT
  ! ms2 = MS2_DEFAULT
  ! mg2 = mg2_default
  ! mu3 = Mu3_DEFAULT
  ! md3 = Md3_DEFAULT
  ! ms3 = MS3_DEFAULT
  ! mg3 = mg3_default

  !basis enumeration

  !!!!!!!############### basis enumeration for three particle ##############!!!!!!!!!!!!!!!!

  allocate(s1dat1(dimtot1), s2dat1(dimtot1), s3dat1(dimtot1), n1dat1(dimtot1), m1dat1(dimtot1)&
    & ,n2dat1(dimtot1), m2dat1(dimtot1), n3dat1(dimtot1), m3dat1(dimtot1), k1dat1(dimtot1)&
    & ,k2dat1(dimtot1), k3dat1(dimtot1), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array'
  end if

  call mpstatesfill3p(nmax1, mj, kt, s1dat1, s2dat1, s3dat1, n1dat1, m1dat1, n2dat1, m2dat1, n3dat1, m3dat1, k1dat1, &
    & k2dat1, k3dat1, nmpstate1)

  if(nmpstate1 .ne. dimtot1)then              !check the basis enumeration and dimension of states
    Print *, 'get the wrong dimension', dimtot1, nmpstate1
  end if

  !!!!!!!################# basis enumeration for four particle ###############!!!!!!!!!!!!!!!!!!

  allocate(s1dat2(dimtot2), s2dat2(dimtot2), s3dat2(dimtot2), s4dat2(dimtot2), &
    & k1dat2(dimtot2), k2dat2(dimtot2), k3dat2(dimtot2), k4dat2(dimtot2), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array1'
  end if
  allocate(n1dat2(dimtot2), m1dat2(dimtot2), n2dat2(dimtot2), m2dat2(dimtot2), n3dat2(dimtot2)&
    & ,m3dat2(dimtot2), n4dat2(dimtot2), m4dat2(dimtot2), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array2'
  end if

  call mpstatesfill4p(nmax2, mj, kt, s1dat2, s2dat2, s3dat2, s4dat2, n1dat2, m1dat2, n2dat2, m2dat2, n3dat2, m3dat2, &
    & n4dat2, m4dat2, k1dat2, k2dat2, k3dat2, k4dat2, nmpstate2)

  if(nmpstate2 .ne. dimtot2)then              !check the basis enumeration and dimension of states
    Print *, 'get the wrong dimension', dimtot2, nmpstate2
  end if

  !!!!!!!!!!!!!!!! basis enumeration for five particle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(s1dat3(dimtot3), s2dat3(dimtot3), s3dat3(dimtot3), s4dat3(dimtot3), s5dat3(dimtot3), &
    & k1dat3(dimtot3), k2dat3(dimtot3), k3dat3(dimtot3), k4dat3(dimtot3), k5dat3(dimtot3), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array1'
  end if
  allocate(n1dat3(dimtot3), m1dat3(dimtot3), n2dat3(dimtot3), m2dat3(dimtot3), n3dat3(dimtot3)&
    & ,m3dat3(dimtot3), n4dat3(dimtot3), m4dat3(dimtot3), n5dat3(dimtot3), m5dat3(dimtot3), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array2'
  end if

  call mpstatesfill5p(nmax3, mj, kt, s1dat3, s2dat3, s3dat3, s4dat3, s5dat3, n1dat3, m1dat3, n2dat3, m2dat3, n3dat3, &
    & m3dat3, n4dat3, m4dat3, n5dat3, m5dat3, k1dat3, k2dat3, k3dat3, k4dat3, k5dat3, nmpstate3)

  if(nmpstate3 .ne. dimtot3)then              !check the basis enumeration and dimension of states
    Print *, 'get the wrong dimension', dimtot3, nmpstate3
  end if

  !!!!!!!!!!!!!!################# the hamiltonian matrix element of the three particle Fock sector ###############!!!!!!!!!!!!!!!!!!!!

  allocate(i_nzk(dimtot1+1), j_nzk(dimtot1*100), hamiltoniankineticvalue(dimtot1*100)&
    & ,i_nztsw(dimtot1+1), j_nztsw(dimtot1*100), hamiltonianinteractswvalue(dimtot1*100)&
    & ,i_nzlsw(dimtot1+1), j_nzlsw(dimtot1*300), hamiltonianinteractlongivalue(dimtot1*300), stat = error)

  if(error .ne. 0) then
    print *, 'cannot allocate hamiltonian3p'
  end if

  call hamiltoniankinetic3p(nmax2, Mj, Kt, i_nzk, j_nzk, hamiltoniankineticvalue)
  call hamiltonianinteractsw3p(nmax2, Mj, Kt, i_nztsw, j_nztsw, hamiltonianinteractswvalue)
  call hamiltonianinteractlongi3p(nmax2, Mj, Kt, i_nzlsw, j_nzlsw, hamiltonianinteractlongivalue)

  If(method .eq. 1) then

    allocate(i_nz_oge(dimtot1+1), j_nz_oge(dimtot1*2000), hamiltonianinteractogevalue(dimtot1*2000), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate OGE'
    end if

    call hamiltonianinteractOGE3p(nmax2, Mj, Kt, i_nz_oge, j_nz_oge, hamiltonianinteractogevalue)

  endif

  If(method .eq. 2) then

    allocate(i_nz_inst(dimtot1+1), j_nz_inst(dimtot1*1000), hamiltonianinteractinstvalue(dimtot1*1000), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate inst'
    end if

    call hamiltonianInstantaneous3p(nmax2, Mj, Kt, i_nz_inst, j_nz_inst, hamiltonianinteractinstvalue)

  endif

  !!!!!!!!!!!!!################## the hamiltonian matrix element of the four particle Fock sector ####################!!!!!!!!!!!!!!!!!

  If(method .eq. 2) then

    allocate(i_nzk2(2*dimtot2+1), j_nzk2(2*dimtot2*100), hamiltoniankineticvalue2(2*dimtot2*100), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate kinetic 4p'
    end if

    call hamiltoniankinetic4p(nmax2, Mj, Kt, i_nzk2, j_nzk2, hamiltoniankineticvalue2)

    allocate(i_nz_vcvqtqg(dimtot1+1), j_nz_vcvqtqg(dimtot1*1000), hamiltonianinteractvcvqtqg(dimtot1*1000), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate Vertex q to qg'
    end if

    call hamiltonianinteractqtoqg(nmax2, Mj, Kt, i_nz_vcvqtqg, j_nz_vcvqtqg, hamiltonianinteractvcvqtqg)

    allocate(i_nz_vcvgtqq(2*dimtot2+1), j_nz_vcvgtqq(2*dimtot2*1000), hamiltonianinteractvcvgtqq(2*dimtot2*100), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate Vertex g to qq'
    end if

    !call hamiltonianinteractgtoqqbar(nmax2, Mj, Kt, i_nz_vcvgtqq, j_nz_vcvgtqq, hamiltonianinteractvcvgtqq)

  endif

  !!!!!!!!!!!!!################## the hamiltonian matrix element of the five particle Fock sector ####################!!!!!!!!!!!!!!!!!

  allocate(i_nzk3(3*dimtot3+1), j_nzk3(3*dimtot3*100), hamiltoniankineticvalue3(3*dimtot3*100), stat = error)

  if(error .ne. 0) then
    print *, 'cannot allocate kinetic 5p'
  end if

  !call hamiltoniankinetic5p(nmax2, Mj, Kt, i_nzk3, j_nzk3, hamiltoniankineticvalue3)

  allocate(i_nztsw3(3*dimtot3+1), j_nztsw3(3*dimtot3*100), hamiltonianinteractswvalue3(3*dimtot3*100), stat = error)

  if(error .ne. 0) then
    print *, 'cannot allocate interaction sw 5p'
  end if

  !call hamiltonianinteractsw5p(nmax2, Mj, Kt, i_nztsw3, j_nztsw3, hamiltonianinteractswvalue3)

  allocate(i_nzlsw3(3*dimtot3+1), j_nzlsw3(3*dimtot3*1000), hamiltonianinteractlongivalue3(3*dimtot3*1000), stat = error)

  if(error .ne. 0) then
    print *, 'cannot allocate interaction sw 5p'
  end if

  !call hamiltonianinteractlongi5p(nmax2, Mj, Kt, i_nzlsw3, j_nzlsw3, hamiltonianinteractlongivalue3)

  If(method .eq. 1) then

    allocate(i_nz_oge3(3*dimtot3+1), j_nz_oge3(3*dimtot3*1000), hamiltonianinteractogevalue3(3*dimtot3*1000), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate interaction OGE 5p'
    end if

    !call hamiltonianinteractOGE5p(nmax2, Mj, Kt, i_nz_oge3, j_nz_oge3, hamiltonianinteractogevalue3)

    allocate(i_nz_oge3t5(dimtot1+1), j_nz_oge3t5(dimtot1*1000), hamiltonianinteractogevalue3t5(dimtot1*1000), stat = error)

    if(error .ne. 0) then
      print *, 'cannot allocate interaction OGE 3t5'
    end if

    !call hamiltonianinteractOGE3to5(nmax2, Mj, Kt, i_nz_oge3t5, j_nz_oge3t5, hamiltonianinteractogevalue3t5)

  endif

  allocate(i_nz_inst3t5(dimtot1+1), j_nz_inst3t5(dimtot1*1000), hamiltonianinteractinst3t5(dimtot1*1000), stat = error)

  if(error .ne. 0) then
    print *, 'cannot allocate interaction OGE 3t5'
  end if

  !call hamiltonianinstantaneous3to5(nmax2, Mj, Kt, method, i_nz_inst3t5, j_nz_inst3t5, hamiltonianinteractinst3t5)

  ! call cpu_time(hamiltonian_time)
  Print*,"hamiltonian finish"

  !!!!!!!!!!!!!!!!!!!!!!!!!!############ diagonalize the hamiltonian matrix ##################!!!!!!!!!!!!!!!!!!!!!!!!!!!

  verbosesw = 0
  renorm_flag = 0
  If(method .eq. 1) then
    dimtot = dimtot1+3*dimtot3
  else if(method .eq. 2) then
    !dimtot = dimtot1+2*dimtot2+3*dimtot3
    dimtot = dimtot1+2*dimtot2
    Print*, "dimtot = dimtot1+2*dimtot2   ----->:"
  else
    ! dimtot = dimtot1+3*dimtot3
    dimtot = dimtot1
  endif

  allocate(initial_vector(dimtot), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array'
  end if

  Print*, "dimtot=",dimtot, "||   dimtot1=",dimtot1, "||   dimtot2=",dimtot2, "||   dimtot3=",dimtot3

  call diag_arpack(nestate, dimtot, method, renorm_flag, initial_vector, verbosesw, ev, evector)

  normalization = 0.D0

  do test2 = 1, nestate
    !test2 = 4
    print*," "
    Print*,"eigenvalue"," || ", "eigenvector"," || "," No.nestate=",test2
    Print*,Sqrt(ev(test2)*(dble(Kt)+0.5D0))
    normalization = 0.D0
    do test = 1, dimtot1
      If(m1dat1(test)+m2dat1(test)+m3dat1(test).eq.0) then
        normalization(test2)=normalization(test2)+evector(test2, test)**2
      endif
    enddo
    print*,"normal1=",normalization(test2), "m1+m2+m3 = 0; normal1 = vector_dim1**2"
    normalization = 0.D0
    do test = 1, dimtot1
      If(abs(m1dat1(test)+m2dat1(test)+m3dat1(test)).eq.1) then
        normalization(test2)=normalization(test2)+evector(test2, test)**2
      endif
    enddo
    print*,"normal1=",normalization(test2), "|m1+m2+m3|=1; normal1 = vector_dim1**2"
    normalization = 0.D0
    do test = 1, dimtot1
      If(abs(m1dat1(test)+m2dat1(test)+m3dat1(test)).eq.2) then
        normalization(test2)=normalization(test2)+evector(test2, test)**2
      endif
    enddo
    print*,"normal1=",normalization(test2), "|m1+m2+m3|=2; normal1 = vector_dim1**2"
    normalization = 0.D0
    do test = 1, dimtot2
      If(m1dat2(test)+m2dat2(test)+m3dat2(test)+m4dat2(test).eq.0) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+test)**2&
          & +evector(test2, dimtot1+dimtot2+test)**2
      endif
    enddo
    print*,"normal2=",normalization(test2), "m1+m2+m3+m4 = 0; normal2 = vector_dim2color1**2+vector_dim2color2**2"
    normalization = 0.D0
    do test = 1, dimtot2
      If(abs(m1dat2(test)+m2dat2(test)+m3dat2(test)+m4dat2(test)).eq.1) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+test)**2&
          & +evector(test2, dimtot1+dimtot2+test)**2
        !print*,normalization(test2), evector(test2, test), test
      endif
    enddo
    print*,"normal2=",normalization(test2), "|m1+m2+m3+m4|=1; normal2 = vector_dim2color1**2+vector_dim2color2**2"
    normalization = 0.D0
    do test = 1, dimtot2
      If(abs(m1dat2(test)+m2dat2(test)+m3dat2(test)+m4dat2(test)).eq.2) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+test)**2&
          & +evector(test2, dimtot1+dimtot2+test)**2
      endif
    enddo
    print*,"normal2=",normalization(test2), "|m1+m2+m3+m4|=2; normal2 = vector_dim2color1**2+vector_dim2color2**2"
    normalization = 0.D0
    do test = 1, dimtot2
      If(abs(m1dat2(test)+m2dat2(test)+m3dat2(test)+m4dat2(test)).eq.3) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+test)**2&
          & +evector(test2, dimtot1+dimtot2+test)**2
      endif
    enddo
    print*,"normal2=",normalization(test2), "|m1+m2+m3+m4|=3; normal2 = vector_dim2color1**2+vector_dim2color2**2"
    normalization = 0.D0
    do test = 1, dimtot3
      If(abs(m1dat3(test)+m2dat3(test)+m3dat3(test)+m4dat3(test)+m5dat3(test)).eq.0) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+2*dimtot2+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+dimtot3+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+2*dimtot3+test)**2
      endif
    enddo
    print*,"normal3=",normalization(test2), "m1+m2+m3+m4+m5 = 0; normal2 = vector_dim3c1**2+vector_dim3c2**2+vector_dim3c3"
    normalization = 0.D0
    do test = 1, dimtot3
      If(abs(m1dat3(test)+m2dat3(test)+m3dat3(test)+m4dat3(test)+m5dat3(test)).eq.1) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+2*dimtot2+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+dimtot3+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+2*dimtot3+test)**2
      endif
    enddo
    print*,"normal3=",normalization(test2), "|m1+m2+m3+m4+m5|=1; normal2 = vector_dim3c1**2+vector_dim3c2**2+vector_dim3c3"
    normalization = 0.D0
    do test = 1, dimtot3
      If(abs(m1dat3(test)+m2dat3(test)+m3dat3(test)+m4dat3(test)+m5dat3(test)).eq.2) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+2*dimtot2+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+dimtot3+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+2*dimtot3+test)**2
      endif
    enddo
    print*,"normal3=",normalization(test2), "|m1+m2+m3+m4+m5|=2; normal2 = vector_dim3c1**2+vector_dim3c2**2+vector_dim3c3"
    normalization = 0.D0
    do test = 1, dimtot3
      If(abs(m1dat3(test)+m2dat3(test)+m3dat3(test)+m4dat3(test)+m5dat3(test)).eq.3) then
        normalization(test2)=normalization(test2)+evector(test2, dimtot1+2*dimtot2+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+dimtot3+test)**2&
          & +evector(test2, dimtot1+2*dimtot2+2*dimtot3+test)**2
      endif
    enddo
    print*,"normal3=",normalization(test2), "|m1+m2+m3+m4+m5|=3; normal2 = vector_dim3c1**2+vector_dim3c2**2+vector_dim3c3"
  enddo

  ! do test = 1, dimtot1
  !   If(abs(m1dat1(test)+m2dat1(test)+m3dat1(test)).eq.0) then
  !     Print*,evector(seqeigen, test), test, s1dat1(test), s2dat1(test), s3dat1(test)
  !   endif
  ! enddo

  deallocate(i_nzk, j_nzk, hamiltoniankineticvalue)
  deallocate(i_nztsw, j_nztsw, hamiltonianinteractswvalue)
  deallocate(i_nzlsw, j_nzlsw, hamiltonianinteractlongivalue)
  deallocate(i_nzk3, j_nzk3, hamiltoniankineticvalue3)
  deallocate(i_nztsw3, j_nztsw3, hamiltonianinteractswvalue3)
  deallocate(i_nzlsw3, j_nzlsw3, hamiltonianinteractlongivalue3)
  deallocate(i_nz_inst3t5, j_nz_inst3t5, hamiltonianinteractinst3t5)
  If(method .eq. 1) then
    deallocate(i_nz_oge, j_nz_oge, hamiltonianinteractogevalue)
    deallocate(i_nz_oge3, j_nz_oge3, hamiltonianinteractogevalue3)
    deallocate(i_nz_oge3t5, j_nz_oge3t5, hamiltonianinteractogevalue3t5)
  else if(method .eq. 2) then
    deallocate(i_nzk2, j_nzk2, hamiltoniankineticvalue2)
    deallocate(i_nz_vcvqtqg, j_nz_vcvqtqg, hamiltonianinteractvcvqtqg)
    deallocate(i_nz_vcvgtqq, j_nz_vcvgtqq, hamiltonianinteractvcvgtqq)
    deallocate(i_nz_inst, j_nz_inst, hamiltonianinteractinstvalue)
  endif

  return
end subroutine wfproduction

