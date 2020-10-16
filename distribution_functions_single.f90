! subroutine distribution_functions_GPD(nmax, Mj, Kt, mu, md, ms, mg, mp, nestate, selectquark, deltamax, lag)
!     use numbers
!     use basis_info

!     implicit none

!     integer:: nmax, Mj, Kt, selectquark, nestate
!     double precision:: deltamax, mu, md, ms, lag, mg, mp

!     double precision, dimension(:), allocatable:: evector, ev_loc
!     double precision:: ev
!     double precision, dimension(:,:), allocatable:: evector_loc!,distributionfunctions
!     double precision, dimension(:,:), allocatable:: distributionfunctions1, distributionfunctions2
!     !double precision, dimension(:,:), allocatable:: distributionfunctionskt1, distributionfunctionskt2
!     double precision, dimension(:,:), allocatable:: distributionfunctions_GPD
!     integer:: functionsdimension, snm, loopnumber1, loopnumber_loc
!     integer, dimension(9):: sums
!     double precision:: xcutoff, nmcutoff, massoffset, dx, dy, step
!     integer:: error, i, j, dimtot1, dimtot2, delta
!     double precision:: normalization
!     integer:: locstatesnum, selectev, numeigen
!     integer:: test


!     massoffset = 0.0D-2
!     nmcutoff = 0.0D0
!     xcutoff = 0.0D0
!     step = 0.1

!     sums(1)=0
!     call snmlen2(Nmax, Mj, snm)
!     call dddnmlen2(Nmax, Mj, sums(2))
!     call uddnmlen2(Nmax, Mj, sums(3))
!     call dudnmlen2(Nmax, Mj, sums(4))
!     call uudnmlen2(Nmax, Mj, sums(5))
!     call ddunmlen2(Nmax, Mj, sums(6))
!     call udunmlen2(Nmax, Mj, sums(7))
!     call duunmlen2(Nmax, Mj, sums(8))
!     call uuunmlen2(Nmax, Mj, sums(9))

! !!$$$$$$$$ which quark can be calculate by using different parameter selectquark $$$$$$$$$!!!!!
! !!$$$$$$$$ when the selectquark = 1, calculating the d quark of proton in the script $$$$!!!!!
! !!$$$$$$$$ when the selectquark = 2, calculating the u quark of proton in the script $$$$!!!!!
! !!$$$$$$$$ when the selectquark = 3, calculating the s quark of proton in the script $$$$!!!!!

! !!$$$$$$ decide the dimension of distribution function on each process $$$$$!!!!!!!!!!!!!

!         functionsdimension = Kt*(floor(deltamax/step)+1)

!         Print*,"the dimensions of matrix of distribution functions = ",functionsdimension

!         allocate(distributionfunctions1(7, kt), distributionfunctions2(7, kt), stat = error)
!             if(error .ne. 0) then
!                 print *, 'cannot allocate basis array'
!             end if

!         allocate(distributionfunctions_GPD(7, functionsdimension), stat = error)
!             if(error .ne. 0) then
!                 print *, 'cannot allocate basis array'
!             end if

! !!$$$$$$ wave function production $$$$$$$$$$$$$$$$$$$!!!!!!!!

!         call dimtotal(nmax, nmax, Mj, Kt, dimtot1)

!         call seperation(dimtot1)

!         allocate(ev_loc(nestate), evector_loc(nloc, nestate), stat = error)
!             if(error .ne. 0) then
!                 print *, 'cannot allocate basis array'
!             end if


!         call wfproduction(nmax, nmax, Mj, Kt, massoffset, nmcutoff, xcutoff, lag, mu, md, ms, mg, mp, nestate, &
!             & ev_loc, evector_loc)


! !!$$$$$$$$ wave function and basis transfer $$$$$$$$$$$$$$$$!!!!!!!!!

!         allocate(evector(nloc), stat = error)
!             if(error .ne. 0) then
!                 print*, "cannot allocate basis array"
!             endif

!         numeigen = seqeigen
!         do selectev = 1, nloc
!             evector(selectev)=evector_loc(selectev, numeigen)
!             ev = ev_loc(numeigen)
!         enddo
!         call MPI_Type_Extent(MPI_Double_precision, mpstates_size, ierr)
!         mpstates_window_size = nloc*mpstates_size
!         !print*,mpstates_window_size

!         call MPI_Win_create(evector, mpstates_window_size, mpstates_size, MPI_INFO_NULL, comm, evector_win, ierr)
!         call MPI_Win_fence(0, evector_win, ierr)

!         call mpi_barrier(comm, ierr)

!         if(myid .eq. root_process) then
!             print*,"window for vector has been created"
!         endif

!         allocate(readk1(nloc+1), readk2(nloc+1), readk3(nloc+1), readn1(nloc+1), &
!             &readm1(nloc+1), readn2(nloc+1), readm2(nloc+1), readn3(nloc+1), readm3(nloc+1), readvector(nloc+1), stat = error)
!             if(error .ne. 0) then
!                 print *, 'cannot allocate basis array in longitudinal confine potential'
!             end if
!         allocate(read1k1(nloc+1), read1k2(nloc+1), read1k3(nloc+1), read1n1(nloc+1), &
!             &read1m1(nloc+1), read1n2(nloc+1), read1m2(nloc+1), read1n3(nloc+1), read1m3(nloc+1), read1vector(nloc+1), stat = error)
!             if(error .ne. 0) then
!                 print *, 'cannot allocate basis array in longitudinal confine potential'
!             end if

! !!!$$$$$$$$$$$$$$$$ calculating distribution function $$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!

!         selectquark = 1
!         loopnumber_loc = 0

!         do delta = 0, floor(deltamax/step)

!             dy = 0.0D0
!             dx = delta*step
!             !dx = 1.8

!             call GPD_for_d_quark(nmax, Mj, Kt, sums, snm, nestate, dimtot1, evector, ev, dx, dy, distributionfunctions1)

!             do i = 1, kt

!                 loopnumber_loc = loopnumber_loc+1

!                 distributionfunctions_GPD(1, loopnumber_loc)=distributionfunctions1(1, i)
!                 distributionfunctions_GPD(2, loopnumber_loc)=distributionfunctions1(2, i)
!                 distributionfunctions_GPD(3, loopnumber_loc)=distributionfunctions1(3, i)!*distributionfunctions1(1, i)
!                 distributionfunctions_GPD(4, loopnumber_loc)=distributionfunctions1(4, i)
!                 distributionfunctions_GPD(5, loopnumber_loc)=distributionfunctions1(5, i)
!                 distributionfunctions_GPD(6, loopnumber_loc)=distributionfunctions1(6, i)
!                 distributionfunctions_GPD(7, loopnumber_loc)=distributionfunctions1(7, i)

!             enddo

!         enddo

!         !call MPI_Reduce(loopnumber_loc, loopnumber1, 1, MPI_INTEGER, MPI_SUM, root_process, comm, ierr)

!         !If(myid .eq. 0) then
!         !    print*,"loopnumber=",loopnumber_loc, loopnumber1
!         !endif

!         call output(nmax, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_GPD)

! !!!!!!!!test the normalization for d quark. and the normalization should be equal to 1 !!!!!!!!

!         !do j = 1, kt-3

!         !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

!         !enddo

!         !Print*,normalization

!     !else if (selectquark .eq. 2) then
!     selectquark = 2
!     loopnumber_loc = 0

!         do delta = 0, floor(deltamax/step)

!             dy = 0.0D0
!             dx = delta*step

!             call GPD_for_u_quark(nmax, Mj, Kt, sums, snm, nestate, dimtot1, evector, ev, dx, dy, distributionfunctions1)

!             do i = 1, kt

!                 loopnumber_loc = loopnumber_loc+1

!                 distributionfunctions_GPD(1, loopnumber_loc)=distributionfunctions1(1, i)
!                 distributionfunctions_GPD(2, loopnumber_loc)=distributionfunctions1(2, i)
!                 distributionfunctions_GPD(3, loopnumber_loc)=distributionfunctions1(3, i)!*distributionfunctions1(1, i)
!                 distributionfunctions_GPD(4, loopnumber_loc)=distributionfunctions1(4, i)
!                 distributionfunctions_GPD(5, loopnumber_loc)=distributionfunctions1(5, i)
!                 distributionfunctions_GPD(6, loopnumber_loc)=distributionfunctions1(6, i)
!                 distributionfunctions_GPD(7, loopnumber_loc)=distributionfunctions1(7, i)

!             enddo

!         enddo

!         call output(nmax, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_GPD)

!         !!!!!!!!!!!test the normalization for u quark. And the normalization should be equal to 2 !!!!!!!!

!         !do j = 1, kt-1

!         !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

!         !enddo


!     selectquark = 3
!     loopnumber_loc = 0

!         do delta = 0, floor(deltamax/step)

!             dy = 0.0D0
!             dx = delta*step

!             call GPD_for_s_quark(nmax, Mj, Kt, sums, snm, nestate, dimtot1, evector, ev, dx, dy, distributionfunctions1)

!             do i = 1, kt

!                 loopnumber_loc = loopnumber_loc+1

!                 distributionfunctions_GPD(1, loopnumber_loc)=distributionfunctions1(1, i)
!                 distributionfunctions_GPD(2, loopnumber_loc)=distributionfunctions1(2, i)
!                 distributionfunctions_GPD(3, loopnumber_loc)=distributionfunctions1(3, i)!*distributionfunctions1(1, i)
!                 distributionfunctions_GPD(4, loopnumber_loc)=distributionfunctions1(4, i)
!                 distributionfunctions_GPD(5, loopnumber_loc)=distributionfunctions1(5, i)
!                 distributionfunctions_GPD(6, loopnumber_loc)=distributionfunctions1(6, i)
!                 distributionfunctions_GPD(7, loopnumber_loc)=distributionfunctions1(7, i)

!             enddo

!         enddo

!         call output(nmax, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_GPD)

!         !!!!!!!!!!!test the normalization for u quark. And the normalization should be equal to 2 !!!!!!!!

!         !do j = 1, kt-1

!         !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

!         !enddo

!         !Print*,normalization


!         call MPI_Win_Free(s1_win, ierr)
!         call MPI_Win_Free(s2_win, ierr)
!         call MPI_Win_Free(s3_win, ierr)
!         call MPI_Win_Free(n1_win, ierr)
!         call MPI_Win_Free(m1_win, ierr)
!         call MPI_Win_Free(n2_win, ierr)
!         call MPI_Win_Free(m2_win, ierr)
!         call MPI_Win_Free(n3_win, ierr)
!         call MPI_Win_Free(m3_win, ierr)
!         call MPI_Win_Free(k1_win, ierr)
!         call MPI_Win_Free(k2_win, ierr)
!         call MPI_Win_Free(k3_win, ierr)
!         call MPI_Win_Free(evector_win, ierr)
!         deallocate(s1dat, s2dat, s3dat, n1dat, m1dat, n2dat, m2dat, n3dat, m3dat, k1dat, k2dat, k3dat, distributionfunctions1, &
!             & distributionfunctions_GPD, ev_loc, evector)
!         deallocate(readk1, readk2, readk3, readn1, readm1, readn2, readm2, readn3, readm3, readvector)
!         deallocate(read1k1, read1k2, read1k3, read1n1, read1m1, read1n2, read1m2, read1n3, read1m3, read1vector)

!     !endif

!     return
! end subroutine distribution_functions_GPD

subroutine distribution_functions_FF(nmax1, nmax2, nmax3, Mj, Kt, nestate, selectquark, deltamax, lag, method)
  use numbers
  use basis_info
  use colorfactor

  implicit none

  integer:: nmax1, nmax2, nmax3, Mj, Kt, selectquark, nestate, method
  double precision:: deltamax, mu, md, ms, lag, mg, mp

  double precision, dimension(:), allocatable:: ev1, evector
  double precision, dimension(:,:), allocatable:: evector1
  double precision, dimension(:,:), allocatable:: distributionfunctions1, distributionfunctions2, distributionfunctions3
  !double precision, dimension(:,:), allocatable:: distributionfunctionskt1, distributionfunctionskt2
  double precision, dimension(:,:), allocatable:: distributionfunctions_FF
  integer:: functionsdimension, snm, loopnumber1, loopnumber_loc
  integer:: sums(9), SumS5p(33), SumS4p(17)
  double precision:: xcutoff, nmcutoff, massoffset, dx, dy, step
  integer:: error, i, j, delta
  double precision:: normalization, formfactor_D, formfactor_P
  double precision:: formfactor_Ds, formfactor_Dp, formfactor_Dd
  integer:: locstatesnum, selectev, numeigen
  integer:: test


  massoffset = 0.0D0
  nmcutoff = 0.0D0
  xcutoff = 0.0D0
  step = 0.1D0

  !!!!!$$$$$$$$$ get basis information $$$$$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!
  !!!!!$$$$$$$$$ for the first Fock sector $$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!
  sums(1)=0
  call snmlen2(Nmax1, Mj, snm)
  call dddnmlen2(Nmax1, Mj, sums(2))
  call uddnmlen2(Nmax1, Mj, sums(3))
  call dudnmlen2(Nmax1, Mj, sums(4))
  call uudnmlen2(Nmax1, Mj, sums(5))
  call ddunmlen2(Nmax1, Mj, sums(6))
  call udunmlen2(Nmax1, Mj, sums(7))
  call duunmlen2(Nmax1, Mj, sums(8))
  call uuunmlen2(Nmax1, Mj, sums(9))

  !!!!!!!!!$$$ for the five particle Fock sector $$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!

  call sumspin5p(Nmax3, Mj, SumS5p)
  call sumspin4p(Nmax2, Mj, SumS4p)


  !!$$$$$$$$ which quark can be calculate by using different parameter selectquark $$$$$$$$$!!!!!
  !!$$$$$$$$ when the selectquark = 1, calculating the d quark of proton in the script $$$$!!!!!
  !!$$$$$$$$ when the selectquark = 2, calculating the u quark of proton in the script $$$$!!!!!

  !if (selectquark .eq. 1) then


  !!$$$$$$ decide the dimension of distribution function on each process $$$$$!!!!!!!!!!!!!

  functionsdimension=(floor(deltamax/step)+1)

  Print*,"the dimensions of matrix of distribution functions = ",functionsdimension

  allocate(distributionfunctions1(7, kt), distributionfunctions2(7, kt), &
    & distributionfunctions3(7, kt), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array'
  end if

  allocate(distributionfunctions_FF(7, functionsdimension), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array'
  end if

  !!$$$$$$ wave function production $$$$$$$$$$$$$$$$$$$!!!!!!!!

  call dimtotal3p(nmax1, Mj, Kt, dimtot1)
  call dimtotal4p(nmax2, Mj, Kt, dimtot2)
  call dimtotal5p(nmax3, Mj, Kt, dimtot3)
  call ColorMatrix()

  Print*,"method=",method

  If(method .eq. 1) then
    dimtot = dimtot1+3*dimtot3
  else if(method .eq. 2) then
    dimtot = dimtot1+2*dimtot2+3*dimtot3
  else
    dimtot = dimtot1+3*dimtot3
  endif

  allocate(ev1(nestate), evector1(nestate, dimtot), evector(dimtot), stat = error)
  if(error .ne. 0) then
    print *, 'cannot allocate basis array'
  end if



  call wfproduction(nmax1, nmax2, nmax3, Mj, Kt, method, massoffset, nmcutoff, xcutoff, nestate, ev1, evector1)


  !!$$$$$$$$ calculate distribution functions $$$$$$$$$$$$$$$$!!!!!!!!!

  loopnumber_loc = 0
  numeigen = seqeigen

  do selectev = 1, dimtot
    evector(selectev)=evector1(numeigen, selectev)
  enddo

  selectquark = 1
  loopnumber_loc = 0

  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    call GPD_for_d_quark1(nmax1, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctions1)
    If (method .eq. 1) then
      call GPD_for_d_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)
    else If(method .eq. 2) then
      call GPD_for_d_quark2(nmax3, Mj, Kt, SumS4p, SumS4p(17), nestate, evector, dx, dy, distributionfunctions2)
      !call GPD_for_d_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)
    endif

    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0
    loopnumber_loc = loopnumber_loc+1

    do i = 1, kt

      If(method .eq. 1) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions3(4, i)

      else if(method .eq. 2) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions2(3, i)&
          & +1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions2(4, i)&
          & +1.0D0/Kt*distributionfunctions3(4, i)

      endif
      ! formfactor_Ds = formfactor_Ds+1.0D0/Kt*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/Kt*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/Kt*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions1(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  do test = 1, loopnumber_loc

    print*, "Form Factor for d quark", distributionfunctions_FF(2, test), "FF_Dirac=", distributionfunctions_FF(3, test), &
      & "FF_Pauli=", distributionfunctions_FF(4, test)

  enddo

  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!!test the normalization for d quark. and the normalization should be equal to 1 !!!!!!!!

  !do j = 1, kt-3

  !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

  !enddo

  !Print*,normalization


  !else if (selectquark .eq. 2) then
  selectquark = 2
  loopnumber_loc = 0


  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    !Print*,"test"

    call GPD_for_u_quark1(nmax1, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctions1)

    If (method .eq. 1) then
      call GPD_for_u_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)
    else if (method .eq. 2) then

      call GPD_for_u_quark2(nmax3, Mj, Kt, SumS4p, SumS4p(17), nestate, evector, dx, dy, distributionfunctions2)
      !call GPD_for_u_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)

    endif

    loopnumber_loc = loopnumber_loc+1
    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0

    do i = 1, kt

      If(method .eq. 1) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions3(4, i)

      else if(method .eq. 2) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions2(3, i)&
          & +1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions2(4, i)&
          & +1.0D0/Kt*distributionfunctions3(4, i)

      endif
      ! formfactor_Ds = formfactor_Ds+1.0D0/dble(Kt)*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/dble(Kt)*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/dble(Kt)*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions1(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  do test = 1, loopnumber_loc

    print*, "Form Factor for u quark", distributionfunctions_FF(2, test), "FF_Dirac=",distributionfunctions_FF(3, test), &
      & "FF_Pauli",distributionfunctions_FF(4, test)

  enddo

  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!!!!!test the normalization for u quark. And the normalization should be equal to 2 !!!!!!!!

  !do j = 1, kt-1

  !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

  !enddo

  !Print*,normalization

  selectquark = 3
  loopnumber_loc = 0

  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    call GPD_for_s_quark1(nmax1, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctions1)
    If(method .eq. 1) then
      call GPD_for_s_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)
    else if(method .eq. 2) then
      call GPD_for_s_quark2(nmax3, Mj, Kt, SumS4p, SumS4p(17), nestate, evector, dx, dy, distributionfunctions2)
      !call GPD_for_s_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)
    endif

    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0
    loopnumber_loc = loopnumber_loc+1

    do i = 1, kt

      If(method .eq. 1) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions3(4, i)

      else if(method .eq. 2) then

        formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions1(3, i)+1.0D0/Kt*distributionfunctions2(3, i)&
          & +1.0D0/Kt*distributionfunctions3(3, i)
        formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions1(4, i)+1.0D0/Kt*distributionfunctions2(4, i)&
          & +1.0D0/Kt*distributionfunctions3(4, i)

      endif
      ! formfactor_Ds = formfactor_Ds+1.0D0/Kt*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/Kt*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/Kt*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions1(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  do test = 1, loopnumber_loc

    print*, "Form Factor for s quark", distributionfunctions_FF(2, test), "FF_Dirac=", distributionfunctions_FF(3, test), &
      & "FF_Pauli=",distributionfunctions_FF(4, test)

  enddo


  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!!!!$$$$$$$$$ calculate the form factor for gluon $$$$$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!!!!

  selectquark = 5
  loopnumber_loc = 0

  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    call GPD_for_gluon(nmax3, Mj, Kt, SumS4p, SumS4p(17), nestate, evector, dx, dy, distributionfunctions2)

    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0
    loopnumber_loc = loopnumber_loc+1

    do i = 1, kt

      formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions2(3, i)
      formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions2(4, i)
      ! formfactor_Ds = formfactor_Ds+1.0D0/Kt*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/Kt*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/Kt*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions2(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  print*, " "
  print*, "calculate the form factor for gluon!"
  do test = 1, loopnumber_loc

    print*, "Form Factor for d quark", distributionfunctions_FF(2, test), distributionfunctions_FF(3, test), &
      & distributionfunctions_FF(4, test)

  enddo

  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!!!!$$$$$$$$$ calculate the form factor for sea quark $$$$$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!!!!

  selectquark = 5
  loopnumber_loc = 0

  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    !call GPD_for_sea1_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)

    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0
    loopnumber_loc = loopnumber_loc+1

    do i = 1, kt

      formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions3(3, i)
      formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions3(4, i)
      ! formfactor_Ds = formfactor_Ds+1.0D0/Kt*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/Kt*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/Kt*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions3(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  print*, " "
  print*, "calculate the form factor for sea quark!"
  do test = 1, loopnumber_loc

    print*, "Form Factor for d quark", distributionfunctions_FF(2, test), distributionfunctions_FF(3, test), &
      & distributionfunctions_FF(4, test)

  enddo

  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!$$$$$$$$$$$$$$$$$$ calculate the form factor for anti sea quark $$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!!!!

  selectquark = 6
  loopnumber_loc = 0

  do delta = 0, floor(deltamax/step)

    dy = 0.0D0
    dx = delta*step

    !call GPD_for_sea2_quark3(nmax3, Mj, Kt, SumS5p, SumS5p(33), nestate, evector, dx, dy, method, distributionfunctions3)

    formfactor_D = 0.D0
    formfactor_P = 0.D0
    formfactor_Ds = 0.D0
    formfactor_Dp = 0.D0
    formfactor_Dd = 0.D0
    loopnumber_loc = loopnumber_loc+1

    do i = 1, kt

      formfactor_D = formfactor_D+1.0D0/Kt*distributionfunctions3(3, i)
      formfactor_P = formfactor_P+1.0D0/Kt*distributionfunctions3(4, i)
      ! formfactor_Ds = formfactor_Ds+1.0D0/Kt*distributionfunctions1(5, i)
      ! formfactor_Dp = formfactor_Dp+1.0D0/Kt*distributionfunctions1(6, i)
      ! formfactor_Dd = formfactor_Dd+1.0D0/Kt*distributionfunctions1(7, i)

    enddo

    distributionfunctions_FF(1, loopnumber_loc)=0.D0
    distributionfunctions_FF(2, loopnumber_loc)=distributionfunctions3(2, 1)
    distributionfunctions_FF(3, loopnumber_loc)=formfactor_D
    distributionfunctions_FF(4, loopnumber_loc)=formfactor_P
    ! distributionfunctions_FF(5, loopnumber_loc)=formfactor_Ds
    ! distributionfunctions_FF(6, loopnumber_loc)=formfactor_Dp
    ! distributionfunctions_FF(7, loopnumber_loc)=formfactor_Dd

  enddo

  print*, " "
  print*, "calculate the form factor for anti-sea quark!"
  do test = 1, loopnumber_loc

    print*, "Form Factor for d quark", distributionfunctions_FF(2, test), distributionfunctions_FF(3, test), &
      & distributionfunctions_FF(4, test)

  enddo

  !call output(nmax1, nmax2, nmax3, Mj, Kt, selectquark, loopnumber_loc, distributionfunctions_FF)

  !!!!!!!!test the normalization for d quark. and the normalization should be equal to 1 !!!!!!!!

  !do j = 1, kt-3

  !    normalization = normalization+(distributionfunctions(1, j+1)-distributionfunctions(1, j))*distributionfunctions(3, j)

  !enddo

  !Print*,normalization

  deallocate(s1dat1, s2dat1, s3dat1, n1dat1, m1dat1, n2dat1, m2dat1, n3dat1, m3dat1, k1dat1, k2dat1, k3dat1)
  deallocate(s1dat2, s2dat2, s3dat2, s4dat2, n1dat2, m1dat2, n2dat2, m2dat2, n3dat2, m3dat2, n4dat2, m4dat2, k1dat2, &
    & k2dat2, k3dat2, k4dat2)
  deallocate(s1dat3, s2dat3, s3dat3, s4dat3, s5dat3, n1dat3, m1dat3, n2dat3, m2dat3, n3dat3, m3dat3, n4dat3, m4dat3, &
    & n5dat3, m5dat3, k1dat3, k2dat3, k3dat3, k4dat3, k5dat3)
  deallocate(distributionfunctions1, distributionfunctions2, distributionfunctions_FF, ev1, evector, evector1)

  !endif

  return
end subroutine distribution_functions_FF


subroutine GPD_for_d_quark1(nmax1, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax1, Mj, Kt, nestate, snm
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision:: dx, dy, bq, deltamax
  double precision, dimension(7, kt):: distributionfunctionsvalue
  double precision:: distributionfunctionsvalue_E, dxcorrect
  double precision:: distribution_H_s_loc, distribution_H_p_loc, distribution_H, distribution_E_loc, distribution_E
  double precision:: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d, distribution_H_loc

  integer, dimension(Kt+1):: kmatrix
  integer:: npu1, mpu1, npd, mpd, npu2, mpu2, pu1, pd, pu2, nku1, mku1, nkd, mkd, nku2, mku2, ku1, kd, ku2
  integer:: accumulate, k_d, k_u, spin, totalspin, initialstate, finalstate
  double precision:: fraction, bp, bk, factor, olaps, eigenvaluek, eigenvaluep
  integer:: summaxp, summaxk, sum1, sum2, i, error, test
  integer:: uplimit, finalindex, focksector
  integer:: initialmin, initialmax


  kmatrix(1)=0
  fraction = 0.0D0
  accumulate = 0
  factor = 0.D0
  dxcorrect = 1.0D-5
  bq = B_DEFAULT

  do i = 0, kt-1
    accumulate = accumulate+Kt-i
    kmatrix(i+2)=accumulate*snm
  enddo


  do k_d = 0, kt-1

    fraction=(dble(k_d)+0.5)/(dble(Kt)+0.5)

    distributionfunctionsvalue(3, k_d+1)=0.D0
    distributionfunctionsvalue(4, k_d+1)=0.D0

    distribution_H_s = 0.D0
    distribution_H_p = 0.D0
    distribution_H_d = 0.D0
    distribution_H = 0.D0
    distribution_E = 0.D0

    do k_u = 0, Kt-k_d-1

      do spin = 1, 8

        initialmin = kmatrix(k_d+1)+k_u*snm+sums(spin)+1
        initialmax = kmatrix(k_d+1)+k_u*snm+sums(spin+1)

        do initialstate = initialmin, initialmax

          eigenvaluep = evector(initialstate)

          if (abs(eigenvaluep).gt.10.0D0**(-13)) then

            npu1 = n1dat1(initialstate)
            mpu1 = m1dat1(initialstate)
            npd = n2dat1(initialstate)
            mpd = m2dat1(initialstate)
            npu2 = n3dat1(initialstate)
            mpu2 = m3dat1(initialstate)
            pu1 = k1dat1(initialstate)
            pd  =k2dat1(initialstate)
            pu2 = k3dat1(initialstate)

          else
            cycle
          endif

          if (k_d .eq. pd .and. k_u .eq. pu1) then

            ! call Talmisquare(nmax, npu1, mpu1, npu2, mpu2, npd, mpd, dble(pu1)+0.5, dble(pu2)+0.5, &
            !             & dble(pd)+0.5, summaxp, talmip)

            do finalstate = kmatrix(k_d+1)+k_u*snm+sums(spin)+1, kmatrix(k_d+1)+k_u*snm+sums(spin+1)

              eigenvaluek = evector(finalstate)

              if (abs(eigenvaluek).gt.10.0D0**(-13)) then

                nku1 = n1dat1(finalstate)
                mku1 = m1dat1(finalstate)
                nkd = n2dat1(finalstate)
                mkd = m2dat1(finalstate)
                nku2 = n3dat1(finalstate)
                mku2 = m3dat1(finalstate)
                ku1 = k1dat1(finalstate)
                kd  =k2dat1(finalstate)
                ku2 = k3dat1(finalstate)

              else
                cycle
              endif

              factor = olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
                & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
                & -(dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx), -(dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
                & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
                & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
                & -(dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx), -(dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
                & olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
                & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
                & (1.0D0-fraction)*(dx), (1.0D0-fraction)*dy)

              distributionfunctionsvalue(1, k_d+1)=fraction
              distributionfunctionsvalue(2, k_d+1)=dx**2+dy**2
              distribution_H = distribution_H+eigenvaluep*eigenvaluek*factor*kt

              ! print*,fraction, eigenvaluep*eigenvaluek*factor*kt, eigenvaluep, eigenvaluek, factor, &
              ! &  initialstate, finalstate

            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            uplimit = sums(10-spin)-sums(9-spin)

            do finalindex = 1, uplimit

              finalstate = kmatrix(k_d+1)+k_u*snm+sums(9-spin)+finalindex

              eigenvaluek = evector(finalstate)
              nku1 = n1dat1(finalstate)
              mku1 = m1dat1(finalstate)
              nkd = n2dat1(finalstate)
              mkd = m2dat1(finalstate)
              nku2 = n3dat1(finalstate)
              mku2 = m3dat1(finalstate)
              ku1 = k1dat1(finalstate)
              kd  =k2dat1(finalstate)
              ku2 = k3dat1(finalstate)


              factor = olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
                & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
                & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
                & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
                & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
                & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
                & olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
                & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
                & -(1.0D0-fraction)*(dx+dxcorrect), -(1.0D0-fraction)*dy) &
                & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2/(dx+dxcorrect)

              distribution_E = distribution_E+eigenvaluep*eigenvaluek*factor*dble(kt)

              ! If(kd .eq. 2) then
              !     Print*,istribution_E, initialstate, finalstate, fraction, mpu1+mpu2+mpd, mku1+mku2+mkd, &
              !         & eigenvaluep, eigenvaluek, factor
              ! endif

            enddo
          endif
        enddo
      enddo
    enddo

    distributionfunctionsvalue(1, k_d+1)=fraction
    distributionfunctionsvalue(2, k_d+1)=dx**2+dy**2
    distributionfunctionsvalue(3, k_d+1)=distribution_H!(distribution_H_s+distribution_H_p+distribution_H_d)!*fraction
    distributionfunctionsvalue(4, k_d+1)=distribution_E
    !distributionfunctionsvalue(5, k_d+1)=distribution_H_s
    !distributionfunctionsvalue(6, k_d+1)=distribution_H_p
    !distributionfunctionsvalue(7, k_d+1)=distribution_H_d

  enddo
  Print*, " "
  Print*,"GPD_for_d_quark"
  Print*, " Fraction","     ||     ","dx**2+dy**2","     ||     ","distribution_H","     ||     ","distribution_E"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_d_quark1

subroutine GPD_for_d_quark2(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(2):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot2-snm, snm

    do spin = 1, 16

      do initialstate = k+sums(spin)+1, k+sums(spin+1)
        !Print*,initialstate

        eigenvaluep(1)=evector(dimtot1+initialstate)
        eigenvaluep(2)=evector(dimtot1+dimtot2+initialstate)

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10)) then

          np1 = n1dat2(initialstate)
          mp1 = m1dat2(initialstate)
          np2 = n2dat2(initialstate)
          mp2 = m2dat2(initialstate)
          np3 = n3dat2(initialstate)
          mp3 = m3dat2(initialstate)
          np4 = n4dat2(initialstate)
          mp4 = m4dat2(initialstate)
          kp1 = k1dat2(initialstate)
          kp2 = k2dat2(initialstate)
          kp3 = k3dat2(initialstate)
          kp4 = k4dat2(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp2+1)=distribution_H(kp2+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)
          !Print*,distribution_H(kp2+1), kp2


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(18-spin)-sums(17-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(17-spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4)*2.0D0/(dx+dxcorrect)


          distribution_E(kp2+1)=distribution_E(kp2+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_d_quark2"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_d_quark2

subroutine GPD_for_d_quark3(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, method, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate, method
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot3-snm, snm

    do spin = 1, 32

      do initialstate = k+sums(spin)+1, k+sums(spin+1)
        !Print*,initialstate

        If (method .eq. 1) then

          eigenvaluep(1)=evector(dimtot1+initialstate)
          eigenvaluep(2)=evector(dimtot1+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot3+initialstate)
        else if (method .eq. 2) then
          eigenvaluep(1)=evector(dimtot1+2*dimtot2+initialstate)
          eigenvaluep(2)=evector(dimtot1+2*dimtot2+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot2+2*dimtot3+initialstate)
        endif

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10).or.abs(eigenvaluep(3)).gt.10.0D0**(-10)) then

          np1 = n1dat3(initialstate)
          mp1 = m1dat3(initialstate)
          np2 = n2dat3(initialstate)
          mp2 = m2dat3(initialstate)
          np3 = n3dat3(initialstate)
          mp3 = m3dat3(initialstate)
          np4 = n4dat3(initialstate)
          mp4 = m4dat3(initialstate)
          np5 = n5dat3(initialstate)
          mp5 = m5dat3(initialstate)
          kp1 = k1dat3(initialstate)
          kp2 = k2dat3(initialstate)
          kp3 = k3dat3(initialstate)
          kp4 = k4dat3(initialstate)
          kp5 = k5dat3(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          If(method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if(method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp2+1)=distribution_H(kp2+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)

          !Print*,distribution_H(kp2+1), kp2


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(34-spin)-sums(33-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(33-spin)+finalindex

          If(method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if(method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, -mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4+mk5+1)*2.0D0/(dx+dxcorrect)


          distribution_E(kp2+1)=distribution_E(kp2+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_d_quark3"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_d_quark3

subroutine GPD_for_u_quark1(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision:: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: npu1, mpu1, npd, mpd, npu2, mpu2, pu1, pu2, pd
  integer:: nku1, mku1, nkd, mkd, nku2, mku2, ku1, ku2, kd
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  double precision, dimension(:,:), allocatable:: talmip1, talmik1, talmip2, talmik2
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error

  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT
  ! Print*,"test"

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot1-snm, snm

    do spin = 1, 8

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        eigenvaluep = evector(initialstate)

        if(abs(eigenvaluep).gt.10.0D0**(-10)) then

          npu1 = n1dat1(initialstate)
          mpu1 = m1dat1(initialstate)
          npd = n2dat1(initialstate)
          mpd = m2dat1(initialstate)
          npu2 = n3dat1(initialstate)
          mpu2 = m3dat1(initialstate)
          pu1 = k1dat1(initialstate)
          pd  =k2dat1(initialstate)
          pu2 = k3dat1(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          eigenvaluek = evector(finalstate)

          if(abs(eigenvaluek).gt.10.0D0**(-10)) then

            nku1 = n1dat1(finalstate)
            mku1 = m1dat1(finalstate)
            nkd = n2dat1(finalstate)
            mkd = m2dat1(finalstate)
            nku2 = n3dat1(finalstate)
            mku2 = m3dat1(finalstate)
            ku1 = k1dat1(finalstate)
            kd  =k2dat1(finalstate)
            ku2 = k3dat1(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(pu1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
            & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
            & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
            & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
            & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
            & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(pu1+1)=distribution_H(pu1+1)+eigenvaluep*eigenvaluek*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(10-spin)-sums(9-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(9-spin)+finalindex

          eigenvaluek = evector(finalstate)

          if(abs(eigenvaluek).gt.10.0D0**(-10)) then

            nku1 = n1dat1(finalstate)
            mku1 = m1dat1(finalstate)
            nkd = n2dat1(finalstate)
            mkd = m2dat1(finalstate)
            nku2 = n3dat1(finalstate)
            mku2 = m3dat1(finalstate)
            ku1 = k1dat1(finalstate)
            kd  =k2dat1(finalstate)
            ku2 = k3dat1(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(pu1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
            & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
            & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
            & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
            & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
            & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          distribution_E(pu1+1)=distribution_E(pu1+1)+eigenvaluep*eigenvaluek*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*," "
  Print*,"GPD_for_u_quark"
  Print*, " Fraction","     ||     ","dx**2+dy**2","     ||     ","distribution_H","     ||     ","distribution_E"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_u_quark1

subroutine GPD_for_u_quark2(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(2):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot2-snm, snm

    do spin = 1, 16

      do initialstate = k+sums(spin)+1, k+sums(spin+1)


        eigenvaluep(1)=evector(dimtot1+initialstate)
        eigenvaluep(2)=evector(dimtot1+dimtot2+initialstate)


        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10)) then

          np1 = n1dat2(initialstate)
          mp1 = m1dat2(initialstate)
          np2 = n2dat2(initialstate)
          mp2 = m2dat2(initialstate)
          np3 = n3dat2(initialstate)
          mp3 = m3dat2(initialstate)
          np4 = n4dat2(initialstate)
          mp4 = m4dat2(initialstate)
          kp1 = k1dat2(initialstate)
          kp2 = k2dat2(initialstate)
          kp3 = k3dat2(initialstate)
          kp4 = k4dat2(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex


          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)


          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp1+1)=distribution_H(kp1+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(18-spin)-sums(17-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(17-spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4)*2.0D0/(dx+dxcorrect)


          distribution_E(kp1+1)=distribution_E(kp1+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_u_quark2"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_u_quark2

subroutine GPD_for_u_quark3(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, method, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate, method
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot3-snm, snm

    do spin = 1, 32

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        If (method .eq. 1) then

          eigenvaluep(1)=evector(dimtot1+initialstate)
          eigenvaluep(2)=evector(dimtot1+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot3+initialstate)

        else if(method .eq. 2) then

          eigenvaluep(1)=evector(dimtot1+2*dimtot2+initialstate)
          eigenvaluep(2)=evector(dimtot1+2*dimtot2+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot2+2*dimtot3+initialstate)

        endif

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10).or.abs(eigenvaluep(3)).gt.10.0D0**(-10)) then

          np1 = n1dat3(initialstate)
          mp1 = m1dat3(initialstate)
          np2 = n2dat3(initialstate)
          mp2 = m2dat3(initialstate)
          np3 = n3dat3(initialstate)
          mp3 = m3dat3(initialstate)
          np4 = n4dat3(initialstate)
          mp4 = m4dat3(initialstate)
          np5 = n5dat3(initialstate)
          mp5 = m5dat3(initialstate)
          kp1 = k1dat3(initialstate)
          kp2 = k2dat3(initialstate)
          kp3 = k3dat3(initialstate)
          kp4 = k4dat3(initialstate)
          kp5 = k5dat3(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif


          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp1+1)=distribution_H(kp1+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(34-spin)-sums(33-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(33-spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp1)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, -mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4+mk5+1)*2.0D0/(dx+dxcorrect)


          distribution_E(kp1+1)=distribution_E(kp1+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_u_quark3"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_u_quark3

subroutine GPD_for_s_quark1(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision:: ev
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H, distribution_E
  double precision:: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: npu1, mpu1, npd, mpd, npu2, mpu2, pu1, pu2, pd
  integer:: nku1, mku1, nkd, mkd, nku2, mku2, ku1, ku2, kd
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  double precision, dimension(:,:), allocatable:: talmip1, talmik1, talmip2, talmik2
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt


    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot1-snm, snm


    do spin = 1, 8

      do initialstate = k+sums(spin)+1, k+sums(spin+1)


        eigenvaluep = evector(initialstate)

        if(abs(eigenvaluep).gt.10.0D0**(-10)) then

          npu1 = n1dat1(initialstate)
          mpu1 = m1dat1(initialstate)
          npd = n2dat1(initialstate)
          mpd = m2dat1(initialstate)
          npu2 = n3dat1(initialstate)
          mpu2 = m3dat1(initialstate)
          pu1 = k1dat1(initialstate)
          pd  =k2dat1(initialstate)
          pu2 = k3dat1(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex


          eigenvaluek = evector(finalstate)

          if(abs(eigenvaluek).gt.10.0D0**(-10)) then

            nku1 = n1dat1(finalstate)
            mku1 = m1dat1(finalstate)
            nkd = n2dat1(finalstate)
            mkd = m2dat1(finalstate)
            nku2 = n3dat1(finalstate)
            mku2 = m3dat1(finalstate)
            ku1 = k1dat1(finalstate)
            kd  =k2dat1(finalstate)
            ku2 = k3dat1(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          ! fractionx=(dble(pu1)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !     & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !     & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !     & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !     & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !     & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
          !     & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !     & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !     & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          ! distribution_H_loc(pu1+1)=distribution_H_loc(pu1+1)+eigenvaluep*eigenvaluek*factor*(kt)


          fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
            & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
            & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
            & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
            & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
            & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          distribution_H(pu2+1)=distribution_H(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(10-spin)-sums(9-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(9-spin)+finalindex

          eigenvaluek = evector(finalstate)

          if(abs(eigenvaluek).gt.10.0D0**(-10)) then

            nku1 = n1dat1(finalstate)
            mku1 = m1dat1(finalstate)
            nkd = n2dat1(finalstate)
            mkd = m2dat1(finalstate)
            nku2 = n3dat1(finalstate)
            mku2 = m3dat1(finalstate)
            ku1 = k1dat1(finalstate)
            kd  =k2dat1(finalstate)
            ku2 = k3dat1(finalstate)


          else
            cycle
          endif



          ! fractionx=(dble(pu1)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !     & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !     & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !     & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !     & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !     & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
          !     & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !     & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !     & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !     & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          ! distribution_E_loc(pu1+1)=distribution_E_loc(pu1+1)+eigenvaluep*eigenvaluek*factor*(kt)


          fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
            & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
            & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
            & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
            & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
            & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          distribution_E(pu2+1)=distribution_E(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)
          ! If(pu2 .eq. 3) then
          !     Print*, distribution_E(pu2+1), pu2, initialstate, finalstate, mpu1+mpu2+mpd, mku1+mku2+mkd, &
          !         & eigenvaluep, eigenvaluek, factor
          ! endif

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo
  print*," "
  Print*,"GPD_for_s_quark"
  Print*, " Fraction","     ||     ","dx**2+dy**2","     ||     ","distribution_H","     ||     ","distribution_E"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_s_quark1

subroutine GPD_for_s_quark2(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot2-snm, snm

    do spin = 1, 16

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        eigenvaluep(1)=evector(dimtot1+initialstate)
        eigenvaluep(2)=evector(dimtot1+dimtot2+initialstate)

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10)) then

          ! if(abs(eigenvaluep).gt.10.0D0**(-10)) then

          np1 = n1dat2(initialstate)
          mp1 = m1dat2(initialstate)
          np2 = n2dat2(initialstate)
          mp2 = m2dat2(initialstate)
          np3 = n3dat2(initialstate)
          mp3 = m3dat2(initialstate)
          np4 = n4dat2(initialstate)
          mp4 = m4dat2(initialstate)
          kp1 = k1dat2(initialstate)
          kp2 = k2dat2(initialstate)
          kp3 = k3dat2(initialstate)
          kp4 = k4dat2(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp3)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp3+1)=distribution_H(kp3+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(18-spin)-sums(17-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(17-spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp3)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt))), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4))/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4))/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4)*2.0D0/(dx+dxcorrect)


          distribution_E(kp3+1)=distribution_E(kp3+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_s_quark2"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_s_quark2

subroutine GPD_for_s_quark3(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, method, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate, method
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot3-snm, snm

    do spin = 1, 32

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        If (method .eq. 1) then

          eigenvaluep(1)=evector(dimtot1+initialstate)
          eigenvaluep(2)=evector(dimtot1+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot3+initialstate)

        else if (method .eq. 2) then

          eigenvaluep(1)=evector(dimtot1+2*dimtot2+initialstate)
          eigenvaluep(2)=evector(dimtot1+2*dimtot2+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot2+2*dimtot3+initialstate)

        endif

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10).or.abs(eigenvaluep(3)).gt.10.0D0**(-10)) then

          ! if(abs(eigenvaluep).gt.10.0D0**(-10)) then

          np1 = n1dat3(initialstate)
          mp1 = m1dat3(initialstate)
          np2 = n2dat3(initialstate)
          mp2 = m2dat3(initialstate)
          np3 = n3dat3(initialstate)
          mp3 = m3dat3(initialstate)
          np4 = n4dat3(initialstate)
          mp4 = m4dat3(initialstate)
          np5 = n5dat3(initialstate)
          mp5 = m5dat3(initialstate)
          kp1 = k1dat3(initialstate)
          kp2 = k2dat3(initialstate)
          kp3 = k3dat3(initialstate)
          kp4 = k4dat3(initialstate)
          kp5 = k5dat3(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp3)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp3+1)=distribution_H(kp3+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(34-spin)-sums(33-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(33-spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp3)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, -mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4+mk5+1)*2.0D0/(dx+dxcorrect)


          distribution_E(kp3+1)=distribution_E(kp3+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_s_quark3"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_s_quark3

subroutine GPD_for_gluon(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(2):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot2-snm, snm

    do spin = 1, 16

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        eigenvaluep(1)=evector(dimtot1+initialstate)
        eigenvaluep(2)=evector(dimtot1+dimtot2+initialstate)

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10)) then

          ! if(abs(eigenvaluep).gt.10.0D0**(-10)) then

          np1 = n1dat2(initialstate)
          mp1 = m1dat2(initialstate)
          np2 = n2dat2(initialstate)
          mp2 = m2dat2(initialstate)
          np3 = n3dat2(initialstate)
          mp3 = m3dat2(initialstate)
          np4 = n4dat2(initialstate)
          mp4 = m4dat2(initialstate)
          kp1 = k1dat2(initialstate)
          kp2 = k2dat2(initialstate)
          kp3 = k3dat2(initialstate)
          kp4 = k4dat2(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp4))/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, mk4, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp4)=distribution_H(kp4)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(18-spin)-sums(17-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(17-spin)+finalindex

          eigenvaluek(1)=evector(dimtot1+finalstate)
          eigenvaluek(2)=evector(dimtot1+dimtot2+finalstate)

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10)) then

            nk1 = n1dat2(finalstate)
            mk1 = m1dat2(finalstate)
            nk2 = n2dat2(finalstate)
            mk2 = m2dat2(finalstate)
            nk3 = n3dat2(finalstate)
            mk3 = m3dat2(finalstate)
            nk4 = n4dat2(finalstate)
            mk4 = m4dat2(finalstate)
            kk1 = k1dat2(finalstate)
            kk2 = k2dat2(finalstate)
            kk3 = k3dat2(finalstate)
            kk4 = k4dat2(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp4))/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4))/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4))/(dble(Kt)+0.5)), nk4, -mk4, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4)*2.0D0/(dx+dxcorrect)


          distribution_E(kp4)=distribution_E(kp4)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u))/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo
  print*," "
  Print*,"GPD_for_gluon"
  Print*, " Fraction","     ||     ","dx**2+dy**2","     ||     ","distribution_H","     ||     ","distribution_E"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_gluon

subroutine GPD_for_sea1_quark3(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, method, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate, method
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error


  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot3-snm, snm

    do spin = 1, 32

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        If (method .eq. 1) then

          eigenvaluep(1)=evector(dimtot1+initialstate)
          eigenvaluep(2)=evector(dimtot1+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot3+initialstate)

        else if (method .eq. 2) then

          eigenvaluep(1)=evector(dimtot1+2*dimtot2+initialstate)
          eigenvaluep(2)=evector(dimtot1+2*dimtot2+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot2+2*dimtot3+initialstate)

        endif

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10).or.abs(eigenvaluep(3)).gt.10.0D0**(-10)) then

          np1 = n1dat3(initialstate)
          mp1 = m1dat3(initialstate)
          np2 = n2dat3(initialstate)
          mp2 = m2dat3(initialstate)
          np3 = n3dat3(initialstate)
          mp3 = m3dat3(initialstate)
          np4 = n4dat3(initialstate)
          mp4 = m4dat3(initialstate)
          np5 = n5dat3(initialstate)
          mp5 = m5dat3(initialstate)
          kp1 = k1dat3(initialstate)
          kp2 = k2dat3(initialstate)
          kp3 = k3dat3(initialstate)
          kp4 = k4dat3(initialstate)
          kp5 = k5dat3(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp4)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, mk4, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp4+1)=distribution_H(kp4+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(34-spin)-sums(33-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(33-spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp4)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, -mk5, &
            & (dble(kp5)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp5)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, -mk4, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4+mk5+1)*2.0D0/(dx+dxcorrect)


          distribution_E(kp4+1)=distribution_E(kp4+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_sea1_quark"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_sea1_quark3

subroutine GPD_for_sea2_quark3(nmax, Mj, Kt, sums, snm, nestate, evector, dx, dy, method, distributionfunctionsvalue)
  use numbers
  use basis_info

  implicit none

  integer:: nmax, Mj, Kt, snm, nestate, method
  double precision:: bq, dx, dy
  integer, dimension(*):: sums
  double precision, dimension(dimtot):: evector
  double precision, dimension(7, Kt):: distributionfunctionsvalue
  double precision, dimension(kt):: distribution_H_loc, distribution_H, distribution_E_loc, distribution_E
  double precision, dimension(3):: eigenvaluep, eigenvaluek
  double precision, dimension(kt):: distribution_H_s_loc, distribution_H_p_loc
  double precision, dimension(kt):: distribution_H_s, distribution_H_p, distribution_H_d_loc, distribution_H_d

  double precision:: test1, dxcorrect, dycorrect

  double precision:: fractionx, bp, bk, factor, olaps
  integer:: np1, mp1, np2, mp2, np3, mp3, np4, mp4, np5, mp5, kp1, kp2, kp3, kp4, kp5
  integer:: nk1, mk1, nk2, mk2, nk3, mk3, nk4, mk4, nk5, mk5, kk1, kk2, kk3, kk4, kk5
  integer:: k_u, k, spin, initialstate, finalstate, finalindex, finalseperate, uplimit
  integer:: summaxp1, summaxk1, summaxp2, summaxk2, looptmp, looptmk, talmitotnumber
  integer:: test, error

  fractionx = 0.D0
  dxcorrect = 10.0D-10
  dxcorrect = 10.0D-10
  bq = B_DEFAULT

  do k_u = 1, Kt

    distribution_H(k_u)=0.D0
    distribution_H_s(k_u)=0.D0
    distribution_H_p(k_u)=0.D0
    distribution_H_d(k_u)=0.D0
    distribution_E(k_u)=0.D0

  enddo

  do k = 0, dimtot3-snm, snm

    do spin = 1, 32

      do initialstate = k+sums(spin)+1, k+sums(spin+1)

        If (method .eq. 1) then

          eigenvaluep(1)=evector(dimtot1+initialstate)
          eigenvaluep(2)=evector(dimtot1+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot3+initialstate)

        else if (method .eq. 2) then

          eigenvaluep(1)=evector(dimtot1+2*dimtot2+initialstate)
          eigenvaluep(2)=evector(dimtot1+2*dimtot2+dimtot3+initialstate)
          eigenvaluep(3)=evector(dimtot1+2*dimtot2+2*dimtot3+initialstate)

        endif

        if(abs(eigenvaluep(1)).gt.10.0D0**(-10).or.abs(eigenvaluep(2))&
          &.gt.10.0D0**(-10).or.abs(eigenvaluep(3)).gt.10.0D0**(-10)) then

          np1 = n1dat3(initialstate)
          mp1 = m1dat3(initialstate)
          np2 = n2dat3(initialstate)
          mp2 = m2dat3(initialstate)
          np3 = n3dat3(initialstate)
          mp3 = m3dat3(initialstate)
          np4 = n4dat3(initialstate)
          mp4 = m4dat3(initialstate)
          np5 = n5dat3(initialstate)
          mp5 = m5dat3(initialstate)
          kp1 = k1dat3(initialstate)
          kp2 = k2dat3(initialstate)
          kp3 = k3dat3(initialstate)
          kp4 = k4dat3(initialstate)
          kp5 = k5dat3(initialstate)

          ! call Talmisquare(nmax, npd, mpd, npu2, mpu2, npu1, mpu1, dble(pd)+0.5, dble(pu2)+0.5, &
          !     & dble(pu1)+0.5, summaxp1, talmip1)

          ! call Talmisquare(nmax, npd, mpd, npu1, mpu1, npu2, mpu2, dble(pd)+0.5, dble(pu1)+0.5, &
          !     & dble(pu2)+0.5, summaxp2, talmip2)

        else
          cycle
        endif

        !!!!!!########### calculating the GPD-H of proton   #############!!!!!!!!!!!

        uplimit = sums(spin+1)-sums(spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            ! call Talmisquare(nmax, nkd, mkd, nku2, mku2, nku1, mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !      & dble(ku1)+0.5, summaxk1, talmik1)
            ! call Talmisquare(nmax, nkd, mkd, nku1, mku1, nku2, mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !      & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif


          fractionx=(dble(kp5)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, mk5, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)

          distribution_H(kp5+1)=distribution_H(kp5+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy)


          !     distribution_H_loc(pu2+1)=distribution_H_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

        !!!!!!!!##############  calculating the GPD-E of proton   ##########!!!!!!!!!!!!!!!!!!!


        uplimit = sums(34-spin)-sums(33-spin)


        do finalindex = 1, uplimit


          finalstate = k+sums(33-spin)+finalindex

          If (method .eq. 1) then

            eigenvaluek(1)=evector(dimtot1+finalstate)
            eigenvaluek(2)=evector(dimtot1+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot3+finalstate)

          else if (method .eq. 2) then

            eigenvaluek(1)=evector(dimtot1+2*dimtot2+finalstate)
            eigenvaluek(2)=evector(dimtot1+2*dimtot2+dimtot3+finalstate)
            eigenvaluek(3)=evector(dimtot1+2*dimtot2+2*dimtot3+finalstate)

          endif

          if(abs(eigenvaluek(1)).gt.10.0D0**(-10).or.abs(eigenvaluek(2))&
            &.gt.10.0D0**(-10).or.abs(eigenvaluek(3)).gt.10.0D0**(-10)) then

            nk1 = n1dat3(finalstate)
            mk1 = m1dat3(finalstate)
            nk2 = n2dat3(finalstate)
            mk2 = m2dat3(finalstate)
            nk3 = n3dat3(finalstate)
            mk3 = m3dat3(finalstate)
            nk4 = n4dat3(finalstate)
            mk4 = m4dat3(finalstate)
            nk5 = n5dat3(finalstate)
            mk5 = m5dat3(finalstate)
            kk1 = k1dat3(finalstate)
            kk2 = k2dat3(finalstate)
            kk3 = k3dat3(finalstate)
            kk4 = k4dat3(finalstate)
            kk5 = k5dat3(finalstate)


            !call Talmisquare(nmax, nkd, -mkd, nku2, -mku2, nku1, -mku1, dble(kd)+0.5, dble(ku2)+0.5, &
            !    & dble(ku1)+0.5, summaxk1, talmik1)

            !call Talmisquare(nmax, nkd, -mkd, nku1, -mku1, nku2, -mku2, dble(kd)+0.5, dble(ku1)+0.5, &
            !    & dble(ku2)+0.5, summaxk2, talmik2)

          else
            cycle
          endif



          fractionx=(dble(kp5)+0.5)/(dble(Kt)+0.5)
          factor = olaps(bq*sqrt((dble(kp2)+0.5)/(dble(Kt)+0.5)), np2, mp2, &
            & bq*sqrt((dble(kk2)+0.5)/(dble(Kt)+0.5)), nk2, -mk2, &
            & (dble(kp2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp2)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp1)+0.5)/(dble(Kt)+0.5)), np1, mp1, &
            & bq*sqrt((dble(kk1)+0.5)/(dble(Kt)+0.5)), nk1, -mk1, &
            & (dble(kp1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp1)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp3)+0.5)/(dble(Kt)+0.5)), np3, mp3, &
            & bq*sqrt((dble(kk3)+0.5)/(dble(Kt)+0.5)), nk3, -mk3, &
            & (dble(kp3)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp3)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp4)+0.5)/(dble(Kt)+0.5)), np4, mp4, &
            & bq*sqrt((dble(kk4)+0.5)/(dble(Kt)+0.5)), nk4, -mk4, &
            & (dble(kp4)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(kp4)+0.5)/(dble(Kt)+0.5)*dy)*&
            & olaps(bq*sqrt((dble(kp5)+0.5)/(dble(Kt)+0.5)), np5, mp5, &
            & bq*sqrt((dble(kk5)+0.5)/(dble(Kt)+0.5)), nk5, -mk5, &
            & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
            & *(0.D0-1.0D0)**(mk1+mk2+mk3+mk4+mk5+1)*2.0D0/(dx+dxcorrect)


          distribution_E(kp5+1)=distribution_E(kp5+1)+eigenvaluep(1)*eigenvaluek(1)*factor*(kt)&
            &+eigenvaluep(2)*eigenvaluek(2)*factor*(kt)+eigenvaluep(3)*eigenvaluek(3)*factor*(kt)


          ! fractionx=(dble(pu2)+0.5)/(dble(Kt)+0.5)
          ! factor = olaps(bq*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)), npd, mpd, &
          !         & bq*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)), nkd, -mkd, &
          !         & (dble(pd)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pd)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)), npu1, mpu1, &
          !         & bq*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)), nku1, -mku1, &
          !         & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect), (dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
          !         & olaps(bq*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)), npu2, mpu2, &
          !         & bq*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)), nku2, -mku2, &
          !         & -(1.0D0-fractionx)*(dx+dxcorrect), -(1.0D0-fractionx)*dy) &
          !         & *(0.D0-1.0D0)**(mku1+mku2+mkd+1)*2.0D0/(dx+dxcorrect)


          !     distribution_E_loc(pu2+1)=distribution_E_loc(pu2+1)+eigenvaluep*eigenvaluek*factor*(kt)

        enddo

      enddo

    enddo

  enddo

  do k_u = 1, Kt

    distributionfunctionsvalue(1, k_u)=(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(2, k_u)=dx**2+dy**2
    distributionfunctionsvalue(3, k_u)=distribution_H(k_u)!*(dble(k_u-1)+0.5)/(dble(Kt)+0.5)
    distributionfunctionsvalue(4, k_u)=distribution_E(k_u)

  enddo

  Print*,"GPD_for_sea2_quark"
  do test = 1, kt
    Print*, distributionfunctionsvalue(1, test), distributionfunctionsvalue(2, test), distributionfunctionsvalue(3, test), &
      & distributionfunctionsvalue(4, test)!,distributionfunctionsvalue(5, test), distributionfunctionsvalue(6, test), &
    !& distributionfunctionsvalue(7, test)
  enddo

  return
end subroutine GPD_for_sea2_quark3

subroutine TalmiSquare(nmax, n1, m1, n2, m2, n3, m3, kp1half, kp2half, kp3half, j, TalmisquareValue)

  implicit none

  integer:: nmax, n1, m1, n2, m2, n3, m3
  double precision:: kp1half, kp2half, kp3half
  double precision, dimension(7, *):: TalmiSquareValue

  integer:: mscale, nscale
  integer:: capnf, capmf, nf, mf, capnm, capmm, nm, mm
  double precision:: TMC
  integer:: error, i, j

  mscale = nmax-3
  nscale = floor(dble(mscale)/2.0D0)
  j = 0

  do capnf = 0, nscale
    do capmf = -mscale, mscale
      do capnm = 0, nscale
        do capmm = -mscale, mscale

          mm = m1+m2-capmm
          nm = n1+n2-capnm+(abs(m1)+abs(m2)-abs(capmm)-abs(mm))/2
          mf = capmm+m3-capmf
          nf = capnm+n3-capnf+(abs(capmm)+abs(m3)-abs(capmf)-abs(mf))/2

          if(nm .ge. 0 .and. nf .ge. 0) then

            j = j+1

            TalmiSquareValue(1, j)=capnf
            TalmiSquareValue(2, j)=capmf
            TalmiSquareValue(3, j)=nf
            TalmiSquareValue(4, j)=mf
            TalmiSquareValue(5, j)=nm
            TalmiSquareValue(6, j)=mm

            TalmiSquareValue(7, j)=TMC(capnm, capmm, nm, mm, n1, m1, n2, m2, sqrt(kp2half/kp1half))*TMC(capnf, capmf, nf, mf, &
              & n3, m3, capnm, capmm, sqrt((kp1half+kp2half)/kp3half))

          endif

        enddo
      enddo
    enddo
  enddo

  return
end subroutine TalmiSquare

subroutine average(loopnumber1, distributionfunctionskt, distributionfunctionskt_1, distributionfunctions)

  implicit none

  integer:: loopnumber1
  double precision, dimension(3, *):: distributionfunctionskt, distributionfunctionskt_1
  double precision, dimension(3, *):: distributionfunctions

  integer:: i, test

  do i = 1, loopnumber1
    distributionfunctions(1, i)=(distributionfunctionskt(1, i)+distributionfunctionskt_1(1, i))/2.D0
    distributionfunctions(2, i)=distributionfunctionskt(2, i)
    distributionfunctions(3, i)=(distributionfunctionskt(3, i)+distributionfunctionskt_1(3, i))/2.0D0
  enddo

  !do test = 1, loopnumber1
  !    Print*,distributionfunctions(1, test), distributionfunctions(2, test), distributionfunctions(3, test)
  !enddo

  return
end subroutine average
