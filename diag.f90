
       subroutine diag_arpack(nestate,dim,method,renorm_flag,initial_vector,verbosesw,ev,evector)

         ! c-----------------------------------------------------------------------
         ! c
         ! c     %------------------------------------------------------%
         ! c     | Storage Declarations:                                |
         ! c     |                                                      |
         ! c     | The maximum dimensions for all arrays are            |
         ! c     | set here to accommodate a problem size of            |
         ! c     | N .le. MAXN                                          |
         ! c     |                                                      |
         ! c     | NEV is the number of eigenvalues requested.          |
         ! c     |     See specifications for ARPACK usage below.       |
         ! c     |                                                      |
         ! c     | NCV is the largest number of basis vectors that will |
         ! c     |     be used in the Implicitly Restarted Arnoldi      |
         ! c     |     Process.  Work per major iteration is            |
         ! c     |     proportional to N*NCV*NCV.                       |
         ! c     |                                                      |
         ! c     | You must set:                                        |
         ! c     |                                                      |
         ! c     | MAXN:   Maximum dimension of the A allowed.          |
         ! c     | MAXNEV: Maximum NEV allowed.                         |
         ! c     | MAXNCV: Maximum NCV allowed.                         |
         ! c     %------------------------------------------------------%
         ! c
         implicit none
         integer  verbosesw,method
         integer  ii,ev_index
         integer  maxn, maxnev, maxncv, ldv, nestate
         integer  error,i
         integer  renorm_flag
         ! c      parameter       (maxn=256, maxnev=10, maxncv=25,
         ! c     $                 ldv=maxn )
         ! parameter       (maxnev=1,maxncv=50)
         ! c      parameter  (nz_index=37911)
         ! parameter  (nz=10000000,nz2=10000000)

         integer  dim
         double precision, dimension(nestate) :: ev
         double precision, dimension(nestate,dim) :: evector
         double precision, dimension(dim) :: initial_vector
         ! double precision  V_anomm,normtest,snorm
         ! double precision  V_anomm2
         integer  msaupd, msaup2, msaitr, logfil, ndigit

         ! c
         ! c     %--------------%
         ! c     | Local Arrays |
         ! c     %--------------%
         ! c
         double precision,dimension(:), allocatable :: workl
         double precision,dimension(:,:), allocatable :: d
         double precision, dimension(:), allocatable :: resid,ax
         double precision, dimension(:), allocatable :: workd
         double precision, dimension(:,:), allocatable :: v
         logical,dimension(:), allocatable :: select
         integer          iparam(11), ipntr(11)
         ! c
         ! c     %---------------%
         ! c     | Local Scalars |
         ! c     %---------------%
         ! c
         character        bmat*1, which*2
         integer          ido, n, nev, ncv, lworkl, info, ierr, j, ishfts, maxitr, mode1, nconv
         logical          rvec
	!integer            rvec
         double precision  tol, sigma
         ! c
         ! c     %------------%
         ! c     | Parameters |
         ! c     %------------%
         ! c
         double precision  zero
         parameter        (zero = 0.0D0)
         ! c
         ! c     %-----------------------------%
         ! c     | BLAS & LAPACK routines used |
         ! c     %-----------------------------%
         ! c
         double precision dnrm2
         ! c
         ! c     %--------------------%
         ! c     | Intrinsic function |
         ! c     %--------------------%
         ! c
         maxn=dim
         maxnev=nestate
         maxncv=min(maxnev*5,maxn)
	!maxnev=20
	!maxncv=85
         ldv=maxn

         allocate(resid(maxn),ax(maxn),workd(3*maxn),d(maxncv,2),workl(maxncv*(maxncv+8)),v(maxn,maxncv),&
            & select(maxncv),stat=error)
         if(error.ne.0) then
            print *, 'cannot allocate work array in diag',maxn,maxncv
         end if

!!$         allocate(resid(maxn),stat=error)
!!$         if(error.ne.0) then
!!$            print *, 'cannot allocate array'
!!$         end if
!!$         allocate(ax(maxn),stat=error)
!!$         if(error.ne.0) then
!!$            print *, 'cannot allocate array'
!!$         end if
!!$         allocate(workd(3*maxn),stat=error)
!!$         if(error.ne.0) then
!!$            print *, 'cannot allocate array'
!!$         end if
!!$         allocate(v(ldv,maxncv),stat=error)
!!$         if(error.ne.0) then
!!$            print *, 'cannot allocate array'
!!$         end if
         ! c
         ! c     %-----------------------%
         ! c     | Executable Statements |
         ! c     %-----------------------%
         ! c
         ! c     %-------------------------------------------------%
         ! c     | The following include statement and assignments |
         ! c     | initiate trace output from the internal         |
         ! c     | actions of ARPACK.  See debug.doc in the        |
         ! c     | DOCUMENTS directory for usage.  Initially, the  |
         ! c     | most useful information will be a breakdown of  |
         ! c     | time spent in the various stages of computation |
         ! c     | given by setting msaupd = 1.                    |
         ! c     %-------------------------------------------------%
         ! c
         ! c      include 'debug.h'
         ndigit = -6
	 !ndigit=-3
         logfil = 6
         ! c      msgets = 0
         msaitr = 0
         ! c      msapps = 0
         msaupd = 1
         msaup2 = 0
         ! c      mseigt = 0
         ! c      mseupd = 0

         ! c     %-----------------------------------------------%
         ! c     |                                               |
         ! c     | Specifications for ARPACK usage are set       |
         ! c     | below:                                        |
         ! c     |                                               |
         ! c     |    1) NEV = 4  asks for 4 eigenvalues to be   |
         ! c     |       computed.                               |
         ! c     |                                               |
         ! c     |    2) NCV = 20 sets the length of the Arnoldi |
         ! c     |       factorization                           |
         ! c     |                                               |
         ! c     |    3) This is a standard problem              |
         ! c     |         (indicated by bmat  = 'I')            |
         ! c     |                                               |
         ! c     |    4) Ask for the NEV eigenvalues of          |
         ! c     |       largest magnitude                       |
         ! c     |         (indicated by which = 'LM')           |
         ! c     |       See documentation in DSAUPD for the     |
         ! c     |       other options SM, LA, SA, LI, SI.       |
         ! c     |                                               |
         ! c     | Note: NEV and NCV must satisfy the following  |
         ! c     | conditions:                                   |
         ! c     |              NEV <= MAXNEV                    |
         ! c     |          NEV + 1 <= NCV <= MAXNCV             |
         ! c     %-----------------------------------------------%
         ! c
         ! c     %-------------------------------------------------%
         ! c     | The following sets dimensions for this problem. |
         ! c     %-------------------------------------------------%
         n = maxn
         nev = maxnev
         ! c  NCV must be greater than NEV and less than or equal to N.
         ! if(n.ge.60) then
         !    ncv   = maxncv
         ! else
         !    ncv   = 5
         ! end if
         ncv=maxncv
         bmat  = 'I'
         which = 'SA'

         ! c      nz_size=nz_index
         ! c
         if ( n .gt. maxn ) then
            print *, ' ERROR with _SSIMP: N is greater than MAXN '
            go to 9000
         else if ( nev .gt. maxnev ) then
            print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
            go to 9000
         else if ( ncv .gt. maxncv ) then
            print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
            go to 9000
         end if
         ! c
         ! c     %-----------------------------------------------------%
         ! c     |                                                     |
         ! c     | Specification of stopping rules and initial         |
         ! c     | conditions before calling DSAUPD                    |
         ! c     |                                                     |
         ! c     | TOL  determines the stopping criterion.             |
         ! c     |                                                     |
         ! c     |      Expect                                         |
         ! c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
         ! c     |               computed   true                       |
         ! c     |                                                     |
         ! c     |      If TOL .le. 0,  then TOL <- macheps            |
         ! c     |           (machine precision) is used.              |
         ! c     |                                                     |
         ! c     | IDO  is the REVERSE COMMUNICATION parameter         |
         ! c     |      used to specify actions to be taken on return  |
         ! c     |      from DSAUPD. (See usage below.)                |
         ! c     |                                                     |
         ! c     |      It MUST initially be set to 0 before the first |
         ! c     |      call to DSAUPD.                                |
         ! c     |                                                     |
         ! c     | INFO on entry specifies starting vector information |
         ! c     |      and on return indicates error codes            |
         ! c     |                                                     |
         ! c     |      Initially, setting INFO=0 indicates that a     |
         ! c     |      random starting vector is requested to         |
         ! c     |      start the ARNOLDI iteration.  Setting INFO to  |
         ! c     |      a nonzero value on the initial call is used    |
         ! c     |      if you want to specify your own starting       |
         ! c     |      vector (This vector must be placed in RESID.)  |
         ! c     |                                                     |
         ! c     | The work array WORKL is used in DSAUPD as           |
         ! c     | workspace.  Its dimension LWORKL is set as          |
         ! c     | illustrated below.                                  |
         ! c     |                                                     |
         ! c     %-----------------------------------------------------%
         ! c
         lworkl = ncv*(ncv+8)
         tol = 1.0D-20
         info = renorm_flag
         ido = 0
         ! c
         ! c     %---------------------------------------------------%
         ! c     | Specification of Algorithm Mode:                  |
         ! c     |                                                   |
         ! c     | This program uses the exact shift strategy        |
         ! c     | (indicated by setting PARAM(1) = 1).              |
         ! c     | IPARAM(3) specifies the maximum number of Arnoldi |
         ! c     | iterations allowed.  Mode 1 of DSAUPD is used     |
         ! c     | (IPARAM(7) = 1). All these options can be changed |
         ! c     | by the user. For details see the documentation in |
         ! c     | DSAUPD.                                           |
         ! c     %---------------------------------------------------%
         ! c
         ishfts = 1
         maxitr = 300000
         mode1 = 1
         ! c
         iparam(1) = ishfts
         ! c
         iparam(3) = maxitr
         ! c
         iparam(7) = mode1

         ! c     %------------------------------------------------%
         ! c     | M A I N   L O O P (Reverse communication loop) |
         ! c     %------------------------------------------------%
         ! c
         resid=initial_vector
         i=0
10       continue
         ! c
         ! c        %---------------------------------------------%
         ! c        | Repeatedly call the routine DSAUPD and take |
         ! c        | actions indicated by parameter IDO until    |
         ! c        | either convergence is indicated or maxitr   |
         ! c        | has been exceeded.                          |
         ! c        %---------------------------------------------%
         ! c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
              ncv, v, ldv, iparam, ipntr, workd, workl,&
              lworkl, info )
         ! c
         if (ido .eq. -1 .or. ido .eq. 1) then
	 !if (ido.ne.99) then
            ! c
            ! c           %--------------------------------------%
            ! c           | Perform matrix vector multiplication |
            ! c           |              y <--- OP*x             |
            ! c           | The user should supply his/her own   |
            ! c           | matrix vector multiplication routine |
            ! c           | here that takes workd(ipntr(1)) as   |
            ! c           | the input, and return the result to  |
            ! c           | workd(ipntr(2)).                     |
            ! c           %--------------------------------------%
            ! c
		!print*,workd(ipntr(1))
            call av (workd(ipntr(1)), workd(ipntr(2)), method)
            !print*,ipntr(1),ipntr(2)
            i=i+1
            !Print*,i
            ! c
            ! c           %-----------------------------------------%
            ! c           | L O O P   B A C K to call DSAUPD again. |
            ! c           %-----------------------------------------%
            ! c
            go to 10
            ! c
         end if
         ! c
         ! c     %----------------------------------------%
         ! c     | Either we have convergence or there is |
         ! c     | an error.                              |
         ! c     %----------------------------------------%
         ! c
         if ( info .lt. 0 ) then
            ! c
            ! c        %--------------------------%
            ! c        | Error message. Check the |
            ! c        | documentation in DSAUPD. |
            ! c        %--------------------------%
            ! c
            print *, ' '
            print *, ' Error with _saupd, info = ', info
            print *, ' Check documentation in _saupd '
            print *, ' '
            ! c
         else
            ! c
            ! c        %-------------------------------------------%
            ! c        | No fatal errors occurred.                 |
            ! c        | Post-Process using DSEUPD.                |
            ! c        |                                           |
            ! c        | Computed eigenvalues may be extracted.    |
            ! c        |                                           |
            ! c        | Eigenvectors may be also computed now if  |
            ! c        | desired.  (indicated by rvec = .true.)    |
            ! c        |                                           |
            ! c        | The routine DSEUPD now called to do this  |
            ! c        | post processing (Other modes may require  |
            ! c        | more complicated post processing than     |
            ! c        | mode1.)                                   |
            ! c        |                                           |
            ! c        %-------------------------------------------%
            ! c
            rvec = .true.
	    !rvec = 1
            ! c
            call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
                 bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                 iparam, ipntr, workd, workl, lworkl, ierr )
            ! c
            ! c         %----------------------------------------------%
            ! c         | Eigenvalues are returned in the first column |
            ! c         | of the two dimensional array D and the       |
            ! c         | corresponding eigenvectors are returned in   |
            ! c         | the first NCONV (=IPARAM(5)) columns of the  |
            ! c         | two dimensional array V if requested.        |
            ! c         | Otherwise, an orthogonal basis for the       |
            ! c         | invariant subspace corresponding to the      |
            ! c         | eigenvalues in D is returned in V.           |
            ! c         %----------------------------------------------%
            ! c
            if ( ierr .ne. 0) then
               ! c
               ! c            %------------------------------------%
               ! c            | Error condition:                   |
               ! c            | Check the documentation of DSEUPD. |
               ! c            %------------------------------------%
               ! c
               print *, ' '
               print *, ' Error with _seupd, info = ', ierr
               print *, ' Check the documentation of _seupd. '
               print *, ' '
               ! c
            else
               ! c
               nconv =  iparam(5)
               do j=1, nconv
                  ! c
                  ! c               %---------------------------%
                  ! c               | Compute the residual norm |
                  ! c               |                           |
                  ! c               |   ||  A*x - lambda*x ||   |
                  ! c               |                           |
                  ! c               | for the NCONV accurately  |
                  ! c               | computed eigenvalues and  |
                  ! c               | eigenvectors.  (iparam(5) |
                  ! c               | indicates how many are    |
                  ! c               | accurate to the requested |
                  ! c               | tolerance)                |
                  ! c               %---------------------------%
                  ! c
                  call av(v(1,j), ax, method)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,2) = dnrm2(n, ax, 1)
                  d(j,2) = d(j,2) / abs(d(j,1))

                  ! c                write(*,*) 'j',v(3,j)
                  ! c
               end do
               ! c
               ! c            %-----------------------------%
               ! c            | Display computed residuals. |
               ! c            %-----------------------------%
               ! c
!!$               call dmout(6, nconv, 2, d, maxncv, -6,'Ritz values and relative residuals')

               ! c             write(*,*) 'E_1',d(1,1),d(1,2)

               ! lowest=d(1,1)


               do ev_index=1,maxnev

                  ev(ev_index)=d(ev_index,1)
                  do ii = 1, n
                    evector(ev_index,ii)=v(ii,ev_index)
                  end do

               end do

               ! V_anomm=0.0D0
               ! V_anomm2=0.0D0
               ! normtest=0.0D0
               ! snorm=0.0D0

               ! do ii=1, 1
               !    snorm=snorm+(v(ii,1))**2.D0
               ! end do

               ! V_anomm=sqrt(V_anomm)
               ! normtest=sqrt(normtest)
!!$               write(*,*) 'normtest',normtest


            end if
            ! c
            ! c         %-------------------------------------------%
            ! c         | Print additional convergence information. |
            ! c         %-------------------------------------------%
            ! c
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit',&
                    ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
            ! c
!!$            print *, ' '
            if(verbosesw.eq.1) then
               print *, '## _SSIMP '
               print *, '## ====== '
               print *, '## '
               print *, '## Size of the matrix is ', n
               print *, '## The number of Ritz values requested is ', nev
               print *, '## The number of Arnoldi vectors generated (NCV) is ', ncv
               print *, '## What portion of the spectrum: ', which
               print *, '## The number of converged Ritz values is ', nconv
               print *, '## The number of Implicit Arnoldi update iterations taken is ', iparam(3)
               print *, '## The number of OP*x is ', iparam(9)
               print *, '## The convergence criterion is ', tol
               print *, '## '
            end if
            ! c
         end if
         ! c
         ! c     %---------------------------%
         ! c     | Done with program dssimp. |
         ! c     %---------------------------%
         ! c



9000     continue

         deallocate(resid,ax,workd,d,workl,v,select)

!!$         deallocate(resid)
!!$         deallocate(ax)
!!$         deallocate(workd)
!!$         deallocate(v)
         return
       end subroutine diag_arpack


       ! c#######################################################################
       ! c        OMA VIRITELMA
       ! c#######################################################################


       subroutine av (v, w, method)

         use basis_info
         use hamiMatrix

         implicit none
         double precision :: v(*), w(dimtot)
         integer :: i,j,method
         integer :: dim

         w=0.D0

         If(method.eq.1) then
            dim=dimtot1
         else if(method.eq.2) then
            dim=dimtot1+2*dimtot2
            !dim=0
         else
            dim=dimtot1
            !dim=0
         endif

         !!!!!!!!!!!!!!! the hamiltonian matrix multiply the vector for the first Fock sector !!!!!!!!!!!!!!!!!!!!!!!
         do i=1,dimtot1
            do j=i_nzk(i),i_nzk(i+1)-1
               w(i)=w(i)+hamiltoniankineticvalue(j)*v(j_nzk(j))
            end do
         end do

         do i=1,dimtot1
            do j=i_nztsw(i),i_nztsw(i+1)-1
               w(i)=w(i)+hamiltonianinteractswvalue(j)*v(j_nztsw(j))
            end do
         end do

         do i=1,dimtot1
            do j=i_nzlsw(i),i_nzlsw(i+1)-1
               w(i)=w(i)+hamiltonianinteractlongivalue(j)*v(j_nzlsw(j))
            end do
         end do



         ! !!!!!!!!!!!!!!! the hamiltonian matrix multiply the vector for the Five particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!
         do i=1,3*dimtot3
         !do i=1,dimtot3
            do j=i_nzk3(i),i_nzk3(i+1)-1
               w(dim+i)=w(dim+i)+hamiltoniankineticvalue3(j)*v(dim+j_nzk3(j))
            end do
         end do

         do i=1,3*dimtot3
            do j=i_nztsw3(i),i_nztsw3(i+1)-1
               w(dim+i)=w(dim+i)+hamiltonianinteractswvalue3(j)*v(dim+j_nztsw3(j))
            end do
         end do

         do i=1,3*dimtot3
            do j=i_nzlsw3(i),i_nzlsw3(i+1)-1
               w(dim+i)=w(dim+i)+hamiltonianinteractlongivalue3(j)*v(dim+j_nzlsw3(j))
            end do
         end do

         !!!!!!!!!!!!!!!! the instantaneous interaction of three particle to five particle multiply the vector !!!!!!!!!!!!!!!!

         do i=1,dimtot1
            do j=i_nz_inst3t5(i),i_nz_inst3t5(i+1)-1
               w(i)=w(i)+hamiltonianinteractinst3t5(j)*v(dim+j_nz_inst3t5(j))
               w(dim+j_nz_inst3t5(j))=w(dim+j_nz_inst3t5(j))+hamiltonianinteractinst3t5(j)*v(i)
            end do
         end do

         If(method.eq.1) then ! method 1 corresponding to the effective interaction

            do i=1,dimtot1
               do j=i_nz_oge(i),i_nz_oge(i+1)-1
                  w(i)=w(i)+hamiltonianinteractogevalue(j)*v(j_nz_oge(j))
               end do
            end do

            do i=1,3*dimtot3
               do j=i_nz_oge3(i),i_nz_oge3(i+1)-1
                  w(dimtot1+i)=w(dimtot1+i)+hamiltonianinteractogevalue3(j)*v(dimtot1+j_nz_oge3(j))
                  !w(i)=w(i)+hamiltonianinteractogevalue3(j)*v(j_nz_oge3(j))
               end do
            end do

            do i=1,dimtot1
               do j=i_nz_oge3t5(i),i_nz_oge3t5(i+1)-1
                  w(i)=w(i)+hamiltonianinteractogevalue3t5(j)*v(dimtot1+j_nz_oge3t5(j))
                  w(dimtot1+j_nz_oge3t5(j))=w(dimtot1+j_nz_oge3t5(j))+hamiltonianinteractogevalue3t5(j)*v(i)
               end do
            end do

         endif

         If(method.eq.2) then ! method2 corresponding to the fundamental interaction

            do i=1,dimtot1
               do j=i_nz_inst(i),i_nz_inst(i+1)-1
                  w(i)=w(i)+hamiltonianinteractinstvalue(j)*v(j_nz_inst(j))
               end do
            end do

            do i=1,2*dimtot2
               do j=i_nzk2(i),i_nzk2(i+1)-1
                  w(dimtot1+i)=w(dimtot1+i)+hamiltoniankineticvalue2(j)*v(dimtot1+j_nzk2(j))
               end do
            end do

            do i=1,dimtot1
               do j=i_nz_vcvqtqg(i),i_nz_vcvqtqg(i+1)-1
                  w(i)=w(i)+hamiltonianinteractvcvqtqg(j)*v(dimtot1+j_nz_vcvqtqg(j))
                  w(dimtot1+j_nz_vcvqtqg(j))=w(dimtot1+j_nz_vcvqtqg(j))+hamiltonianinteractvcvqtqg(j)*v(i)
               end do
            end do

            do i=1,2*dimtot2
               do j=i_nz_vcvgtqq(i),i_nz_vcvgtqq(i+1)-1
                  w(dimtot1+i)=w(dimtot1+i)+hamiltonianinteractvcvgtqq(j)*v(dim+j_nz_vcvgtqq(j))
                  w(dim+j_nz_vcvgtqq(j))=w(dim+j_nz_vcvgtqq(j))+hamiltonianinteractvcvgtqq(j)*v(dimtot1+i)
               end do
            end do

         endif

         return
       end subroutine av
