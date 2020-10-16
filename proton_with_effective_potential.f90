program proton_with_effective_potential
  use, intrinsic:: iso_fortran_env
  use basis_info
  use numbers
  use colorfactor

  implicit none
  ! c==============time parameter of system_clock======
  real(kind = 4):: udsgTimeStart, udsgTimeEnd  !To record the truncation of time.
  ! c==============time parameter of system_clock======

  integer:: nmax1, nmax2, nmax3  ! For three Fock Sector
  integer:: nmax, Mj, Kt  ! input parameter
  double precision:: lag, deltamax  ! input parameter
  integer:: nestate  ! nestate is the number of eigenvalue
  double precision:: b  ! harmonic oscillator scale

  ! integer:: KroneckerDelta  ! external function
  ! double precision:: ifactoreval, adota, ifactor, adotb, olaps  ! external function
  double precision, dimension(:), allocatable:: ev1, ev2  ! defined the variable to calculate the mirror parity
  double precision, dimension(:,:), allocatable:: evector1, evector2  !defined the variable to calculate the mirror parity
  integer:: selectquark
  integer:: setparameter     !swatch parameter of inputting parameter
  double precision:: massoffset, nmcutoff, xcutoff  ! the parameter correspond to the mass renormalization
  integer:: error, method
  !the various about reading from the terminal
  integer:: nargs, ii, nmpstate, i
  character(len = 30), dimension(30):: arg
  double precision:: test
  double precision, external:: TMC
  double precision:: deltam_in(2), deltam_out(2)


  call cpu_time(udsgTimeStart)

  !Mj = 0
  !b = ME_DEFAULT*G_DEFAULT**2/4.0D0/PI
  b = B_DEFAULT
  lag = 10.1D0
  nestate = 1
  selectquark = 1
  massoffset = 0.D0
  nmcutoff = 0.D0
  xcutoff = 0.D0
  select_distribution_function = 0
  setparameter = 0
  method = -1
  seqeigen = 1
  coupling3 = coupling_default3
  coupling5 = coupling_default5
  coupling3t5 = coupling_default3t5

  nargs = command_argument_count()
  do ii = 1, nargs
    call get_command_argument(ii, arg(ii))
  end do
  if(nargs .lt. 2) then
    print *, "Too few parameters."
    print *, "Examples: ./pos.out 8 10  for nmax = 8, ktot = 10, Mj = 0"
    print *, "Examples: ./pos.out 8 10 1 for nmax = 8, ktot = 10, Mj = 1"
    stop
  end if

  do ii = 1, nargs
    select case (ii)
    case(1)
      read(arg(ii), *) nmax  !harmonic oscillator truncation
    case(2)
      read(arg(ii), *) Kt    !longitudinal momentum truncation
    case(3)
      read(arg(ii), *) Mj    !the total spin of proton at z direction
    case default

      if(arg(ii).eq."-lag") then
        setparameter = 1
        cycle
      endif
      if(arg(ii).eq."-numeigenvalue") then
        setparameter = 2
        cycle
      endif
      if(arg(ii).eq."-singlequark") then
        select_distribution_function = 1
        cycle
      endif
      if(arg(ii).eq."-wavefunction") then
        select_distribution_function = 2
        cycle
      endif

      if(arg(ii).eq."-selectmethod") then
        setparameter = 6
        cycle
      endif

      if(arg(ii).eq."-seqeigen") then
        setparameter = 7
        cycle
      endif
      if(arg(ii).eq."-coupling3p") then
        setparameter = 8
        cycle
      endif
      if(arg(ii).eq."-coupling5p") then
        setparameter = 9
        cycle
      endif
      if(arg(ii).eq."-coupling3t5p") then
        setparameter = 10
        cycle
      endif

      !!!!!!!!!!! set parameter of GPD. !!!!!!!!!!!!!

      if(arg(ii).eq."-FF") then
        !setparameter = 4
        select_distribution_function = 3
        cycle
      endif
      if(arg(ii).eq."-selectquark") then
        setparameter = 3
        cycle
      endif
      if(arg(ii).eq."-deltamax") then
        setparameter = 5
        cycle
      endif
    end select

    select case (setparameter)
    case(1)
      read(arg(ii), *) lag
    case(2)
      read(arg(ii), *) nestate
    case(3)
      if(arg(ii).eq."-u") then
        selectquark = 2
      else if (arg(ii).eq."-d") then
        selectquark = 1
      endif
    case(4)
      select_distribution_function = 2
    case(5)
      read(arg(ii), *) deltamax
    case(6)
      read(arg(ii), *) method
    case(7)
      read(arg(ii), *) seqeigen
    case(8)
      read(arg(ii), *) coupling3
    case(9)
      read(arg(ii), *) coupling5
    case(10)
      read(arg(ii), *) coupling3t5
    case default
    endselect

  end do

  nmax1 = nmax
  nmax2 = nmax
  nmax3 = nmax
  Print *,"Truncation parameter for three Fock Sector"
  print *, "Nmax1=",nmax1, " || "," Nmax2=",nmax2, " || ","Nmax3=",nmax3

  select case (select_distribution_function)

  case(1)

    call renormalization(nmax1, nmax2, nmax3, Mj, Kt, method, 1, nestate, deltam_in, deltam_out)

  case(2)

    call dimtotal3p(nmax1, Mj, Kt, dimtot1)
    call dimtotal4p(nmax2, Mj, Kt, dimtot2)
    call dimtotal5p(nmax3, Mj, Kt, dimtot3)
    call ColorMatrix()

    If(method .eq. 1) then
      dimtot = dimtot1+3*dimtot3
    else if(method .eq. 2) then
      dimtot = dimtot1+2*dimtot2+3*dimtot3
    else
      dimtot = dimtot1+3*dimtot3
    endif

    allocate(ev1(nestate), evector1(nestate, dimtot), stat = error)
    if(error .ne. 0) then
      print *, 'cannot allocate basis array'
    end if

    call wfproduction(nmax1, nmax2, nmax3, Mj, Kt, method, massoffset, nmcutoff, xcutoff, nestate, ev1, evector1)
    !call wfproduction(nmax1, nmax2, -Mj, Kt, massoffset, nmcutoff, xcutoff, lag, b, mu, md, ms, mg, nestate, dimtot, ev2, evector2)

    !call mparity(nmax1, nmax2, Mj, Kt, nestate, dimtot, ev1, evector1, ev2, evector2)
    !call mparity1(nmax1, nmax2, Mj, Kt, nestate, dimtot, ev1, evector1, ev2, evector2)

  case(3)

    !Print*,"test"

    call distribution_functions_FF(nmax1, nmax2, nmax3, Mj, Kt, nestate, selectquark, deltamax, lag, method)


  case default
  endselect

  call cpu_time(udsgTimeEnd)
  print*," "
  write(*,*) "TimeUsed:", nint((udsgTimeEnd-udsgTimeStart)/60.0D0), "mins"
  write(*,*) "It will show 0 mins when the process is less than 1 min!"

end program proton_with_effective_potential
