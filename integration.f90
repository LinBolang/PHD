program integration

    implicit none

    integer :: row,col,nmax,Kt,model
    double precision :: mass1(3),mass2(5),massg(2),b
    double precision,dimension(:),allocatable :: integral3,integral5,integral3to5
    integer :: error,ii,nargs
    integer :: integral_dim3, integral_dim5, integral_dim3to5
    character(len=30),dimension(30) :: arg

    mass1(1)=0.2    !First Fock Sector quark mass
    mass1(2)=0.2    !First Fock Sector quark mass
    mass1(3)=0.2    !First Fock Sector quark mass

    massg(1)=0.05   !First Fock Sector gluon mass

    mass2(1)=0.15   !Second Fock Sector quark mass
    mass2(2)=0.15   !Second Fock Sector quark mass
    mass2(3)=0.15   !Second Fock Sector quark mass
    mass2(4)=0.15   !Second Fock Sector quark mass
    mass2(5)=0.15   !Second Fock Sector quark mass

    massg(2)=0.05   !Second Fock Sector gluon mass

    b=0.6

    model=4

    nargs=command_argument_count()
    do ii=1, nargs
        call get_command_argument(ii,arg(ii))
    end do

    do ii=1,nargs
        select case (ii)
        case(1)
            read(arg(ii),*) nmax  !harmonic oscillator truncation
        case(2)
            read(arg(ii),*) Kt    !longitudinal momentum truncation
        case(3)
            read(arg(ii),*) model   !1 for produce 3 particles Fock sector integration list, 2 for produce 5 particles sector integration list, 
                                    !3 for produce 3 particles to 5 particles integration list, 4 for produce all
        case default
        end select
    enddo

    If(model.eq.1) then

        row=(nmax)**2
        col=kt*(kt+1)*(2*kt+1)/6

        integral_dim3=row*col*3

        Print*, integral_dim3

        allocate(integral3(integral_dim3),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

        call IEDintegration3p(nmax,Kt,mass1,massg(1),b,integral3)

        call output3p(nmax,Kt,mass1,massg(1),b,integral_dim3,integral3)

        deallocate(integral3)

    else if(model.eq.2) then

        row=(nmax-2)**2
        col=kt*(kt+1)*(2*kt+1)/6

        integral_dim5=row*col*10

        Print*, integral_dim5

        allocate(integral5(integral_dim5),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

        call IEDintegration5p(nmax,Kt,mass2,massg(2),b,integral5)

        call output5p(nmax,Kt,mass2,massg(2),b,integral_dim5,integral5)

        deallocate(integral5)

    else if(model.eq.3) then

        row=(nmax-4)**2
        col=Kt*(kt-1)*(kt+1)/6

        integral_dim3to5=row*col*3

        allocate(integral3to5(integral_dim3to5),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

        print*,integral_dim3to5

        call IEDintegration3to5(nmax,Kt,mass1,mass2,b,integral3to5)

        call output3to5(nmax,Kt,mass1,mass2,b,integral_dim3to5,integral3to5)

        deallocate(integral3to5)

    else

        integral_dim3=(nmax)**2*kt*(kt+1)*(2*kt+1)*3/6
        integral_dim5=(nmax-2)**2*kt*(kt+1)*(2*kt+1)*10/6
        integral_dim3to5=(nmax-4)**2*Kt*(kt-1)*(kt+1)*3/6

        print*,integral_dim3,integral_dim5,integral_dim3to5

        allocate(integral3(integral_dim3),integral5(integral_dim5),integral3to5(integral_dim3to5),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

        call IEDintegration3p(nmax,Kt,mass1,massg(1),b,integral3)

        call output3p(nmax,Kt,mass1,massg(1),b,integral_dim3,integral3)

        call IEDintegration5p(nmax,Kt,mass2,massg(2),b,integral5)

        call output5p(nmax,Kt,mass2,massg(2),b,integral_dim5,integral5)

        call IEDintegration3to5(nmax,Kt,mass1,mass2,b,integral3to5)

        call output3to5(nmax,Kt,mass1,mass2,b,integral_dim3to5,integral3to5)

        deallocate(integral3,integral5,integral3to5)

    endif

end program integration

subroutine IEDintegration3p(nmax2,Kt,mass1,massg,b,integration)

    implicit none

    integer :: nmax2,Kt
    double precision :: mass1(3),massg,b
    double precision, dimension(*) :: integration

    double precision :: x1,x2,x1prime,x2prime
    integer :: kp1,kp2,kk1,kk2,i,j,capNprime,LagN,selectmass
    double precision :: c1,c2,delta(3)

    integer :: nvec,flags,method,nregions,neval,info
    integer,external :: integrand3p
    double precision :: relerr,abserr,para(5),Integral(1),Error(1),prob(1)
    integer :: mineval,maxeval

    i=0
    flags=0
    ! mineval=0.0D0
    ! maxeval=0.99999999D0
    mineval=0
    maxeval=100000000
    method=13
    relerr=1.0D-14
    abserr=1.0D-14



    !print*,"test"
    do capNprime=0,nmax2-1
        do LagN=0,nmax2-1
            do kp1=0,Kt-1
                do kp2=0,Kt-1-kp1
    ! capNprime=1
    ! LagN=0
    ! kp1=0
    ! kp2=0
    ! kk1=0
                    x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
                    x2=(dble(kp2)+0.5)/(dble(Kt)+0.5)
                    do kk1=0,kp1+kp2

                        kk2=kp1+kp2-kk1
    
                        x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
                        x2prime=(dble(kk2)+0.5)/(dble(Kt)+0.5)
                        c1=(sqrt(x1prime*x2)-sqrt(x1*x2prime))**2/(x1+x2)
                        c2=(sqrt(x1prime*x2)+sqrt(x1*x2prime))**2/(x1+x2)
                        delta(1)=(x1-x1prime)**2*(mass1(1)**2/x1/x1prime+mass1(2)**2/x2/x2prime)+2*massg**2
                        delta(2)=(x1-x1prime)**2*(mass1(1)**2/x1/x1prime+mass1(3)**2/x2/x2prime)+2*massg**2
                        delta(3)=(x1-x1prime)**2*(mass1(2)**2/x1/x1prime+mass1(3)**2/x2/x2prime)+2*massg**2

                        do selectmass=1,3
                            i=i+1

                            para(1)=dble(capNprime)
                            para(2)=dble(LagN)
                            para(3)=c1
                            para(4)=c2
                            para(5)=delta(selectmass)/b**2

                            call cuhre(2,1, integrand3p, para, 1, relerr, abserr, flags, mineval, maxeval,&
                                & method, "", -1, nregions, neval, info,  Integral, Error, prob)

                            If(info.eq.1) then
                                print*,Integral(1),error(1),info,i,j
                            endif

                            integration(i)=Integral(1)

                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


   return
end subroutine IEDintegration3p

integer function integrand3p(ndim,x,ncomp,f,para)
    implicit none
    integer :: ndim,ncomp
    double precision :: x(ndim),f(ncomp),para(5)
    double precision,external :: laguerre

    f(1)=exp(0-(x(1)+x(2)-2*x(1)*x(2))/2/(1-x(1))/(1-x(2)))*laguerre(int(para(1)),0,x(1)/(1-x(1)))* &
        & laguerre(int(para(2)),0,x(2)/(1-x(2)))/((para(3)*x(1)+para(4)*x(2)-(para(3)+para(4))*x(1)*x(2))*(1-x(1))* &
        & (1-x(2))+(1-x(1))**2*(1-x(2))**2*para(5))

    integrand3p=0

end function integrand3p

subroutine output3p(Nmax,Kt,mass1,mg,b,dim,integration)

    implicit none

    integer :: Nmax,kt,dim
    double precision,dimension(*) :: integration
    double precision :: mass1(3),mg,b
    character(len=4) :: nmaxchar,ktchar
    character(len=128) :: outputfile
    character(len=16) :: mgchar,muchar,mdchar,mschar,bvalue
    integer :: i,j

    write(nmaxchar,'(I4)') nmax
    write(ktchar,'(I4)')   Kt
    write(mgchar,'(F16.4)')   mg
    write(muchar,'(F16.4)')   mass1(1)
    write(mdchar,'(F16.4)')   mass1(2)
    write(mschar,'(F16.4)')   mass1(3)
    write(bvalue,'(F16.4)')   b

    outputfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mu_"//&
        & trim(adjustl(muchar))//"_md_"//trim(adjustl(mdchar))//"_ms_"//trim(adjustl(mschar))//"_mg_"//&
        & trim(adjustl(mgchar))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (36,file=outputfile,status="replace")

    write(36,*) Nmax,Kt

    do i=1,dim
            write(36,*) integration(i)
    enddo

    close(36)

end subroutine output3p

!!!!!!!!!!!!!!############### one gluon exchange interaction for 5 particle Fock sector ##########################!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IEDintegration5p(nmax2,Kt,mass2,massg,b,integration)

    implicit none

    integer :: nmax2,Kt
    double precision :: mass2(5),massg,b
    double precision, dimension(*) :: integration

    double precision :: x1,x2,x1prime,x2prime
    integer :: kp1,kp2,kk1,kk2,i,j,capNprime,LagN,selectmass
    double precision :: c1,c2,delta(10)

    integer :: nvec,flags,method,nregions,neval,info
    integer,external :: integrand
    double precision :: relerr,abserr,para(5),Integral(1),Error(1),prob(1)
    ! double precision :: mineval,maxeval
    integer :: mineval,maxeval

    i=0
    flags=0
    ! mineval=0.0D0
    ! maxeval=0.999D0
    mineval=0
    maxeval=100000000
    method=13
    relerr=1.0D-13
    abserr=1.0D-13
    info=-999

    !print*,"test"
    do capNprime=0,nmax2-3
        do LagN=0,nmax2-3
            do kp1=0,Kt-1
                do kp2=0,Kt-1-kp1
    ! capNprime=0
    ! LagN=2
    ! kp1=0
    ! kp2=0
    ! kk1=0
                    x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
                    x2=(dble(kp2)+0.5)/(dble(Kt)+0.5)
                    do kk1=0,kp1+kp2
                        kk2=kp1+kp2-kk1
    
                        x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
                        x2prime=(dble(kk2)+0.5)/(dble(Kt)+0.5)
                        c1=(sqrt(x1prime*x2)-sqrt(x1*x2prime))**2/(x1+x2)
                        c2=(sqrt(x1prime*x2)+sqrt(x1*x2prime))**2/(x1+x2)
                        delta(1)=(x1-x1prime)**2*(mass2(1)**2/x1/x1prime+mass2(2)**2/x2/x2prime)+2*massg**2
                        delta(2)=(x1-x1prime)**2*(mass2(1)**2/x1/x1prime+mass2(3)**2/x2/x2prime)+2*massg**2
                        delta(3)=(x1-x1prime)**2*(mass2(1)**2/x1/x1prime+mass2(4)**2/x2/x2prime)+2*massg**2
                        delta(4)=(x1-x1prime)**2*(mass2(1)**2/x1/x1prime+mass2(5)**2/x2/x2prime)+2*massg**2
                        delta(5)=(x1-x1prime)**2*(mass2(2)**2/x1/x1prime+mass2(3)**2/x2/x2prime)+2*massg**2
                        delta(6)=(x1-x1prime)**2*(mass2(2)**2/x1/x1prime+mass2(4)**2/x2/x2prime)+2*massg**2
                        delta(7)=(x1-x1prime)**2*(mass2(2)**2/x1/x1prime+mass2(5)**2/x2/x2prime)+2*massg**2
                        delta(8)=(x1-x1prime)**2*(mass2(3)**2/x1/x1prime+mass2(4)**2/x2/x2prime)+2*massg**2
                        delta(9)=(x1-x1prime)**2*(mass2(3)**2/x1/x1prime+mass2(5)**2/x2/x2prime)+2*massg**2
                        delta(10)=(x1-x1prime)**2*(mass2(4)**2/x1/x1prime+mass2(5)**2/x2/x2prime)+2*massg**2

                        do selectmass=1,10
                            i=i+1

                            para(1)=dble(capNprime)
                            para(2)=dble(LagN)
                            para(3)=c1
                            para(4)=c2
                            para(5)=delta(selectmass)/b**2

                            call cuhre(2,1, integrand, para, 1, relerr, abserr, flags, mineval, maxeval,&
                                        & method, "", -1, nregions, neval, info,  Integral, Error, prob)

                            integration(i)=Integral(1)
                            If(info.eq.1) then
                                Print*,info,Integral(1),Error(1),c1,c2,para(5)
                            endif

                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


   return
end subroutine IEDintegration5p

 integer function integrand(ndim,x,ncomp,f,para)
!integer function integrand(ndim,x,ncomp,f)
    implicit none
    integer :: ndim,ncomp
    double precision :: x(ndim),f(ncomp),para(5)
    double precision,external :: laguerre

    
    f(1)=exp(0-(x(1)+x(2)-2*x(1)*x(2))/2/(1-x(1))/(1-x(2)))*laguerre(int(para(1)),0,x(1)/(1-x(1)))* &
        & laguerre(int(para(2)),0,x(2)/(1-x(2)))/((para(3)*x(1)+para(4)*x(2)-(para(3)+para(4))*x(1)*x(2))*(1-x(1))* &
        & (1-x(2))+(1-x(1))**2*(1-x(2))**2*para(5))
    !f(1)=exp(0-(x(1)+x(2))/2.0D0)!*laguerre(int(para(1)),0,x(1))*laguerre(int(para(2)),0,x(2))/(para(3)*x(1)+para(4)*x(2)+para(5))
     !f(1)=1.0D0

    integrand=0

end function integrand

subroutine output5p(Nmax,Kt,mass2,mg,b,dim,integration)

    implicit none

    integer :: Nmax,kt,dim
    double precision,dimension(*) :: integration
    double precision :: mass2(5),mg,b
    character(len=4) :: nmaxchar,ktchar
    character(len=128) :: outputfile
    character(len=16) :: mgchar,muchar,mdchar,mschar,msea1char,msea2char,bvalue
    integer :: i,j

    write(nmaxchar,'(I4)') nmax
    write(ktchar,'(I4)')   Kt
    write(mgchar,'(F16.4)')   mg
    write(muchar,'(F16.4)')   mass2(1)
    write(mdchar,'(F16.4)')   mass2(2)
    write(mschar,'(F16.4)')   mass2(3)
    write(msea1char,'(F16.4)')   mass2(4)
    write(msea2char,'(F16.4)')   mass2(5)
    write(bvalue,'(F16.4)')   b

    outputfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mu_"//&
        & trim(adjustl(muchar))//"_md_"//trim(adjustl(mdchar))//"_ms_"//trim(adjustl(mschar))//"_msea1_"//&
        & trim(adjustl(msea1char))//"_msea2_"//trim(adjustl(msea2char))//"_mg_"//&
        & trim(adjustl(mgchar))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (40,file=outputfile,status="replace")

    write(40,*) Nmax,Kt

    do i=1,dim
            write(40,*) integration(i)
    enddo

    close(40)

end subroutine output5p

!!!!!!!!!!!!!!############### one gluon exchange interaction for 3 particle Fock sector to 5 particle Fock sector #################!!!!!!!!!!!!!

subroutine IEDintegration3to5(nmax2,Kt,mass1,mass2,b,integration3to5)

    implicit none

    integer :: nmax2,Kt
    double precision :: mass1(3),mass2(5),b
    double precision, dimension(*) :: integration3to5

    double precision :: x1,x1prime,x2,x3
    integer :: kp1,kk1,kk2,kk3,i,j,n1,n2,selectmass
    double precision :: c1,c2,delta(3)

    integer :: nvec,flags,method,nregions,neval,info
    integer,external :: integrand3to5
    double precision :: relerr,abserr,para(5),Integral(1),Error(1),prob(1)
    !double precision :: mineval,maxeval
    integer :: mineval,maxeval

    i=0
    flags=0
    mineval=0
    maxeval=100000000
    ! mineval=0.D0
    ! maxeval=0.999999999D0
    method=13
    relerr=1.0D-13
    abserr=1.0D-13
    info=-999

    !print*,"test"
    do n1=0,nmax2-5
        do n2=0,nmax2-5
            do kp1=1,Kt-1

                x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
                do kk1=0,kp1-1
                    do kk2=0,kp1-kk1-1

                        kk3=kp1-kk1-kk2-1
    
                        x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
                        x2=(dble(kk2)+0.5)/(dble(Kt)+0.5)
                        x3=(dble(kk3)+0.5)/(dble(Kt)+0.5)
                        c1=-(x2+x3)
                        delta(1)=mass1(1)**2/x1-mass2(1)**2/x1prime-mass2(4)**2/x2-mass2(5)**2/x3
                        delta(2)=mass1(2)**2/x1-mass2(2)**2/x1prime-mass2(4)**2/x2-mass2(5)**2/x3
                        delta(3)=mass1(3)**2/x1-mass2(3)**2/x1prime-mass2(4)**2/x2-mass2(5)**2/x3
                        do selectmass=1,3

                            i=i+1

                            para(1)=dble(n1)
                            para(2)=dble(n2)
                            para(3)=c1
                            para(4)=delta(selectmass)/b**2

                            call cuhre(2,1, integrand3to5, para, 1, relerr, abserr, flags, mineval, maxeval,&
                                & method, "", -1, nregions, neval, info,  Integral, Error, prob)

                            If(info.eq.1) then
                                print*,Integral(1),error(1),info,para(1),para(2),para(3),para(4)
                            endif

                            integration3to5(i)=Integral(1)

                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


   return
end subroutine IEDintegration3to5

integer function integrand3to5(ndim,x,ncomp,f,para)
    implicit none
    integer :: ndim,ncomp
    double precision :: x(ndim),f(ncomp),para(4)
    double precision,external :: laguerre
    double precision :: PI=acos(-1.D0)

    !f(1)=exp(-(x(1)+x(2))/2)!*laguerre(int(para(1)),0,x(1))* &
    !     & laguerre(int(para(2)),0,x(2))/para(3)/(x(1)+x(2)-para(4))/4/PI
    f(1)=exp(0-(x(1)+x(2)-2*x(1)*x(2))/2/(1-x(1))/(1-x(2)))*laguerre(int(para(1)),0,x(1)/(1-x(1)))* &
         & laguerre(int(para(2)),0,x(2)/(1-x(2)))/para(3)/((x(1)+x(2)-2*x(1)*x(2))/(1-x(1))/(1-x(2))-para(4))&
         & /(1-x(1))**2/(1-x(2))**2/4/PI

    integrand3to5=0

end function integrand3to5

subroutine output3to5(Nmax,Kt,mass1,mass2,b,dim,integration)

    implicit none

    integer :: Nmax,kt,dim
    double precision,dimension(*) :: integration
    double precision :: mass1(3),mass2(5),b
    character(len=4) :: nmaxchar,ktchar
    character(len=128) :: outputfile
    character(len=16) :: mu1char,md1char,ms1char,bvalue
    character(len=16) :: mu2char,md2char,ms2char,msea1char,msea2char
    integer :: i,j

    write(nmaxchar,'(I4)') nmax
    write(ktchar,'(I4)')   Kt
    write(mu1char,'(F16.4)')   mass1(1)
    write(md1char,'(F16.4)')   mass1(2)
    write(ms1char,'(F16.4)')   mass1(3)
    write(mu2char,'(F16.4)')   mass2(1)
    write(md2char,'(F16.4)')   mass2(2)
    write(ms2char,'(F16.4)')   mass2(3)
    write(msea1char,'(F16.4)')   mass2(4)
    write(msea2char,'(F16.4)')   mass2(5)
    write(bvalue,'(F16.4)')   b

    outputfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mass1_"//&
        & trim(adjustl(mu1char))//"_"//trim(adjustl(md1char))//"_"//trim(adjustl(ms1char))//"_mass2_"//&
        & trim(adjustl(mu2char))//"_"//trim(adjustl(md2char))//"_"//trim(adjustl(ms2char))//"_"//&
        & trim(adjustl(msea1char))//"_"//trim(adjustl(msea2char))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (45,file=outputfile,status="replace")

    write(45,*) Nmax,Kt

    do i=1,dim
            write(45,*) integration(i)
    enddo

    close(45)

end subroutine output3to5


!!!!!!!!!!!################# Laguerre polynomials ###################!!!!!!!!!!!!!!!!!!!!!!!!

double precision function laguerre(n,m,p)
implicit none
integer :: n,m,error
double precision :: p,L
double precision ::laguerre0m,laguerre1m
double precision :: laguerre0,laguerre1
double precision,dimension(:),allocatable :: l_part
integer :: k,j

allocate(l_part(n+1),stat=error)
if(error.ne.0) then
   print *, 'cannot allocate array'
end if

laguerre0m=laguerre0(m,p)
laguerre1m=laguerre1(m,p)

if(n.eq.0)then
  Laguerre=laguerre0m
else if(n.eq.1)then
  laguerre=laguerre1m
else
  l_part(1) = laguerre0m
  l_part(2) = laguerre1m
  do k=3,n+1
    l_part(k) = ((2*k-3+m-p)*l_part(k-1)-(k-2+m)*l_part(k-2))/(k-1)
  end do
  laguerre=l_part(n+1)
end if



deallocate(l_part)
return
end function

double precision function laguerre0(m,p)
implicit none
integer :: m
double precision :: p


laguerre0=1

return
end function laguerre0

double precision function laguerre1(m,p)
implicit none
integer :: m
double precision :: p


laguerre1= 1 + m - p

return
end function laguerre1
