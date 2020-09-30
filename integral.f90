!!!!! this part is the function of One-gluon exchange for two particles in one Fock sector !!!!!!!!!!!!!!!!!!!
subroutine searchintegration3p(Nmax,Kt,CapN,LagN,kp1,kp2,kk1,order)
    implicit none

    integer :: Nmax,Kt,CapN,LagN,kp1,kp2,kk1
    integer, dimension(2) :: order

    order(1)=CapN*(Nmax)*Kt*(kt+1)/2+LagN*Kt*(kt+1)/2+(Kt*(kt+1)-((Kt-kp1)*(Kt-kp1+1)))/2+kp2+1
    order(2)=kk1+1

    return
end subroutine searchintegration3p

subroutine IEDintegration3p(nmax2,Kt,mass1,mass2,mass3,mg,b)
    use basis_info

    implicit none

    integer :: nmax2,Kt
    double precision :: mass1,mass2,mass3,mg,b,integral
    character(len=128) :: readfile
    character(len=4) :: nmaxchar,ktchar
    integer :: nmaxtest,kttest,i,j,k,imax,imin,jmax,kmax,h
    character(len=16) :: mgchar,muchar,mdchar,mschar,bvalue
    integer :: initialvalue,listtest

    write(nmaxchar,'(I4)') nmax2
    write(ktchar,'(I4)')   Kt
    write(mgchar,'(F16.4)')   mg
    write(muchar,'(F16.4)')   mass1
    write(mdchar,'(F16.4)')   mass2
    write(mschar,'(F16.4)')   mass3
    write(bvalue,'(F16.4)')   b

    readfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mu_"//&
        & trim(adjustl(muchar))//"_md_"//trim(adjustl(mdchar))//"_ms_"//trim(adjustl(mschar))//"_mg_"//&
        & trim(adjustl(mgchar))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (36,file=readfile)

    read(36,*) nmaxtest,kttest

    if(nmaxtest.eq.nmax2.and.kttest.eq.Kt) then
        !print*,"right"
        continue
    else
        print*,"get the wrong Nmax and Kt",nmax2,Kt
        stop
    endif

    kmax=(nmax2)**2!*Kt**2
    imax=0
    listtest=0

    do h=1,kmax
        imin = 1+imax
        imax = kt*(kt+1)/2+imin-1
        initialvalue=1
        jmax=0
        do i=imin,imax
            if(jmax+1.le.kt) then
                jmax = jmax+1
            else
                initialvalue=initialvalue+1
                jmax = initialvalue
            endif
            do j=1,jmax
                do k=1,3
                    read(36,*) integral
                    integration(i,j,k)=integral
                    listtest=listtest+1
                    
                enddo
            enddo
        enddo
    enddo


    close (36)

    return
end subroutine IEDintegration3p

double precision function IED(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,k)
    use numbers
    use basis_info

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,k
    double precision :: mass1,mass2,massg

    double precision :: x1,x2,x1prime,x2prime,constant,summation
    integer :: lambda,m,capNprime,capN,mu,nu,lagN,order(2)
    double precision :: talmi3
    double precision,external ::TMC
    integer,external :: KroneckerDelta

    x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
    x2=(dble(kp2)+0.5)/(dble(Kt)+0.5)
    x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
    x2prime=(dble(kk2)+0.5)/(dble(Kt)+0.5)
    lambda=mp1+mp2
    summation=0.0D0
    constant=0-KroneckerDelta(kp1+kp2,kk1+kk2)*sqrt(x1*x2*x1prime*x2prime)/2/PI/(x1+x2)
    order(1)=-1
    order(2)=-1

    if(mp1+mp2.eq.mk1+mk2) then
        do m=-(Nmax-3),(Nmax-3)
            do capNprime=0,Nmax
                do capN=0,Nmax

                    mu=np1+np2+(abs(mp1)+abs(mp2)-abs(m)-abs(lambda-m))/2
                    nu=nk1+nk2+(abs(mk1)+abs(mk2)-abs(m)-abs(lambda-m))/2
                    lagN=mu+nu-2*capN+abs(m)-capNprime

                    if(mu-capN.ge.0.and.nu-capN.ge.0.and.lagN.ge.0) then

                        call searchintegration3p(Nmax,Kt,CapNprime,LagN,kp1,kp2,kk1,order)

                        talmi3=TMC(capN,lambda-m,mu-capN,m,np1,mp1,np2,mp2,sqrt(x2/x1))* &
                            & TMC(capN,lambda-m,nu-capN,m,nk1,mk1,nk2,mk2,sqrt(x2prime/x1prime))* &
                            & TMC(capNprime,0,LagN,0,mu-capN,m,nu-capN,-m,1.0D0)

                        summation=summation+talmi3*integration(order(1),order(2),k)
                        ! print*,talmi3,integration(order(1),order(2),k),order(1),order(2)

                    endif
                enddo
            enddo
        enddo

    endif

    IED=constant*summation

    return
end function IED

double precision function IED2momentum(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,b,i,j,k)

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,i,j,k
    double precision :: mass1,mass2,massg,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED
    integer,external :: KroneckerDelta,unitstep
    double precision :: firstterm,secondterm,thirdterm,fourthterm
    integer :: h1


    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kp2)+0.5
    longimomentum(3)=dble(kk1)+0.5
    longimomentum(4)=dble(kk2)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=np2
    trans(2,2)=mp2
    trans(3,1)=nk1
    trans(3,2)=-mk1
    trans(4,1)=nk2
    trans(4,2)=-mk2

    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    am(i)=-1
    am(j)=1
    constant=b**2*sqrt(longimomentum(i)*longimomentum(j)/ktot**2)

    If(trans(i,2).le.0.and.trans(j,2).ge.0) then

        an(i)=-1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)))* &
            & IED(Nmax,Kt,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),np2+an(2)+KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(i,3),mk1-am(3),nk2+an(4)+KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)))*& 
            & IED(Nmax,Kt,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),np2+an(2)+KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(j,3),mk1-am(3),nk2+an(4)+KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(i,1)-1)*unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)))*&
            & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).le.0.and.trans(j,2).lt.0) then

        an(i)=-1;
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+1))* &
            & IED(Nmax,Kt,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),np2+an(2)+KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(i,3),mk1-am(3),nk2+an(4)+KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)))*& 
            & IED(Nmax,Kt,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),np2+an(2)-KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(j,3),mk1-am(3),nk2+an(4)-KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)+1))*&
            & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).gt.0.and.trans(j,2).ge.0) then

        an(i)=1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)))* &
            & IED(Nmax,Kt,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),np2+an(2)-KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(i,3),mk1-am(3),nk2+an(4)-KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)+1))*& 
            & IED(Nmax,Kt,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),np2+an(2)+KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(j,3),mk1-am(3),nk2+an(4)+KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+1)*dble(trans(j,1)))*&
            & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).gt.0.and.trans(j,2).lt.0) then

        an(i)=1
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+1))* &
            & IED(Nmax,Kt,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),np2+an(2)-KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(i,3),mk1-am(3),nk2+an(4)-KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)+1))*& 
            & IED(Nmax,Kt,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),np2+an(2)-KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(j,3),mk1-am(3),nk2+an(4)-KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=sqrt(dble(trans(i,1)+1)*dble(trans(j,1)+1))*&
            & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    endif

    IED2momentum=constant*(firstterm+secondterm+thirdterm+fourthterm)

end function IED2momentum

double precision function IED1momentum(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,b,i,j,k)

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,i,j,k
    double precision :: mass1,mass2,massg,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED
    integer,external :: KroneckerDelta,unitstep
    integer :: h1
    double precision :: firstterm,secondterm

        !print*,i,j

    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kp2)+0.5
    longimomentum(3)=dble(kk1)+0.5
    longimomentum(4)=dble(kk2)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=np2
    trans(2,2)=mp2
    trans(3,1)=nk1
    trans(3,2)=-mk1 
    trans(4,1)=nk2
    trans(4,2)=-mk2
    constant=b*sqrt(longimomentum(i)/ktot)
    secondterm=0.0D0
    firstterm=0.0D0

    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    if(j.eq.1.and.trans(i,2).le.0) then

        am(i)=-1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.1.and.trans(i,2).gt.0) then

        am(i)=-1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.2.and.trans(i,2).ge.0) then

        am(i)=1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.2.and.trans(i,2).lt.0) then

        am(i)=1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    end if

    IED1momentum=constant*(firstterm+secondterm)

    return
end function IED1momentum

!!!!!!!!!!!!!!!!!!!!!!!!!########################read integration list for five particle Fock sector #######################!!!!!!!!!!!!!!!!!

subroutine searchintegration5p(Nmax,Kt,CapN,LagN,kp1,kp2,kk1,order)
    implicit none

    integer :: Nmax,Kt,CapN,LagN,kp1,kp2,kk1
    integer, dimension(2) :: order

    order(1)=CapN*(Nmax-2)*Kt*(kt+1)/2+LagN*Kt*(kt+1)/2+(Kt*(kt+1)-((Kt-kp1)*(Kt-kp1+1)))/2+kp2+1
    order(2)=kk1+1

    return
end subroutine searchintegration5p

subroutine IEDintegration5p(nmax2,Kt,mass1,mass2,mass3,mass4,mass5,mg,b)
    use basis_info

    implicit none

    integer :: nmax2,Kt
    double precision :: mass1,mass2,mass3,mass4,mass5,mg,b,integral
    character(len=128) :: readfile
    character(len=4) :: nmaxchar,ktchar
    integer :: nmaxtest,kttest,i,j,k,imax,imin,jmax,kmax,h
    character(len=16) :: mgchar,muchar,mdchar,mschar,msea1char,msea2char,bvalue
    integer :: initialvalue,listtest

    write(nmaxchar,'(I4)') nmax2
    write(ktchar,'(I4)')   Kt
    write(mgchar,'(F16.4)')   mg
    write(muchar,'(F16.4)')   mass1
    write(mdchar,'(F16.4)')   mass2
    write(mschar,'(F16.4)')   mass3
    write(msea1char,'(F16.4)')   mass4
    write(msea2char,'(F16.4)')   mass5
    write(bvalue,'(F16.4)')   b

    readfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mu_"//&
        & trim(adjustl(muchar))//"_md_"//trim(adjustl(mdchar))//"_ms_"//trim(adjustl(mschar))//"_msea1_"//&
        & trim(adjustl(msea1char))//"_msea2_"//trim(adjustl(msea2char))//"_mg_"//&
        & trim(adjustl(mgchar))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (36,file=readfile)

    read(36,*) nmaxtest,kttest

    if(nmaxtest.eq.nmax2.and.kttest.eq.Kt) then
        !print*,"right"
        continue
    else
        print*,"get the wrong Nmax and Kt",nmax2,Kt
        stop
    endif

    kmax=(nmax2-2)**2!*Kt**2
    imax=0
    listtest=0

    do h=1,kmax
        imin = 1+imax
        imax = kt*(kt+1)/2+imin-1
        initialvalue=1
        jmax=0
        do i=imin,imax
            if(jmax+1.le.kt) then
                jmax = jmax+1
            else
                initialvalue=initialvalue+1
                jmax = initialvalue
            endif
            do j=1,jmax
                do k=1,10
                    read(36,*) integral
                    integration(i,j,k)=integral
                    listtest=listtest+1
                    
                enddo
            enddo
        enddo
    enddo


    close (36)

    return
end subroutine IEDintegration5p
    

double precision function IED5p(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,k)
    use numbers
    use basis_info

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,k
    double precision :: mass1,mass2,massg

    double precision :: x1,x2,x1prime,x2prime,constant,summation
    integer :: lambda,m,capNprime,capN,mu,nu,lagN,order(2)
    double precision :: talmi3
    double precision,external ::TMC
    integer,external :: KroneckerDelta

    x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
    x2=(dble(kp2)+0.5)/(dble(Kt)+0.5)
    x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
    x2prime=(dble(kk2)+0.5)/(dble(Kt)+0.5)
    lambda=mp1+mp2
    summation=0.0D0
    constant=0-KroneckerDelta(kp1+kp2,kk1+kk2)*sqrt(x1*x2*x1prime*x2prime)/2/PI/(x1+x2)
    order(1)=-1
    order(2)=-1

    if(mp1+mp2.eq.mk1+mk2) then
        do m=-(Nmax-5),(Nmax-5)
            do capNprime=0,Nmax
                do capN=0,Nmax

                    mu=np1+np2+(abs(mp1)+abs(mp2)-abs(m)-abs(lambda-m))/2
                    nu=nk1+nk2+(abs(mk1)+abs(mk2)-abs(m)-abs(lambda-m))/2
                    lagN=mu+nu-2*capN+abs(m)-capNprime

                    if(mu-capN.ge.0.and.nu-capN.ge.0.and.lagN.ge.0) then

                        call searchintegration5p(Nmax,Kt,CapNprime,LagN,kp1,kp2,kk1,order)

                        talmi3=TMC(capN,lambda-m,mu-capN,m,np1,mp1,np2,mp2,sqrt(x2/x1))* &
                            & TMC(capN,lambda-m,nu-capN,m,nk1,mk1,nk2,mk2,sqrt(x2prime/x1prime))* &
                            & TMC(capNprime,0,LagN,0,mu-capN,m,nu-capN,-m,1.0D0)

                        summation=summation+talmi3*integration(order(1),order(2),k)

                        !Print*,CapNprime,LagN,order(1),order(2),integration(order(1),order(2),k)

                    endif
                enddo
            enddo
        enddo

    endif

    IED5p=constant*summation

    return
end function IED5p

double precision function IED2momentum5p(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,b,i,j,k)

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,i,j,k
    double precision :: mass1,mass2,massg,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED5p
    integer,external :: KroneckerDelta,unitstep
    double precision :: firstterm,secondterm,thirdterm,fourthterm
    integer :: h1


    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kp2)+0.5
    longimomentum(3)=dble(kk1)+0.5
    longimomentum(4)=dble(kk2)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=np2
    trans(2,2)=mp2
    trans(3,1)=nk1
    trans(3,2)=-mk1
    trans(4,1)=nk2
    trans(4,2)=-mk2

    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    am(i)=-1
    am(j)=1
    constant=b**2*sqrt(longimomentum(i)*longimomentum(j)/ktot**2)

    If(trans(i,2).le.0.and.trans(j,2).ge.0) then

        an(i)=-1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)))* &
            & IED5p(Nmax,Kt,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),np2+an(2)+KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(i,3),mk1-am(3),nk2+an(4)+KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)))*& 
            & IED5p(Nmax,Kt,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),np2+an(2)+KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(j,3),mk1-am(3),nk2+an(4)+KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(i,1)-1)*unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)))*&
            & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).le.0.and.trans(j,2).lt.0) then

        an(i)=-1;
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+1))* &
            & IED5p(Nmax,Kt,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),np2+an(2)+KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(i,3),mk1-am(3),nk2+an(4)+KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)))*& 
            & IED5p(Nmax,Kt,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),np2+an(2)-KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(j,3),mk1-am(3),nk2+an(4)-KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)+1))*&
            & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).gt.0.and.trans(j,2).ge.0) then

        an(i)=1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)))* &
            & IED5p(Nmax,Kt,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),np2+an(2)-KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(i,3),mk1-am(3),nk2+an(4)-KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)+1))*& 
            & IED5p(Nmax,Kt,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),np2+an(2)+KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)+KroneckerDelta(j,3),mk1-am(3),nk2+an(4)+KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+1)*dble(trans(j,1)))*&
            & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    else if(trans(i,2).gt.0.and.trans(j,2).lt.0) then

        an(i)=1
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+1))* &
            & IED5p(Nmax,Kt,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),np2+an(2)-KroneckerDelta(i,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(i,3),mk1-am(3),nk2+an(4)-KroneckerDelta(i,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)+1))*& 
            & IED5p(Nmax,Kt,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),np2+an(2)-KroneckerDelta(j,2),mp2+am(2),&
            & nk1+an(3)-KroneckerDelta(j,3),mk1-am(3),nk2+an(4)-KroneckerDelta(j,4),mk2-am(4),kp1,kp2,kk1,kk2,mass1,&
            & mass2,massg,k)

        fourthterm=sqrt(dble(trans(i,1)+1)*dble(trans(j,1)+1))*&
            & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2,kk1,kk2,&
            & mass1,mass2,massg,k)

    endif

    IED2momentum5p=constant*(firstterm+secondterm+thirdterm+fourthterm)

end function IED2momentum5p

double precision function IED1momentum5p(Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mass1,mass2,massg,b,i,j,k)

    implicit none

    integer :: Nmax,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,i,j,k
    double precision :: mass1,mass2,massg,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED5p
    integer,external :: KroneckerDelta,unitstep
    integer :: h1
    double precision :: firstterm,secondterm

        !print*,i,j

    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kp2)+0.5
    longimomentum(3)=dble(kk1)+0.5
    longimomentum(4)=dble(kk2)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=np2
    trans(2,2)=mp2
    trans(3,1)=nk1
    trans(3,2)=-mk1 
    trans(4,1)=nk2
    trans(4,2)=-mk2
    constant=b*sqrt(longimomentum(i)/ktot)
    secondterm=0.0D0
    firstterm=0.0D0

    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    if(j.eq.1.and.trans(i,2).le.0) then

        am(i)=-1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.1.and.trans(i,2).gt.0) then

        am(i)=-1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.2.and.trans(i,2).ge.0) then

        am(i)=1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    else if(j.eq.2.and.trans(i,2).lt.0) then

        am(i)=1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED5p(Nmax,Kt,np1,mp1+am(1),np2,mp2+am(2),nk1,mk1-am(3),nk2,mk2-am(4),kp1,kp2,kk1,kk2,mass1,mass2,massg,k)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED5p(Nmax,Kt,np1+an(1),mp1+am(1),np2+an(2),mp2+am(2),nk1+an(3),mk1-am(3),nk2+an(4),mk2-am(4),kp1,kp2, &
          & kk1,kk2,mass1,mass2,massg,k)

    end if

    IED1momentum5p=constant*(firstterm+secondterm)

    return
end function IED1momentum5p

double precision function OgeOneFockSector5p(Nmax2,Kt,b,mass1,mass2,mg,k,kp1,kp2,kk1,kk2,sp1,sp2,sk1,sk2,np1,mp1,np2,mp2,nk1,mk1,&
    & nk2,mk2)
    implicit none

    integer :: Nmax2,Kt,k,kp1,kp2,kk1,kk2,sp1,sp2,sk1,sk2,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2
    double precision :: b,mass1,mass2,mg,mu,md
    double precision :: V12,kp1half,kk1half,kp2half,kk2half,Pplus
    double precision :: IED2momentum5p,IED1momentum5p,IED5p

    mu=mass1
    md=mass2
    kp1half=dble(kp1)+0.5D0
    kp2half=dble(kp2)+0.5D0
    kk1half=dble(kk1)+0.5D0
    kk2half=dble(kk2)+0.5D0
    Pplus=Kt+0.5
    V12=0.D0

    if(sp1.eq.sk1.and.sp2.eq.sk2) then

        V12=(mu**2*Pplus**2/kp1half/kk1half+md**2*Pplus**2/kp2half/kk2half)*IED5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,&
          & kk1,kk2,mu,md,mg,k)

        if(sp1.eq.1.and.sp2.eq.1) then

        V12=V12+IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,k)*Pplus**2/kk1half/kp1half &
            & +IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,k)*Pplus**2/kk2half/kp2half &
            & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,k)*Pplus**2/kk2half/kp1half &
            & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,k)*Pplus**2/kk1half/kp2half 

        else if(sp1.eq.-1.and.sp2.eq.-1) then

        V12=V12+IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,3,k)*Pplus**2/kk1half/kp1half &
          & +IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,4,k)*Pplus**2/kk2half/kp2half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,3,k)*Pplus**2/kk1half/kp2half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,4,k)*Pplus**2/kk2half/kp1half

        else if(sp1.eq.1.and.sp2.eq.-1) then

        V12=V12+IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,k)*Pplus**2/kk1half/kp1half &
          & +IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,4,k)*Pplus**2/kk2half/kp2half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,4,k)*Pplus**2/kk2half/kk1half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,k)*Pplus**2/kp1half/kp2half

        else if(sp1.eq.-1.and.sp2.eq.1) then

        V12=V12+IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,3,k)*Pplus**2/kk1half/kp1half &
          & +IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,k)*Pplus**2/kk2half/kp2half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,3,k)*Pplus**2/kk2half/kk1half &
          & -IED2momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,k)*Pplus**2/kp1half/kp2half 

        endif

    else if(sp1.eq.sk1.and.sp2.ne.sk2) then

        if(sp1.eq.1.and.sp2.eq.1) then

          V12=md*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,k))*Pplus**2/kp2half/kk2half+ &
            & Pplus**2*(kk2half-kp2half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,k) &
            & /kp1half/kp2half/kk2half)

        else if(sp1.eq.-1.and.sp2.eq.1) then

          V12=md*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,k))*Pplus**2/kp2half/kk2half+ &
            & Pplus**2*(kk2half-kp2half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,k) &
            & /kk1half/kp2half/kk2half)

        else if(sp1.eq.1.and.sp2.eq.-1) then

          V12=md*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,k))*Pplus**2/kp2half/kk2half+ &
            & Pplus**2*(kp2half-kk2half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,k) &
            & /kk1half/kp2half/kk2half)

        else if(sp1.eq.-1.and.sp2.eq.-1) then

          V12=md*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,k))*Pplus**2/kp2half/kk2half+ &
            & Pplus**2*(kp2half-kk2half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,k) &
            & /kp1half/kp2half/kk2half)

        endif

    else if(sp1.ne.sk1.and.sp2.eq.sk2) then

        if(sp1.eq.1.and.sp2.eq.1) then

          V12=mu*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,k))*Pplus**2/kp1half/kk1half+ &
            & Pplus**2*(kk1half-kp1half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,k) &
            & /kp1half/kp2half/kk1half)

        else if(sp1.eq.1.and.sp2.eq.-1) then

          V12=mu*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,k))*Pplus**2/kp1half/kk1half+ &
            & Pplus**2*(kk1half-kp1half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,k) &
            & /kp1half/kk1half/kk2half)

        else if(sp1.eq.-1.and.sp2.eq.1) then

          V12=mu*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,k))*Pplus**2/kp1half/kk1half+ &
            & Pplus**2*(kp1half-kk1half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,k) &
            & /kp1half/kk2half/kk1half)

        else if(sp1.eq.-1.and.sp2.eq.-1) then

          V12=mu*((IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,k)- &
            & IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,k))*Pplus**2/kp1half/kk1half+ &
            & Pplus**2*(kp1half-kk1half)*IED1momentum5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,k) &
            & /kp1half/kp2half/kk1half)

        endif

    else if(sp2.eq.sk1.and.sp1.eq.sk2.and.sp1.ne.sp2) then

        V12=mu*md*Pplus**2*(kk1half-kp1half)*(kk2half-kp2half)/kk1half/kk2half/kp1half/kp2half* &
          & IED5p(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,k)

    else if(sp1.eq.sp2.and.sk1.eq.sk2.and.sp1.ne.sk1) then

        V12=0.D0

    endif

   OgeOneFockSector5p=2*V12

    return
end function OgeOneFockSector5p

!!!!!!!!############# the function of one-gluon exchange for five particle system ###########!!!!!!!!!!!!!!!!!!!!

subroutine searchintegration3t5(Nmax,Kt,n1,n2,kp1,kk1,kk2,kk3,order)
    implicit none

    integer :: Nmax,Kt,n1,n2,kp1,kk1,kk2,kk3
    integer, dimension(2) :: order

    order(1)=n1*(Nmax-4)+n2+1
    order(2)=kp1*(kp1-1)*(kp1+1)/6+kk1*kp1-kk1*(kk1-1)/2+kk2+1

    return
end subroutine searchintegration3t5

subroutine IEDintegration3t5(nmax2,Kt,mass1,mass2,b)
    use basis_info

    implicit none

    integer :: nmax2,Kt
    double precision :: mass1(3),mass2(5),b,integral
    character(len=128) :: readfile
    character(len=4) :: nmaxchar,ktchar
    integer :: n1,n2,kp1,kk1,kk2,k,i,j,nmaxtest,kttest
    character(len=16) :: mu1char,md1char,ms1char,bvalue
    character(len=16) :: mu2char,md2char,ms2char,msea1char,msea2char
    integer :: listtest!,initialvalue

    write(nmaxchar,'(I4)') nmax2
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

    readfile="integration_"//"nmax_"//trim(adjustl(nmaxchar))//"_Kmax_"//trim(adjustl(ktchar))//"_mass1_"//&
        & trim(adjustl(mu1char))//"_"//trim(adjustl(md1char))//"_"//trim(adjustl(ms1char))//"_mass2_"//&
        & trim(adjustl(mu2char))//"_"//trim(adjustl(md2char))//"_"//trim(adjustl(ms2char))//"_"//&
        & trim(adjustl(msea1char))//"_"//trim(adjustl(msea2char))//"_b_"//trim(adjustl(bvalue))//".dat"

    open (36,file=readfile)

    read(36,*) nmaxtest,kttest

    if(nmaxtest.eq.nmax2.and.kttest.eq.Kt) then
        print*,"right"
        continue
    else
        print*,"get the wrong Nmax and Kt",nmax2,Kt
        stop
    endif

    i=0

    do n1=0,nmax2-5
        do n2=0,nmax2-5
            i=i+1
            j=0
            do kp1=1,kt-1
                do kk1=0,kp1-1
                    do kk2=0,kp1-kk1-1
                        j=j+1
                        do k=1,3
                            read(36,*) integral
                            integration3t5(i,j,k)=integral
                            listtest=listtest+1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


    close (36)

    return
end subroutine IEDintegration3t5

double precision function IED3t5(Nmax,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)
    use numbers
    use basis_info

    implicit none

    integer :: Nmax,Kt,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3
    double precision :: mass1,mass1prime,mass2,mass3

    double precision :: x1,x1prime,x2,x3,constant,summation
    integer :: n1,n2,order(2),median,median2,k
    double precision :: talmi3
    double precision,external ::TMC
    integer,external :: KroneckerDelta

    x1=(dble(kp1)+0.5)/(dble(Kt)+0.5)
    x1prime=(dble(kk1)+0.5)/(dble(Kt)+0.5)
    x2=(dble(kk2)+0.5)/(dble(Kt)+0.5)
    x3=(dble(kk3)+0.5)/(dble(Kt)+0.5)
    summation=0.0D0
    constant=KroneckerDelta(kp1,kk1+kk2+kk3+1)*sqrt(x1prime*x2*x3/x1)
    order(1)=-1
    order(2)=-1

    if(np1.ge.0.and.nk1.ge.0.and.nk2.ge.0.and.nk3.ge.0) then
        If(constant.ne.0) then
            do n1=0,Nmax-5
                do n2=0,Nmax-5

                    median=nk2+nk3-n1+(abs(mk2)+abs(mk3)-abs(mk2+mk3))/2
                    median2=2*median+2*nk1+abs(mk2+mk3)+abs(mk1)-2*np1-abs(mp1)-2*n2

                    if(median.ge.0.and.median2.eq.0.and.mp1.eq.mk1+mk2+mk3) then

                        call searchintegration3t5(Nmax,Kt,n1,n2,kp1,kk1,kk2,kk3,order)

                        talmi3=TMC(median,-mk2-mk3,n1,0,nk2,-mk2,nk3,-mk3,sqrt(x3/x2))* &
                            & TMC(np1,-mp1,n2,0,median,-mk2-mk3,nk1,-mk1,sqrt(x1prime/(x2+x3)))

                        summation=summation+talmi3*integration3t5(order(1),order(2),k)

                        ! Print*,integration3t5(order(1),order(2),k),order(1),order(2),n1,n2,kp1,kk1,kk2,kk3
                        ! Print*,integration3t5(1,8,k)

                    endif
                enddo
            enddo
        endif
    endif

    IED3t5=constant*summation

    return
end function IED3t5

double precision function IED2momentum3t5(Nmax,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,&
    & mass3,b,i,j)

    implicit none

    integer :: Nmax,Kt,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,i,j,k
    double precision :: mass1,mass1prime,mass2,mass3,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED3t5
    integer,external :: KroneckerDelta,unitstep
    double precision :: firstterm,secondterm,thirdterm,fourthterm
    integer :: h1


    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kk1)+0.5
    longimomentum(3)=dble(kk2)+0.5
    longimomentum(4)=dble(kk3)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=nk1
    trans(2,2)=-mk1
    trans(3,1)=nk2
    trans(3,2)=-mk2
    trans(4,1)=nk3
    trans(4,2)=-mk3

    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    am(i)=-1
    am(j)=1
    constant=b**2*sqrt(longimomentum(i)*longimomentum(j)/ktot**2)

    If(trans(i,2).le.0.and.trans(j,2).ge.0) then

        an(i)=-1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)))* &
            & IED3t5(Nmax,Kt,k,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),nk1+an(2)+KroneckerDelta(i,2),mk1-am(2),&
            & nk2+an(3)+KroneckerDelta(i,3),mk2-am(3),nk3+an(4)+KroneckerDelta(i,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)))*& 
            & IED3t5(Nmax,Kt,k,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),nk1+an(2)+KroneckerDelta(j,2),mk1-am(2),&
            & nk2+an(3)+KroneckerDelta(j,3),mk2-am(3),nk3+an(4)+KroneckerDelta(j,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        fourthterm=unitstep(trans(i,1)-1)*unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)))*&
            & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1,kk2,kk3,&
            & mass1,mass1prime,mass2,mass3)

    else if(trans(i,2).le.0.and.trans(j,2).lt.0) then

        an(i)=-1;
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2))+1)*dble(trans(j,1)+1))* &
            & IED3t5(Nmax,Kt,k,np1+an(1)+KroneckerDelta(i,1),mp1+am(1),nk1+an(2)+KroneckerDelta(i,2),mk1-am(2),&
            & nk2+an(3)+KroneckerDelta(i,3),mk2-am(3),nk3+an(4)+KroneckerDelta(i,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        thirdterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)))*& 
            & IED3t5(Nmax,Kt,k,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),nk1+an(2)-KroneckerDelta(j,2),mk1-am(2),&
            & nk2+an(3)-KroneckerDelta(j,3),mk2-am(3),nk3+an(4)-KroneckerDelta(j,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        fourthterm=unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1))*dble(trans(j,1)+1))*&
            & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1,kk2,kk3,&
            & mass1,mass1prime,mass2,mass3)

    else if(trans(i,2).gt.0.and.trans(j,2).ge.0) then

        an(i)=1
        an(j)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))+1))* &
            & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)))* &
            & IED3t5(Nmax,Kt,k,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),nk1+an(2)-KroneckerDelta(i,2),mk1-am(2),&
            & nk2+an(3)-KroneckerDelta(i,3),mk2-am(3),nk3+an(4)-KroneckerDelta(i,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2))+1)*dble(trans(i,1)+1))*& 
            & IED3t5(Nmax,Kt,k,np1+an(1)+KroneckerDelta(j,1),mp1+am(1),nk1+an(2)+KroneckerDelta(j,2),mk1-am(2),&
            & nk2+an(3)+KroneckerDelta(j,3),mk2-am(3),nk3+an(4)+KroneckerDelta(j,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        fourthterm=unitstep(trans(j,1)-1)*sqrt(dble(trans(i,1)+1)*dble(trans(j,1)))*&
            & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1,kk2,kk3,&
            & mass1,mass1prime,mass2,mass3)

    else if(trans(i,2).gt.0.and.trans(j,2).lt.0) then

        an(i)=1
        an(j)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+abs(trans(j,2))))* &
            & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-sqrt(dble(trans(i,1)+abs(trans(i,2)))*dble(trans(j,1)+1))* &
            & IED3t5(Nmax,Kt,k,np1+an(1)-KroneckerDelta(i,1),mp1+am(1),nk1+an(2)-KroneckerDelta(i,2),mk1-am(2),&
            & nk2+an(3)-KroneckerDelta(i,3),mk2-am(3),nk3+an(4)-KroneckerDelta(i,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        thirdterm=-sqrt(dble(trans(j,1)+abs(trans(j,2)))*dble(trans(i,1)+1))*& 
            & IED3t5(Nmax,Kt,k,np1+an(1)-KroneckerDelta(j,1),mp1+am(1),nk1+an(2)-KroneckerDelta(j,2),mk1-am(2),&
            & nk2+an(3)-KroneckerDelta(j,3),mk2-am(3),nk3+an(4)-KroneckerDelta(j,4),mk3-am(4),kp1,kk1,kk2,kk3,mass1,&
            & mass1prime,mass2,mass3)

        fourthterm=sqrt(dble(trans(i,1)+1)*dble(trans(j,1)+1))*&
            & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1,kk2,kk3,&
            & mass1,mass1prime,mass2,mass3)

    endif

    IED2momentum3t5=constant*(firstterm+secondterm+thirdterm+fourthterm)

end function IED2momentum3t5

double precision function IED1momentum3t5(Nmax,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,&
    & mass2,mass3,b,i,j)

    implicit none

    integer :: Nmax,Kt,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,i,j,k
    double precision :: mass1,mass1prime,mass2,mass3,b

    double precision :: ktot,longimomentum(4),constant
    integer :: trans(4,2),am(4),an(4)
    double precision,external :: IED3t5
    integer,external :: KroneckerDelta,unitstep
    integer :: h1
    double precision :: firstterm,secondterm

        !print*,i,j

    ktot=dble(kt)+0.5
    longimomentum(1)=dble(kp1)+0.5
    longimomentum(2)=dble(kk1)+0.5
    longimomentum(3)=dble(kk2)+0.5
    longimomentum(4)=dble(kk3)+0.5
    trans(1,1)=np1
    trans(1,2)=mp1
    trans(2,1)=nk1
    trans(2,2)=-mk1
    trans(3,1)=nk2
    trans(3,2)=-mk2 
    trans(4,1)=nk3
    trans(4,2)=-mk3
    constant=b*sqrt(longimomentum(i)/ktot)
    secondterm=0.0D0
    firstterm=0.0D0
    do h1=1,4
        am(h1)=0
        an(h1)=0
    enddo

    if(j.eq.1.and.trans(i,2).le.0) then

        am(i)=-1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1, &
          & kk2,kk3,mass1,mass1prime,mass2,mass3)

    else if(j.eq.1.and.trans(i,2).gt.0) then

        am(i)=-1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1, &
          & kk2,kk3,mass1,mass1prime,mass2,mass3)

    else if(j.eq.2.and.trans(i,2).ge.0) then

        am(i)=1
        an(i)=-1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))+1))* &
          & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-unitstep(trans(i,1)-1)*sqrt(dble(trans(i,1)))* &
          & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1, &
          & kk2,kk3,mass1,mass1prime,mass2,mass3)

    else if(j.eq.2.and.trans(i,2).lt.0) then

        am(i)=1
        an(i)=1

        firstterm=sqrt(dble(trans(i,1)+abs(trans(i,2))))* &
          & IED3t5(Nmax,Kt,k,np1,mp1+am(1),nk1,mk1-am(2),nk2,mk2-am(3),nk3,mk3-am(4),kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        secondterm=-sqrt(dble(trans(i,1)+1))* &
          & IED3t5(Nmax,Kt,k,np1+an(1),mp1+am(1),nk1+an(2),mk1-am(2),nk2+an(3),mk2-am(3),nk3+an(4),mk3-am(4),kp1,kk1, &
          & kk2,kk3,mass1,mass1prime,mass2,mass3)

    end if

    IED1momentum3t5=constant*(firstterm+secondterm)

    return
end function IED1momentum3t5

subroutine OgeThreeToFive(Nmax2,Kt,b,k,mass1,mass1prime,mass2,mass3,kp1,sp1,np1,mp1,kk1,kk2,kk3,sk1,sk2,sk3,nk1,&
    &  mk1,nk2,mk2,nk3,mk3,Vertex)
    implicit none

    integer :: Nmax2,Kt,kp1,sp1,np1,mp1,kk1,kk2,kk3,sk1,sk2,sk3,nk1,mk1,nk2,mk2,nk3,mk3,k
    double precision :: b,mass1,mass1prime,mass2,mass3
    double precision :: V12,kp1half,kk1half,kk2half,kk3half,Pplus,Vertex
    double precision,external :: IED2momentum3t5,IED1momentum3t5,IED3t5

    kp1half=dble(kp1)+0.5D0
    kk1half=dble(kk1)+0.5D0
    kk2half=dble(kk2)+0.5D0
    kk3half=dble(kk3)+0.5D0
    Pplus=Kt+0.5
    V12=0.D0

    if(sp1.eq.sk1.and.sk2.ne.sk3) then

        V12=(mass1*mass1prime*Pplus**2/kp1half/kk1half-mass2*mass3*Pplus**2/kk2half/kk3half)*IED3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2 &
            & ,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3)

        if(sk2.eq.-1.and.sp1.eq.1) then

          V12=V12+IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
            & *Pplus**2/kk1half/kp1half &
            & +IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,3)&
            & *Pplus**2/kk2half/kk3half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,3)&
            & *Pplus**2/kk2half/kk1half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,1)&
            & *Pplus**2/kp1half/kk3half 

        else if(sk2.eq.1.and.sp1.eq.1) then

          V12=V12+IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
            & *Pplus**2/kk1half/kp1half &
            & +IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,4)&
            & *Pplus**2/kk2half/kk3half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,4)&
            & *Pplus**2/kk1half/kk3half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,1)&
            & *Pplus**2/kk2half/kp1half

        else if(sk2.eq.1.and.sp1.eq.-1) then

          V12=V12+IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,4)&
            & *Pplus**2/kk2half/kk3half &
            & +IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
            & *Pplus**2/kk1half/kp1half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,4)&
            & *Pplus**2/kk3half/kp1half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,2)&
            & *Pplus**2/kk1half/kk2half

        else if(sk2.eq.-1.and.sp1.eq.-1) then

          V12=V12+IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,3)&
            & *Pplus**2/kk2half/kk3half &
            & +IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
            & *Pplus**2/kk1half/kp1half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,3)&
            & *Pplus**2/kk2half/kp1half &
            & -IED2momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,2)&
            & *Pplus**2/kk1half/kk3half 

        endif

    else if(sp1.eq.sk1.and.sk2.eq.sk3) then

        if(sp1.eq.1.and.sk2.eq.-1) then

          V12=mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,2)&
            & *Pplus**2/kk2half/kk3half- &
            & mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
            & *Pplus**2/kp1half/kk3half+ &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,2)&
            & *Pplus**2/kk3half/kk2half+ &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
            & *Pplus**2/kp1half/kk2half

        else if(sp1.eq.1.and.sk2.eq.1) then

          V12=mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
            & *Pplus**2/kk1half/kk3half- &
            & mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,1)&
            & *Pplus**2/kk2half/kk3half+ &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
            & *Pplus**2/kk1half/kk2half- &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,1)&
            & *Pplus**2/kk3half/kk2half

        else if(sp1.eq.-1.and.sk2.eq.1) then

          V12=mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,1)&
            & *Pplus**2/kp1half/kk3half- &
            & mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,1)&
            & *Pplus**2/kk2half/kk3half+ &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,1)&
            & *Pplus**2/kp1half/kk2half- &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,1)&
            & *Pplus**2/kk3half/kk2half

        else if(sp1.eq.-1.and.sk2.eq.-1) then

          V12=mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,2)&
            & *Pplus**2/kk2half/kk3half- &
            & mass3*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,2)&
            & *Pplus**2/kk1half/kk3half+ &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,2)&
            & *Pplus**2/kk3half/kk2half- &
            & mass2*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,2)&
            & *Pplus**2/kk1half/kk2half

        endif

    else if(sp1.ne.sk1.and.sk2.ne.sk3) then

        if(sp1.eq.1.and.sk2.eq.-1) then

        V12=mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,2)&
        & *Pplus**2/kp1half/kk2half- &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,2)&
        & *Pplus**2/kk1half/kk2half+ &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
        & *Pplus**2/kp1half/kk1half- &
        & mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,2)&
        & *Pplus**2/kp1half/kk1half

        else if(sp1.eq.1.and.sk2.eq.1) then

        V12=mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,2)&
        & *Pplus**2/kp1half/kk3half- &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,2)&
        & *Pplus**2/kk1half/kk3half+ &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,2)&
        & *Pplus**2/kp1half/kk1half- &
        & mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,2)&
        & *Pplus**2/kp1half/kk1half

        else if(sp1.eq.-1.and.sk2.eq.1) then


        V12=mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
        & *Pplus**2/kp1half/kk1half- &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,1)&
        & *Pplus**2/kp1half/kk1half+ &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,1)&
        & *Pplus**2/kk1half/kk2half- &
        & mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,3,1)&
        & *Pplus**2/kp1half/kk2half

        else if(sp1.eq.-1.and.sk2.eq.-1) then

        V12=mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,2,1)&
        & *Pplus**2/kp1half/kk1half- &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,1,1)&
        & *Pplus**2/kp1half/kk1half+ &
        & mass1prime*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,1)&
        & *Pplus**2/kk1half/kk3half- &
        & mass1*IED1momentum3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,mass1,mass1prime,mass2,mass3,b,4,1)&
        & *Pplus**2/kp1half/kk3half

        endif

    else if(sp1.ne.sk1.and.sp1.eq.sk2.and.sp1.eq.sk3) then

        V12=(mass1prime*mass3*Pplus**2/kk1half/kk3half+mass1prime*mass2*Pplus**2/kk1half/kk2half-mass1*mass3*Pplus**2/kp1half &
          & /kk3half-mass1*mass2*Pplus**2/kp1half/kk2half)*IED3t5(nmax2,Kt,k,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3,kp1,kk1,kk2,kk3,&
          & mass1,mass1prime,mass2,mass3)

    else if(sp1.ne.sk1.and.sk1.eq.sk2.and.sk1.eq.sk3) then

        V12=0.D0

    endif

    Vertex=2*V12

    !Print*,OgeThreeToFive

    return
end subroutine OgeThreeToFive

!!!!!!!################ Laguerre polynomials #########!!!!!!!!!!
double precision function laguerre(n,m,p)
implicit none
integer :: n,m,error
double precision :: p
double precision ::laguerre0m,laguerre1m
double precision :: laguerre0,laguerre1
double precision,dimension(:),allocatable :: l_part
integer :: k

allocate(l_part(n+1),stat=error)
if(error.ne.0) then
   print *, 'cannot allocate array'
end if

laguerre0m=laguerre0()
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

double precision function laguerre0()
implicit none

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