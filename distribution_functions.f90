subroutine distribution_functions_GPD(Nmax1,Nmax2,Nmax3,Mj,Kt,b,nestate,selectquark,deltamax,lag)
    use numbers
    use basis_info

    implicit none

    integer :: Nmax,Mj,Kt,selectquark,nestate
    double precision :: deltamax,b,lag
    double precision :: mu1,md1,mg1,ms1,mu2,md2,mg2,ms2,mu2,md2,mg2,ms2

    double precision, dimension(:),allocatable :: ev1,ev2
    double precision, dimension(:,:),allocatable :: distributionfunctions,evector1,evector2
    double precision, dimension(:,:),allocatable :: distributionfunctions1,distributionfunctions2
    double precision, dimension(:,:),allocatable :: distributionfunctionskt1,distributionfunctionskt2
    integer :: functionsdimension,snm,loopnumber1,loopnumber2
    integer, dimension(9) :: sums
    double precision :: xcutoff,nmcutoff,massoffset,dx,dy,step
    integer :: error,i,j,dimtot1,dimtot2,delta
    double precision :: normalization


    massoffset=0.0D-2
    nmcutoff=0.0D0
    xcutoff=0.0D0
    step=0.1

    sums(1)=0
    call snmlen2(Nmax,Mj,snm)
    call dddnmlen2(Nmax,Mj,sums(2))
    call uddnmlen2(Nmax,Mj,sums(3))
    call dudnmlen2(Nmax,Mj,sums(4))
    call uudnmlen2(Nmax,Mj,sums(5))
    call ddunmlen2(Nmax,Mj,sums(6))
    call udunmlen2(Nmax,Mj,sums(7))
    call duunmlen2(Nmax,Mj,sums(8))
    call uuunmlen2(Nmax,Mj,sums(9))


    if (selectquark.eq.1) then

        functionsdimension=Kt*(floor(deltamax/step)+1)

        Print*,"the dimensions of matrix of distribution functions = ",functionsdimension

        allocate(distributionfunctions(4,functionsdimension),distributionfunctions1(4,kt),&
            & distributionfunctions2(4,kt),stat=error)
            if(error.ne.0) then
                print *, 'cannot allocate basis array'
            end if

        call dimtotal(nmax,nmax,Mj,Kt,dimtot1)
        allocate(ev1(nestate),evector1(nestate,dimtot1),stat=error)
            if(error.ne.0) then
                print *, 'cannot allocate basis array'
            end if
        call wfproduction(nmax1,nmax2,nmax3,Mj,Kt,massoffset,nmcutoff,xcutoff,lag,b,nestate,dimtot1,ev1,evector1)
!	Print*,"test1"

        loopnumber1=0

        do delta=0,floor(deltamax/step)

            dy=0.0D0
            dx=delta*step

            call GPD_for_d_quark(nmax,Mj,Kt,b,sums,snm,nestate,dimtot1,evector1,dx,dy,distributionfunctions1)

            do i=1,kt

                loopnumber1=loopnumber1+1

                distributionfunctions(1,loopnumber1)=distributionfunctions1(1,i)
                distributionfunctions(2,loopnumber1)=distributionfunctions1(2,i)
                distributionfunctions(3,loopnumber1)=distributionfunctions1(3,i)
                distributionfunctions(4,loopnumber1)=distributionfunctions1(4,i)

            enddo

        enddo

        deallocate(s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat,distributionfunctions1)

        call output(nmax,Mj,Kt,selectquark,loopnumber1,distributionfunctions)

        !!!!!!!!test the normalization for d quark. and the normalization should be equal to 1 !!!!!!!!

        !do j=1,kt-3

        !    normalization=normalization+(distributionfunctions(1,j+1)-distributionfunctions(1,j))*distributionfunctions(3,j)

        !enddo

        !Print*,normalization
    
        deallocate(distributionfunctions,ev1,evector1)

    else if (selectquark.eq.2) then

        functionsdimension=(Kt)*(floor(deltamax/step)+1)

        Print*,"the dimensions of matrix of distribution functions = ",functionsdimension

         allocate(distributionfunctions(3,functionsdimension),distributionfunctions1(3,kt),stat=error)
            if(error.ne.0) then
                print *, 'cannot allocate basis array'
            end if

        call dimtotal(nmax,nmax,Mj,Kt,dimtot1)

        allocate(ev1(nestate),evector1(nestate,dimtot1),stat=error)
            if(error.ne.0) then
                print *, 'cannot allocate basis array'
            end if
        call wfproduction(nmax,nmax,Mj,Kt,massoffset,nmcutoff,xcutoff,lag,b,mu,md,mg,nestate,dimtot1,ev1,evector1)

       


        loopnumber1=0
        normalization=0.D0

        do delta=0,floor(deltamax/step)

            dy=0.0D0
            dx=delta*step

            call GPD_for_u_quark(nmax,Mj,Kt,b,sums,snm,nestate,dimtot1,evector1,dx,dy,distributionfunctions1)

            do i=1,kt

                loopnumber1=loopnumber1+1

                distributionfunctions(1,loopnumber1)=distributionfunctions1(1,i)
                distributionfunctions(2,loopnumber1)=distributionfunctions1(2,i)
                distributionfunctions(3,loopnumber1)=distributionfunctions1(3,i)

            enddo

        enddo

        call output(nmax,Mj,Kt,selectquark,loopnumber1,distributionfunctions)

        !!!!!!!!!!!test the normalization for u quark. And the normalization should be equal to 2 !!!!!!!!

        !do j=1,kt-1

        !    normalization=normalization+(distributionfunctions(1,j+1)-distributionfunctions(1,j))*distributionfunctions(3,j)

        !enddo

        !Print*,normalization

        deallocate(s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat,distributionfunctions1,&
            & distributionfunctions,ev1,evector1)

    endif

    return
end subroutine distribution_functions_GPD


subroutine GPD_for_d_quark(nmax,Mj,Kt,b,sums,snm,nestate,dimtot,evector,dx,dy,distributionfunctionsvalue)
    use numbers
    use basis_info

    implicit none

    integer :: nmax,Mj,Kt,nestate,dimtot,snm
    integer, dimension(*) :: sums
    double precision,dimension(nestate,dimtot) :: evector
    double precision :: dx,dy,b,deltamax
    double precision,dimension(4,kt) :: distributionfunctionsvalue
    double precision :: distributionfunctionsvalue_E,dxcorrect,mq

    integer, dimension(Kt+1):: kmatrix
    integer :: npu1,mpu1,npd,mpd,npu2,mpu2,pu1,pd,pu2,nku1,mku1,nkd,mkd,nku2,mku2,ku1,kd,ku2
    integer :: accumulate,k_d,k_u,spin,initialstate,finalstate,talmitotnumber
    double precision :: fraction,bp,bk,factor,olaps
    double precision, dimension(:,:), allocatable :: talmip,talmik
    integer :: summaxp,summaxk,sum1,sum2,i,error,test
    integer :: uplimit,finalindex

    talmitotnumber=(2*(nmax-3)+1)**2*(floor(dble(nmax-3)/2.0D0)+1)**2

    allocate(talmip(7,talmitotnumber),talmik(7,talmitotnumber),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

    kmatrix(1)=0
    fraction=0.0D0
    accumulate=0
    factor=0.D0
    dxcorrect=1.0D-10
    mq=Md_DEFAULT

    do i=0,kt-1
        accumulate=accumulate+Kt-i
        kmatrix(i+2)=accumulate*snm
    enddo


    do k_d=0,kt-1

        fraction=(dble(k_d)+0.5)/(dble(Kt)+0.5)

        distributionfunctionsvalue(3,k_d+1)=0.D0
        distributionfunctionsvalue(4,k_d+1)=0.D0

        distributionfunctionsvalue_E=0.D0

        do k_u=0,Kt-k_d-1

            do spin=1,8

                do initialstate=kmatrix(k_d+1)+k_u*snm+sums(spin)+1,kmatrix(k_d+1)+k_u*snm+sums(spin+1)

                    npu1=n1dat(initialstate)
                    mpu1=m1dat(initialstate)
                    npd =n2dat(initialstate)
                    mpd =m2dat(initialstate)
                    npu2=n3dat(initialstate)
                    mpu2=m3dat(initialstate)
                    pu1 =k1dat(initialstate)
                    pd  =k2dat(initialstate)
                    pu2 =k3dat(initialstate)

                    if (k_d.eq.pd.and.k_u.eq.pu1) then

                        do finalstate=kmatrix(k_d+1)+k_u*snm+sums(spin)+1,kmatrix(k_d+1)+k_u*snm+sums(spin+1)

                            nku1=n1dat(finalstate)
                            mku1=m1dat(finalstate)
                            nkd =n2dat(finalstate)
                            mkd =m2dat(finalstate)
                            nku2=n3dat(finalstate)
                            mku2=m3dat(finalstate)
                            ku1 =k1dat(finalstate)
                            kd  =k2dat(finalstate)
                            ku2 =k3dat(finalstate)

                            if (k_d.eq.kd.and.k_u.eq.ku1) then

                                call Talmisquare(nmax,npu1,mpu1,npu2,mpu2,npd,mpd,dble(pu1)+0.5,dble(pu2)+0.5,&
                                    & dble(pd)+0.5,summaxp,talmip)
                                call Talmisquare(nmax,nku1,mku1,nku2,mku2,nkd,mkd,dble(ku1)+0.5,dble(ku2)+0.5,&
                                    & dble(kd)+0.5,summaxk,talmik)

                                do sum1=1,summaxp

                                    do sum2=1,summaxk

                                        if(talmip(1,sum1).eq.talmik(1,sum2).and.talmip(2,sum1).eq.talmik(2,sum2)) then

                                            if(talmip(5,sum1).eq.talmik(5,sum2).and.talmip(6,sum1).eq.talmik(6,sum2)) then


                                                bp=b*sqrt(dble(pu1+pu2+1)*(dble(pd)+0.5))/(dble(Kt)+0.5)
                                                bk=b*sqrt(dble(ku1+ku2+1)*(dble(kd)+0.5))/(dble(Kt)+0.5)
                                                factor=olaps(bp,nint(talmip(3,sum1)),nint(talmip(4,sum1)),bk,&
                                                    & nint(talmik(3,sum2)),nint(talmik(4,sum2)),(1.0D0-fraction)*dx,&
                                                    & (1.0D0-fraction)*dy)

                                                distributionfunctionsvalue(1,k_d+1)=fraction
                                                distributionfunctionsvalue(2,k_d+1)=dx**2+dy**2
                                                distributionfunctionsvalue(3,k_d+1)=distributionfunctionsvalue(3,k_d+1)&
                                                    & +talmip(7,sum1)*evector(nestate,initialstate)*talmik(7,sum2)*&
                                                    & evector(nestate,finalstate)*factor*kt!*fraction
                                                !distributionfunctionsvalue(3,k_d+1)=&
                                                !    & talmip(7,sum1)*evector(nestate,initialstate)*talmik(7,sum2)*&
                                                !    & evector(nestate,finalstate)*factor*kt!*fraction
                                                !Print*,distributionfunctionsvalue(3,k_d+1),k_d,talmip(7,sum1),evector(nestate,&
						!	& initialstate),factor
                                            endif

                                        endif

                                    enddo

                                enddo

                            endif

                        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        uplimit=sums(10-spin)-sums(9-spin)

                        do finalindex=1,uplimit

                            finalstate=kmatrix(k_d+1)+k_u*snm+sums(9-spin)+finalindex

                            nku1=n1dat(finalstate)
                            mku1=m1dat(finalstate)
                            nkd =n2dat(finalstate)
                            mkd =m2dat(finalstate)
                            nku2=n3dat(finalstate)
                            mku2=m3dat(finalstate)
                            ku1 =k1dat(finalstate)
                            kd  =k2dat(finalstate)
                            ku2 =k3dat(finalstate)

                            !call Talmisquare(nmax,nku1,-mku1,nku2,-mku2,nkd,-mkd,dble(ku1)+0.5,dble(ku2)+0.5,&
                            !    & dble(kd)+0.5,summaxk,talmik)

                            !do sum1=1,summaxp

                                    !do sum2=1,summaxk

                                        !if(talmip(1,sum1).eq.talmik(1,sum2).and.talmip(2,sum1).eq.talmik(2,sum2)) then

                                            !if(talmip(5,sum1).eq.talmik(5,sum2).and.talmip(6,sum1).eq.talmik(6,sum2)) then


                                                !bp=b*sqrt(dble(pu1+pu2+1)*(dble(pd)+0.5))/(dble(Kt)+0.5)
                                                !bk=b*sqrt(dble(ku1+ku2+1)*(dble(kd)+0.5))/(dble(Kt)+0.5)

                                                !factor=olaps(bp,nint(talmip(3,sum1)),nint(talmip(4,sum1)),bk,&
                                                !    & nint(talmik(3,sum2)),nint(talmik(4,sum2)),(1.0D0-fraction)*(dx+dxcorrect),&
                                                !    & (1.0D0-fraction)*dy)*(0.D0-1.0D0)**(mku1+mku2+mkd)*2*mq/&
                                                !    & (dx+dxcorrect)
                                                factor=olaps(b*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)),npu1,mpu1, &
                                                    & b*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)),nku1,-mku1, &
                                            & (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect),(dble(pu1)+0.5)/(dble(Kt)+0.5)*dy)*&
                                                    & olaps(b*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)),npu2,mpu2, &
                                                    & b*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)),nku2,-mku2, &
                                            & (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect),(dble(pu2)+0.5)/(dble(Kt)+0.5)*dy)*&
                                                    olaps(b*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)),npd,mpd, &
                                                    & b*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)),nkd,-mkd, &
                                                    & -(1.0D0-fraction)*(dx+dxcorrect),-(1.0D0-fraction)*dy) &
                                                    & *(0.D0-1.0D0)**(mku1+mku2+mkd)*2*mq/(dx+dxcorrect)

                                                !distributionfunctionsvalue_E=distributionfunctionsvalue_E+&
                                                !    & talmip(7,sum1)*evector(nestate,initialstate)*talmik(7,sum2)*&
                                                !    & evector(nestate,finalstate)*factor*(kt)!*fraction
                                                distributionfunctionsvalue_E=evector(nestate,initialstate)*&
                                                    & evector(nestate,finalstate)*factor*(kt)

                                                distributionfunctionsvalue(4,k_d+1)=distributionfunctionsvalue(4,k_d+1)&
                                                    & +distributionfunctionsvalue_E


                                                !Print*,distributionfunctionsvalue_E,initialstate,finalstate,factor,&
                                                !    & evector(nestate,initialstate),evector(nestate,finalstate)
                                                !Print*,factor,olaps(b*sqrt((dble(pu1)+0.5)/(dble(Kt)+0.5)),npu1,mpu1, &
                                                !    & b*sqrt((dble(ku1)+0.5)/(dble(Kt)+0.5)),nku1,-mku1, &
                                            !& (dble(pu1)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect),(dble(pu1)+0.5)/(dble(Kt)+0.5)*dy),&
                                            !    & olaps(b*sqrt((dble(pu2)+0.5)/(dble(Kt)+0.5)),npu2,mpu2, &
                                            !        & b*sqrt((dble(ku2)+0.5)/(dble(Kt)+0.5)),nku2,-mku2, &
                                            !& (dble(pu2)+0.5)/(dble(Kt)+0.5)*(dx+dxcorrect),(dble(pu2)+0.5)/(dble(Kt)+0.5)*dy),&
                                            !    & olaps(b*sqrt((dble(pd)+0.5)/(dble(Kt)+0.5)),npd,mpd, &
                                            !        & b*sqrt((dble(kd)+0.5)/(dble(Kt)+0.5)),nkd,-mkd, &
                                            !        & -(1.0D0-fraction)*(dx+dxcorrect),-(1.0D0-fraction)*dy),&
                                            !    & initialstate,finalstate
                                            !endif
                                        !endif
                                    !enddo
                                !enddo


                        enddo

                    endif

                enddo

            enddo

        enddo

        

    enddo



    Print*,"test"
    do test=1,kt
        Print*, distributionfunctionsvalue(1,test),distributionfunctionsvalue(2,test),distributionfunctionsvalue(3,test), &
        & distributionfunctionsvalue(4,test)
    enddo


    return
end subroutine GPD_for_d_quark

subroutine GPD_for_u_quark(nmax,Mj,Kt,b,sums,snm,nestate,dimtot,evector,dx,dy,distributionfunctionsvalue)
    use numbers
    use basis_info

    implicit none

    integer :: nmax,Mj,Kt,snm,nestate,dimtot
    double precision :: b,dx,dy
    integer, dimension(*) :: sums
    double precision,dimension(nestate,dimtot) :: evector

    double precision,dimension(3,Kt) :: distributionfunctionsvalue

    double precision :: fraction,bp,bk,factor,olaps
    integer :: npu1,mpu1,npd,mpd,npu2,mpu2,pu1,pu2,pd
    integer :: nku1,mku1,nkd,mkd,nku2,mku2,ku1,ku2,kd
    integer ::k_u,k,spin,initialstate,finalstate
    double precision, dimension(:,:), allocatable :: talmip,talmik
    integer :: summaxp,summaxk,looptmp,looptmk,talmitotnumber
    integer :: test,error

    talmitotnumber=(2*(nmax-3)+1)**2*(floor(dble(nmax-3)/2.0D0)+1)**2

    allocate(talmip(7,talmitotnumber),talmik(7,talmitotnumber),stat=error)
        if(error.ne.0) then
            print *, 'cannot allocate basis array'
        end if

    fraction=0.D0

    do k_u = 0 , Kt-1

        fraction=(k_u+0.5)/(Kt+0.5)

        distributionfunctionsvalue(3,k_u+1)=0.D0

        do k=0,dimtot-snm,snm

            do spin=1,8

                do initialstate=k+sums(spin)+1,k+sums(spin+1)

                    npu1=n1dat(initialstate)
                    mpu1=m1dat(initialstate)
                    npd =n2dat(initialstate)
                    mpd =m2dat(initialstate)
                    npu2=n3dat(initialstate)
                    mpu2=m3dat(initialstate)
                    pu1 =k1dat(initialstate)
                    pd  =k2dat(initialstate)
                    pu2 =k3dat(initialstate)

                    do finalstate=k+sums(spin)+1,k+sums(spin+1)

                        nku1=n1dat(finalstate)
                        mku1=m1dat(finalstate)
                        nkd =n2dat(finalstate)
                        mkd =m2dat(finalstate)
                        nku2=n3dat(finalstate)
                        mku2=m3dat(finalstate)
                        ku1 =k1dat(finalstate)
                        kd  =k2dat(finalstate)
                        ku2 =k3dat(finalstate)

                        if (pu1.eq.k_u) then

                            call Talmisquare(nmax,npd,mpd,npu2,mpu2,npu1,mpu1,dble(pd)+0.5,dble(pu2)+0.5,&
                                & dble(pu1)+0.5,summaxp,talmip)
                            call Talmisquare(nmax,nkd,mkd,nku2,mku2,nku1,mku1,dble(kd)+0.5,dble(ku2)+0.5,&
                                & dble(ku1)+0.5,summaxk,talmik)

                        endif

                        !if (pu2.eq.k_u) then

                        !    call Talmisquare(nmax,npd,mpd,npu1,mpu1,npu2,mpu2,dble(pd)+0.5,dble(pu1)+0.5,&
                        !        & dble(pu2)+0.5,summaxp,talmip)
                        !    call Talmisquare(nmax,nkd,mkd,nku1,mku1,nku2,mku2,dble(kd)+0.5,dble(ku1)+0.5,&
                        !        & dble(ku2)+0.5,summaxk,talmik)

                        !endif

                        do looptmp=1,summaxp

                            do looptmk=1,summaxk

                                if(talmip(1,looptmp).eq.talmik(1,looptmk).and.talmip(2,looptmp).eq.talmik(2,looptmk).and.&
                                    & talmip(5,looptmp).eq.talmik(5,looptmk).and.talmip(6,looptmp).eq.talmik(6,looptmk)) then

                                    if (pu1.eq.k_u) then

                                        bp=b*sqrt(dble(pu2+pd+1)*(dble(pu1)+0.5))/(dble(Kt)+0.5)
                                        bk=b*sqrt(dble(ku2+kd+1)*(dble(ku1)+0.5))/(dble(Kt)+0.5)

                                        factor=olaps(bp,nint(talmip(3,looptmp)),nint(talmip(4,looptmp)),bk,&
                                            & nint(talmik(3,looptmk)),nint(talmik(4,looptmk)),(1.0D0-fraction)*dx,&
                                            & (1.0D0-fraction)*dy)

                                    !else if (pu2.eq.k_u) then

                                    !    bp=b*sqrt(dble(pu1+pd+1)*(dble(pu2)+0.5))/(dble(Kt)+0.5)
                                    !    bk=b*sqrt(dble(ku1+kd+1)*(dble(ku2)+0.5))/(dble(Kt)+0.5)

                                    !    factor=olaps(bp,nint(talmip(3,looptmp)),nint(talmip(4,looptmp)),bk,&
                                    !        & nint(talmik(3,looptmk)),nint(talmik(4,looptmk)),(1.0D0-fraction)*dx,&
                                    !        & (1.0D0-fraction)*dy)

                                    else

                                        factor=0.D0

                                    endif

                                    distributionfunctionsvalue(1,k_u+1)=fraction
                                    distributionfunctionsvalue(2,k_u+1)=dx**2+dy**2
                                    distributionfunctionsvalue(3,k_u+1)=distributionfunctionsvalue(3,k_u+1)&
                                        & +talmip(7,looptmp)*evector(nestate,initialstate)*talmik(7,looptmk)*&
                                        & evector(nestate,finalstate)*factor*(kt)!*fraction

                                endif

                            enddo

                        enddo

                    enddo

                enddo

            enddo

        enddo

    enddo

    Print*,"test"
    do test=1,kt
        Print*, distributionfunctionsvalue(1,test),distributionfunctionsvalue(2,test),distributionfunctionsvalue(3,test)
    enddo

    return
end subroutine GPD_for_u_quark


subroutine TalmiSquare(nmax,n1,m1,n2,m2,n3,m3,kp1half,kp2half,kp3half,j,TalmisquareValue)

    implicit none

    integer :: nmax,n1,m1,n2,m2,n3,m3
    double precision :: kp1half,kp2half,kp3half
    double precision, dimension(7,*) :: TalmiSquareValue

    integer :: mscale,nscale
    integer :: capnf,capmf,nf,mf,capnm,capmm,nm,mm
    double precision :: TMC
    integer :: error,i,j

    mscale=nmax-3
    nscale=floor(dble(mscale)/2.0D0)
    j=0

    do capnf=0,nscale
        do capmf=-mscale,mscale
            do capnm=0,nscale
                do capmm=-mscale,mscale

                    mm=m1+m2-capmm
                    nm=n1+n2-capnm+(abs(m1)+abs(m2)-abs(capmm)-abs(mm))/2
                    mf=capmm+m3-capmf
                    nf=capnm+n3-capnf+(abs(capmm)+abs(m3)-abs(capmf)-abs(mf))/2

                    if(nm.ge.0.and.nf.ge.0) then

                        j=j+1

                        TalmiSquareValue(1,j)=capnf
                        TalmiSquareValue(2,j)=capmf
                        TalmiSquareValue(3,j)=nf
                        TalmiSquareValue(4,j)=mf
                        TalmiSquareValue(5,j)=nm
                        TalmiSquareValue(6,j)=mm

                        TalmiSquareValue(7,j)=TMC(capnm,capmm,nm,mm,n1,m1,n2,m2,sqrt(kp2half/kp1half))*TMC(capnf,capmf,nf,mf,&
                            & n3,m3,capnm,capmm,sqrt((kp1half+kp2half)/kp3half))

                    endif

                enddo
            enddo
        enddo
    enddo

   return
end subroutine TalmiSquare

subroutine average(loopnumber1,distributionfunctionskt,distributionfunctionskt_1,distributionfunctions)

    implicit none

    integer :: loopnumber1
    double precision, dimension(3,*) :: distributionfunctionskt,distributionfunctionskt_1
    double precision, dimension(3,*) :: distributionfunctions

    integer :: i,test

    do i=1,loopnumber1
        distributionfunctions(1,i)=(distributionfunctionskt(1,i)+distributionfunctionskt_1(1,i))/2.D0
        distributionfunctions(2,i)=distributionfunctionskt(2,i)
        distributionfunctions(3,i)=(distributionfunctionskt(3,i)+distributionfunctionskt_1(3,i))/2.0D0
    enddo

    !do test=1,loopnumber1
    !    Print*,distributionfunctions(1,test),distributionfunctions(2,test),distributionfunctions(3,test)
    !enddo

    return
end subroutine average
