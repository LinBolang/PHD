! In this script, we mainly calculate the Hamiltonian of proton. And subroutine of hamiltonian call by wfproduction.f90
! In this script, we need to use a few functions which are defined in function.f90
! The first part is the hamiltonian for three particle Fock sector where include Kinetic, transverse confining potential, longitudinal confining potential, and one-gluon exchange potential.

  subroutine hamiltoniankinetic3p(nmax2,Mj,Kt,i_nzk,j_nzk,hamiltoniankineticvalue)
    use numbers
    use basis_info
    use colorfactor
    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mu,md,lag,ms
    integer,dimension(*) :: i_nzk,j_nzk
    double precision, dimension(*) :: hamiltoniankineticvalue

    integer :: nz,jlist(100)
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,nk2,nk3,mk1,mk2,mk3
    integer :: i,j,i1,i2,i3,i4,j1,j2,j3,j4,j5,j6
    complex*16 ::imageunit,fourierphase
    integer :: fourierphasen,fourierphasem

    double precision :: kp1half,kp2half,kp3half,Pplus
    double precision :: Ememe,Ep1p1,Ep3p3,Ep2p2,Ep1p2,Ep1p3,Ep3p2,shift
    double precision :: kinetic1,kinetic2,kinetic3,hamitot,hamicm,hamirel
    complex*16 :: lagrangeterm
    integer :: fidelta
    double precision :: adotaq,adotbq
    integer :: test

    mu=Mu1_DEFAULT
    md=Md1_DEFAULT
    ms=Ms1_DEFAULT
    b =B_DEFAULT
    lag=LAGRANGEMULTIPLIER

    nz=1
    imageunit=(0.D0,1.D0)
  
    do i=1,dimtot1

      i_nzk(i)=nz

      kp1=k1dat1(i)
      kp2=k2dat1(i)
      kp3=k3dat1(i)
      sp1=s1dat1(i)
      sp2=s2dat1(i)
      sp3=s3dat1(i)
      np1=n1dat1(i)
      mp1=m1dat1(i)
      np2=n2dat1(i)
      mp2=m2dat1(i)
      np3=n3dat1(i)
      mp3=m3dat1(i)

      do i1=-1,1
        do i2=-1,1
          do i3=-1,1,2
            j1=(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i1,mp1+i3,np2+i2&
              &,mp2-i3,np3,mp3,jlist(j1))
            j2=18+(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2+i1,mp2+i3,&
              &np3+i2,mp3-i3,jlist(j2))
            j3=36+(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i1,mp1+i3,np2,mp2,&
              & np3+i2,mp3-i3,jlist(j3))
          enddo
        enddo
      enddo

      do i4=-1,1,2
        j4=54+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i4,mp1,np2,mp2,np3,mp3,jlist(j4))
        j5=56+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2+i4,mp2,np3,mp3,jlist(j5))
        j6=58+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3+i4,mp3,jlist(j6))
      enddo

      call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3,mp3,jlist(61))

  		do j=1,61

        if(jlist(j).lt.0) then
          cycle
        endif

        kk1=k1dat1(jlist(j))
        kk2=k2dat1(jlist(j))
        kk3=k3dat1(jlist(j))
        sk1=s1dat1(jlist(j))
        sk2=s2dat1(jlist(j))
        sk3=s3dat1(jlist(j))
        nk1=n1dat1(jlist(j))
        mk1=m1dat1(jlist(j))
        nk2=n2dat1(jlist(j))
        mk2=m2dat1(jlist(j))
        nk3=n3dat1(jlist(j))
        mk3=m3dat1(jlist(j))


        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        Pplus   = kp1half+kp2half+kp3half
        fourierphasen = -np1-np2-np3+nk1+nk2+nk3
        Fourierphasem = -abs(mp1)-abs(mp2)-abs(mp3)+abs(mk1)+abs(mk2)+abs(mk3)
        fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem      

        Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
          & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
          & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Ep3p3=adotaq(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)

        Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)*fidelta(np3,mp3&
          & ,nk3,mk3,sp3,sk3,kp3,kk3)

        Ep3p2=adotbq(np3,mp3,nk3,mk3,np2,mp2,nk2,mk2,sp3,sk3,sp2,sk2,kp3,kk3,kp2,kk2)*fidelta(np1,mp1,&
          & nk1,mk1,sp1,sk1,kp1,kk1)

        Ep1p3=adotbq(np1,mp1,nk1,mk1,np3,mp3,nk3,mk3,sp1,sk1,sp3,sk3,kp1,kk1,kp3,kk3)*fidelta(np2,mp2,nk2&
          & ,mk2,sp2,sk2,kp2,kk2)

        shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        kinetic1=(mu**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
        kinetic2=(md**2*Ememe+b**2*(kp2half/Pplus)*Ep2p2)/kp2half
        kinetic3=(ms**2*Ememe+b**2*(kp3half/Pplus)*Ep3p3)/kp3half
            

        hamitot=kinetic2+kinetic3+kinetic1
        hamicm =b**2*(Ep1p1*kp1half/Pplus+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)+Ep3p3*kp3half/Pplus&
          & + 2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+ 2*Ep3p2*sqrt(kp2half*kp3half/Pplus**2)+Ep2p2&
          & *kp2half/Pplus)/Pplus
        hamirel=hamitot-hamicm
        lagrangeterm=lag*b**2*((Ep1p1*kp1half/Pplus+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
          & +Ep3p3*kp3half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep3p2*sqrt(kp2half*kp3half/&
          & Pplus**2)+Ep2p2*kp2half/Pplus)+kappa4*(Ep1p1*kp1half/Pplus+2*Ep1p3*sqrt(kp1half*kp3half/&
          & Pplus**2)+Ep3p3*kp3half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep3p2*sqrt(kp2half&
          & *kp3half/Pplus**2)+Ep2p2*kp2half/Pplus)*fourierphase-shift)/Pplus

            


        if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

          j_nzk(nz)=jlist(j)
          hamiltoniankineticvalue(nz)=hamitot!hamirel+Real(lagrangeterm)
          nz=nz+1

        endif
      enddo
    enddo

    i_nzk(dimtot1+1)=nz
    !print*,"test",ms,mu,md

    Print*, "kinetichami_3p"
    ! do test=1,nz
    !  Print*, hamiltoniankineticvalue(test), i_nzk(test),j_nzk(test),test
    ! enddo

    return
  end subroutine hamiltoniankinetic3p

  subroutine hamiltonianinteractsw3p(nmax2,Mj,Kt,i_nztsw,j_nztsw,hamiltonianinteractswvalue)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b
    integer, dimension(*) :: i_nztsw,j_nztsw
    double precision, dimension(*) :: hamiltonianinteractswvalue

    integer :: nz,jlist(100)
    double precision :: hamisingleparticle
    complex*16 :: hamiCM,hamiinteract
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,nk2,nk3,mk1,mk2,mk3
    double precision :: Es1s1,Es2s2,Es3s3
    double precision ::Epu1pu1,Epu2pu2,Epdpd,Ememe,Epu1pu2,Epu1pd,Epu2pd,shifts
    double precision ::adotaq,adotbq,adotas
    integer ::fidelta,fourierphasen,fourierphasem
    complex*16 :: imageunit=(0.D0,1.D0),fourierphase
    double precision :: kp1half,kp2half,kp3half,Pplus

    !integer :: test
    integer :: i,j,i1,i2,i3,i4,j1,j2,j3,j4,j5,j6

    nz=1
    hamiinteract=(0.D0,0.D0)
    hamisingleparticle=(0.D0,0.D0)
    hamiCM=(0.D0,0.D0)
    b=B_DEFAULT

    do i=1,dimtot1

      i_nztsw(i)=nz

      kp1=k1dat1(i)
      kp2=k2dat1(i)
      kp3=k3dat1(i)
      sp1=s1dat1(i)
      sp2=s2dat1(i)
      sp3=s3dat1(i)
      np1=n1dat1(i)
      mp1=m1dat1(i)
      np2=n2dat1(i)
      mp2=m2dat1(i)
      np3=n3dat1(i)
      mp3=m3dat1(i)

      do i1=-1,1
        do i2=-1,1
          do i3=-1,1,2
            j1=(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i1,mp1+i3,np2+i2&
                  &,mp2-i3,np3,mp3,jlist(j1))
            j2=18+(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2+i1,mp2+i3,&
                &np3+i2,mp3-i3,jlist(j2))
            j3=36+(i1+1)*6+(i2+1)*2+(i3+1)/2+1
            call search(nmax2,nmax2,Mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i1,mp1+i3,np2,mp2,&
                  & np3+i2,mp3-i3,jlist(j3))
          enddo
        enddo
      enddo

      do i4=-1,1,2
        j4=54+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1+i4,mp1,np2,mp2,np3,mp3,jlist(j4))
        j5=56+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2+i4,mp2,np3,mp3,jlist(j5))
        j6=58+(i4+1)/2+1
        call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3+i4,mp3,jlist(j6))
      enddo
  
      call search(nmax2,nmax2,mj,Kt,kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3,mp3,jlist(61))


      do j=1,61

        if(jlist(j).lt.0) then
          cycle
        endif

        kk1=k1dat1(jlist(j))
        kk2=k2dat1(jlist(j))
        kk3=k3dat1(jlist(j))
        sk1=s1dat1(jlist(j))
        sk2=s2dat1(jlist(j))
        sk3=s3dat1(jlist(j))
        nk1=n1dat1(jlist(j))
        mk1=m1dat1(jlist(j))
        nk2=n2dat1(jlist(j))
        mk2=m2dat1(jlist(j))
        nk3=n3dat1(jlist(j))
        mk3=m3dat1(jlist(j))

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        Pplus=kp1half+kp2half+kp3half

        fourierphasen = -np1-np2-np3+nk1+nk2+nk3
        Fourierphasem = -abs(mp1)-abs(mp2)-abs(mp3)+abs(mk1)+abs(mk2)+abs(mk3)
        fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem



        Es1s1=adotas(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Es2s2=adotas(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Es3s3=adotas(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
          & fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)

        Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
          & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Epu1pu1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
          & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Epu2pu2=adotaq(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)

        Epdpd=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

        Epu1pd=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)*fidelta(np3,mp3&
          & ,nk3,mk3,sp3,sk3,kp3,kk3)

        Epu2pd=adotbq(np3,mp3,nk3,mk3,np2,mp2,nk2,mk2,sp3,sk3,sp2,sk2,kp3,kk3,kp2,kk2)*fidelta(np1,mp1,&
          & nk1,mk1,sp1,sk1,kp1,kk1)

        Epu1pu2=adotbq(np1,mp1,nk1,mk1,np3,mp3,nk3,mk3,sp1,sk1,sp3,sk3,kp1,kk1,kp3,kk3)*fidelta(np2,mp2,nk2&
          & ,mk2,sp2,sk2,kp2,kk2)

        shifts=2*(Sqrt(kappat1)*b**2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)/Pplus

        hamisingleparticle=kappat1*b**2*(Es1s1+Es2s2+Es3s3)/Pplus
        hamiCM=kappat1*b**2*((Epu1pu1*kp1half/Pplus+2*Epu1pu2*sqrt(kp1half*kp3half/Pplus**2)+Epu2pu2*&
          & kp3half/Pplus+2*Epu1pd*Sqrt(kp1half*kp2half/Pplus**2)+2*Epu2pd*sqrt(kp3half*kp2half/&
          & Pplus**2)+Epdpd*kp2half/Pplus)*fourierphase)/Pplus
        hamiinteract=hamisingleparticle-hamiCM-shifts
            

        if (Abs(Real(hamiinteract)).gt.epsilon(Real(hamiinteract))) then
          if(Abs(Aimag(hamiinteract)).lt.epsilon(Aimag(hamiinteract))) then

            j_nztsw(nz)=jlist(j)

            hamiltonianinteractswvalue(nz)=Real(hamiinteract)
  
            nz=nz+1
          endif
        endif
      enddo
    enddo
    i_nztsw(dimtot1+1)=nz

    Print*, "transhami_3p"
    !do test=1,nz
    !  Print*, hamiltonianinteractswvalue(test),i_nztsw(test),j_nztsw(test)
    !enddo


    return
  end subroutine hamiltonianinteractsw3p

  subroutine hamiltonianinteractlongi3p(nmax2,Mj,Kt,i_nzlsw,j_nzlsw,hamiltonianinteractlongivalue)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1,mass2,mass3
    integer, dimension(*) :: i_nzlsw,j_nzlsw
    double precision, dimension(*) :: hamiltonianinteractlongivalue

    integer :: i,j,nz,fidelta,KroneckerDelta
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,nk2,nk3,mk1,mk2,mk3
    double precision :: kp1half,kp2half,kp3half,Pplus,kk1half,kk2half,kk3half
    double precision :: V12,V13,V23,hamiinteract,shift
    double precision :: longipartial_2,longitrans!,longipartial1
    integer :: test
    !double precision :: test_longi


    b=B_DEFAULT
    mass1=Mu2_DEFAULT
    mass2=Md2_DEFAULT
    mass3=MS2_DEFAULT

    nz=1

    do i=1,dimtot1

      i_nzlsw(i)=nz

      kp1=k1dat1(i)
      kp2=k2dat1(i)
      kp3=k3dat1(i)
      sp1=s1dat1(i)
      sp2=s2dat1(i)
      sp3=s3dat1(i)
      np1=n1dat1(i)
      mp1=m1dat1(i)
      np2=n2dat1(i)
      mp2=m2dat1(i)
      np3=n3dat1(i)
      mp3=m3dat1(i)

      do j=1,dimtot1

        kk1=k1dat1(j)
        kk2=k2dat1(j)
        kk3=k3dat1(j)
        sk1=s1dat1(j)
        sk2=s2dat1(j)
        sk3=s3dat1(j)
        nk1=n1dat1(j)
        mk1=m1dat1(j)
        nk2=n2dat1(j)
        mk2=m2dat1(j)
        nk3=n3dat1(j)
        mk3=m3dat1(j)

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        kk1half=dble(kk1)+0.5D0
        kk2half=dble(kk2)+0.5D0
        kk3half=dble(kk3)+0.5D0
        Pplus  =kp1half+kp2half+kp3half

        !V12=longipartial5(Kt,kp1,kp2,kk1,kk2,kp3)*longitrans(Nmax2,Pplus,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,&
        V12=longipartial_2(Kt,kp1,kp2,kk1,kk2)*longitrans(Nmax2,Pplus,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,&
          & kp1half,kp2half,kk1half,kk2half)*KroneckerDelta(mp1+mp2,mk1+mk2)*KroneckerDelta(sp1,sk1)*&
          & KroneckerDelta(sp2,sk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)/(mass1+mass2)**2

        !V13=longipartial5(Kt,kp1,kp3,kk1,kk3,kp2)*longitrans(Nmax2,Pplus,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,&
        V13=longipartial_2(Kt,kp1,kp3,kk1,kk3)*longitrans(Nmax2,Pplus,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,&
          & kp1half,kp3half,kk1half,kk3half)*KroneckerDelta(mp1+mp3,mk1+mk3)*KroneckerDelta(sp1,sk1)*&
          & KroneckerDelta(sp3,sk3)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)/(mass1+mass3)**2

        !V23=longipartial5(Kt,kp2,kp3,kk2,kk3,kp1)*longitrans(Nmax2,Pplus,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,&
        V23=longipartial_2(Kt,kp2,kp3,kk2,kk3)*longitrans(Nmax2,Pplus,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,&
          & kp2half,kp3half,kk2half,kk3half)*KroneckerDelta(mp2+mp3,mk2+mk3)*KroneckerDelta(sp2,sk2)*&
          & KroneckerDelta(sp3,sk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)/(mass2+mass3)**2

        shift=sqrt(kappal1)*b**2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
          & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)/Pplus

        hamiinteract=(V12+V13+V23)*(0.D0-kappal1*b**4)/Pplus!-shift
        ! hamiinteract=(V12)*(0.D0-kappal1*b**4)/Pplus

        if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
          j_nzlsw(nz)=j

          hamiltonianinteractlongivalue(nz)=hamiinteract

          nz=nz+1

        endif

      enddo
    enddo

    i_nzlsw(dimtot1+1)=nz

    Print*, "longihami_3p"
    ! do test=1,nz
    !   Print*,hamiltonianinteractlongivalue(test),i_nzlsw(test),j_nzlsw(test),test
    ! enddo

    return
  end subroutine hamiltonianinteractlongi3p

  subroutine hamiltonianinteractOGE3p(nmax2,Mj,Kt,i_nz_oge,j_nz_oge,hamiltonianinteractogevalue)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1,mass2,mass3,mass4
    integer,dimension(*) :: i_nz_oge,j_nz_oge
    double precision,dimension(*) :: hamiltonianinteractogevalue

    integer :: dim1,dim2,error,nz,i,j
    double precision :: kp1half,kp2half,kp3half,kk1half,kk2half,kk3half,Pplus
    double precision,external :: IED,IED1momentum,IED2momentum
    integer,external :: fidelta
    double precision :: V12,V13,V23,couplings,hamiinteract
    double precision :: mu,md,ms,mg
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,mk1,nk2,mk2,nk3,mk3
    !double precision :: testintegration
    integer :: testvalue,test

    dim1=(nmax2+1)**2*Kt**2
    dim2=Kt**2
    b=B_DEFAULT
    mass1=Mu2_DEFAULT
    mass2=Md2_DEFAULT
    mass3=MS2_DEFAULT
    mass4=mg2_default
    couplings=coupling3


    allocate(integration(dim1,dim2,3),stat=error)
    if(error.ne.0) then
        print *, 'cannot allocate basis array'
    end if

    call IEDintegration3p(nmax2,Kt,mass1,mass2,mass3,mass4,b)
    !call dimtotal(nmax1,nmax2,Mj,Kt,dimtot)
    !call cpu_time(integrate_time)

    !Print*,"oge integration matrix : read finish" , integrate_time-start_time
    nz=1

    do i=1,dimtot1
    !do i=1199,1199
    !do i=9,9
      i_nz_oge(i)=nz
      
        kp1=k1dat1(i)
        kp2=k2dat1(i)
        kp3=k3dat1(i)
        sp1=s1dat1(i)
        sp2=s2dat1(i)
        sp3=s3dat1(i)
        np1=n1dat1(i)
        mp1=m1dat1(i)
        np2=n2dat1(i)
        mp2=m2dat1(i)
        np3=n3dat1(i)
        mp3=m3dat1(i)

      do j=1,dimtot1
      !do j=1199,1199
      !do j=11,11
          kk1=k1dat1(j)
          kk2=k2dat1(j)
          kk3=k3dat1(j)
          sk1=s1dat1(j)
          sk2=s2dat1(j)
          sk3=s3dat1(j)
          nk1=n1dat1(j)
          mk1=m1dat1(j)
          nk2=n2dat1(j)
          mk2=m2dat1(j)
          nk3=n3dat1(j)
          mk3=m3dat1(j)

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        kk1half=dble(kk1)+0.5D0
        kk2half=dble(kk2)+0.5D0
        kk3half=dble(kk3)+0.5D0
        Pplus  =kp1half+kp2half+kp3half

        mu=mass1
        md=mass2
        ms=mass3
        mg=mass4

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 2 particle #####################!!!!!!!!!!!
        testvalue=fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
        V12=0.D0
        ! Print*,sp1,sp2,sk1,sk2

        If(testvalue.ne.0) then
          if(sp1.eq.sk1.and.sp2.eq.sk2) then

          V12=(mu**2*Pplus**2/kp1half/kk1half+md**2*Pplus**2/kp2half/kk2half)*IED(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,&
              & kk1,kk2,mu,md,mg,1)
          ! V12=IED(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,&
          !     & kk1,kk2,mu,md,mg,1)

          ! Print*,V12

            if(sp1.eq.1.and.sp2.eq.1) then

          V12=V12+IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,1)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,1)*Pplus**2/kk2half/kp2half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,1)*Pplus**2/kk2half/kp1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,1)*Pplus**2/kk1half/kp2half 

            else if(sp1.eq.-1.and.sp2.eq.-1) then

          V12=V12+IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,3,1)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,4,1)*Pplus**2/kk2half/kp2half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,3,1)*Pplus**2/kk1half/kp2half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,4,1)*Pplus**2/kk2half/kp1half

            else if(sp1.eq.1.and.sp2.eq.-1) then

          V12=V12+IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,1)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,4,1)*Pplus**2/kk2half/kp2half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,4,1)*Pplus**2/kk2half/kk1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,1)*Pplus**2/kp1half/kp2half

            else if(sp1.eq.-1.and.sp2.eq.1) then

          V12=V12+IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,3,1)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,1)*Pplus**2/kk2half/kp2half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,3,1)*Pplus**2/kk2half/kk1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,1)*Pplus**2/kp1half/kp2half 

            endif

          else if(sp1.eq.sk1.and.sp2.ne.sk2) then

            if(sp1.eq.1.and.sp2.eq.1) then

            V12=md*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,1))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kk2half-kp2half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,1) &
              & /kp1half/kp2half/kk2half)

            else if(sp1.eq.-1.and.sp2.eq.1) then

            V12=md*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,1))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kk2half-kp2half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,1) &
              & /kk1half/kp2half/kk2half)

            else if(sp1.eq.1.and.sp2.eq.-1) then

            V12=md*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,1))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kp2half-kk2half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,1) &
              & /kk1half/kp2half/kk2half)

            else if(sp1.eq.-1.and.sp2.eq.-1) then

            V12=md*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,1))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kp2half-kk2half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,1) &
              & /kp1half/kp2half/kk2half)

            endif

          else if(sp1.ne.sk1.and.sp2.eq.sk2) then

            if(sp1.eq.1.and.sp2.eq.1) then

            V12=mu*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,1))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kk1half-kp1half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,2,1) &
              & /kp1half/kp2half/kk1half)

            else if(sp1.eq.1.and.sp2.eq.-1) then

            V12=mu*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,2,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,2,1))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kk1half-kp1half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,2,1) &
              & /kp1half/kk1half/kk2half)

            else if(sp1.eq.-1.and.sp2.eq.1) then

            V12=mu*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,1))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kp1half-kk1half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,4,1,1) &
              & /kp1half/kk2half/kk1half)

            else if(sp1.eq.-1.and.sp2.eq.-1) then

            V12=mu*((IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,3,1,1)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,1,1,1))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kp1half-kk1half)*IED1momentum(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,b,2,1,1) &
              & /kp1half/kp2half/kk1half)

            endif

          else if(sp2.eq.sk1.and.sp1.eq.sk2.and.sp1.ne.sp2) then

          V12=mu*md*Pplus**2*(kk1half-kp1half)*(kk2half-kp2half)/kk1half/kk2half/kp1half/kp2half* &
            & IED(nmax2,Kt,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1,kp2,kk1,kk2,mu,md,mg,1)

          else if(sp1.eq.sp2.and.sk1.eq.sk2.and.sp1.ne.sk1) then

            V12=0.D0

          endif

        endif
        ! Print*,V12

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 3 particle #####################!!!!!!!!!!!

        testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)
        V13=0.D0

        If(testvalue.ne.0) then

          if(sp1.eq.sk1.and.sp3.eq.sk3) then

          V13=(mu**2*Pplus**2/kp1half/kk1half+ms**2*Pplus**2/kp3half/kk3half)*IED(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,&
              & kk1,kk3,mu,ms,mg,2)

            if(sp1.eq.1.and.sp3.eq.1) then

          V13=V13+IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,1,2)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,2,2)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,1,2)*Pplus**2/kk3half/kp1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,2,2)*Pplus**2/kk1half/kp3half 

            else if(sp1.eq.-1.and.sp3.eq.-1) then

          V13=V13+IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,3,2)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,4,2)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,3,2)*Pplus**2/kk1half/kp3half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,4,2)*Pplus**2/kk3half/kp1half

            else if(sp1.eq.1.and.sp3.eq.-1) then

          V13=V13+IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,1,2)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,4,2)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,4,2)*Pplus**2/kk3half/kk1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,1,2)*Pplus**2/kp1half/kp3half

            else if(sp1.eq.-1.and.sp3.eq.1) then

          V13=V13+IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,3,2)*Pplus**2/kk1half/kp1half &
              & +IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,2,2)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,3,2)*Pplus**2/kk3half/kk1half &
              & -IED2momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,2,2)*Pplus**2/kp1half/kp3half 

            endif

          else if(sp1.eq.sk1.and.sp3.ne.sk3) then

            if(sp1.eq.1.and.sp3.eq.1) then

            V13=ms*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,2,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,2,2))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kk3half-kp3half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,2,2) &
              & /kp1half/kp3half/kk3half)

            else if(sp1.eq.-1.and.sp3.eq.1) then

            V13=ms*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,2,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,2,2))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kk3half-kp3half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,2,2) &
              & /kk1half/kp3half/kk3half)

            else if(sp1.eq.1.and.sp3.eq.-1) then

            V13=ms*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,1,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,1,2))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kp3half-kk3half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,1,2) &
              & /kk1half/kp3half/kk3half)

            else if(sp1.eq.-1.and.sp3.eq.-1) then

            V13=ms*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,1,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,1,2))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kp3half-kk3half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,1,2) &
              & /kp1half/kp3half/kk3half)

            endif

          else if(sp1.ne.sk1.and.sp3.eq.sk3) then

            if(sp1.eq.1.and.sp3.eq.1) then

            V13=mu*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,2,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,2,2))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kk1half-kp1half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,2,2) &
              & /kp1half/kp3half/kk1half)

            else if(sp1.eq.1.and.sp3.eq.-1) then

            V13=mu*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,2,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,2,2))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kk1half-kp1half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,2,2) &
              & /kp1half/kk1half/kk3half)

            else if(sp1.eq.-1.and.sp3.eq.1) then

            V13=mu*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,1,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,1,2))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kp1half-kk1half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,4,1,2) &
              & /kp1half/kk3half/kk1half)

            else if(sp1.eq.-1.and.sp3.eq.-1) then

            V13=mu*((IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,3,1,2)- &
              & IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,1,1,2))*Pplus**2/kp1half/kk1half+ &
              & Pplus**2*(kp1half-kk1half)*IED1momentum(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,b,2,1,2) &
              & /kp1half/kp3half/kk1half)

            endif

          else if(sp3.eq.sk1.and.sp1.eq.sk3.and.sp1.ne.sp3) then

          V13=mu*ms*Pplus**2*(kk1half-kp1half)*(kk3half-kp3half)/kk1half/kk3half/kp1half/kp3half* &
            & IED(nmax2,Kt,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,kp1,kp3,kk1,kk3,mu,ms,mg,2)

          else if(sp1.eq.sp3.and.sk1.eq.sk3.and.sp1.ne.sk1) then

            V13=0.D0

          endif

        endif

        !!!!!!!!!!########### calculate the one gluon exchange for 2 and 3 particle #####################!!!!!!!!!!!
        testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)
        V23=0.D0

        If(testvalue.ne.0) then

          if(sp2.eq.sk2.and.sp3.eq.sk3) then

          V23=(md**2*Pplus**2/kp2half/kk2half+ms**2*Pplus**2/kp3half/kk3half)*IED(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,&
              & kk2,kk3,ms,md,mg,3)

            if(sp2.eq.1.and.sp3.eq.1) then

          V23=V23+IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,1,3)*Pplus**2/kk2half/kp2half &
              & +IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,2,3)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,1,3)*Pplus**2/kk3half/kp2half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,2,3)*Pplus**2/kk2half/kp3half 

            else if(sp2.eq.-1.and.sp3.eq.-1) then

          V23=V23+IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,3,3)*Pplus**2/kk2half/kp2half &
              & +IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,4,3)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,3,3)*Pplus**2/kk2half/kp3half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,4,3)*Pplus**2/kk3half/kp2half

            else if(sp2.eq.1.and.sp3.eq.-1) then

          V23=V23+IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,1,3)*Pplus**2/kk2half/kp2half &
              & +IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,4,3)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,4,3)*Pplus**2/kk3half/kk2half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,1,3)*Pplus**2/kp2half/kp3half

            else if(sp2.eq.-1.and.sp3.eq.1) then

          V23=V23+IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,3,3)*Pplus**2/kk2half/kp2half &
              & +IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,2,3)*Pplus**2/kk3half/kp3half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,3,3)*Pplus**2/kk3half/kk2half &
              & -IED2momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,2,3)*Pplus**2/kp2half/kp3half 

            endif

          else if(sp2.eq.sk2.and.sp3.ne.sk3) then

            if(sp2.eq.1.and.sp3.eq.1) then

          V23=ms*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,2,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,2,3))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kk3half-kp3half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,2,3) &
              & /kp2half/kp3half/kk3half)

            else if(sp2.eq.-1.and.sp3.eq.1) then

          V23=ms*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,2,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,2,3))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kk3half-kp3half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,2,3) &
              & /kk2half/kp3half/kk3half)

            else if(sp2.eq.1.and.sp3.eq.-1) then

          V23=ms*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,1,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,1,3))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kp3half-kk3half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,1,3) &
              & /kk2half/kp3half/kk3half)

            else if(sp2.eq.-1.and.sp3.eq.-1) then

          V23=ms*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,1,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,1,3))*Pplus**2/kp3half/kk3half+ &
              & Pplus**2*(kp3half-kk3half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,1,3) &
              & /kp2half/kp3half/kk3half)

            endif

          else if(sp2.ne.sk2.and.sp3.eq.sk3) then

            if(sp2.eq.1.and.sp3.eq.1) then

            V23=md*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,2,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,2,3))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kk2half-kp2half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,2,3) &
              & /kp2half/kp3half/kk2half)

            else if(sp2.eq.1.and.sp3.eq.-1) then

            V23=md*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,2,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,2,3))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kk2half-kp2half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,2,3) &
              & /kp2half/kk2half/kk3half)

            else if(sp2.eq.-1.and.sp3.eq.1) then

            V23=md*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,1,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,1,3))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kp2half-kk2half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,4,1,3) &
              & /kp2half/kk3half/kk2half)

            else if(sp2.eq.-1.and.sp3.eq.-1) then

            V23=md*((IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,3,1,3)- &
              & IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,1,1,3))*Pplus**2/kp2half/kk2half+ &
              & Pplus**2*(kp2half-kk2half)*IED1momentum(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,b,2,1,3) &
              & /kp2half/kp3half/kk2half)

            endif

          else if(sp3.eq.sk2.and.sp2.eq.sk3.and.sp2.ne.sp3) then

            V23=md*ms*Pplus**2*(kk2half-kp2half)*(kk3half-kp3half)/kk2half/kk3half/kp2half/kp3half* &
              & IED(nmax2,Kt,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,kp2,kp3,kk2,kk3,ms,md,mg,3)

          else if(sp2.eq.sp3.and.sk2.eq.sk3.and.sp2.ne.sk2) then

            V23=0.D0

          endif

        endif

        hamiinteract=4.0D0*couplings/3.0D0/Pplus**2*(V12+V13+V23)
        ! hamiinteract=4.0D0*couplings/3.0D0/Pplus**2*(V12)

        if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
          j_nz_oge(nz)=j

          hamiltonianinteractogevalue(nz)=hamiinteract

          nz=nz+1

        endif


      enddo
    enddo

    i_nz_oge(dimtot1+1)=nz



    Print*, "one_gluon_exchange_3p"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue(test),i_nz_oge(test),j_nz_oge(test),test
    ! enddo


    deallocate(integration)

    return
  end subroutine hamiltonianinteractOGE3p

  subroutine hamiltonianInstantaneous3p(nmax2,Mj,Kt,i_nz_inst,j_nz_inst,hamiltonianinteractinstvalue)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,binst
    integer,dimension(*) :: i_nz_inst,j_nz_inst
    double precision,dimension(*) :: hamiltonianinteractinstvalue

    integer :: dim1,dim2,nz,i,j
    double precision :: kp1half,kp2half,kp3half,kk1half,kk2half,kk3half,Pplus
    ! double precision,external :: IED,IED1momentum,IED2momentum
    integer,external :: fidelta
    double precision :: V12,V13,V23,couplings,hamiinteract
    !double precision :: mu,md,ms,mg
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,mp1,np2,mp2,np3,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,mk1,nk2,mk2,nk3,mk3
    !double precision :: testintegration
    integer :: testvalue,test

    dim1=(nmax2+1)**2*Kt**2
    dim2=Kt**2
    b=B_DEFAULT
    binst=binst_default
    couplings=coupling3
    !couplings=0

    nz=1

    do i=1,dimtot1
    ! do i=7,7

      i_nz_inst(i)=nz
      
        kp1=k1dat1(i)
        kp2=k2dat1(i)
        kp3=k3dat1(i)
        sp1=s1dat1(i)
        sp2=s2dat1(i)
        sp3=s3dat1(i)
        np1=n1dat1(i)
        mp1=m1dat1(i)
        np2=n2dat1(i)
        mp2=m2dat1(i)
        np3=n3dat1(i)
        mp3=m3dat1(i)

      do j=1,dimtot1
      ! do j=491,491

          kk1=k1dat1(j)
          kk2=k2dat1(j)
          kk3=k3dat1(j)
          sk1=s1dat1(j)
          sk2=s2dat1(j)
          sk3=s3dat1(j)
          nk1=n1dat1(j)
          mk1=m1dat1(j)
          nk2=n2dat1(j)
          mk2=m2dat1(j)
          nk3=n3dat1(j)
          mk3=m3dat1(j)

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        kk1half=dble(kk1)+0.5D0
        kk2half=dble(kk2)+0.5D0
        kk3half=dble(kk3)+0.5D0
        Pplus  =kp1half+kp2half+kp3half

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 2 particle #####################!!!!!!!!!!!

        testvalue=fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
        V12=0.D0

        If(testvalue.ne.0.and.sp1.eq.sk1.and.sp2.eq.sk2) then
          
          call Instantaneous3p(Nmax2,Kt,b,binst,kp1,kp2,kk1,kk2,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,V12)

        endif


        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 3 particle #####################!!!!!!!!!!!

        testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)
        V13=0.D0

        If(testvalue.ne.0.and.sp1.eq.sk1.and.sp3.eq.sk3) then

          call Instantaneous3p(Nmax2,Kt,b,binst,kp1,kp3,kk1,kk3,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,V13)

        endif


        !!!!!!!!!!########### calculate the one gluon exchange for 2 and 3 particle #####################!!!!!!!!!!!

        testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)
        V23=0.D0

        If(testvalue.ne.0.and.sp2.eq.sk2.and.sp3.eq.sk3) then

          call Instantaneous3p(Nmax2,Kt,b,binst,kp2,kp3,kk2,kk3,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,V23)

        endif

        hamiinteract=-couplings**2*2/3/Pplus*(V12+V13+V23)

        if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
          j_nz_inst(nz)=j

          hamiltonianinteractinstvalue(nz)=hamiinteract

          nz=nz+1

        endif


      enddo
    enddo

    i_nz_inst(dimtot1+1)=nz



    Print*, "Instantaneous_3p"
    ! do test=1,nz
    !   Print*,hamiltonianinteractinstvalue(test),i_nz_inst(test),j_nz_inst(test),test
    ! enddo

    return
  end subroutine hamiltonianInstantaneous3p

! The second part is the hamiltonian for the four particle (qqqg) Fock sector which include kinetic term and q->qg interaction,g->qqbar interaction.

subroutine hamiltoniankinetic4p(nmax2,Mj,Kt,i_nzk2,j_nzk2,hamiltoniankineticvalue2)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,lag,mass1,mass2,mass3,mass4
    integer,dimension(*) :: i_nzk2,j_nzk2
    double precision, dimension(*) :: hamiltoniankineticvalue2

    integer :: nz,color
    integer :: kp1,kp2,kp3,kp4,sp1,sp2,sp3,sp4,np1,np2,np3,np4,mp1,mp2,mp3,mp4
    integer :: kk1,kk2,kk3,kk4,sk1,sk2,sk3,sk4,nk1,nk2,nk3,nk4,mk1,mk2,mk3,mk4
    integer :: k,spin,initialstate,finalstate,SumS(17)
    complex*16 ::imageunit,fourierphase
    integer :: fourierphasen,fourierphasem

    double precision :: kp1half,kp2half,kp3half,kp4half,Pplus
    double precision :: Ememe,Ep1p1,Ep2p2,Ep3p3,Ep4p4,Ep1p2,Ep1p3,Ep1p4,Ep2p3,Ep2p4,Ep3p4,shift
    double precision :: kinetic1,kinetic2,kinetic3,kinetic4,hamitot,hamicm,hamirel
    complex*16 :: lagrangeterm
    integer :: fidelta
    double precision :: adotaq,adotbq
    ! integer :: i,j,k

    !call dimtotal(nmax2,Mj,Kt,dimtot1)

    b=B_DEFAULT
    lag=LAGRANGEMULTIPLIER
    mass1=Mu2_DEFAULT
    mass2=Md2_DEFAULT
    mass3=MS2_DEFAULT
    mass4=Mg2_DEFAULT


    call sumspin4p(Nmax2,Mj,SumS)

    nz=1
    imageunit=(0.D0,1.D0)

    do color=1,2
  
      do k=0,dimtot2-SumS(17),SumS(17)
        do Spin=1,16
          do initialstate=k+SumS(Spin)+1,k+SumS(Spin+1)

            i_nzk2((color-1)*dimtot2+initialstate)=nz

            kp1=k1dat2(initialstate)
            kp2=k2dat2(initialstate)
            kp3=k3dat2(initialstate)
            kp4=k4dat2(initialstate)
            sp1=s1dat2(initialstate)
            sp2=s2dat2(initialstate)
            sp3=s3dat2(initialstate)
            sp4=s4dat2(initialstate)
            np1=n1dat2(initialstate)
            mp1=m1dat2(initialstate)
            np2=n2dat2(initialstate)
            mp2=m2dat2(initialstate)
            np3=n3dat2(initialstate)
            mp3=m3dat2(initialstate)
            np4=n4dat2(initialstate)
            mp4=m4dat2(initialstate)

            kp1half=dble(kp1)+0.5D0
            kp2half=dble(kp2)+0.5D0
            kp3half=dble(kp3)+0.5D0
            kp4half=dble(kp4)
            Pplus  =dble(Kt)+0.5D0


            do finalstate=k+SumS(Spin)+1,k+SumS(Spin+1)


              kk1=k1dat2(finalstate)
              kk2=k2dat2(finalstate)
              kk3=k3dat2(finalstate)
              kk4=k4dat2(finalstate)
              sk1=s1dat2(finalstate)
              sk2=s2dat2(finalstate)
              sk3=s3dat2(finalstate)
              sk4=s4dat2(finalstate)
              nk1=n1dat2(finalstate)
              mk1=m1dat2(finalstate)
              nk2=n2dat2(finalstate)
              mk2=m2dat2(finalstate)
              nk3=n3dat2(finalstate)
              mk3=m3dat2(finalstate)
              nk4=n4dat2(finalstate)
              mk4=m4dat2(finalstate)


        
              fourierphasen = -np1-np2-np3-np4+nk1+nk2+nk3+nk4
              Fourierphasem = -abs(mp1)-abs(mp2)-abs(mp3)-abs(mp4)+abs(mk1)+abs(mk2)+abs(mk3)+abs(mk4)
              fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem      

              Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep3p3=adotaq(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep4p4=adotaq(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

              Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep1p3=adotbq(np1,mp1,nk1,mk1,np3,mp3,nk3,mk3,sp1,sk1,sp3,sk3,kp1,kk1,kp3,kk3)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep1p4=adotbq(np1,mp1,nk1,mk1,np4,mp4,nk4,mk4,sp1,sk1,sp4,sk4,kp1,kk1,kp4,kk4)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

              Ep2p3=adotbq(np3,mp3,nk3,mk3,np2,mp2,nk2,mk2,sp3,sk3,sp2,sk2,kp3,kk3,kp2,kk2)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep2p4=adotbq(np2,mp2,nk2,mk2,np4,mp4,nk4,mk4,sp2,sk2,sp4,sk4,kp2,kk2,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

              Ep3p4=adotbq(np3,mp3,nk3,mk3,np4,mp4,nk4,mk4,sp3,sk3,sp4,sk4,kp3,kk3,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)

              shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              kinetic1=(mass1**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
              kinetic2=(mass2**2*Ememe+b**2*(kp2half/Pplus)*Ep2p2)/kp2half
              kinetic3=(mass3**2*Ememe+b**2*(kp3half/Pplus)*Ep3p3)/kp3half
              kinetic4=(mass4**2*Ememe+b**2*(kp4half/Pplus)*Ep4p4)/kp4half
            
              hamitot=kinetic1+kinetic2+kinetic3+kinetic4

              hamicm =b**2*(Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus &
                & +2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)&
                & +2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2))/Pplus

              hamirel=hamitot-hamicm

              lagrangeterm=lag*b**2*((Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)&
                & +2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2))&
                & +kappa4*(Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)&
                & +2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2))*fourierphase-shift)/Pplus


              if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

                j_nzk2(nz)=(color-1)*dimtot2+finalstate
                hamiltoniankineticvalue2(nz)=hamirel+Real(lagrangeterm)
                nz=nz+1

              endif
            enddo
          enddo
        enddo
      enddo
    enddo

    ! do i=1,dimtot2
    !   i_nzk2(dimtot2+i)=nz-1+i_nzk2(i)
    !   do j=i_nzk2(i),i_nzk2(i+1)-1
    !     j_nzk2(nz-1+j)=dimtot2+j_nzk2(j)
    !     hamiltoniankineticvalue2(nz-1+j)=hamiltoniankineticvalue2(j)
    !   enddo
    ! enddo

    i_nzk2(2*dimtot2+1)=nz

    Print*, "kinetichami_4p"
    !do test=1,nz
    !  Print*, hamiltoniankineticvalue(test), i_nzk(test),j_nzk(test),test
    !enddo

    return
  end subroutine hamiltoniankinetic4p

subroutine hamiltonianinteractqtoqg(nmax2,Mj,Kt,i_nz_vcvqtqg,j_nz_vcvqtqg,hamiltonianinteractvcvqtqg)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1(3),mass2(4)
    integer,dimension(*) :: i_nz_vcvqtqg,j_nz_vcvqtqg
    double precision,dimension(*) :: hamiltonianinteractvcvqtqg

    integer :: nz,initialstate,finalstate,color
    double precision :: Pplus,CFvalue
    integer,external :: fidelta
    double precision :: V11,V22,V33,couplings,hamiinteract
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,kk4,sk1,sk2,sk3,sk4,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4
    integer :: testvalue,test

    b=B_DEFAULT
    ! mass1(1)=Mu1_DEFAULT
    ! mass1(2)=Md1_DEFAULT
    ! mass1(3)=MS1_DEFAULT
    ! mass2(1)=Mu2_DEFAULT
    ! mass2(2)=Md2_DEFAULT
    ! mass2(3)=MS2_DEFAULT
    ! mass2(4)=Mg2_DEFAULT
    mass1(1)=Mu4_DEFAULT
    mass1(2)=Md4_DEFAULT
    mass1(3)=MS4_DEFAULT
    mass2(1)=Mu4_DEFAULT
    mass2(2)=Md4_DEFAULT
    mass2(3)=MS4_DEFAULT
    mass2(4)=Mg2_DEFAULT

    couplings=coupling3
    nz=1
    Pplus=dble(Kt)+0.5D0

    do initialstate=1,dimtot1
    !do initialstate=3652,3652

      i_nz_vcvqtqg(initialstate)=nz
      
      kp1=k1dat1(initialstate)
      kp2=k2dat1(initialstate)
      kp3=k3dat1(initialstate)
      sp1=s1dat1(initialstate)
      sp2=s2dat1(initialstate)
      sp3=s3dat1(initialstate)
      np1=n1dat1(initialstate)
      mp1=m1dat1(initialstate)
      np2=n2dat1(initialstate)
      mp2=m2dat1(initialstate)
      np3=n3dat1(initialstate)
      mp3=m3dat1(initialstate)

      do color=1,2

       do finalstate=1,dimtot2
      !do finalstate=303,303

          kk1=k1dat2(finalstate)
          kk2=k2dat2(finalstate)
          kk3=k3dat2(finalstate)
          kk4=k4dat2(finalstate)
          sk1=s1dat2(finalstate)
          sk2=s2dat2(finalstate)
          sk3=s3dat2(finalstate)
          sk4=s4dat2(finalstate)
          nk1=n1dat2(finalstate)
          mk1=m1dat2(finalstate)
          nk2=n2dat2(finalstate)
          mk2=m2dat2(finalstate)
          nk3=n3dat2(finalstate)
          mk3=m3dat2(finalstate)
          nk4=n4dat2(finalstate)
          mk4=m4dat2(finalstate)

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Four particle Fock Sector quark 1 and gluon #####################!!!!!!!!!!!
        
          testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactorqtqg(color,1)

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            call VCVinteractionqtqg(Nmax2,Kt,b,mass1(1),mass2(1),kp1,sp1,np1,mp1,kk1,kk4,sk1,sk4,nk1,mk1,nk4,mk4,V11)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif


        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 2 to Four particle Fock Sector quark 2 and gluon #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactorqtqg(color,2)

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            call VCVinteractionqtqg(Nmax2,Kt,b,mass1(2),mass2(2),kp2,sp2,np2,mp2,kk2,kk4,sk2,sk4,nk2,mk2,nk4,mk4,V22)
            V22=CFvalue*V22
          else
            V22=0.D0
          endif

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 3 to Four particle Fock Sector quark 3 and gluon #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)
          CFvalue=colorfactorqtqg(color,3)

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            call VCVinteractionqtqg(Nmax2,Kt,b,mass1(3),mass2(3),kp3,sp3,np3,mp3,kk3,kk4,sk3,sk4,nk3,mk3,nk4,mk4,V33)
            V33=CFvalue*V33
          else
            V33=0.D0
          endif

          hamiinteract=couplings/Pplus*(V11+V22+V33)
          !Print*,v11,v22,v33

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
            j_nz_vcvqtqg(nz)=(color-1)*dimtot2+finalstate

            hamiltonianinteractvcvqtqg(nz)=hamiinteract

            nz=nz+1

          endif
        enddo
      enddo
    enddo

    i_nz_vcvqtqg(dimtot1+1)=nz



    Print*, "Vertex_qtqg"
    ! do test=1,nz
    !  Print*,hamiltonianinteractvcvqtqg(test),i_nz_vcvqtqg(test),j_nz_vcvqtqg(test),test
    ! enddo

    return
  end subroutine hamiltonianinteractqtoqg

subroutine hamiltonianinteractgtoqqbar(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1(4),mass2(5)
    integer,dimension(*) :: i_nz_vcvgtqq,j_nz_vcvgtqq
    double precision,dimension(*) :: hamiltonianinteractvcvgtqq

    integer :: nz,initialstate,finalstate,color1,color2
    double precision :: Pplus,CFvalue
    integer,external :: fidelta
    double precision :: V11,couplings,hamiinteract
    integer :: kp1,kp2,kp3,kp4,sp1,sp2,sp3,sp4,np1,np2,np3,mp1,mp2,mp3,np4,mp4
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4,nk5,mk5
    integer :: testvalue,test

    b=B_DEFAULT
    mass1(1)=Mu2_DEFAULT
    mass1(2)=Md2_DEFAULT
    mass1(3)=MS2_DEFAULT
    mass1(4)=Mg2_DEFAULT
    mass2(1)=Mu3_DEFAULT
    mass2(2)=Md3_DEFAULT
    mass2(3)=MS3_DEFAULT
    mass2(4)=Md3_DEFAULT
    mass2(5)=MS3_DEFAULT

    couplings=coupling5
    nz=1
    Pplus=dble(Kt)+0.5D0

    do color1=1,2
      do initialstate=1,dimtot2
    !do initialstate=21,21

        i_nz_vcvgtqq((color1-1)*dimtot2+initialstate)=nz
      
        kp1=k1dat2(initialstate)
        kp2=k2dat2(initialstate)
        kp3=k3dat2(initialstate)
        kp4=k4dat2(initialstate)
        sp1=s1dat2(initialstate)
        sp2=s2dat2(initialstate)
        sp3=s3dat2(initialstate)
        sp4=s4dat2(initialstate)
        np1=n1dat2(initialstate)
        mp1=m1dat2(initialstate)
        np2=n2dat2(initialstate)
        mp2=m2dat2(initialstate)
        np3=n3dat2(initialstate)
        mp3=m3dat2(initialstate)
        np4=n4dat2(initialstate)
        mp4=m4dat2(initialstate)

        do color2=1,3
          do finalstate=1,dimtot3
      !do finalstate=4,4

            kk1=k1dat3(finalstate)
            kk2=k2dat3(finalstate)
            kk3=k3dat3(finalstate)
            kk4=k4dat3(finalstate)
            kk5=k5dat3(finalstate)
            sk1=s1dat3(finalstate)
            sk2=s2dat3(finalstate)
            sk3=s3dat3(finalstate)
            sk4=s4dat3(finalstate)
            sk5=s5dat3(finalstate)
            nk1=n1dat3(finalstate)
            mk1=m1dat3(finalstate)
            nk2=n2dat3(finalstate)
            mk2=m2dat3(finalstate)
            nk3=n3dat3(finalstate)
            mk3=m3dat3(finalstate)
            nk4=n4dat3(finalstate)
            mk4=m4dat3(finalstate)
            nk5=n5dat3(finalstate)
            mk5=m5dat3(finalstate)

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Four particle Fock Sector quark 1 and gluon #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
              & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
            CFvalue=colorfactorgtqq(color1,color2)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              call VCVinteractiongtqq(Nmax2,Kt,b,mass2(4),mass2(5),kp4,sp4,np4,mp4,kk4,kk5,sk4,sk5,nk4,mk4,nk5,mk5,V11)
              V11=V11*CFvalue
            else
              V11=0.D0
            endif

            hamiinteract=couplings/Pplus*(V11)

            if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
              j_nz_vcvgtqq(nz)=(color2-1)*dimtot3+finalstate

              hamiltonianinteractvcvgtqq(nz)=hamiinteract

              nz=nz+1

            endif
          enddo
        enddo
      enddo
    enddo

    i_nz_vcvgtqq(2*dimtot2+1)=nz



    Print*, "Vertex_gtqqbar"
    ! do test=1,nz
    !  Print*,hamiltonianinteractvcvgtqq(test),i_nz_vcvgtqq(test),j_nz_vcvgtqq(test),test
    ! enddo

    return
  end subroutine hamiltonianinteractgtoqqbar

!!!!!!##### the third part is the hamiltonian for five particle Fock sector (qqqqqbar) which include the kinetic term, transverse confining potential, longitudinal confining potential, OGE interaction, q->qqqbar interaction, instaneous interaction

subroutine hamiltoniankinetic5p(nmax2,Mj,Kt,i_nzk3,j_nzk3,hamiltoniankineticvalue3)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,lag,mass1,mass2,mass3,mass4,mass5
    integer,dimension(*) :: i_nzk3,j_nzk3
    double precision, dimension(*) :: hamiltoniankineticvalue3

    integer :: nz,color
    integer :: kp1,kp2,kp3,kp4,kp5,sp1,sp2,sp3,sp4,sp5,np1,np2,np3,np4,np5,mp1,mp2,mp3,mp4,mp5
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,nk4,nk5,mk1,mk2,mk3,mk4,mk5
    integer :: k,spin,initialstate,finalstate,SumS(33)
    complex*16 ::imageunit,fourierphase
    integer :: fourierphasen,fourierphasem

    double precision :: kp1half,kp2half,kp3half,kp4half,kp5half,Pplus
    double precision :: Ememe,Ep1p1,Ep3p3,Ep2p2,Ep4p4,Ep5p5,Ep1p2,Ep1p3,Ep1p4,Ep1p5,Ep2p3,Ep2p4,Ep2p5,&
      & Ep3p4,Ep3p5,Ep4p5,shift
    double precision :: kinetic1,kinetic2,kinetic3,kinetic4,kinetic5,hamitot,hamicm,hamirel
    complex*16 :: lagrangeterm
    integer :: fidelta
    double precision :: adotaq,adotbq
    integer :: test

    ! set parameter
    b=B_DEFAULT
    lag=LAGRANGEMULTIPLIER
    mass1=Mu3_DEFAULT
    mass2=Md3_DEFAULT
    mass3=MS3_DEFAULT
    mass4=Mu3_DEFAULT
    mass5=Md3_DEFAULT

    call sumspin5p(Nmax2,Mj,SumS)

    nz=1
    imageunit=(0.D0,1.D0)

    do color=1,3
    ! color=1
  
      do k=0,dimtot3-SumS(33),SumS(33)
        do Spin=1,32
          do initialstate=k+SumS(Spin)+1,k+SumS(Spin+1)
          ! do initialstate=101,101

            i_nzk3((color-1)*dimtot3+initialstate)=nz

            kp1=k1dat3(initialstate)
            kp2=k2dat3(initialstate)
            kp3=k3dat3(initialstate)
            kp4=k4dat3(initialstate)
            kp5=k5dat3(initialstate)
            sp1=s1dat3(initialstate)
            sp2=s2dat3(initialstate)
            sp3=s3dat3(initialstate)
            sp4=s4dat3(initialstate)
            sp5=s5dat3(initialstate)
            np1=n1dat3(initialstate)
            mp1=m1dat3(initialstate)
            np2=n2dat3(initialstate)
            mp2=m2dat3(initialstate)
            np3=n3dat3(initialstate)
            mp3=m3dat3(initialstate)
            np4=n4dat3(initialstate)
            mp4=m4dat3(initialstate)
            np5=n5dat3(initialstate)
            mp5=m5dat3(initialstate)

            kp1half=dble(kp1)+0.5D0
            kp2half=dble(kp2)+0.5D0
            kp3half=dble(kp3)+0.5D0
            kp4half=dble(kp4)+0.5D0
            kp5half=dble(kp5)+0.5D0
            Pplus  =dble(Kt)+0.5D0

            do finalstate=k+SumS(Spin)+1,k+SumS(Spin+1)
            !do finalstate=1,dimtot3


              kk1=k1dat3(finalstate)
              kk2=k2dat3(finalstate)
              kk3=k3dat3(finalstate)
              kk4=k4dat3(finalstate)
              kk5=k5dat3(finalstate)
              sk1=s1dat3(finalstate)
              sk2=s2dat3(finalstate)
              sk3=s3dat3(finalstate)
              sk4=s4dat3(finalstate)
              sk5=s5dat3(finalstate)
              nk1=n1dat3(finalstate)
              mk1=m1dat3(finalstate)
              nk2=n2dat3(finalstate)
              mk2=m2dat3(finalstate)
              nk3=n3dat3(finalstate)
              mk3=m3dat3(finalstate)
              nk4=n4dat3(finalstate)
              mk4=m4dat3(finalstate)
              nk5=n5dat3(finalstate)
              mk5=m5dat3(finalstate)

              fourierphasen = -np1-np2-np3-np4-np5+nk1+nk2+nk3+nk4+nk5
              Fourierphasem = -abs(mp1)-abs(mp2)-abs(mp3)-abs(mp4)-abs(mp5)+abs(mk1)+abs(mk2)+&
                              &abs(mk3)+abs(mk4)+abs(mk5)
              fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem      

              Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep3p3=adotaq(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep4p4=adotaq(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep5p5=adotaq(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p3=adotbq(np1,mp1,nk1,mk1,np3,mp3,nk3,mk3,sp1,sk1,sp3,sk3,kp1,kk1,kp3,kk3)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p4=adotbq(np1,mp1,nk1,mk1,np4,mp4,nk4,mk4,sp1,sk1,sp4,sk4,kp1,kk1,kp4,kk4)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p5=adotbq(np1,mp1,nk1,mk1,np5,mp5,nk5,mk5,sp1,sk1,sp5,sk5,kp1,kk1,kp5,kk5)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep2p3=adotbq(np2,mp2,nk2,mk2,np3,mp3,nk3,mk3,sp2,sk2,sp3,sk3,kp2,kk2,kp3,kk3)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p4=adotbq(np2,mp2,nk2,mk2,np4,mp4,nk4,mk4,sp2,sk2,sp4,sk4,kp2,kk2,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p5=adotbq(np2,mp2,nk2,mk2,np5,mp5,nk5,mk5,sp2,sk2,sp5,sk5,kp2,kk2,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep3p4=adotbq(np3,mp3,nk3,mk3,np4,mp4,nk4,mk4,sp3,sk3,sp4,sk4,kp3,kk3,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep3p5=adotbq(np3,mp3,nk3,mk3,np5,mp5,nk5,mk5,sp3,sk3,sp5,sk5,kp3,kk3,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep4p5=adotbq(np4,mp4,nk4,mk4,np5,mp5,nk5,mk5,sp4,sk4,sp5,sk5,kp4,kk4,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

              shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              kinetic1=(mass1**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
              kinetic2=(mass2**2*Ememe+b**2*(kp2half/Pplus)*Ep2p2)/kp2half
              kinetic3=(mass3**2*Ememe+b**2*(kp3half/Pplus)*Ep3p3)/kp3half
              kinetic4=(mass4**2*Ememe+b**2*(kp4half/Pplus)*Ep4p4)/kp4half
              kinetic5=(mass5**2*Ememe+b**2*(kp5half/Pplus)*Ep5p5)/kp5half

              hamitot=kinetic1+kinetic2+kinetic3+kinetic4+kinetic5

              hamicm=b**2*(Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +Ep5p5*kp5half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep1p5*sqrt(kp1half*kp5half/Pplus**2) &
                & +2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)+2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2) &
                & +2*Ep2p5*sqrt(kp2half*kp5half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2) &
                & +2*Ep3p5*sqrt(kp3half*kp5half/Pplus**2)+2*Ep4p5*sqrt(kp4half*kp5half/Pplus**2))/Pplus

              hamirel=hamitot-hamicm

              lagrangeterm=lag*b**2*((Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +Ep5p5*kp5half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep1p5*sqrt(kp1half*kp5half/Pplus**2) &
                & +2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)+2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2) &
                & +2*Ep2p5*sqrt(kp2half*kp5half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2) &
                & +2*Ep3p5*sqrt(kp3half*kp5half/Pplus**2)+2*Ep4p5*sqrt(kp4half*kp5half/Pplus**2))&
                & +kappa4*(Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +Ep5p5*kp5half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep1p5*sqrt(kp1half*kp5half/Pplus**2) &
                & +2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)+2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2) &
                & +2*Ep2p5*sqrt(kp2half*kp5half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2) &
                & +2*Ep3p5*sqrt(kp3half*kp5half/Pplus**2)+2*Ep4p5*sqrt(kp4half*kp5half/Pplus**2))*fourierphase-shift)/Pplus

            


              if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

                j_nzk3(nz)=(color-1)*dimtot3+finalstate
                hamiltoniankineticvalue3(nz)=hamirel+Real(lagrangeterm)
                nz=nz+1

              endif
            enddo
          enddo
        enddo
      enddo
    enddo

    i_nzk3(3*dimtot3+1)=nz
    ! i_nzk3(dimtot3+1)=nz

    Print*, "kinetichami_5p"
    ! do test=1,nz
    !   Print*, hamiltoniankineticvalue3(test), i_nzk3(test),j_nzk3(test),test
    ! enddo

    return
  end subroutine hamiltoniankinetic5p

  subroutine hamiltonianinteractsw5p(nmax2,Mj,Kt,i_nztsw3,j_nztsw3,hamiltonianinteractswvalue3)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b
    integer, dimension(*) :: i_nztsw3,j_nztsw3
    double precision, dimension(*) :: hamiltonianinteractswvalue3

    integer :: nz,color
    double precision :: hamisingleparticle
    complex*16 :: hamiCM,hamiinteract,hamicancel
    integer :: kp1,kp2,kp3,kp4,kp5,sp1,sp2,sp3,sp4,sp5,np1,np2,np3,np4,np5,mp1,mp2,mp3,mp4,mp5
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,nk4,nk5,mk1,mk2,mk3,mk4,mk5
    integer :: k,spin,initialstate,finalstate,SumS(33)
    
    double precision :: Es1s1,Es2s2,Es3s3,Es4s4,Es5s5
    double precision :: Ememe,Ep1p1,Ep3p3,Ep2p2,Ep4p4,Ep5p5,Ep1p2,Ep1p3,Ep1p4,Ep1p5,Ep2p3,Ep2p4,Ep2p5,&
      & Ep3p4,Ep3p5,Ep4p5,shifts
    double precision ::adotaq,adotbq,adotas
    integer ::fidelta,fourierphasen,fourierphasem
    complex*16 :: imageunit=(0.D0,1.D0),fourierphase
    double precision :: kp1half,kp2half,kp3half,kp4half,kp5half,Pplus
    integer :: test


    nz=1
    b=B_DEFAULT
    hamiinteract=(0.D0,0.D0)
    hamisingleparticle=(0.D0,0.D0)
    hamiCM=(0.D0,0.D0)

    call sumspin5p(Nmax2,Mj,SumS)

    do color=1,3
      do k=0,dimtot3-SumS(33),SumS(33)
        do Spin=1,32
          do initialstate=k+SumS(Spin)+1,k+SumS(Spin+1)

            i_nztsw3((color-1)*dimtot3+initialstate)=nz

            kp1=k1dat3(initialstate)
            kp2=k2dat3(initialstate)
            kp3=k3dat3(initialstate)
            kp4=k4dat3(initialstate)
            kp5=k5dat3(initialstate)
            sp1=s1dat3(initialstate)
            sp2=s2dat3(initialstate)
            sp3=s3dat3(initialstate)
            sp4=s4dat3(initialstate)
            sp5=s5dat3(initialstate)
            np1=n1dat3(initialstate)
            mp1=m1dat3(initialstate)
            np2=n2dat3(initialstate)
            mp2=m2dat3(initialstate)
            np3=n3dat3(initialstate)
            mp3=m3dat3(initialstate)
            np4=n4dat3(initialstate)
            mp4=m4dat3(initialstate)
            np5=n5dat3(initialstate)
            mp5=m5dat3(initialstate)

            kp1half=dble(kp1)+0.5D0
            kp2half=dble(kp2)+0.5D0
            kp3half=dble(kp3)+0.5D0
            kp4half=dble(kp4)+0.5D0
            kp5half=dble(kp5)+0.5D0
            Pplus  =dble(kt)+0.5D0

            do finalstate=k+SumS(Spin)+1,k+SumS(Spin+1)


              kk1=k1dat3(finalstate)
              kk2=k2dat3(finalstate)
              kk3=k3dat3(finalstate)
              kk4=k4dat3(finalstate)
              kk5=k5dat3(finalstate)
              sk1=s1dat3(finalstate)
              sk2=s2dat3(finalstate)
              sk3=s3dat3(finalstate)
              sk4=s4dat3(finalstate)
              sk5=s5dat3(finalstate)
              nk1=n1dat3(finalstate)
              mk1=m1dat3(finalstate)
              nk2=n2dat3(finalstate)
              mk2=m2dat3(finalstate)
              nk3=n3dat3(finalstate)
              mk3=m3dat3(finalstate)
              nk4=n4dat3(finalstate)
              mk4=m4dat3(finalstate)
              nk5=n5dat3(finalstate)
              mk5=m5dat3(finalstate)

              fourierphasen = -np1-np2-np3-np4-np5+nk1+nk2+nk3+nk4+nk5
              Fourierphasem = -abs(mp1)-abs(mp2)-abs(mp3)-abs(mp4)-abs(mp5)+abs(mk1)+abs(mk2)+&
                              &abs(mk3)+abs(mk4)+abs(mk5)
              fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem



              Es1s1=adotas(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Es2s2=adotas(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Es3s3=adotas(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Es4s4=adotas(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Es5s5=adotas(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep3p3=adotaq(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep4p4=adotaq(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep5p5=adotaq(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p3=adotbq(np1,mp1,nk1,mk1,np3,mp3,nk3,mk3,sp1,sk1,sp3,sk3,kp1,kk1,kp3,kk3)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p4=adotbq(np1,mp1,nk1,mk1,np4,mp4,nk4,mk4,sp1,sk1,sp4,sk4,kp1,kk1,kp4,kk4)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep1p5=adotbq(np1,mp1,nk1,mk1,np5,mp5,nk5,mk5,sp1,sk1,sp5,sk5,kp1,kk1,kp5,kk5)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep2p3=adotbq(np3,mp3,nk3,mk3,np2,mp2,nk2,mk2,sp3,sk3,sp2,sk2,kp3,kk3,kp2,kk2)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p4=adotbq(np2,mp2,nk2,mk2,np4,mp4,nk4,mk4,sp2,sk2,sp4,sk4,kp2,kk2,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep2p5=adotbq(np2,mp2,nk2,mk2,np5,mp5,nk5,mk5,sp2,sk2,sp5,sk5,kp2,kk2,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep3p4=adotbq(np3,mp3,nk3,mk3,np4,mp4,nk4,mk4,sp3,sk3,sp4,sk4,kp3,kk3,kp4,kk4)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)

              Ep3p5=adotbq(np3,mp3,nk3,mk3,np5,mp5,nk5,mk5,sp3,sk3,sp5,sk5,kp3,kk3,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)

              Ep4p5=adotbq(np4,mp4,nk4,mk4,np5,mp5,nk5,mk5,sp4,sk4,sp5,sk5,kp4,kk4,kp5,kk5)&
                & *fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

              shifts=2*(Sqrt(kappat3)*b**2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
                & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/Pplus

              hamisingleparticle=kappat3*b**2*(Es1s1+Es2s2+Es3s3+Es4s4+Es5s5)/Pplus

              hamiCM=kappat3*b**2*((Ep1p1*kp1half/Pplus+Ep2p2*kp2half/Pplus+Ep3p3*kp3half/Pplus+Ep4p4*kp4half/Pplus  &
                & +Ep5p5*kp5half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+2*Ep1p3*sqrt(kp1half*kp3half/Pplus**2)&
                & +2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2)+2*Ep1p5*sqrt(kp1half*kp5half/Pplus**2) &
                & +2*Ep2p3*sqrt(kp2half*kp3half/Pplus**2)+2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2) &
                & +2*Ep2p5*sqrt(kp2half*kp5half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2) &
                & +2*Ep3p5*sqrt(kp3half*kp5half/Pplus**2)+2*Ep4p5*sqrt(kp4half*kp5half/Pplus**2))*fourierphase)/Pplus

              If(color.eq.1) then
                hamicancel=kappat3*b**2*((Ep1p1*(kp4half+kp5half)/Pplus+Ep2p2*(kp4half+kp5half)/Pplus &
                  & +Ep3p3*(kp4half+kp5half)/Pplus+Ep4p4*(kp1half+kp2half+kp3half)/Pplus &
                  & +Ep5p5*(kp1half+kp2half+kp3half)/Pplus+2*Ep1p4*sqrt(kp1half*kp4half/Pplus**2) &
                  & +2*Ep1p5*sqrt(kp1half*kp5half/Pplus**2)+2*Ep2p4*sqrt(kp2half*kp4half/Pplus**2) &
                  & +2*Ep2p5*sqrt(kp2half*kp5half/Pplus**2)+2*Ep3p4*sqrt(kp3half*kp4half/Pplus**2) &
                  & +2*Ep3p5*sqrt(kp3half*kp5half/Pplus**2))*fourierphase)/Pplus
              else
                hamicancel=0.D0
              endif

              hamiinteract=hamisingleparticle-hamiCM-hamicancel-shifts
            

              if (Abs(Real(hamiinteract)).gt.epsilon(Real(hamiinteract))) then
                if(Abs(Aimag(hamiinteract)).lt.epsilon(Aimag(hamiinteract))) then

                  j_nztsw3(nz)=(color-1)*dimtot3+finalstate

                  hamiltonianinteractswvalue3(nz)=Real(hamiinteract)
  
                  nz=nz+1
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
    i_nztsw3(3*dimtot3+1)=nz

    Print*, "transhami_5p"
    ! do test=1,nz
    !   Print*, hamiltonianinteractswvalue3(test),i_nztsw3(test),j_nztsw3(test)
    ! enddo


    return
  end subroutine hamiltonianinteractsw5p

  subroutine hamiltonianinteractlongi5p(nmax2,Mj,Kt,i_nzlsw3,j_nzlsw3,hamiltonianinteractlongivalue3)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1,mass2,mass3,mass4,mass5
    integer, dimension(*) :: i_nzlsw3,j_nzlsw3
    double precision, dimension(*) :: hamiltonianinteractlongivalue3

    integer :: nz,fidelta,KroneckerDelta,color
    integer :: kp1,kp2,kp3,kp4,kp5,sp1,sp2,sp3,sp4,sp5,np1,np2,np3,mp1,mp2,mp3,np4,mp4,np5,mp5
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4,nk5,mk5
    double precision :: kp1half,kp2half,kp3half,kp4half,kp5half,Pplus,kk1half,kk2half,kk3half,&
      & kk4half,kk5half
    double precision :: V12,V13,V14,V15,V23,V24,V25,V34,V35,V45,hamiinteract
    double precision :: longipartial5_2,longitrans,longipartial_2,shift
    integer :: initialstate,finalstate
    ! integer :: test


    b=B_DEFAULT
    mass1=Mu3_DEFAULT
    mass2=Md3_DEFAULT
    mass3=MS3_DEFAULT
    mass4=Mu3_DEFAULT
    mass5=Md3_DEFAULT

    nz=1

    do color=1,3
      do initialstate=1,dimtot3
    !do initialstate=511,511

        i_nzlsw3((color-1)*dimtot3+initialstate)=nz

        kp1=k1dat3(initialstate)
        kp2=k2dat3(initialstate)
        kp3=k3dat3(initialstate)
        kp4=k4dat3(initialstate)
        kp5=k5dat3(initialstate)
        sp1=s1dat3(initialstate)
        sp2=s2dat3(initialstate)
        sp3=s3dat3(initialstate)
        sp4=s4dat3(initialstate)
        sp5=s5dat3(initialstate)
        np1=n1dat3(initialstate)
        mp1=m1dat3(initialstate)
        np2=n2dat3(initialstate)
        mp2=m2dat3(initialstate)
        np3=n3dat3(initialstate)
        mp3=m3dat3(initialstate)
        np4=n4dat3(initialstate)
        mp4=m4dat3(initialstate)
        np5=n5dat3(initialstate)
        mp5=m5dat3(initialstate)

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        kp4half=dble(kp4)+0.5D0
        kp5half=dble(kp5)+0.5D0
        Pplus  =kp1half+kp2half+kp3half+kp4half+kp5half

        do finalstate=1,dimtot3
      !do finalstate=511,511

          kk1=k1dat3(finalstate)
          kk2=k2dat3(finalstate)
          kk3=k3dat3(finalstate)
          kk4=k4dat3(finalstate)
          kk5=k5dat3(finalstate)
          sk1=s1dat3(finalstate)
          sk2=s2dat3(finalstate)
          sk3=s3dat3(finalstate)
          sk4=s4dat3(finalstate)
          sk5=s5dat3(finalstate)
          nk1=n1dat3(finalstate)
          mk1=m1dat3(finalstate)
          nk2=n2dat3(finalstate)
          mk2=m2dat3(finalstate)
          nk3=n3dat3(finalstate)
          mk3=m3dat3(finalstate)
          nk4=n4dat3(finalstate)
          mk4=m4dat3(finalstate)
          nk5=n5dat3(finalstate)
          mk5=m5dat3(finalstate)

        
          kk1half=dble(kk1)+0.5D0
          kk2half=dble(kk2)+0.5D0
          kk3half=dble(kk3)+0.5D0
          kk4half=dble(kk4)+0.5D0
          kk5half=dble(kk5)+0.5D0

          V12=longipartial_2(Kt,kp1,kp2,kk1,kk2)*longitrans(Nmax2,Pplus,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,&
            & kp1half,kp2half,kk1half,kk2half)*KroneckerDelta(mp1+mp2,mk1+mk2)*KroneckerDelta(sp1,sk1)*&
            & KroneckerDelta(sp2,sk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*&
            & fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass1+mass2)**2

          V13=longipartial_2(Kt,kp1,kp3,kk1,kk3)*longitrans(Nmax2,Pplus,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3,&
            & kp1half,kp3half,kk1half,kk3half)*KroneckerDelta(mp1+mp3,mk1+mk3)*KroneckerDelta(sp1,sk1)*&
            & KroneckerDelta(sp3,sk3)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
            & fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass1+mass3)**2

          V23=longipartial_2(Kt,kp2,kp3,kk2,kk3)*longitrans(Nmax2,Pplus,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3,&
            & kp2half,kp3half,kk2half,kk3half)*KroneckerDelta(mp2+mp3,mk2+mk3)*KroneckerDelta(sp2,sk2)*&
            & KroneckerDelta(sp3,sk3)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass2+mass3)**2

          V45=longipartial_2(Kt,kp4,kp5,kk4,kk5)*longitrans(Nmax2,Pplus,np4,mp4,np5,mp5,nk4,mk4,nk5,mk5,&
            & kp4half,kp5half,kk4half,kk5half)*KroneckerDelta(mp4+mp5,mk4+mk5)*KroneckerDelta(sp4,sk4)*&
            & KroneckerDelta(sp5,sk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)/(mass4+mass5)**2

        If(color.ne.1) then

          V14=longipartial_2(Kt,kp1,kp4,kk1,kk4)*longitrans(Nmax2,Pplus,np1,mp1,np4,mp4,nk1,mk1,nk4,mk4,&
            & kp1half,kp4half,kk1half,kk4half)*KroneckerDelta(mp1+mp4,mk1+mk4)*KroneckerDelta(sp1,sk1)*&
            & KroneckerDelta(sp4,sk4)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
            & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass1+mass4)**2

          V15=longipartial_2(Kt,kp1,kp5,kk1,kk5)*longitrans(Nmax2,Pplus,np1,mp1,np5,mp5,nk1,mk1,nk5,mk5,&
            & kp1half,kp5half,kk1half,kk5half)*KroneckerDelta(mp1+mp5,mk1+mk5)*KroneckerDelta(sp1,sk1)*&
            & KroneckerDelta(sp5,sk5)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
            & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)/(mass1+mass5)**2

          V24=longipartial_2(Kt,kp2,kp4,kk2,kk4)*longitrans(Nmax2,Pplus,np2,mp2,np4,mp4,nk2,mk2,nk4,mk4,&
            & kp2half,kp4half,kk2half,kk4half)*KroneckerDelta(mp2+mp4,mk2+mk4)*KroneckerDelta(sp2,sk2)*&
            & KroneckerDelta(sp4,sk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass2+mass4)**2

          V25=longipartial_2(Kt,kp2,kp5,kk2,kk5)*longitrans(Nmax2,Pplus,np2,mp2,np5,mp5,nk2,mk2,nk5,mk5,&
            & kp2half,kp5half,kk2half,kk5half)*KroneckerDelta(mp2+mp5,mk2+mk5)*KroneckerDelta(sp2,sk2)*&
            & KroneckerDelta(sp5,sk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)/(mass2+mass5)**2

          V34=longipartial_2(Kt,kp3,kp4,kk3,kk4)*longitrans(Nmax2,Pplus,np3,mp3,np4,mp4,nk3,mk3,nk4,mk4,&
            & kp3half,kp4half,kk3half,kk4half)*KroneckerDelta(mp3+mp4,mk3+mk4)*KroneckerDelta(sp3,sk3)*&
            & KroneckerDelta(sp4,sk4)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/(mass3+mass4)**2

          V35=longipartial_2(Kt,kp3,kp5,kk3,kk5)*longitrans(Nmax2,Pplus,np3,mp3,np5,mp5,nk3,mk3,nk5,mk5,&
            & kp3half,kp5half,kk3half,kk5half)*KroneckerDelta(mp3+mp5,mk3+mk5)*KroneckerDelta(sp3,sk3)*&
            & KroneckerDelta(sp5,sk5)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*&
            & fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)/(mass3+mass5)**2

        else
          V14=0.D0
          V15=0.D0
          V24=0.D0
          V25=0.D0
          V34=0.D0
          V35=0.D0
        endif

        shift=(Sqrt(kappal3)*b**2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
          & *fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)&
          & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)*fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)/Pplus


          hamiinteract=(V12+V13+V14+V15+V23+V24+V25+V34+V35+V45)*(0.D0-kappal3*b**4)/Pplus-shift
        ! hamiinteract=V12*(0.D0-kappal3*b**4)/Pplus

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
            j_nzlsw3(nz)=(color-1)*dimtot3+finalstate

            hamiltonianinteractlongivalue3(nz)=hamiinteract

            nz=nz+1

          endif
        enddo
      enddo
    enddo

    i_nzlsw3(3*dimtot3+1)=nz

    Print*, "longihami_5p"
    ! do test=1,nz
    !   Print*,hamiltonianinteractlongivalue3(test),i_nzlsw3(test),j_nzlsw3(test),test
    ! enddo

    return
  end subroutine hamiltonianinteractlongi5p

  subroutine hamiltonianinteractOGE5p(nmax2,Mj,Kt,i_nz_oge3,j_nz_oge3,hamiltonianinteractogevalue3)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1,mass2,mass3,mass4,mass5,massg
    integer,dimension(*) :: i_nz_oge3,j_nz_oge3
    double precision,dimension(*) :: hamiltonianinteractogevalue3

    integer :: dim1,dim2,error,nz,initialstate,finalstate,color1,color2
    double precision :: kp1half,kp2half,kp3half,kp4half,kp5half,kk1half,kk2half,kk3half,kk4half,kk5half,Pplus
    double precision :: CFvalue
    double precision,external :: OgeOneFockSector5p
    integer,external :: fidelta
    double precision :: V12,V13,V14,V15,V23,V24,V25,V34,V35,V45,couplings,hamiinteract
    double precision :: mg!,mu,md,ms
    integer :: kp1,kp2,kp3,kp4,kp5,sp1,sp2,sp3,sp4,sp5,np1,np2,np3,mp1,mp2,mp3,np4,mp4,np5,mp5
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4,nk5,mk5
    !double precision :: testintegration
    integer :: testvalue
    integer :: test

    dim1=(nmax2+1)**2*Kt**2
    dim2=Kt**2
    b=B_DEFAULT
    mass1=Mu4_DEFAULT
    mass2=Md4_DEFAULT
    mass3=MS4_DEFAULT
    mass4=Mu4_DEFAULT
    mass5=Md4_DEFAULT
    massg=mg4_default
    couplings=coupling5


    allocate(integration(dim1,dim2,10),stat=error)
    if(error.ne.0) then
        print *, 'cannot allocate basis array'
    end if

    call IEDintegration5p(nmax2,Kt,mass1,mass2,mass3,mass4,mass5,massg,b)
    !call dimtotal(nmax1,nmax2,Mj,Kt,dimtot)
    !call cpu_time(integrate_time)

    !Print*,"oge integration matrix : read finish" , integrate_time-start_time
    nz=1

    do color1=1,3
      do initialstate=1,dimtot3

        i_nz_oge3((color1-1)*dimtot3+initialstate)=nz
      
        kp1=k1dat3(initialstate)
        kp2=k2dat3(initialstate)
        kp3=k3dat3(initialstate)
        kp4=k4dat3(initialstate)
        kp5=k5dat3(initialstate)
        sp1=s1dat3(initialstate)
        sp2=s2dat3(initialstate)
        sp3=s3dat3(initialstate)
        sp4=s4dat3(initialstate)
        sp5=s5dat3(initialstate)
        np1=n1dat3(initialstate)
        mp1=m1dat3(initialstate)
        np2=n2dat3(initialstate)
        mp2=m2dat3(initialstate)
        np3=n3dat3(initialstate)
        mp3=m3dat3(initialstate)
        np4=n4dat3(initialstate)
        mp4=m4dat3(initialstate)
        np5=n5dat3(initialstate)
        mp5=m5dat3(initialstate)

        kp1half=dble(kp1)+0.5D0
        kp2half=dble(kp2)+0.5D0
        kp3half=dble(kp3)+0.5D0
        kp4half=dble(kp4)+0.5D0
        kp5half=dble(kp5)+0.5D0
        Pplus  =kp1half+kp2half+kp3half+kp4half+kp5half

        do color2=1,3
          do finalstate=1,dimtot3

            kk1=k1dat3(finalstate)
            kk2=k2dat3(finalstate)
            kk3=k3dat3(finalstate)
            kk4=k4dat3(finalstate)
            kk5=k5dat3(finalstate)
            sk1=s1dat3(finalstate)
            sk2=s2dat3(finalstate)
            sk3=s3dat3(finalstate)
            sk4=s4dat3(finalstate)
            sk5=s5dat3(finalstate)
            nk1=n1dat3(finalstate)
            mk1=m1dat3(finalstate)
            nk2=n2dat3(finalstate)
            mk2=m2dat3(finalstate)
            nk3=n3dat3(finalstate)
            mk3=m3dat3(finalstate)
            nk4=n4dat3(finalstate)
            mk4=m4dat3(finalstate)
            nk5=n5dat3(finalstate)
            mk5=m5dat3(finalstate)

            kk1half=dble(kk1)+0.5D0
            kk2half=dble(kk2)+0.5D0
            kk3half=dble(kk3)+0.5D0
            kk4half=dble(kk4)+0.5D0
            kk5half=dble(kk5)+0.5D0
            mg=massg

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 2 particle #####################!!!!!!!!!!!
            testvalue=fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,1)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V12=OgeOneFockSector5p(Nmax2,Kt,b,mass1,mass2,mg,1,kp1,kp2,kk1,kk2,sp1,sp2,sk1,sk2,np1,mp1,np2,mp2,nk1,&
                & mk1,nk2,mk2)
              V12=CFvalue*V12
            else
              V12=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 3 particle #####################!!!!!!!!!!!

            testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,2)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V13=OgeOneFockSector5p(Nmax2,Kt,b,mass1,mass3,mg,2,kp1,kp3,kk1,kk3,sp1,sp3,sk1,sk3,np1,mp1,np3,mp3,nk1,mk1,nk3,mk3)
              V13=CFvalue*V13
            else
              V13=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 4 particle #####################!!!!!!!!!!!

            testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,3)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V14=OgeOneFockSector5p(Nmax2,Kt,b,mass1,mass4,mg,3,kp1,kp4,kk1,kk4,sp1,sp4,sk1,sk4,np1,mp1,np4,mp4,nk1,mk1,nk4,mk4)
              V14=CFvalue*V14
            else
              V14=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 1 and 5 particle #####################!!!!!!!!!!!

            testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3) &
              & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)
            CFvalue=colorfactor1F(color1,color2,4)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V15=OgeOneFockSector5p(Nmax2,Kt,b,mass1,mass5,mg,4,kp1,kp5,kk1,kk5,sp1,sp5,sk1,sk5,np1,mp1,np5,mp5,nk1,mk1,nk5,mk5)
              V15=CFvalue*V15
            else
              V15=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 2 and 3 particle #####################!!!!!!!!!!!

            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,5)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V23=OgeOneFockSector5p(Nmax2,Kt,b,mass2,mass3,mg,5,kp2,kp3,kk2,kk3,sp2,sp3,sk2,sk3,np2,mp2,np3,mp3,nk2,mk2,nk3,mk3)
              V23=CFvalue*V23
            else
              V23=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 2 and 4 particle #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,6)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V24=OgeOneFockSector5p(Nmax2,Kt,b,mass2,mass4,mg,6,kp2,kp4,kk2,kk4,sp2,sp4,sk2,sk4,np2,mp2,np4,mp4,nk2,mk2,nk4,mk4)
              V24=CFvalue*V24
            else
              V24=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 2 and 5 particle #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3) &
              & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)
            CFvalue=colorfactor1F(color1,color2,7)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V25=OgeOneFockSector5p(Nmax2,Kt,b,mass2,mass5,mg,7,kp2,kp5,kk2,kk5,sp2,sp5,sk2,sk5,np2,mp2,np5,mp5,nk2,mk2,nk5,mk5)
              V25=CFvalue*V25
            else
              V25=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 3 and 4 particle #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2) &
              & *fidelta(np5,mp5,nk5,mk5,sp5,sk5,kp5,kk5)
            CFvalue=colorfactor1F(color1,color2,8)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V34=OgeOneFockSector5p(Nmax2,Kt,b,mass3,mass4,mg,8,kp3,kp4,kk3,kk4,sp3,sp4,sk3,sk4,np3,mp3,np4,mp4,nk3,mk3,nk4,mk4)
              V34=CFvalue*V34
            else
              V34=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 3 and 5 particle #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2) &
              & *fidelta(np4,mp4,nk4,mk4,sp4,sk4,kp4,kk4)
            CFvalue=colorfactor1F(color1,color2,9)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V35=OgeOneFockSector5p(Nmax2,Kt,b,mass3,mass5,mg,9,kp3,kp5,kk3,kk5,sp3,sp5,sk3,sk5,np3,mp3,np5,mp5,nk3,mk3,nk5,mk5)
              V35=CFvalue*V35
            else
              V35=0.D0
            endif

        !!!!!!!!!!########### calculate the one gluon exchange for 4 and 5 particle #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2) &
              & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
            CFvalue=colorfactor1F(color1,color2,10)
            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              V45=OgeOneFockSector5p(Nmax2,Kt,b,mass4,mass5,mg,10,kp4,kp5,kk4,kk5,sp4,sp5,sk4,sk5,np4,mp4,np5,mp5,nk4,mk4,nk5,mk5)
              V45=CFvalue*V45
            else
              V45=0.D0
            endif


            hamiinteract=couplings/Pplus**2*(-V12-V13-V14+V15-V23-V24+V25-V34+V35+V45)

            if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
              j_nz_oge3(nz)=(color2-1)*dimtot3+finalstate

              hamiltonianinteractogevalue3(nz)=hamiinteract

              nz=nz+1

            endif
          enddo
        enddo
      enddo
    enddo

    i_nz_oge3(3*dimtot3+1)=nz



    Print*, "one_gluon_exchange_5p"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue3(test),i_nz_oge3(test),j_nz_oge3(test),test
    ! enddo


    deallocate(integration)

    return
  end subroutine hamiltonianinteractOGE5p

  subroutine hamiltonianinteractOGE3to5(nmax2,Mj,Kt,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass1(3),mass2(5)
    integer,dimension(*) :: i_nz_oge3t5,j_nz_oge3t5
    double precision,dimension(*) :: hamiltonianinteractogevalue3t5

    integer :: dim1,dim2,error,nz,initialstate,finalstate,color
    double precision :: Pplus,CFvalue
    integer,external :: fidelta
    double precision :: V11,V22,V33,couplings,hamiinteract
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4,nk5,mk5
    !double precision :: testintegration
    integer :: testvalue,test

    dim1=(nmax2-4)**2
    dim2=Kt*(kt-1)*(kt+1)/6
    b=B_DEFAULT
    mass1(1)=Mu2_DEFAULT
    mass1(2)=Md2_DEFAULT
    mass1(3)=MS2_DEFAULT
    mass2(1)=Mu4_DEFAULT
    mass2(2)=Md4_DEFAULT
    mass2(3)=MS4_DEFAULT
    mass2(4)=Md4_DEFAULT
    mass2(5)=MS4_DEFAULT

    ! mass1(1)=M3t5_DEFAULT
    ! mass1(2)=M3t5_DEFAULT
    ! mass1(3)=M3t5_DEFAULT
    ! mass2(1)=M3t5_DEFAULT
    ! mass2(2)=M3t5_DEFAULT
    ! mass2(3)=M3t5_DEFAULT
    ! mass2(4)=M3t5_DEFAULT
    ! mass2(5)=M3t5_DEFAULT

    couplings=coupling3t5

    allocate(integration3t5(dim1,dim2,3),stat=error)
    if(error.ne.0) then
        print *, 'cannot allocate basis array'
    end if

    call IEDintegration3t5(nmax2,Kt,mass1,mass2,b)

    nz=1

    do initialstate=1,dimtot1
    ! do initialstate=12,12

      i_nz_oge3t5(initialstate)=nz
      
      kp1=k1dat1(initialstate)
      kp2=k2dat1(initialstate)
      kp3=k3dat1(initialstate)
      sp1=s1dat1(initialstate)
      sp2=s2dat1(initialstate)
      sp3=s3dat1(initialstate)
      np1=n1dat1(initialstate)
      mp1=m1dat1(initialstate)
      np2=n2dat1(initialstate)
      mp2=m2dat1(initialstate)
      np3=n3dat1(initialstate)
      mp3=m3dat1(initialstate)

      Pplus=Kt+0.5

      do color=1,3
        do finalstate=1,dimtot3
      ! do finalstate=263,263

          kk1=k1dat3(finalstate)
          kk2=k2dat3(finalstate)
          kk3=k3dat3(finalstate)
          kk4=k4dat3(finalstate)
          kk5=k5dat3(finalstate)
          sk1=s1dat3(finalstate)
          sk2=s2dat3(finalstate)
          sk3=s3dat3(finalstate)
          sk4=s4dat3(finalstate)
          sk5=s5dat3(finalstate)
          nk1=n1dat3(finalstate)
          mk1=m1dat3(finalstate)
          nk2=n2dat3(finalstate)
          mk2=m2dat3(finalstate)
          nk3=n3dat3(finalstate)
          mk3=m3dat3(finalstate)
          nk4=n4dat3(finalstate)
          mk4=m4dat3(finalstate)
          nk5=n5dat3(finalstate)
          mk5=m5dat3(finalstate)

        !!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Five particle Fock Sector quark 1 and sea quark #####################!!!!!!!!!!!
        
          testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactor3t5(color,1)
          V11=0.D0

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
 
            Call OgeThreeToFive(Nmax2,Kt,b,1,mass1(1),mass2(1),mass2(4),mass2(5),kp1,sp1,np1,mp1,kk1,kk4,kk5,sk1,sk4,sk5,nk1,&
              & mk1,nk4,mk4,nk5,mk5,V11)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 2 to Five particle Fock Sector quark 2 and sea quark #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactor3t5(color,2)
          V22=0.D0

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then

            Call OgeThreeToFive(Nmax2,Kt,b,2,mass1(2),mass2(2),mass2(4),mass2(5),kp2,sp2,np2,mp2,kk2,kk4,kk5,sk2,sk4,sk5,nk2,&
              & mk2,nk4,mk4,nk5,mk5,V22)
            V22=CFvalue*V22
          else
            V22=0.D0
          endif

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 2 to Five particle Fock Sector quark 2 and sea quark #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)
          CFvalue=colorfactor3t5(color,3)
          V33=0.D0

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then

            Call OgeThreeToFive(Nmax2,Kt,b,3,mass1(3),mass2(3),mass2(4),mass2(5),kp3,sp3,np3,mp3,kk3,kk4,kk5,sk3,sk4,sk5,nk3,&
              & mk3,nk4,mk4,nk5,mk5,V33)
            V33=CFvalue*V33

          else
            V33=0.D0
          endif

          hamiinteract=-couplings/Pplus**2*(V11+V22+V33)
        !hamiinteract=V11

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
            j_nz_oge3t5(nz)=(color-1)*dimtot3+finalstate

            hamiltonianinteractogevalue3t5(nz)=hamiinteract

            nz=nz+1

          endif
        enddo
      enddo
    enddo

    i_nz_oge3t5(dimtot1+1)=nz



    Print*, "one_gluon_exchange_3t5"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue3t5(test),i_nz_oge3t5(test),j_nz_oge3t5(test),test
    ! enddo


    deallocate(integration3t5)

    return
  end subroutine hamiltonianinteractOGE3to5

  subroutine hamiltonianinstantaneous3to5(nmax2,Mj,Kt,method,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt,method
    double precision :: b,binst,mass1(3),mass2(5)
    integer,dimension(*) :: i_nz_oge3t5,j_nz_oge3t5
    double precision,dimension(*) :: hamiltonianinteractogevalue3t5

    integer :: dim1,dim2,nz,initialstate,finalstate,color
    double precision :: Pplus,CFvalue
    double precision,external :: Instantaneous3t5
    integer,external :: fidelta
    double precision :: V11,V22,V33,couplings,hamiinteract
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,kk4,kk5,sk1,sk2,sk3,sk4,sk5,nk1,nk2,nk3,mk1,mk2,mk3,nk4,mk4,nk5,mk5
    !double precision :: testintegration
    integer :: testvalue,test

    dim1=(nmax2+1)**2*Kt**2
    dim2=Kt**2
    b=B_DEFAULT
    binst=binst_default
    mass1(1)=Mu1_DEFAULT
    mass1(2)=Md1_DEFAULT
    mass1(3)=MS1_DEFAULT
    mass2(1)=Mu3_DEFAULT
    mass2(2)=Md3_DEFAULT
    mass2(3)=MS3_DEFAULT
    mass2(4)=Md3_DEFAULT
    mass2(5)=MS3_DEFAULT

    If(method.eq.1) then

      couplings=coupling3t5*4*PI

    else if (method.eq.2) then

      couplings=coupling3t5**2

    end if

    nz=1
    Pplus=dble(Kt)+0.5D0

    do initialstate=1,dimtot1
    ! do initialstate=11,11

      i_nz_oge3t5(initialstate)=nz
      
      kp1=k1dat1(initialstate)
      kp2=k2dat1(initialstate)
      kp3=k3dat1(initialstate)
      sp1=s1dat1(initialstate)
      sp2=s2dat1(initialstate)
      sp3=s3dat1(initialstate)
      np1=n1dat1(initialstate)
      mp1=m1dat1(initialstate)
      np2=n2dat1(initialstate)
      mp2=m2dat1(initialstate)
      np3=n3dat1(initialstate)
      mp3=m3dat1(initialstate)

      do color=1,3
        do finalstate=1,dimtot3
      ! do finalstate=19,19

          kk1=k1dat3(finalstate)
          kk2=k2dat3(finalstate)
          kk3=k3dat3(finalstate)
          kk4=k4dat3(finalstate)
          kk5=k5dat3(finalstate)
          sk1=s1dat3(finalstate)
          sk2=s2dat3(finalstate)
          sk3=s3dat3(finalstate)
          sk4=s4dat3(finalstate)
          sk5=s5dat3(finalstate)
          nk1=n1dat3(finalstate)
          mk1=m1dat3(finalstate)
          nk2=n2dat3(finalstate)
          mk2=m2dat3(finalstate)
          nk3=n3dat3(finalstate)
          mk3=m3dat3(finalstate)
          nk4=n4dat3(finalstate)
          mk4=m4dat3(finalstate)
          nk5=n5dat3(finalstate)
          mk5=m5dat3(finalstate)

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Five particle Fock Sector quark 1 and sea quark #####################!!!!!!!!!!!
        
          testvalue=fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactor3t5(color,1)
          V11=0.D0

          If(testvalue.ne.0.and.sp1.eq.sk1.and.sk4.eq.-sk5.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            V11=Instantaneous3t5(Nmax2,Kt,b,binst,kp1,kk1,kk4,kk5,np1,mp1,nk1,mk1,nk4,mk5,nk4,mk5)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif


        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 2 to Five particle Fock Sector quark 2 and sea quark #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
          CFvalue=colorfactor3t5(color,2)
          V22=0.D0

          If(testvalue.ne.0.and.sp2.eq.sk2.and.sk4.eq.-sk5.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            V22=Instantaneous3t5(Nmax2,Kt,b,binst,kp2,kk2,kk4,kk5,np2,mp2,nk2,mk2,nk4,mk5,nk4,mk5)
            V22=CFvalue*V22
          else
            V22=0.D0
          endif

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 2 to Five particle Fock Sector quark 2 and sea quark #####################!!!!!!!!!!!

          testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)
          CFvalue=colorfactor3t5(color,3)
          V33=0.D0

          If(testvalue.ne.0.and.sp3.eq.sk3.and.sk4.eq.-sk5.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            V33=Instantaneous3t5(Nmax2,Kt,b,binst,kp3,kk3,kk4,kk5,np3,mp3,nk3,mk3,nk4,mk5,nk4,mk5)
            V33=CFvalue*V33
          else
            V33=0.D0
          endif

          hamiinteract=-couplings/Pplus*(V11+V22+V33)

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then

            j_nz_oge3t5(nz)=(color-1)*dimtot3+finalstate

            hamiltonianinteractogevalue3t5(nz)=hamiinteract

            nz=nz+1

          endif
        enddo
      enddo
    enddo

    i_nz_oge3t5(dimtot1+1)=nz



    Print*, "Instantaneous3t5"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue3t5(test),i_nz_oge3t5(test),j_nz_oge3t5(test),test
    ! enddo


    return
  end subroutine hamiltonianinstantaneous3to5
