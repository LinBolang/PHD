subroutine hamiltoniankineticSinglequark(nmax2,Mj,Kt,method,renormass,mg,i_nzk,j_nzk,hamiltoniankineticvalue)
    use numbers
    use basis_info
    use colorfactor
    implicit none

    integer :: nmax2,Mj,Kt,method
    double precision :: b,mu,md,lag,ms,renormass,mg
    integer,dimension(*) :: i_nzk,j_nzk
    double precision, dimension(*) :: hamiltoniankineticvalue

    integer :: nz,jlist(100)
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,nk2,nk3,mk1,mk2,mk3
    integer :: i,j,i1,i2,i3,i4,j1,j2,j3,j4,j5,j6
    complex*16 ::imageunit,fourierphase
    integer :: fourierphasen,fourierphasem
    integer :: sumspin(5),color,k,spin

    double precision :: kp1half,kp2half,kp3half,Pplus
    double precision :: Ememe,Ep1p1,Ep3p3,Ep2p2,Ep1p2,Ep1p3,Ep3p2,shift
    double precision :: kinetic1,kinetic2,kinetic3,hamitot,hamicm,hamirel
    complex*16 :: lagrangeterm
    integer :: fidelta
    double precision :: adotaq,adotbq
    integer :: test

    mu=Mu1_DEFAULT+renormass
    md=Md1_DEFAULT+renormass
    ms=Ms1_DEFAULT+renormass
    b =B_DEFAULT
    lag=LAGRANGEMULTIPLIER

    nz=1
    imageunit=(0.D0,1.D0)

!!!!!!!!! kinetic term for one particle systerm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    i=1
    i_nzk(i)=nz
    kp1=k1dat(i)
    kp2=k2dat(i)
    kp3=k3dat(i)
    sp1=s1dat(i)
    sp2=s2dat(i)
    sp3=s3dat(i)
    np1=n1dat(i)
    mp1=m1dat(i)
    np2=n2dat(i)
    mp2=m2dat(i)
    np3=n3dat(i)
    mp3=m3dat(i)
    
    kk1=k1dat(i)
    kk2=k2dat(i)
    kk3=k3dat(i)
    sk1=s1dat(i)
    sk2=s2dat(i)
    sk3=s3dat(i)
    nk1=n1dat(i)
    mk1=m1dat(i)
    nk2=n2dat(i)
    mk2=m2dat(i)
    nk3=n3dat(i)
    mk3=m3dat(i)

    Pplus=dble(Kt)+0.5D0
    kp1half=dble(kp1)+0.5D0
    fourierphasen = -np1+nk1
    Fourierphasem = -abs(mp1)+abs(mk1)
    fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem      

    Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
      & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
    Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
      & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
    shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
          & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

    kinetic1=(mu**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
    hamicm=b**2*(Ep1p1*kp1half/Pplus)/Pplus
    hamirel= kinetic1-hamicm
    lagrangeterm=lag*b**2*(Ep1p1*kp1half/Pplus+kappa4*(Ep1p1*kp1half/Pplus)*fourierphase-shift)/Pplus

    if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

      j_nzk(nz)=i
      hamiltoniankineticvalue(nz)=hamirel+Real(lagrangeterm)
      nz=nz+1

    endif

!!!!!!!!!!!!!!!!!!!!!!!! kinetic term for the two particle system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    call SumS2p(Nmax2,Mj,SumSpin)
    mu=Mu2_DEFAULT
    md=Md2_DEFAULT
    ms=Ms2_DEFAULT

    If(method.eq.2) then
        do k=dimtot1,dimtot1+2*dimtot2-Sumspin(5),Sumspin(5)
          do spin=1,4
            do i=k+SumSpin(spin)+1,k+SumSpin(spin+1)

              i_nzk(i)=nz
              kp1=k1dat(i)
              kp2=k2dat(i)
              kp3=k3dat(i)
              sp1=s1dat(i)
              sp2=s2dat(i)
              sp3=s3dat(i)
              np1=n1dat(i)
              mp1=m1dat(i)
              np2=n2dat(i)
              mp2=m2dat(i)
              np3=n3dat(i)
              mp3=m3dat(i)
              kp1half=dble(kp1)+0.5D0
              kp2half=dble(kp2)
              Pplus  =dble(Kt)+0.5D0

              do j=k+SumSpin(spin)+1,k+SumSpin(spin+1)

                kk1=k1dat(j)
                kk2=k2dat(j)
                kk3=k3dat(j)
                sk1=s1dat(j)
                sk2=s2dat(j)
                sk3=s3dat(j)
                nk1=n1dat(j)
                mk1=m1dat(j)
                nk2=n2dat(j)
                mk2=m2dat(j)
                nk3=n3dat(j)
                mk3=m3dat(j)

                fourierphasen = -np1-np2+nk1+nk2
                Fourierphasem = -abs(mp1)-abs(mp2)+abs(mk1)+abs(mk2)
                fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem

                Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

                kinetic1=(mu**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
                kinetic2=(mg**2*Ememe+b**2*(kp2half/Pplus)*Ep2p2)/kp2half
                hamicm  =b**2*(Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus)/Pplus
                hamirel = kinetic1+kinetic2-hamicm
                lagrangeterm=lag*b**2*((Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus) &
                  & +kappa4*(Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus)*fourierphase &
                  & -shift)/Pplus

                if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

                  j_nzk(nz)=j
                  hamiltoniankineticvalue(nz)=hamirel+Real(lagrangeterm)
                  nz=nz+1

                endif
              enddo
            enddo
          enddo
        enddo
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kinetic term for the three particle system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    mu=Mu3_DEFAULT
    md=Md3_DEFAULT
    ms=Ms3_DEFAULT

    do i=dimtot1+2*dimtot2+1,dimtot1+2*dimtot2+3*dimtot3

      i_nzk(i)=nz

      kp1=k1dat(i)
      kp2=k2dat(i)
      kp3=k3dat(i)
      sp1=s1dat(i)
      sp2=s2dat(i)
      sp3=s3dat(i)
      np1=n1dat(i)
      mp1=m1dat(i)
      np2=n2dat(i)
      mp2=m2dat(i)
      np3=n3dat(i)
      mp3=m3dat(i)

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

        kk1=k1dat(dimtot1+2*dimtot2+jlist(j))
        kk2=k2dat(dimtot1+2*dimtot2+jlist(j))
        kk3=k3dat(dimtot1+2*dimtot2+jlist(j))
        sk1=s1dat(dimtot1+2*dimtot2+jlist(j))
        sk2=s2dat(dimtot1+2*dimtot2+jlist(j))
        sk3=s3dat(dimtot1+2*dimtot2+jlist(j))
        nk1=n1dat(dimtot1+2*dimtot2+jlist(j))
        mk1=m1dat(dimtot1+2*dimtot2+jlist(j))
        nk2=n2dat(dimtot1+2*dimtot2+jlist(j))
        mk2=m2dat(dimtot1+2*dimtot2+jlist(j))
        nk3=n3dat(dimtot1+2*dimtot2+jlist(j))
        mk3=m3dat(dimtot1+2*dimtot2+jlist(j))


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
          hamiltoniankineticvalue(nz)=hamirel+Real(lagrangeterm)
          nz=nz+1

        endif
      enddo
    enddo

    i_nzk(dimtot+1)=nz
    !print*,"test",ms,mu,md

    Print*, "kinetichami_3p"
    ! do test=1,nz
    !  Print*, hamiltoniankineticvalue(test), i_nzk(test),j_nzk(test),test
    ! enddo

    return
end subroutine hamiltoniankineticSinglequark

subroutine hamiltonianinteractqtoqgSingle(nmax2,Mj,Kt,particlenumber,renormass,i_nz_vcvqtqg,j_nz_vcvqtqg,hamiltonianinteractvcvqtqg)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt,particlenumber
    double precision :: b,mass1(3),mass2(4),renormass
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
    mass1(1)=Mu1_DEFAULT+renormass
    mass1(2)=Md1_DEFAULT+renormass
    mass1(3)=MS1_DEFAULT+renormass
    mass2(1)=Mu2_DEFAULT
    mass2(2)=Md2_DEFAULT
    mass2(3)=MS2_DEFAULT
    mass2(4)=Mg2_DEFAULT

    couplings=coupling3
    nz=1
    Pplus=dble(Kt)+0.5D0

    initialstate=1

      i_nz_vcvqtqg(initialstate)=nz
      
      kp1=k1dat(initialstate)
      kp2=k2dat(initialstate)
      kp3=k3dat(initialstate)
      sp1=s1dat(initialstate)
      sp2=s2dat(initialstate)
      sp3=s3dat(initialstate)
      np1=n1dat(initialstate)
      mp1=m1dat(initialstate)
      np2=n2dat(initialstate)
      mp2=m2dat(initialstate)
      np3=n3dat(initialstate)
      mp3=m3dat(initialstate)

        do finalstate=dimtot1+1,dimtot1+2*dimtot2

          kk1=k1dat(finalstate)
          kk2=k2dat(finalstate)
          kk3=k3dat(finalstate)
          sk1=s1dat(finalstate)
          sk2=s2dat(finalstate)
          sk3=s3dat(finalstate)
          nk1=n1dat(finalstate)
          mk1=m1dat(finalstate)
          nk2=n2dat(finalstate)
          mk2=m2dat(finalstate)
          nk3=n3dat(finalstate)
          mk3=m3dat(finalstate)

          testvalue=fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

          If(finalstate.ge.dimtot1+1.and.finalstate.le.dimtot1+dimtot2) then

            CFvalue=colorfactorqtqg(1,particlenumber)

          else if (finalstate.ge.dimtot1+dimtot2+1.and.finalstate.le.dimtot1+2*dimtot2) then

            CFvalue=colorfactorqtqg(2,particlenumber)

          endif


          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            call VCVinteractionqtqg(Nmax2,Kt,b,mass1(particlenumber),mass2(particlenumber),kp1,sp1,np1,mp1,kk1&
              & ,kk4,sk1,sk4,nk1,mk1,nk4,mk4,V11)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif

          hamiinteract=couplings/Pplus*(V11)

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
            j_nz_vcvqtqg(nz)=finalstate

            hamiltonianinteractvcvqtqg(nz)=hamiinteract

            nz=nz+1

          endif
        enddo

    i_nz_vcvqtqg(2)=nz



    Print*, "Vertex_qtqg"
    ! do test=1,nz
    !  Print*,hamiltonianinteractvcvqtqg(test),i_nz_vcvqtqg(test),j_nz_vcvqtqg(test),test
    ! enddo

    return
  end subroutine hamiltonianinteractqtoqgsingle

subroutine hamiltonianinteractgtoqqSingle(nmax2,Mj,Kt,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
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
    mass2(4)=Md3_DEFAULT
    mass2(5)=MS3_DEFAULT

    couplings=coupling5
    nz=1
    Pplus=dble(Kt)+0.5D0

    do initialstate=dimtot1+1,dimtot1+2*dimtot2

      i_nz_vcvgtqq(initialstate)=nz
      
      kp1=k1dat(initialstate)
      kp2=k2dat(initialstate)
      kp3=k3dat(initialstate)
      sp1=s1dat(initialstate)
      sp2=s2dat(initialstate)
      sp3=s3dat(initialstate)
      np1=n1dat(initialstate)
      mp1=m1dat(initialstate)
      np2=n2dat(initialstate)
      mp2=m2dat(initialstate)
      np3=n3dat(initialstate)
      mp3=m3dat(initialstate)

      do finalstate=dimtot1+2*dimtot2+1,dimtot1+2*dimtot2+3*dimtot3

        kk1=k1dat(finalstate)
        kk2=k2dat(finalstate)
        kk3=k3dat(finalstate)
        sk1=s1dat(finalstate)
        sk2=s2dat(finalstate)
        sk3=s3dat(finalstate)
        nk1=n1dat(finalstate)
        mk1=m1dat(finalstate)
        nk2=n2dat(finalstate)
        mk2=m2dat(finalstate)
        nk3=n3dat(finalstate)
        mk3=m3dat(finalstate)


        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Four particle Fock Sector quark 1 and gluon #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)

              If(initialstate.ge.dimtot1+1.and.initialstate.le.dimtot1+dimtot2) then
                color1=1
                if(finalstate.ge.dimtot1+2*dimtot2+1.and.finalstate.le.dimtot1+2*dimtot2+dimtot3) then
                  color2=1
                else if(finalstate.ge.dimtot1+2*dimtot2+dimtot3+1.and.finalstate.le.dimtot1+2*dimtot2+2*dimtot3) then
                  color2=2
                else if(finalstate.ge.dimtot1+2*dimtot2+2*dimtot3+1.and.finalstate.le.dimtot1+2*dimtot2+3*dimtot3) then
                  color2=3
                endif
              else if(initialstate.ge.dimtot1+dimtot2+1.and.initialstate.le.dimtot1+2*dimtot2) then
                color1=2
                if(finalstate.ge.dimtot1+2*dimtot2+1.and.finalstate.le.dimtot1+2*dimtot2+dimtot3) then
                  color2=1
                else if(finalstate.ge.dimtot1+2*dimtot2+dimtot3+1.and.finalstate.le.dimtot1+2*dimtot2+2*dimtot3) then
                  color2=2
                else if(finalstate.ge.dimtot1+2*dimtot2+2*dimtot3+1.and.finalstate.le.dimtot1+2*dimtot2+3*dimtot3) then
                  color2=3
                endif
              endif

            CFvalue=colorfactorgtqq(color1,color2)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              call VCVinteractiongtqq(Nmax2,Kt,b,mass2(4),mass2(5),kp4,sp4,np4,mp4,kk4,kk5,sk4,sk5,nk4,mk4,nk5,mk5,V11)
              V11=V11*CFvalue
            else
              V11=0.D0
            endif

            hamiinteract=couplings/Pplus*(V11)

            if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
              j_nz_vcvgtqq(nz)=finalstate

              hamiltonianinteractvcvgtqq(nz)=hamiinteract

              nz=nz+1

            endif
      enddo
    enddo

    i_nz_vcvgtqq(2*dimtot2+1)=nz



    Print*, "Vertex_gtqqbar"
    ! do test=1,nz
    !  Print*,hamiltonianinteractvcvgtqq(test),i_nz_vcvgtqq(test),j_nz_vcvgtqq(test),test
    ! enddo

    return
end subroutine hamiltonianinteractgtoqqSingle

subroutine hamiltonianinteractOGE1to3(nmax2,Mj,Kt,particlenumber,renormass,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt,particlenumber
    double precision :: b,mass1(3),mass2(5),renormass
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
    mass1(1)=Mu2_DEFAULT+renormass
    mass1(2)=Md2_DEFAULT+renormass
    mass1(3)=MS2_DEFAULT+renormass
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
    initialstate=1

      i_nz_oge3t5(initialstate)=nz
      
      kp1=k1dat(initialstate)
      kp2=k2dat(initialstate)
      kp3=k3dat(initialstate)
      sp1=s1dat(initialstate)
      sp2=s2dat(initialstate)
      sp3=s3dat(initialstate)
      np1=n1dat(initialstate)
      mp1=m1dat(initialstate)
      np2=n2dat(initialstate)
      mp2=m2dat(initialstate)
      np3=n3dat(initialstate)
      mp3=m3dat(initialstate)

      Pplus=Kt+0.5

        do finalstate=dimtot1+1,dimtot1+3*dimtot3

          kk1=k1dat(finalstate)
          kk2=k2dat(finalstate)
          kk3=k3dat(finalstate)
          sk1=s1dat(finalstate)
          sk2=s2dat(finalstate)
          sk3=s3dat(finalstate)
          nk1=n1dat(finalstate)
          mk1=m1dat(finalstate)
          nk2=n2dat(finalstate)
          mk2=m2dat(finalstate)
          nk3=n3dat(finalstate)
          mk3=m3dat(finalstate)

          if(finalstate.ge.dimtot1+1.and.finalstate.le.dimtot1+dimtot3) then
                  color=1
          else if(finalstate.ge.dimtot1+dimtot3+1.and.finalstate.le.dimtot1+2*dimtot3) then
                  color=2
          else if(finalstate.ge.dimtot1+2*dimtot3+1.and.finalstate.le.dimtot1+3*dimtot3) then
                  color=3
          endif

        !!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Five particle Fock Sector quark 1 and sea quark #####################!!!!!!!!!!!
        
          testvalue=1
          CFvalue=colorfactor3t5(color,particlenumber)
          V11=0.D0

          If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
 
            Call OgeThreeToFive(Nmax2,Kt,b,1,mass1(particlenumber),mass2(particlenumber),mass2(4),mass2(5),&
              & kp1,sp1,np1,mp1,kk1,kk4,kk5,sk1,sk4,sk5,nk1,mk1,nk4,mk4,nk5,mk5,V11)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif

          hamiinteract=-couplings/Pplus**2*(V11)
        !hamiinteract=V11

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
            j_nz_oge3t5(nz)=finalstate

            hamiltonianinteractogevalue3t5(nz)=hamiinteract

            nz=nz+1

          endif
        enddo

    i_nz_oge3t5(dimtot1+1)=nz



    Print*, "one_gluon_exchange_3t5"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue3t5(test),i_nz_oge3t5(test),j_nz_oge3t5(test),test
    ! enddo


    deallocate(integration3t5)

    return
  end subroutine hamiltonianinteractOGE1to3

subroutine hamiltonianinstantaneous1to3(nmax2,Mj,Kt,method,particlenumber,i_nz_oge3t5,j_nz_oge3t5,hamiltonianinteractogevalue3t5)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt,method,particlenumber
    double precision :: b,binst
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

    couplings=coupling3t5

    nz=1
    Pplus=dble(Kt)+0.5D0

    initialstate=1

      i_nz_oge3t5(initialstate)=nz
      
      kp1=k1dat(initialstate)
      kp2=k2dat(initialstate)
      kp3=k3dat(initialstate)
      sp1=s1dat(initialstate)
      sp2=s2dat(initialstate)
      sp3=s3dat(initialstate)
      np1=n1dat(initialstate)
      mp1=m1dat(initialstate)
      np2=n2dat(initialstate)
      mp2=m2dat(initialstate)
      np3=n3dat(initialstate)
      mp3=m3dat(initialstate)

      If(method.eq.1) then
        dim1=dimtot1+1
        dim2=dimtot1+3*dimtot3
      else if(method.eq.2) then
        dim1=dimtot1+2*dimtot2+1
        dim2=dimtot1+2*dimtot2+3*dimtot3
      endif

        do finalstate=dim1,dim2

          kk1=k1dat(finalstate)
          kk2=k2dat(finalstate)
          kk3=k3dat(finalstate)
          sk1=s1dat(finalstate)
          sk2=s2dat(finalstate)
          sk3=s3dat(finalstate)
          nk1=n1dat(finalstate)
          mk1=m1dat(finalstate)
          nk2=n2dat(finalstate)
          mk2=m2dat(finalstate)
          nk3=n3dat(finalstate)
          mk3=m3dat(finalstate)

          if(finalstate.ge.dim1.and.finalstate.le.dim1+dimtot3-1) then
                  color=1
          else if(finalstate.ge.dim1+dimtot3.and.finalstate.le.dim1+2*dimtot3-1) then
                  color=2
          else if(finalstate.ge.dim1+2*dimtot3.and.finalstate.le.dim2) then
                  color=3
          endif

        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Five particle Fock Sector quark 1 and sea quark #####################!!!!!!!!!!!
        
          testvalue=1
          CFvalue=colorfactor3t5(color,particlenumber)
          V11=0.D0

          If(testvalue.ne.0.and.sp1.eq.sk1.and.sk4.eq.-sk5.and.abs(CFvalue).gt.epsilon(CFvalue)) then
            V11=Instantaneous3t5(Nmax2,Kt,b,binst,kp1,kk1,kk4,kk5,np1,mp1,nk1,mk1,nk4,mk5,nk4,mk5)
            V11=CFvalue*V11
          else
            V11=0.D0
          endif

          hamiinteract=-couplings/Pplus**2*(V11)

          if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then

            j_nz_oge3t5(nz)=finalstate

            hamiltonianinteractogevalue3t5(nz)=hamiinteract

            nz=nz+1

          endif
        enddo

    i_nz_oge3t5(dimtot1+1)=nz



    Print*, "Instantaneous3t5"
    ! do test=1,nz
    !  Print*,hamiltonianinteractogevalue3t5(test),i_nz_oge3t5(test),j_nz_oge3t5(test),test
    ! enddo


    return
  end subroutine hamiltonianinstantaneous1to3


subroutine hamiltoniankineticSinglephoton(nmax2,Mj,Kt,mg,i_nzk,j_nzk,hamiltoniankineticvalue)
    use numbers
    use basis_info
    use colorfactor
    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mg,lag,mu,md
    integer,dimension(*) :: i_nzk,j_nzk
    double precision, dimension(*) :: hamiltoniankineticvalue

    integer :: nz,jlist(100)
    integer :: kp1,kp2,kp3,sp1,sp2,sp3,np1,np2,np3,mp1,mp2,mp3
    integer :: kk1,kk2,kk3,sk1,sk2,sk3,nk1,nk2,nk3,mk1,mk2,mk3
    integer :: i,j,i1,i2,i3,i4,j1,j2,j3,j4,j5,j6
    complex*16 ::imageunit,fourierphase
    integer :: fourierphasen,fourierphasem
    integer :: sumspin(5),k,spin

    double precision :: kp1half,kp2half,kp3half,Pplus
    double precision :: Ememe,Ep1p1,Ep3p3,Ep2p2,Ep1p2,Ep1p3,Ep3p2,shift
    double precision :: kinetic1,kinetic2,kinetic3,hamitot,hamicm,hamirel
    complex*16 :: lagrangeterm
    integer :: fidelta
    double precision :: adotaq,adotbq
    integer :: test

    mg=mg
    b =B_DEFAULT
    lag=LAGRANGEMULTIPLIER

    nz=1
    imageunit=(0.D0,1.D0)

!!!!!!!!! kinetic term for one particle systerm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    Do i=1,1

      i_nzk(i)=nz
      kp1=k1dat(i)
      kp2=k2dat(i)
      kp3=k3dat(i)
      sp1=s1dat(i)
      sp2=s2dat(i)
      sp3=s3dat(i)
      np1=n1dat(i)
      mp1=m1dat(i)
      np2=n2dat(i)
      mp2=m2dat(i)
      np3=n3dat(i)
      mp3=m3dat(i)
      
      kk1=k1dat(i)
      kk2=k2dat(i)
      kk3=k3dat(i)
      sk1=s1dat(i)
      sk2=s2dat(i)
      sk3=s3dat(i)
      nk1=n1dat(i)
      mk1=m1dat(i)
      nk2=n2dat(i)
      mk2=m2dat(i)
      nk3=n3dat(i)
      mk3=m3dat(i)

      Pplus=dble(Kt)
      kp1half=dble(kp1)
      fourierphasen = -np1+nk1
      Fourierphasem = -abs(mp1)+abs(mk1)
      fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem      

      Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
        & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
      Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
        & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
      shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*&
            & fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

      kinetic1=(mg**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
      hamicm=b**2*(Ep1p1*kp1half/Pplus)/Pplus
      hamirel= kinetic1-hamicm
      lagrangeterm=lag*b**2*(Ep1p1*kp1half/Pplus+kappa4*(Ep1p1*kp1half/Pplus)*fourierphase-shift)/Pplus

      !if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

        j_nzk(nz)=i
        hamiltoniankineticvalue(nz)=hamirel+Real(lagrangeterm)+epsilon(Abs(hamirel+Real(lagrangeterm)))
        nz=nz+1
      !endif
    enddo

!!!!!!!!!!!!!!!!!!!!!!!! kinetic term for the two particle system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    call sums2q(Nmax2,Mj,SumSpin)
    mu=Mu3_DEFAULT
    md=Md3_DEFAULT

        !do k=dimtot1,dimtot1+3*dimtot2-Sumspin(5),Sumspin(5)
        do k=dimtot1,dimtot1+dimtot2-Sumspin(5),Sumspin(5)
          do spin=1,4
            do i=k+SumSpin(spin)+1,k+SumSpin(spin+1)

              i_nzk(i)=nz
              kp1=k1dat(i)
              kp2=k2dat(i)
              kp3=k3dat(i)
              sp1=s1dat(i)
              sp2=s2dat(i)
              sp3=s3dat(i)
              np1=n1dat(i)
              mp1=m1dat(i)
              np2=n2dat(i)
              mp2=m2dat(i)
              np3=n3dat(i)
              mp3=m3dat(i)
              kp1half=dble(kp1)+0.5D0
              kp2half=dble(kp2)+0.5D0
              Pplus  =dble(Kt)

              do j=k+SumSpin(spin)+1,k+SumSpin(spin+1)

                kk1=k1dat(j)
                kk2=k2dat(j)
                kk3=k3dat(j)
                sk1=s1dat(j)
                sk2=s2dat(j)
                sk3=s3dat(j)
                nk1=n1dat(j)
                mk1=m1dat(j)
                nk2=n2dat(j)
                mk2=m2dat(j)
                nk3=n3dat(j)
                mk3=m3dat(j)

                fourierphasen = -np1-np2+nk1+nk2
                Fourierphasem = -abs(mp1)-abs(mp2)+abs(mk1)+abs(mk2)
                fourierphase  = (-1.D0)**fourierphasen*imageunit**fourierphasem

                Ememe=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep1p1=adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep2p2=adotaq(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                Ep1p2=adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)
                shift=2*fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)*fidelta(np2,mp2,nk2,mk2,sp2,sk2,kp2,kk2)&
                  & *fidelta(np3,mp3,nk3,mk3,sp3,sk3,kp3,kk3)

                kinetic1=(mu**2*Ememe+b**2*(kp1half/Pplus)*Ep1p1)/kp1half
                kinetic2=(md**2*Ememe+b**2*(kp2half/Pplus)*Ep2p2)/kp2half
                hamicm  =b**2*(Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus)/Pplus
                hamirel = kinetic1+kinetic2-hamicm
                lagrangeterm=lag*b**2*((Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus) &
                  & +kappa4*(Ep1p1*kp1half/Pplus+2*Ep1p2*sqrt(kp1half*kp2half/Pplus**2)+Ep2p2*kp2half/Pplus)*fourierphase &
                  & -shift)/Pplus

                if (abs(hamirel+Real(lagrangeterm)).gt.epsilon(Abs(hamirel+Real(lagrangeterm)))) then

                  j_nzk(nz)=j
                  hamiltoniankineticvalue(nz)=hamirel+Real(lagrangeterm)
                  nz=nz+1

                endif
              enddo
            enddo
          enddo
        enddo

    i_nzk(dimtot+1)=nz
    !print*,"test",ms,mu,md

    Print*, "kinetichami_3p"
    ! do test=1,nz
    !  Print*, hamiltoniankineticvalue(test), i_nzk(test),j_nzk(test),test
    ! enddo

    return
end subroutine hamiltoniankineticSinglephoton

subroutine hamiltonianinteractgtoqqphoton(nmax2,Mj,Kt,color1,i_nz_vcvgtqq,j_nz_vcvgtqq,hamiltonianinteractvcvgtqq)
    use numbers
    use basis_info
    use colorfactor

    implicit none

    integer :: nmax2,Mj,Kt
    double precision :: b,mass2(5)
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
    mass2(4)=Mu3_DEFAULT
    mass2(5)=Md3_DEFAULT

    couplings=coupling5
    nz=1
    Pplus=dble(Kt)

    initialstate=1
      i_nz_vcvgtqq(initialstate)=nz
      
      kp1=k1dat(initialstate)
      kp2=k2dat(initialstate)
      kp3=k3dat(initialstate)
      sp1=s1dat(initialstate)
      sp2=s2dat(initialstate)
      sp3=s3dat(initialstate)
      np1=n1dat(initialstate)
      mp1=m1dat(initialstate)
      np2=n2dat(initialstate)
      mp2=m2dat(initialstate)
      np3=n3dat(initialstate)
      mp3=m3dat(initialstate)

      do finalstate=dimtot1+1,dimtot1+3*dimtot2
      !do finalstate=4,4

        kk1=k1dat(finalstate)
        kk2=k2dat(finalstate)
        kk3=k3dat(finalstate)
        sk1=s1dat(finalstate)
        sk2=s2dat(finalstate)
        sk3=s3dat(finalstate)
        nk1=n1dat(finalstate)
        mk1=m1dat(finalstate)
        nk2=n2dat(finalstate)
        mk2=m2dat(finalstate)
        nk3=n3dat(finalstate)
        mk3=m3dat(finalstate)


        !!!!!!!!!!########### calculate the one gluon exchange from Three particle Fock Sector quark 1 to Four particle Fock Sector quark 1 and gluon #####################!!!!!!!!!!!
        
            testvalue=fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)

                if(finalstate.ge.dimtot1+1.and.finalstate.le.dimtot1+dimtot2) then
                  color2=1
                else if(finalstate.ge.dimtot1+dimtot2+1.and.finalstate.le.dimtot1+2*dimtot2) then
                  color2=2
                else if(finalstate.ge.dimtot1+2*dimtot2+1.and.finalstate.le.dimtot1+3*dimtot2) then
                  color2=3
                endif

            CFvalue=colorfactorgtqq(color1,color2)

            If(testvalue.ne.0.and.abs(CFvalue).gt.epsilon(CFvalue)) then
              call VCVinteractiongtqq(Nmax2,Kt,b,mass2(4),mass2(5),kp4,sp4,np4,mp4,kk4,kk5,sk4,sk5,nk4,mk4,nk5,mk5,V11)
              V11=V11*CFvalue
            else
              V11=0.D0
            endif

            hamiinteract=couplings/Pplus*(V11)

            if(Abs(hamiinteract).gt.epsilon(Abs(hamiinteract)))then
              j_nz_vcvgtqq(nz)=finalstate

              hamiltonianinteractvcvgtqq(nz)=hamiinteract

              nz=nz+1

            endif
      enddo
    !enddo

    i_nz_vcvgtqq(2*dimtot2+1)=nz



    Print*, "Vertex_gtqqbar"
    ! do test=1,nz
    !  Print*,hamiltonianinteractvcvgtqq(test),i_nz_vcvgtqq(test),j_nz_vcvgtqq(test),test
    ! enddo

    return
end subroutine hamiltonianinteractgtoqqphoton