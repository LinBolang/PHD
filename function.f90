       FUNCTION tetraco(n1,n2,n3,n4)
         IMPLICIT NONE
         INTEGER n1,n2,n3,n4
         double precision tetraco,factlog
         if(min(n1,n2,n3,n4)>=0) then
            tetraco=(exp(factlog(n1+n2+n3+n4)-factlog(n1)&
                 -factlog(n2)-factlog(n3)-factlog(n4)))
         else
            tetraco=0.0D0
         end if
       END FUNCTION tetraco

! In this part we defined a function of calculating binomial coefficient

       FUNCTION binoco(n1,n2)
         IMPLICIT NONE
         INTEGER n1,n2
         double precision binoco,factlog
         if(n1-n2>=0.and.n2>=0) then
            binoco=(exp(factlog(n1)-factlog(n1-n2)-factlog(n2)))
         else

            binoco=0.0D0
         end if
       END FUNCTION binoco

       function factlog(n)
         implicit none
         integer n
         integer i
         integer tmax
         parameter (tmax=100000)
         double precision a(1:tmax+1)
         double precision lnfact,factlog
         save a
         logical,save :: init=.true.

         if(init) then
            lnfact=0.0D0
            do i=1,tmax
               a(i)=lnfact
               lnfact=lnfact+log(dble(i))
            end do
            a(tmax+1)=lnfact
            init=.false.
         end if

         if(n<0) then

            factlog=-100000.0D0
         else
            if(n<=tmax) then
               factlog=a(n+1)
            else
               lnfact=a(tmax+1)
               do i=tmax+1,n
                  lnfact=lnfact+log(dble(i))
               end do
               factlog=lnfact
            end if
         end if
       end function factlog

       subroutine error(n)
         implicit none
         integer n

         if(n.eq.1) then
            write(*,*) 'fact neg arg!'
         else
            write(*,*) 'other error'
         end if
         return
       end subroutine error

! In this part, we defined a few functions like KroneckerDelta
! and those functions mainly call by hamiltonian.f90

       integer function KroneckerDelta(i,j)
           implicit none

           integer i,j

           if(i.eq.j) then
              KroneckerDelta=1
           else
              KroneckerDelta=0
           end if

           return
         end function KroneckerDelta

         double precision function ifactor(n1,m1,n2,m2)
           implicit none

           integer n1,m1,n2,m2,KroneckerDelta

           if(m1.ge.0.and.m2.eq.m1+1) then
              ifactor= -(Sqrt(dble(n1))*KroneckerDelta(-1 + n1,n2)) + Sqrt(dble(1 + m1 + n1))*KroneckerDelta(n1,n2)
           else if(m1.gt.0.and.m2.eq.m1-1) then
              ifactor=  Sqrt(dble(m1 + n1))*KroneckerDelta(n1,n2) - Sqrt(dble(1 + n1))*KroneckerDelta(1 + n1,n2)
           else if(m1.lt.0.and.m2.eq.m1+1) then
              ifactor= Sqrt(dble(-m1 + n1))*KroneckerDelta(n1,n2) - Sqrt(dble(1 + n1))*KroneckerDelta(1 + n1,n2)
           else if(m1.le.0.and.m2.eq.m1-1) then
              ifactor= -(Sqrt(dble(n1))*KroneckerDelta(-1 + n1,n2)) + Sqrt(dble(1 - m1 + n1))*KroneckerDelta(n1,n2)
           else
              ifactor=0
           end if

           return
         end function ifactor

         integer function fidelta(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)    

         ! this function give us normalization between initial and final states

          implicit none

          integer :: np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1
          integer :: KroneckerDelta

          fidelta=KroneckerDelta(np1,nk1)*KroneckerDelta(mp1,mk1)*KroneckerDelta(sp1,sk1)*KroneckerDelta(kp1,kk1)

          return
        end function fidelta

        ! adotaq and adotbq are in momentum space

         double precision function adotaq(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)
           implicit none

           integer :: np1,mp1,nk1,mk1,KroneckerDelta
           integer :: sp1,sk1,kp1,kk1

           adotaq=KroneckerDelta(mk1,mp1)*((1 + 2*np1 + Abs(mp1))*KroneckerDelta(nk1,np1) - KroneckerDelta(1,Max(nk1,np1) -&
                & Min(nk1,np1))*Sqrt(dble(Max(nk1,np1))*(Abs(mp1) + Max(nk1,np1))))*KroneckerDelta(sp1,sk1)&
                &*KroneckerDelta(kp1,kk1)

           return
         end function adotaq

         double precision function adotbq(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)
           implicit none

           integer :: np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sp2,sk1,sk2,kp1,kp2,kk1,kk2
           integer :: KroneckerDelta
           double precision :: ifactor

           adotbq = (ifactor(np1,mp1,nk1,mk1)*ifactor(np2,mp2,nk2,mk2)*(KroneckerDelta(mk1,1 + mp1)*KroneckerDelta(mk2,-1 +&
                & mp2) + KroneckerDelta(mk1,-1 + mp1)*KroneckerDelta(mk2,1 + mp2))*KroneckerDelta(sp1,sk1)&
                &*KroneckerDelta(sp2,sk2)*KroneckerDelta(kp1,kk1)*KroneckerDelta(kp2,kk2))/2.0D0

           return
         end function adotbq

         ! adotas and adotbs are in coordinate space

         double precision function adotas(np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1)
            implicit none

            integer :: np1,mp1,nk1,mk1,sp1,sk1,kp1,kk1
            integer :: KroneckerDelta

            adotas = (KroneckerDelta(np1,nk1)*(2*np1+abs(mp1)+1)+KroneckerDelta(max(np1,nk1)-min(np1,nk1),1)*sqrt(dble(max(np1,nk1)&
              &*(max(np1,nk1)+abs(mp1)))))*KroneckerDelta(mp1,mk1)*KroneckerDelta(sp1,sk1)*KroneckerDelta(kp1,kk1)

            return
          end function adotas

          double precision function adotbs(np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2)
            implicit none

            integer :: np1,mp1,nk1,mk1,np2,mp2,nk2,mk2,sp1,sk1,sp2,sk2,kp1,kk1,kp2,kk2
            integer :: KroneckerDelta,factor
            double precision :: ifactor

            factor = np1+np2-nk1-nk2+(abs(mp1)+abs(mp2)-abs(mk1)-abs(mk2))/2

            adotbs = (0-1)**factor*(KroneckerDelta(mp1+1,mk1)*KroneckerDelta(mp2-1,mk2)+KroneckerDelta(mp1-1,mk1)&
              &*KroneckerDelta(mp2+1,mk2))*ifactor(np1,mp1,nk1,mk1)*ifactor(np2,mp2,nk2,mk2)*KroneckerDelta(sp1,sk1)&
              &*KroneckerDelta(sp2,sk2)*KroneckerDelta(kp1,kk1)*KroneckerDelta(kp2,kk2)

            return
          end function adotbs

          double precision function longipartial5(Kt,kp1,kp2,kk1,kk2,kp3)
            implicit none

            integer :: Kt,kp1,kp2,kp3,kk1,kk2
            integer :: km1,KroneckerDelta
            double precision :: value0,value1
            integer,external :: unitstep

            longipartial5=0.d0

            do km1=0,Kt-kp3

              value0=(dble(KroneckerDelta(kp1+1,km1)-KroneckerDelta(kp1-1,km1))*5.D0/6.D0+dble(KroneckerDelta(kp1-2,km1)-&
                & KroneckerDelta(kp1+2,km1))*5.D0/21.D0+dble(KroneckerDelta(kp1+3,km1)-KroneckerDelta(kp1-3,km1))*5.D0/&
                & 84.D0+dble(KroneckerDelta(kp1-4,km1)-KroneckerDelta(kp1+4,km1))*5.D0/504.D0+dble(KroneckerDelta(kp1+5,&
                & km1)-KroneckerDelta(kp1-5,km1))/1260.D0)*unitstep(Kt-kp3-km1-1)

              value1=dble(KroneckerDelta(Kt-kp3-km1-2,kk2)-KroneckerDelta(Kt-kp3-km1,kk2))*5.D0/6.D0+dble(KroneckerDelta(&
                & Kt-kp3-km1+1,kk2)-KroneckerDelta(Kt-kp3-km1-3,kk2))*5.D0/21.D0+dble(KroneckerDelta(Kt-kp3-km1-4,kk2)-&
                & KroneckerDelta(Kt-kp3-km1+2,kk2))*5.D0/84.D0+dble(KroneckerDelta(Kt-kp3-km1+3,kk2)-KroneckerDelta(&
                & Kt-kp3-km1-5,kk2))*5.D0/504.D0+dble(KroneckerDelta(Kt-kp3-km1-6,kk2)-KroneckerDelta(Kt-kp3-km1+4,kk2))/&
                & 1260.D0

              longipartial5=longipartial5+value0*value1*(dble(km1)+0.5D0)*(dble(Kt-kp3-km1)-0.5D0)

            enddo

            return
          end function longipartial5

          double precision function longipartial5_2(Kt,kp1,kp2,kk1,kk2)
            implicit none

            integer :: Kt,kp1,kp2,kk1,kk2,kmax
            integer :: km1,KroneckerDelta
            double precision :: value0,value1
            integer,external :: unitstep

            longipartial5_2=0.d0
            kmax=kp1+kp2

            do km1=0,kmax

              value0=(dble(KroneckerDelta(kp1+1,km1)-KroneckerDelta(kp1-1,km1))*5.D0/6.D0+dble(KroneckerDelta(kp1-2,km1)-&
                & KroneckerDelta(kp1+2,km1))*5.D0/21.D0+dble(KroneckerDelta(kp1+3,km1)-KroneckerDelta(kp1-3,km1))*5.D0/&
                & 84.D0+dble(KroneckerDelta(kp1-4,km1)-KroneckerDelta(kp1+4,km1))*5.D0/504.D0+dble(KroneckerDelta(kp1+5,&
                & km1)-KroneckerDelta(kp1-5,km1))/1260.D0)*unitstep(kmax-km1)

              value1=dble(KroneckerDelta(kmax-km1-1,kk2)-KroneckerDelta(kmax-km1+1,kk2))*5.D0/6.D0+dble(KroneckerDelta(&
                & kmax-km1+2,kk2)-KroneckerDelta(kmax-km1-2,kk2))*5.D0/21.D0+dble(KroneckerDelta(kmax-km1-3,kk2)-&
                & KroneckerDelta(kmax-km1+3,kk2))*5.D0/84.D0+dble(KroneckerDelta(kmax-km1+4,kk2)-KroneckerDelta(&
                & kmax-km1-4,kk2))*5.D0/504.D0+dble(KroneckerDelta(kmax-km1-5,kk2)-KroneckerDelta(kmax-km1+5,kk2))/&
                & 1260.D0

              longipartial5_2=longipartial5_2+value0*value1*(dble(km1)+0.5D0)*(dble(kmax-km1)+0.5D0)

            enddo

            return
          end function longipartial5_2

          double precision function longipartial(Kt,kp1,kp2,kk1,kk2,kp3)
            implicit none

            integer :: Kt,kp1,kp2,kp3,kk1,kk2
            integer :: km1,KroneckerDelta
            double precision :: value0,value1
            integer,external:: unitstep
            integer,dimension(2) :: SumList
            integer :: i

            longipartial=0.d0
            SumList(1)=2*kp1+1
            SumList(2)=2*kp1-1

            !do km1=0,2*(Kt-kp3-1)
            do i=1,2

              km1=SumList(i)

              value0=dble(KroneckerDelta(2*kp1+1,km1)-KroneckerDelta(2*kp1-1,km1))*unitstep(2*kt-2*kp3-km1-2)*Unitstep(km1)
              ! value0=(dble(KroneckerDelta(kp1+1,km1)-KroneckerDelta(kp1-1,km1))*5.D0/6.D0+dble(KroneckerDelta(kp1-2,km1)-&
              !   & KroneckerDelta(kp1+2,km1))*5.D0/21.D0+dble(KroneckerDelta(kp1+3,km1)-KroneckerDelta(kp1-3,km1))*5.D0/&
              !   & 84.D0+dble(KroneckerDelta(kp1-4,km1)-KroneckerDelta(kp1+4,km1))*5.D0/504.D0+dble(KroneckerDelta(kp1+5,&
              !   & km1)-KroneckerDelta(kp1-5,km1))/1260.D0)*unitstep(Kt-kp3-km1-1)

              value1=dble(KroneckerDelta(2*Kt-2*kp3-2-km1-1,2*kk2)-KroneckerDelta(2*kt-2*kp3-2-km1+1,2*kk2))
              ! value1=dble(KroneckerDelta(Kt-kp3-km1-2,kk2)-KroneckerDelta(Kt-kp3-km1,kk2))*5.D0/6.D0+dble(KroneckerDelta(&
              !   & Kt-kp3-km1+1,kk2)-KroneckerDelta(Kt-kp3-km1-3,kk2))*5.D0/21.D0+dble(KroneckerDelta(Kt-kp3-km1-4,kk2)-&
              !   & KroneckerDelta(Kt-kp3-km1+2,kk2))*5.D0/84.D0+dble(KroneckerDelta(Kt-kp3-km1+3,kk2)-KroneckerDelta(&
              !   & Kt-kp3-km1-5,kk2))*5.D0/504.D0+dble(KroneckerDelta(Kt-kp3-km1-6,kk2)-KroneckerDelta(Kt-kp3-km1+4,kk2))/&
              !   & 1260.D0

              longipartial=longipartial+value0*value1*(dble(km1)/2.0D0+0.5D0)*(dble(2*Kt-2*kp3-2-km1)/2.0D0+0.5D0)

              !Print*,longipartial5,km1,kp1,kk1,kp2,kk2,value0,value1

            enddo

            return
          end function longipartial

          double precision function longipartial_2(Kt,kp1,kp2,kk1,kk2)
            implicit none

            integer :: Kt,kp1,kp2,kk1,kk2,kmax
            integer :: km1,KroneckerDelta
            double precision :: value0,value1
            integer,external:: unitstep
            integer,dimension(2) :: SumList
            integer :: i

            longipartial_2=0.d0
            SumList(1)=2*kp1+1
            SumList(2)=2*kp1-1
            kmax=kp1+kp2

            !do km1=0,2*(Kt-kp3-1)
            do i=1,2

              km1=SumList(i)

              value0=dble(KroneckerDelta(2*kp1+1,km1)-KroneckerDelta(2*kp1-1,km1))*unitstep(2*kmax-km1)*Unitstep(km1)
              ! value0=(dble(KroneckerDelta(kp1+1,km1)-KroneckerDelta(kp1-1,km1))*5.D0/6.D0+dble(KroneckerDelta(kp1-2,km1)-&
              !   & KroneckerDelta(kp1+2,km1))*5.D0/21.D0+dble(KroneckerDelta(kp1+3,km1)-KroneckerDelta(kp1-3,km1))*5.D0/&
              !   & 84.D0+dble(KroneckerDelta(kp1-4,km1)-KroneckerDelta(kp1+4,km1))*5.D0/504.D0+dble(KroneckerDelta(kp1+5,&
              !   & km1)-KroneckerDelta(kp1-5,km1))/1260.D0)*unitstep(Kt-kp3-km1-1)

              value1=dble(KroneckerDelta(2*kmax-km1-1,2*kk2)-KroneckerDelta(2*kmax-km1+1,2*kk2))
              ! value1=dble(KroneckerDelta(Kt-kp3-km1-2,kk2)-KroneckerDelta(Kt-kp3-km1,kk2))*5.D0/6.D0+dble(KroneckerDelta(&
              !   & Kt-kp3-km1+1,kk2)-KroneckerDelta(Kt-kp3-km1-3,kk2))*5.D0/21.D0+dble(KroneckerDelta(Kt-kp3-km1-4,kk2)-&
              !   & KroneckerDelta(Kt-kp3-km1+2,kk2))*5.D0/84.D0+dble(KroneckerDelta(Kt-kp3-km1+3,kk2)-KroneckerDelta(&
              !   & Kt-kp3-km1-5,kk2))*5.D0/504.D0+dble(KroneckerDelta(Kt-kp3-km1-6,kk2)-KroneckerDelta(Kt-kp3-km1+4,kk2))/&
              !   & 1260.D0

              longipartial_2=longipartial_2+value0*value1*(dble(km1)/2.0D0+0.5D0)*(dble(2*kmax-km1)/2.0D0+0.5D0)

              !Print*,longipartial5,km1,kp1,kk1,kp2,kk2,value0,value1

            enddo

            return
          end function longipartial_2

          double precision function longitrans(Nmax,Pplus,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,kp1half,kp2half,kk1half,kk2half)
            use numbers

            implicit none

            integer :: Nmax,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2
            double precision :: Pplus,kp1half,kp2half,kk1half,kk2half

            integer :: factorp,factork,sumn,summ,N,M
            double precision :: olaps,TMC,tandeltap,tandeltak
            double precision :: meddlevalue,b,bp,bk

            !Print*,"test"

            tandeltap=sqrt(kp2half/kp1half)
            tandeltak=sqrt(kk2half/kk1half)
            sumn=floor(dble(Nmax-3)/2.D0)
            summ=Nmax-3
            longitrans=0.D0
            b=B_DEFAULT
            bp=b*sqrt(kp1half*kp2half*(kp1half+kp2half)/Pplus**3)
            bk=b*sqrt(kk1half*kk2half*(kk1half+kk2half)/Pplus**3)

            !Print*,sumn,summ
            If(mp1+mp2.eq.mk1+mk2) then

              do N=0,sumn
                do M=-summ,summ

                  factorp=np1+np2-N+(Abs(mp1)+Abs(mp2)-Abs(M)-Abs(mp1+mp2-M))/2
                  factork=nk1+nk2-N+(Abs(mk1)+Abs(mk2)-Abs(M)-Abs(mk1+mk2-M))/2

                  If(factorp.ge.0.and.factork.ge.0)then

                    meddlevalue=TMC(N,M,factorp,mp1+mp2-M,np1,mp1,np2,mp2,tandeltap)*TMC(N,M,factork,mk1+mk2-M,nk1,mk1,nk2,&
                      & mk2,tandeltak)*olaps(bp,factorp,mp1+mp2-M,bk,factork,mk1+mk2-M,0.D0,0.D0)

                    longitrans=longitrans+meddlevalue

                  endif

                enddo
              enddo

            else

              longitrans=0.D0

            endif

            return
          end function longitrans

          integer function unitstep(x)
              implicit none

              integer :: x

              if(x.ge.0) then
                unitstep=1
              else
                unitstep=0
              endif

              return
            end function unitstep

!!!!!!!!!############## vector coupling vertex for q to qg ###########################!!!!!!!!!!!!!!!!!!!!!!!

subroutine VCVinteractionqtqg(Nmax2,Kt,b,mass1,mass2,kp1,sp1,np1,mp1,kk1,kk2,sk1,sk2,nk1,mk1,nk2,mk2,Vertex)
  implicit none

  integer :: Nmax2,Kt,kp1,sp1,np1,mp1,kk1,kk2,sk1,sk2,nk1,mk1,nk2,mk2
  double precision :: b,mass1,mass2,tandelta,constant
  double precision :: kp1half,kk1half,kk2half,Pplus
  double precision :: PI=acos(-1.D0),Vertex
  integer :: n,m
  double precision,external :: TMC,vertexspinor
          
  kp1half=dble(kp1)+0.5D0
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)
  Pplus  =dble(Kt)+0.5D0

  tandelta=Sqrt(kk2half/kk1half)
  constant=b**2/PI/Sqrt(2*kk2half)
  Vertex=0.D0

  If(kp1.eq.kk1+kk2) then
    m=(sp1-sk1)/2-sk2
    n=nk1+nk2-np1+(abs(mk1)+abs(mk2)-abs(mp1)-abs(m))/2
    If(m+mp1.eq.mk1+mk2.and.n.ge.0.and.abs(m).le.1) then
      Vertex=constant*TMC(np1,mp1,n,m,nk1,mk1,nk2,mk2,tandelta)&
        & *vertexspinor(Kt,n,b,mass1,mass2,kp1,kk1,kk2,sp1,sk1,sk2)
    endif
  endif

  !Print*,vertex,sp1,sk1,sk2,n,m

  return
end subroutine VCVinteractionqtqg

double precision function vertexspinor(Kt,n,b,mass1,mass2,kp1,kk1,kk2,sp1,sk1,sk2)
  implicit none

  integer :: Kt,n,kp1,kk1,kk2,sp1,sk1,sk2
  double precision :: mass1,mass2,b
  double precision :: kp1half,kk1half,kk2half,Pplus
  integer,external :: KroneckerDelta

  kp1half=dble(kp1)+0.5D0
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)
  Pplus  =dble(Kt)+0.5D0

  if(sp1.eq.1.and.sk1.eq.1) then
    vertexspinor=2*Sqrt(dble(n+1))*(KroneckerDelta(sk2,1)+kk1half/kp1half*KroneckerDelta(sk2,-1))*(-1)**(n+1)
  else if(sp1.eq.-1.and.sk1.eq.-1) then
    vertexspinor=2*Sqrt(dble(n+1))*(kk1half/kp1half*KroneckerDelta(sk2,1)+KroneckerDelta(sk2,-1))*(-1)**(n+1)
  else if(sp1.eq.1.and.sk1.eq.-1) then
    vertexspinor=-1/b*(mass1*Pplus/kp1half-mass2*Pplus/kk1half)*Sqrt(kk2half*kk1half/Pplus/kp1half)*&
      & KroneckerDelta(sk2,1)*(-1)**(n+1)
  else if(sp1.eq.-1.and.sk1.eq.1) then
    vertexspinor=1/b*(mass1*Pplus/kp1half-mass2*Pplus/kk1half)*Sqrt(kk2half*kk1half/Pplus/kp1half)*&
      & KroneckerDelta(sk2,-1)*(-1)**(n+1)
  endif

  return
end function vertexspinor

!!!!!!!!!############## vector coupling vertex for g to qqbar ###########################!!!!!!!!!!!!!!!!!!!!!!!

subroutine VCVinteractiongtqq(Nmax2,Kt,b,mass1,mass2,kp1,sp1,np1,mp1,kk1,kk2,sk1,sk2,nk1,mk1,nk2,mk2,Vertex)
  implicit none

  integer :: Nmax2,Kt,kp1,sp1,np1,mp1,kk1,kk2,sk1,sk2,nk1,mk1,nk2,mk2
  double precision :: b,mass1,mass2,tandelta,constant
  double precision :: kp1half,kk1half,kk2half,Pplus
  double precision :: PI=acos(-1.D0),Vertex
  integer :: n,m
  double precision,external :: TMC,vertexspinorgtqq
          
  kp1half=dble(kp1)
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)+0.5D0
  Pplus  =dble(Kt)+0.5D0

  tandelta=Sqrt(kk2half/kk1half)
  constant=b**2/PI/Sqrt(2*kp1half)
  Vertex=0.D0

  If(kp1.eq.kk1+kk2+1) then
    m=sp1-(sk1+sk2)/2
    n=nk1+nk2-np1+(abs(mk1)+abs(mk2)-abs(mp1)-abs(m))/2
    If(m+mp1.eq.mk1+mk2.and.n.ge.0) then
      Vertex=constant*TMC(np1,mp1,n,m,nk1,mk1,nk2,mk2,tandelta)&
        & *vertexspinorgtqq(Kt,n,b,mass1,mass2,kp1,kk1,kk2,sp1,sk1,sk2)
    endif
  endif

  return
end subroutine VCVinteractiongtqq

double precision function vertexspinorgtqq(Kt,n,b,mass1,mass2,kp1,kk1,kk2,sp1,sk1,sk2)
  implicit none

  integer :: Kt,n,kp1,kk1,kk2,sp1,sk1,sk2
  double precision :: mass1,mass2,b
  double precision :: kp1half,kk1half,kk2half,Pplus
  integer,external :: KroneckerDelta

  kp1half=dble(kp1)
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)+0.5D0
  Pplus  =dble(Kt)+0.5D0

  if(sk1.eq.1.and.sk2.eq.-1) then
    vertexspinorgtqq=2*Sqrt(dble(n+1))*(kk1half/kp1half*KroneckerDelta(sp1,1)-kk2half/kp1half*KroneckerDelta(sp1,-1))*(-1)**(n)
  else if(sk1.eq.-1.and.sk2.eq.1) then
    vertexspinorgtqq=2*Sqrt(dble(n+1))*(-kk2half/kp1half*KroneckerDelta(sp1,1)+kk1half/kp1half*KroneckerDelta(sp1,-1))*(-1)**(n)
  else if(sk1.eq.1.and.sk2.eq.1) then
    vertexspinorgtqq=1/b*(mass1*Pplus/kk1half+mass2*Pplus/kk2half)*Sqrt(kk2half*kk1half/Pplus/kp1half)&
      & *KroneckerDelta(sp1,1)*(-1)**(n)
  else if(sp1.eq.-1.and.sk1.eq.-1) then
    vertexspinorgtqq=-1/b*(mass1*Pplus/kk1half+mass2*Pplus/kk2half)*Sqrt(kk2half*kk1half/Pplus/kp1half)&
      & *KroneckerDelta(sp1,-1)*(-1)**(n)
  endif

  return
end function vertexspinorgtqq

!!!!!!!!!!!!!!################# Instantaneous term for three particle Fock sector ###########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Instantaneous3p(Nmax2,Kt,b,binst,kp1,kp2,kk1,kk2,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2,Vertex)
  implicit none

  integer :: Nmax2,Kt,kp1,kp2,kk1,kk2,np1,mp1,np2,mp2,nk1,mk1,nk2,mk2
  double precision :: b,binst,kp1half,kp2half,kk1half,kk2half,Pplus

  double precision,external :: TMC,instolap
  double precision :: constant,PI=acos(-1.D0),Vertex,Inst3p
  integer :: n,nprime,capn

  kp1half=dble(kp1)+0.5D0
  kp2half=dble(kp2)+0.5D0
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)+0.5D0
  Pplus  =dble(Kt)+0.5D0
  Inst3p=0.0D0

  If(kp1.ne.kk1) then
    If(kp1+kp2.eq.kk1+kk2) then
      If(mp1+mp2.eq.mk1+mk2) then

        constant=1/PI*Sqrt(kp1half*kp2half*kk1half*kk2half)/(kp1half+kp2half)/(kp1half-kk1half)**2

        do capn=0,Nmax2-3

          n=np1+np2-capn+(abs(mp1)+abs(mp2)-abs(mp1+mp2))/2
          nprime=nk1+nk2-capn+(abs(mk1)+abs(mk2)-abs(mk1+mk2))/2

          If(n.ge.0.and.nprime.ge.0) then

            Inst3p=Inst3p+TMC(capn,mp1+mp2,n,0,np1,mp1,np2,mp2,Sqrt(kp2half/kp1half)) &
              & *TMC(capn,mk1+mk2,nprime,0,nk1,mk1,nk2,mk2,Sqrt(kk2half/kk1half))*instolap(nmax2,n,nprime,b,&
              & binst,kp1half,kp2half,kk1half,kk2half,Pplus)

          endif
        enddo
      endif
    endif
  endif

  Vertex=constant*Inst3p

  return
end subroutine Instantaneous3p

double precision function instolap(Nmax,n,nprime,b,binst,kp1half,kp2half,kk1half,kk2half,Pplus)
  implicit none

  integer :: n,nprime,Nmax
  double precision :: b,binst,kp1half,kp2half,kk1half,kk2half,Pplus,PI=acos(-1.D0),norm
  double precision,external :: TMC,olaps
  integer :: capn,nr
  double precision :: xk,xp

  !instolap=b**2*(-1)**(n+nprime)
  instolap=0.D0
  xp=kp1half*kp2half/(kp1half+kp2half)/Pplus
  xk=kk1half*kk2half/(kk1half+kk2half)/Pplus
  norm=binst*Sqrt(abs(kp1half-kk1half)/Pplus)/Sqrt(4*Pi)

  do capn=0,Nmax-3

    nr=nprime+n-capn

    if(nr.ge.0) then

      instolap=instolap+TMC(capn,0,nr,0,n,0,nprime,0,1)*b*(-1)**capn/Sqrt(Pi)*olaps(b,nr,0,&
        & binst*Sqrt(abs(kp1half-kk1half)/Pplus),0,0,0.D0,0.D0)*norm

    endif

  enddo

  return
end function instolap

!!!!!!!!!!!!!!!!!############################ instantaneous interaction between three particle system and five particle system ##############!!!!!!!!!!!!!!!!!!!!

double precision function Instantaneous3t5(Nmax2,Kt,b,binst,kp1,kk1,kk2,kk3,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3)
  implicit none

  integer :: Nmax2,Kt,kp1,kk1,kk2,kk3,np1,mp1,nk1,mk1,nk2,mk2,nk3,mk3
  double precision :: b,binst,kp1half,kk1half,kk2half,kk3half,Pplus

  double precision,external :: TMC,instolap3t5
  double precision :: constant,PI=acos(-1.D0)
  integer :: n,nprime,capn

  kp1half=dble(kp1)+0.5D0
  kk1half=dble(kk1)+0.5D0
  kk2half=dble(kk2)+0.5D0
  kk3half=dble(kk3)+0.5D0
  Pplus  =dble(Kt)+0.5D0
  Instantaneous3t5=0.0D0

  If(kp1.eq.kk1+kk2+kk3+1) then
    If(mp1.eq.mk1+mk2+mk3) then

      constant=1/PI*Sqrt(kk1half*kk2half*kk3half/kp1half)/(kp1half-kk1half)**2

      do capn=0,Nmax2-5

        nprime=nk2+nk3-capn+(abs(mk2)+abs(mk3)-abs(mk2+mk3))/2
        n=capn+nk1-np1+(abs(mk2+mk3)+abs(mk1)-abs(mp1))/2

        If(n.ge.0.and.nprime.ge.0) then

          Instantaneous3t5=Instantaneous3t5+TMC(capn,mk2+mk3,nprime,0,nk2,mk2,nk3,mk3,Sqrt(kk3half/kk2half)) &
            & *TMC(np1,mp1,n,0,capn,mk2+mk3,nk1,mk1,Sqrt(kk1half/(kk2half+kk3half)))*instolap3t5(nmax2,n,nprime,b,&
            & binst,kp1half,kk1half,kk2half,kk3half,Pplus)

        endif
      enddo
    endif
  endif

  return
end function Instantaneous3t5

double precision function instolap3t5(Nmax,n,nprime,b,binst,kp1half,kk1half,kk2half,kk3half,Pplus)
  implicit none

  integer :: n,nprime,Nmax
  double precision :: b,binst,kp1half,kk1half,kk2half,kk3half,Pplus,PI=acos(-1.D0),norm
  double precision,external :: olaps

  instolap3t5=0.D0
  norm=binst*Sqrt(abs(kp1half-kk1half)/Pplus)/Sqrt(4*Pi)

  instolap3t5=b*(-1)**nprime/Sqrt(Pi)*olaps(b,n,0,binst*Sqrt(abs(kp1half-kk1half)/Pplus),0,0,0.D0,0.D0)*norm

  return
end function instolap3t5