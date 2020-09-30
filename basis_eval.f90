!!!!!!!!###### dimension calculation and search formula for three particle Fock Sector #####!!!!!!!!!!!!!!!
       subroutine dimtotal3p(nmax2, mj, kt, dim)
         implicit none
         integer nmax1, nmax2, mj, kt, dim

         dim=0
         call dimtotal2(nmax2, mj, kt, dim)
         
         return
       end subroutine dimtotal3p

       subroutine dimtotal2(nmax2, mj, kt, dim2)
         implicit none
         integer nmax2, mj, kt, dim2, eval
         double precision binoco

         call sum2m1m2n1n2n3s1s2s3(nmax2, mj, eval)
         dim2=nint(binoco(Kt+1,2))*eval!nint((binoco(Kt+1,2)-floor((dble(Kt)+1.0d0)/2.0d0))/2)

         return
       end subroutine dimtotal2

       
       subroutine sum2m1m2n1n2n3s1s2s3(nm, mj, eval)
         implicit none
         integer nm, mj, eval, eval0, s1index, s2index, s3index

         eval=0
         do s1index=-1,1,2
            do s2index=-1,1,2
               do s3index=-1,1,2
                  
                  call sum2m1m2n1n2n3(2,nm,mj,s1index,s2index,s3index,eval0)
                  eval=eval+eval0

               end do
            end do
         end do

         return
       end subroutine sum2m1m2n1n2n3s1s2s3


       subroutine sum2m1m2n1n2n3(nmax1, nmax2, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax1, nmax2, nu, nl, mjindex, s1index, s2index, s3index, eval
         double precision s1, s2, s3, mj

         nu=nmax2-3
         nl=nmax1-2-1
         eval=0
         
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0


         if(mod(nu-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nu=nu
         else
            nu=nu-1
         end if

         if(mod(nl-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nl=nl
         else
            nl=nl-1
         end if

         if(nu.ge.nl+1) then
            if(nu.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval - ((2 + nu - nint(abs(mj - s1 - s2 - s3)))*(-40 + 10*mj**2 - 24*nu - 3*nu**2 + (10*s1**2 + 20*s1*s2 &
                    &+ 10 *s2**2 + 20*s1*s3 + 20*s2*s3 + 10*s3**2 - 20*mj*(s1 + s2 + s3)) - 9*(4 + nu)*nint(abs(mj - s1 - s2 -&
                    & s3)) - 18*nint(abs(mj - s1 - s2 - s3))**2)*(24 + 10*nu + nu**2 - 2*(5 + nu)*nint(abs(mj - s1 - s2 - s3)) +&
                    & nint(abs(mj - s1 - s2 - s3))**2))/1920
            end if
            if(nl.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval + ((2 + nl - nint(abs(mj - s1 - s2 - s3)))*(-40 + 10*mj**2 - 24*nl - 3*nl**2 + (10*s1**2 + 20*s1*s2 &
                    &+ 10 *s2**2 + 20*s1*s3 + 20*s2*s3 + 10*s3**2 - 20*mj*(s1 + s2 + s3)) - 9*(4 + nl)*nint(abs(mj - s1 - s2 -&
                    & s3)) - 18*nint(abs(mj - s1 - s2 - s3))**2)*(24 + 10*nl + nl**2 - 2*(5 + nl)*nint(abs(mj - s1 - s2 - s3)) +&
                    & nint(abs(mj - s1 - s2 - s3))**2))/1920
            end if
         end if
         
         return
       end subroutine sum2m1m2n1n2n3
       

       subroutine sum2m1(nmax1, nmax2, n1, n2, n3, mjindex, s1index, s2index, s3index, m2, eval)
         implicit none
         integer nmax1, nmax2, n1, n2, n3, mjindex, s1index, s2index, s3index, m2, eval, nu, nl
         double precision s1, s2, s3, mj

         nu=nmax2-3 - 2*n1 - 2*n2 - 2*n3 - abs(m2)
         nl=nmax1-2-1 - 2*n1 - 2*n2 - 2*n3 - abs(m2)
         eval=0
         
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         
         if(mod(nu-nint(abs(mj-s1-s2-s3-m2)),2).eq.0) then
            nu=nu
         else
            nu=nu-1
         end if

         if(mod(nl-nint(abs(mj-s1-s2-s3-m2)),2).eq.0) then
            nl=nl
         else
            nl=nl-1
         end if

         if(nu.ge.nl+1) then
            if(nu.ge.nint(abs(mj-s1-s2-s3-m2))) then
               eval=eval + nu + 1
            end if
            if(nl.ge.nint(abs(mj-s1-s2-s3-m2))) then
               eval=eval - (nl + 1)
            end if
         end if
         
         return
       end subroutine sum2m1
       

       subroutine sum2m1m2(nmax1, nmax2, n1, n2, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax1, nmax2, n1, n2, n3, mjindex, s1index, s2index, s3index, eval, nu, nl
         double precision s1, s2, s3, mj

         nu=nmax2-3 - 2*n1 - 2*n2 - 2*n3
         nl=nmax1-2-1 - 2*n1 - 2*n2 - 2*n3
         eval=0
         
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         
         if(mod(nu-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nu=nu
         else
            nu=nu-1
         end if

         if(mod(nl-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nl=nl
         else
            nl=nl-1
         end if

         if(nu.ge.nl+1) then
            if(nu.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval + (4 + 3*nu*(2 + nu) - nint((mj - s1 - s2 - s3))**2)/4
            end if
            if(nl.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval - (4 + 3*nl*(2 + nl) - nint((mj - s1 - s2 - s3))**2)/4
            end if
         end if
         
         return
       end subroutine sum2m1m2
       

       subroutine sum2m1m2n1(nmax1, nmax2, n2, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax1, nmax2, n2, n3, mjindex, s1index, s2index, s3index, eval, nu, nl
         double precision s1, s2, s3, mj

         nu=nmax2-3 - 2*n2 - 2*n3
         nl=nmax1-2-1 - 2*n2 - 2*n3
         eval=0
         
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         
         if(mod(nu-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nu=nu
         else
            nu=nu-1
         end if

         if(mod(nl-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nl=nl
         else
            nl=nl-1
         end if

         if(nu.ge.nl+1) then
            if(nu.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval + nint(((2 + Nu - Abs(mj - s1 - s2 - s3))*(4 - mj**2 + 4*Nu + Nu**2 - s1**2 - 2*s1*s2 - s2**2 - 2*s1*s3 &
                    &- 2*s2*s3 - s3**2 + 2*mj*(s1 + s2 + s3) + (2 + Nu)*Abs(mj - s1 - s2 - s3) + Abs(mj - s1 - s2 - s3)**2))/8)
            end if
            if(nl.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval - nint(((2 + Nl - Abs(mj - s1 - s2 - s3))*(4 - mj**2 + 4*Nl + Nl**2 - s1**2 - 2*s1*s2 - s2**2 - 2*s1*s3 &
                    &- 2*s2*s3 - s3**2 + 2*mj*(s1 + s2 + s3) + (2 + Nl)*Abs(mj - s1 - s2 - s3) + Abs(mj - s1 - s2 - s3)**2))/8)
            end if
         end if
         
         return
       end subroutine sum2m1m2n1

       subroutine sum2m1m2n1n2(nmax1, nmax2, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax1, nmax2, n3, mjindex, s1index, s2index, s3index, eval, nu, nl
         double precision s1, s2, s3, mj

         nu=nmax2-3 - 2*n3
         nl=nmax1-2-1 - 2*n3
         eval=0
         
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         
         if(mod(nu-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nu=nu
         else
            nu=nu-1
         end if

         if(mod(nl-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nl=nl
         else
            nl=nl-1
         end if

         if(nu.ge.nl+1) then
            if(nu.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval + nint(((2 + Nu - Abs(mj - s1 - s2 - s3))*(4 + Nu - Abs(mj - s1 - s2 - s3))*(8 - 2*mj**2 + 6*Nu + Nu**2 &
                    &- 2*s1**2 - 4*s1*s2 - 2*s2**2 - 4*s1*s3 - 4*s2*s3 - 2*s3**2 + 4*mj*(s1 + s2 + s3) + 2*(3 + Nu)*Abs(mj - s1 -&
                    & s2 - s3) + 3*Abs(mj - s1 - s2 - s3)**2))/64)
            end if
            if(nl.ge.nint(abs(mj-s1-s2-s3))) then
               eval=eval - nint(((2 + Nl - Abs(mj - s1 - s2 - s3))*(4 + Nl - Abs(mj - s1 - s2 - s3))*(8 - 2*mj**2 + 6*Nl + Nl**2 &
                    &- 2*s1**2 - 4*s1*s2 - 2*s2**2 - 4*s1*s3 - 4*s2*s3 - 2*s3**2 + 4*mj*(s1 + s2 + s3) + 2*(3 + Nl)*Abs(mj - s1 -&
                    & s2 - s3) + 3*Abs(mj - s1 - s2 - s3)**2))/64)
            end if
         end if
         
         return
       end subroutine sum2m1m2n1n2
       
       
       subroutine snmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval

         call sum2m1m2n1n2n3s1s2s3(nm, mj, eval)

         return
       end subroutine snmlen2

       subroutine dddnmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval

         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval)

         return
       end subroutine dddnmlen2


       subroutine uddnmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval

         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval1)
         eval=eval0+eval1

         return
       end subroutine uddnmlen2

       
       subroutine dudnmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval2, eval

         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, -1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval1)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval2)
         
         eval=eval0+eval1+eval2

         return
       end subroutine dudnmlen2


       subroutine uudnmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval2, eval3, eval

         call sum2m1m2n1n2n3(2,nm, mj, 1, 1, -1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, -1, eval1)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval2)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval3)
         
         eval=eval0+eval1+eval2+eval3

         return
       end subroutine uudnmlen2


       subroutine ddunmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval2, eval3, eval4, eval

         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, 1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, 1, 1, -1, eval1)
         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, -1, eval2)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval3)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval4)
         
         eval=eval0+eval1+eval2+eval3+eval4

         return
       end subroutine ddunmlen2


       subroutine udunmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval2, eval3, eval4, eval5, eval

         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, 1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, 1, eval1)
         call sum2m1m2n1n2n3(2,nm, mj, 1, 1, -1, eval2)
         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, -1, eval3)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval4)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval5)
         
         eval=eval0+eval1+eval2+eval3+eval4+eval5

         return
       end subroutine udunmlen2


       subroutine duunmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval0, eval1, eval2, eval3, eval4, eval5, eval6, eval

         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, 1, eval0)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, 1, eval1)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, 1, eval2)
         call sum2m1m2n1n2n3(2,nm, mj, 1, 1, -1, eval3)
         call sum2m1m2n1n2n3(2,nm, mj, -1, 1, -1, eval4)
         call sum2m1m2n1n2n3(2,nm, mj, 1, -1, -1, eval5)
         call sum2m1m2n1n2n3(2,nm, mj, -1, -1, -1, eval6)
         
         eval=eval0+eval1+eval2+eval3+eval4+eval5+eval6

         return
       end subroutine duunmlen2

       
       subroutine uuunmlen2(nm, mj, eval)
         implicit none
         integer nm, mj, eval

         call sum2m1m2n1n2n3s1s2s3(nm, mj, eval)

         return
       end subroutine uuunmlen2


       subroutine ksnmseg2(nm, mj, kt, k1, k2, eval)
         implicit none
         integer nm, mj, kt, k1, k2, eval0, eval

         call snmlen2(nm, mj, eval0)
         eval=eval0*(k1-k2*(1+k2-2*(Kt+1))/2)

         return
       end subroutine ksnmseg2


       subroutine m1m2n1n2n3len2(nmax, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax, nm, n3, mjindex, s1index, s2index, s3index, eval
         double precision s1, s2, s3, mj

         nm=nmax
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         eval=0

         if(mod(nm-3-2*n3-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nm=nm
         else
            nm=nm-1
         end if

         if(nm-3-2*n3.ge.nint(abs(mj-s1-s2-s3))) then
            eval=nint((n3*(15 - 28*(-1 + n3) - 32*(-1 + n3)**2 + 72*(-1 + n3)**3 + 48*(-1 + n3)**4 + 60*(-1 + n3)*Nm - 120*(-1 +&
                 & n3) **2 *Nm - 120*(-1 + n3)**3*Nm - 30*Nm**2 + 60*(-1 + n3)*Nm**2 + 120*(-1 + n3)**2*Nm**2 - 60*(-1 + n3)*Nm&
                 &**3 + 15 *Nm**4 - 10*Mj**2*(-3 + 2*(-1 + n3) + 4*(-1 + n3)**2 - 6*(-1 + n3)*Nm + 3*Nm**2) + 30*s1**2 - 20*(-1 +&
                 & n3)*s1**2 - 40*(-1 + n3)**2*s1**2 + 60*(-1 + n3)*Nm*s1**2 - 30*Nm**2*s1**2 + 60*s1*s2 - 40*(-1 + n3)*s1*s2 -&
                 & 80*(-1 + n3)**2 *s1*s2 + 120*(-1 + n3)*Nm*s1*s2 - 60*Nm**2*s1*s2 + 30*s2**2 - 20*(-1 + n3)*s2**2 - 40*(-1 +&
                 & n3)**2*s2**2 + 60*(-1 + n3) *Nm*s2**2 - 30*Nm**2*s2**2 + 60*s1*s3 - 40*(-1 + n3)*s1*s3 - 80*(-1 + n3)**2*s1*s3&
                 & + 120*(-1 + n3)*Nm*s1*s3 - 60 *Nm**2*s1*s3 + 60*s2*s3 - 40*(-1 + n3)*s2*s3 - 80*(-1 + n3)**2*s2*s3 + 120*(-1 +&
                 & n3)*Nm*s2*s3 - 60*Nm**2*s2 *s3 + 30*s3**2 - 20*(-1 + n3)*s3**2 - 40*(-1 + n3)**2*s3**2 + 60*(-1 + n3)*Nm*s3**2&
                 & - 30*Nm**2*s3**2 + 20*Mj*(-3 + 2*(-1 + n3) + 4*(-1 + n3)**2 - 6*(-1 + n3)*Nm + 3*Nm**2)*(s1 + s2 + s3) - 60*(&
                 &-1 + n3 - Nm)*(-Mj + s1 + s2 + s3) **2 *Abs(Mj - s1 - s2 - s3) - 30*(2 + Mj**2 + s1**2 + s2**2 + 2*s2*s3 + s3&
                 &**2 + 2*s1*(s2 + s3) - 2*Mj*(s1 + s2 + s3)) *Abs(Mj - s1 - s2 - s3)**2 + 60*(-1 + n3 - Nm)*Abs(Mj - s1 - s2 -&
                 & s3)**3 + 45*Abs(Mj - s1 - s2 - s3)**4)) /960)
         end if

         return
       end subroutine m1m2n1n2n3len2

       
       subroutine m1m2n1n2len2(nmax, n2, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax, nm, n2, n3, mjindex, s1index, s2index, s3index, eval
         double precision s1, s2, s3, mj

         nm=nmax
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         eval=0

         if(mod(nm-3-2*n2-2*n3-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nm=nm
         else
            nm=nm-1
         end if

         if(nm-3-2*n2-2*n3.ge.nint(abs(mj-s1-s2-s3))) then
            eval=nint(-(n2*((-n2 - 2*n3 + Nm)*(-1 + Mj**2 - 2*(-1 + n2)**2 - 4*n3 - 4*n3**2 + 2*Nm + 4*n3*Nm - Nm**2 + (-1 + n2) &
                 &*(-4 - 4*n3 + 2*Nm) + s1**2 + 2*s1*s2 + s2**2 + 2*s1*s3 + 2*s2*s3 + s3**2 - 2*Mj*(s1 + s2 + s3)) - (-Mj + s1 +&
                 & s2 + s3)**2*Abs(Mj - s1 - s2 - s3) + Abs(Mj - s1 - s2 - s3)**3))/8)
         end if

         return
       end subroutine m1m2n1n2len2


       subroutine m1m2n1len2(nmax, n1, n2, n3, mjindex, s1index, s2index, s3index, eval)
         implicit none
         integer nmax, nm, n1, n2, n3, mjindex, s1index, s2index, s3index, eval
         double precision s1, s2, s3, mj

         nm=nmax
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0
         eval=0

         if(mod(nm-3-2*n1-2*n2-2*n3-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nm=nm
         else
            nm=nm-1
         end if

         if(nm-3-2*n1-2*n2-2*n3.ge.nint(abs(mj-s1-s2-s3))) then    
            eval=nint((n1*(13 - Mj**2 + 4*(-1 + n1)**2 + 24*n2 + 12*n2**2 + 24*n3 + 24*n2*n3 + 12*n3**2 + 2*(-1 + n1)*(7 + 6*n2 +&
                 & 6*n3 - 3*Nm) - 12*Nm - 12*n2*Nm - 12*n3*Nm + 3*Nm**2 - s1**2 - 2*s1*s2 - s2**2 - 2*s1*s3 - 2*s2*s3 - s3**2 + 2&
                 &*Mj*(s1 + s2 + s3)))/4)
         end if

         
         return
       end subroutine m1m2n1len2


       subroutine m1m2len2(nmax, n1, n2, n3, mjindex, m2, s1index, s2index, s3index, eval)
         implicit none
         integer nmax, nm, n1, n2, n3, mjindex, m2, s1index, s2index, s3index, eval
         double precision s1,s2,s3,mj


         nm=nmax
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0

         eval=0

         if(mod(nm-3-2*n1-2*n2-2*n3-nint(abs(mj-s1-s2-s3)),2).eq.0) then
            nm=nm
         else
            nm=nm-1
         end if

         if(nm-3-2*n1-2*n2-2*n3-abs(m2).ge.nint(abs(mj-s1-s2-s3-m2))) then
            if(m2.lt.1) then            
               eval=nint(((-7 + 2*m2 + Mj - s1 - s2 - s3 - 6*n1 - 6*n2 - 6*n3 + 3*Nm )*(-3 + 2*m2 - Mj +s1 + s2 + s3 - &
                & 2*n1 - 2*n2 - 2*n3 + Nm ))/8)
            else
               eval=nint((-4*m2**2 - 4*m2*(3 + 4*n1 + 4*n2 + 4*n3 - 2*Nm) - (3 + Mj - s1 - s2 - s3 + 2*n1 + 2*n2 + 2*n3 - Nm )*(-7&
                    & + Mj- s1 - s2 - s3 - 6*n1 - 6*n2 - 6*n3 + 3*Nm ))/8)
            end if
         end if
      
         
         return
       end subroutine m1m2len2


       subroutine m1len2(nmax, n1, n2, n3, mjindex, m1, m2, s1index, s2index, s3index, eval)
         implicit none
         integer nmax, nm, n1, n2, n3, mjindex, m1, m2, s1index, s2index, s3index, eval
         double precision s1, s2, s3, mj

         nm=nmax
         s1=dble(s1index)/2.D0
         s2=dble(s2index)/2.D0
         s3=dble(s3index)/2.D0
         mj=dble(mjindex)/2.D0

         eval=0

         if(mod(nm-3-2*n1-2*n2-2*n3-abs(m2)-nint(abs(mj-s1-s2-s3-m2)),2).eq.0) then
            nm=nm
         else
            nm=nm-1
         end if



         if(nm-3-2*n1-2*n2-2*n3-abs(m2)-nint(abs(mj-s1-s2-s3-m2)).ge.0) then
            eval=nint(1 + m1 + ( Nm - 3 - Mj + s1 + s2 + s3 + m2 - 2*n1 - 2*n2 - 2*n3  - Abs(m2))/2)
         end if
         
         return
       end subroutine m1len2
       

       subroutine search2(nm1, nm2, mj, kt, k1, k2, k3, s1index, s2index, s3index, n1, m1, n2, m2, n3, m3, eval)
         implicit none
         integer nm1, nm2, mj, kt, k1, k2, k3, s1index, s2index, s3index, n1, m1, n2, m2, n3, m3
         integer eval0, eval1, eval2, eval3, eval4, eval5, eval6, eval7, eval8, eval9, eval10, eval11, eval12, eval13, eval

         eval0=0
         eval1=0
         eval2=0
         eval3=0
         eval4=0
         eval5=0
         eval6=0
         eval7=0
         eval8=0

         
         call ksnmseg2(nm2,mj,kt,k1,k2,eval1)

         if(s1index.eq.1.and.s2index.eq.-1.and.s3index.eq.-1) then
            call dddnmlen2(nm2,mj,eval2)
         end if
         
         if(s1index.eq.-1.and.s2index.eq.1.and.s3index.eq.-1) then
            call uddnmlen2(nm2,mj,eval3)
         end if

         if(s1index.eq.1.and.s2index.eq.1.and.s3index.eq.-1) then
            call dudnmlen2(nm2,mj,eval4)
         end if

         if(s1index.eq.-1.and.s2index.eq.-1.and.s3index.eq.1) then
            call uudnmlen2(nm2,mj,eval5)
         end if

         if(s1index.eq.1.and.s2index.eq.-1.and.s3index.eq.1) then
            call ddunmlen2(nm2,mj,eval6)
         end if

         if(s1index.eq.-1.and.s2index.eq.1.and.s3index.eq.1) then
            call udunmlen2(nm2,mj,eval7)
         end if

         if(s1index.eq.1.and.s2index.eq.1.and.s3index.eq.1) then
            call duunmlen2(nm2,mj,eval8)
         end if

         call m1m2n1n2n3len2(nm2, n3, mj, s1index, s2index, s3index, eval9)
         call m1m2n1n2len2(nm2, n2, n3, mj, s1index, s2index, s3index, eval10)
         call m1m2n1len2(nm2, n1, n2, n3, mj, s1index, s2index, s3index, eval11)
         call m1m2len2(nm2, n1, n2, n3, mj, m2, s1index, s2index, s3index, eval12)
         call m1len2(nm2, n1, n2, n3, mj, m1, m2, s1index, s2index, s3index, eval13)

         eval=eval0+eval1+eval2+eval3+eval4+eval5+eval6+eval7+eval8+eval9+eval10+eval11+eval12+eval13

          !print *, eval0, eval1,eval2,eval3,eval4,eval5,eval6,eval7,eval8,eval9,eval10,eval11,eval12,eval13

         return
       end subroutine search2


       subroutine search(nm1, nm2, mj, kt, k1, k2, k3, s1index, s2index, s3index, n1, m1, n2, m2, n3, m3, eval)
         implicit none
         integer nm1, nm2, mj, kt, k1, k2, k3, s1index, s2index, s3index, n1, m1, n2, m2, n3, m3, eval
         
         eval=-1

         if(abs(s3index).eq.1.and.k1+k2+k3.eq.kt &
              &-1.and.k1.ge.0.and.k2.ge.0.and.k3.ge.0.and.abs(s1index).eq.1.and.abs(s2index).eq.1.and.2*n1+2*n2+2*n3+abs(m1) &
              &+abs(m2)+abs(m3)+3.le.nm2.and.(mj-s1index-s2index-s3index)/2-m1-m2-m3.eq.0.and.n1.ge.0.and.n2.ge.0.and.n3.ge.0) then
            call search2(nm1, nm2, mj, kt, k1, k2, k3, s1index, s2index, s3index, n1, m1, n2, m2, n3, m3, eval)
            !Print *, "test"
         end if

         return
       end subroutine search

!!!!!!!!!! dimension calculation for qqqg Fock Sector !!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine dimtotal4p(Nmax,Mj,Kt,dim)
        implicit none

    integer :: Nmax,Mj,Kt,dim
    integer :: K
    integer,dimension(17) :: spin4p

    call SumK4p(Kt,K)
    call Sumspin4p(Nmax,Mj,spin4p)

    dim=k*spin4p(17)

        return
    end subroutine dimtotal4p

    subroutine SumK4p(Kt,LenK)
        implicit none
    integer :: Kt,LenK
    integer :: k1,k2,k3
    LenK=0

    do k1=0,kt-2
        do k2=0,kt-2-k1
                LenK=LenK+kt-k1-k2-1
        enddo
    enddo

        return
    end subroutine SumK4p

    subroutine sumspin4p(Nmax,Mj,SumS)
        use basis_info
        implicit none
    integer :: Nmax, Mj
    integer, dimension(17) :: SumS ! the matrix length equal to 2^4+1
    integer, dimension(6)  :: dimNM !record different spin combination corresponding to the number of N&M combination
    integer :: s1index,s2index,s3index,s4index
    integer :: i,m

    i=1
    SumS(1)=0
    call sumnm4p(Nmax,Mj,dimNM)

      do s4index=-1,1,2
        do s3index=-1,1,2
          do s2index=-1,1,2
            do s1index=-1,1,2
                i=i+1
                m=(mj-s1index-s2index-s3index)/2-s4index
                SumS(i)=SumS(i-1)+dimNM(m-(Mj-7)/2)
            enddo
          enddo
        enddo
      enddo

        return
    end subroutine sumspin4p

    subroutine sumnm4p(Nmax,Mj,dim)
        implicit none

    integer :: Nmax,Mj
    integer, dimension(6) :: dim ! the dimension come from 5/2-(-5/2)=6
    integer :: n1,m1,n2,m2,n3,m3,n4,m4
    integer :: i,j,m,sumN

    do j=1,6

        i=0
        m=(Mj-5)/2+j-1
        dim(j)=0

        do n4=0,floor(dble(nmax-4)/2.0D0)
          do n3=0,floor(dble(nmax-4)/2.0D0)-n4
            do n2=0,floor(dble(nmax-4)/2.0D0)-n4-n3
              do n1=0,floor(dble(nmax-4)/2.0D0)-n4-n3-n2
                do m3=-floor(dble(nmax-4-2*n1-2*n2-2*n3-2*n4-m)/2.0D0)&
                    &,ceiling(dble(nmax-4-2*n1-2*n2-2*n3-2*n4+m)/2.0D0)
                  do m2=-floor(dble(nmax-4-2*n1-2*n2-2*n3-2*n4-m+m3)/2.0D0)&
                      &,ceiling(dble(nmax-4-2*n1-2*n2-2*n3-2*n4+m-m3)/2.0D0)
                    do m1=-floor(dble(nmax-4-2*n1-2*n2-2*n3-2*n4-m+m3+m2)/2.0D0)&
                        &,ceiling(dble(nmax-4-2*n1-2*n2-2*n3-2*n4+m-m3-m2)/2.0D0)

                        m4=m-m1-m2-m3
                        SumN=2*(n1+n2+n3+n4)+abs(m1)+abs(m2)+abs(m3)+abs(m4)+4

                        if(SumN.le.Nmax) then
                            i=i+1
                        endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

      dim(j)=i

    enddo


        return
    end subroutine sumnm4p

!!!!!!!!!! dimension calculation for qqqqqbar Fock Sector !!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dimtotal5p(Nmax,Mj,Kt,dim)
        implicit none
    integer :: Nmax,Mj,Kt,dim
    integer :: K5p
    integer, dimension(33) :: spin5p

    call SumK5p(Kt,K5p)
    call Sumspin5p(Nmax,Mj,spin5p)

    dim=K5p*spin5p(33)

        return
    end subroutine dimtotal5p

    subroutine SumK5p(Kt,LenK)
        implicit none
    integer :: Kt,LenK
    integer :: k1,k2,k3
    LenK=0

    do k1=0,kt-2
        do k2=0,kt-2-k1
            do k3=0,kt-2-k2-k1
                LenK=LenK+kt-k1-k2-k3-1
            enddo
        enddo
    enddo

        return
    end subroutine SumK5p

    subroutine sumspin5p(Nmax,Mj,SumS)
        implicit none
    integer :: Nmax, Mj
    integer, dimension(33) :: SumS ! the matrix length equal to 2^5+1
    integer, dimension(6)  :: dimNM !record different spin combination corresponding to the number of N&M combination
    integer :: s1index,s2index,s3index,s4index,s5index
    integer :: i,m

    i=1
    SumS(1)=0
    call sumnm5p(Nmax,Mj,dimNM)

    do s5index=-1,1,2
      do s4index=-1,1,2
        do s3index=-1,1,2
          do s2index=-1,1,2
            do s1index=-1,1,2
                i=i+1
                m=(mj-s1index-s2index-s3index-s4index-s5index)/2
                SumS(i)=SumS(i-1)+dimNM(m-(Mj-7)/2)
            enddo
          enddo
        enddo
      enddo
    enddo

        return
    end subroutine sumspin5p

    subroutine sumnm5p(Nmax,Mj,dim)
        implicit none

    integer :: Nmax,Mj
    integer, dimension(6) :: dim ! the dimension come from 5/2-(-5/2)=6
    integer :: n1,m1,n2,m2,n3,m3,n4,m4,n5,m5
    integer :: i,j,m,sumN

    do j=1,6

        i=0
        m=(Mj-7)/2+j
        dim(j)=0

      do n5=0,floor(dble(nmax-5)/2)
        do n4=0,floor(dble(nmax-5)/2)-n5
          do n3=0,floor(dble(nmax-5)/2)-n5-n4
            do n2=0,floor(dble(nmax-5)/2)-n5-n4-n3
              do n1=0,floor(dble(nmax-5)/2)-n5-n4-n3-n2
                do m4=-floor(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5-m)/2)&
                  &,ceiling(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5+m)/2)
                  do m3=-floor(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m4)/2)&
                    &,ceiling(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m4)/2)
                    do m2=-floor(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m4+m3)/2)&
                      &,ceiling(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m4-m3)/2)
                      do m1=-floor(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m4+m3+m2)/2)&
                        &,ceiling(dble(nmax-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m4-m3-m2)/2)

                        m5=m-m1-m2-m3-m4
                        SumN=2*(n1+n2+n3+n4+n5)+abs(m1)+abs(m2)+abs(m3)+abs(m4)+abs(m5)+5

                        if(SumN.le.Nmax) then
                            i=i+1
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      dim(j)=i

    enddo


        return
    end subroutine sumnm5p