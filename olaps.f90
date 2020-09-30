       function olaps(b1,n1,m1,b2,n2,m2,dx,dy)
         implicit none
         integer n1,m1,n2,m2,kp1,kp2,particleid
         integer i,j,k,l,m
         integer nmax,olaps_cache_max,k_spstates,olaps_cache_max_k_spstates
         integer error
         double precision b1,b2,dx,dy,norm,K_tot,jp1
         double precision factlog
         double precision dx1,dy1,dx2,dy2
         
         complex(8) eval0,eval1,eval2

         double precision olaps


         call dolaps(1.0D0/b1,n1,m1,1.0D0/b2,n2,m2,dx,dy,eval0)



         norm=(2.D0*(-1)**(n1 + n2)*b1*b2*sqrt(exp(factlog(n1) + factlog(n2) + factlog(n1+abs(m1)) + factlog(n2+abs(m2)))))&
            &/(b1 **2 + b2**2)  


         olaps=norm*dble(eval0)

         return
       end function olaps
       
       subroutine dolaps(b1,n1,m1,b2,n2,m2,dx,dy,eval)
         implicit none
         integer p,q,r,s,t,u,v,rmin,rmax,w,wmin,wmax,vmin,vmax,z,zmin,zmax,n1,m1,n2,m2
         double precision lx, ly, l, alpha, cos2theta, sin2theta, costheta, sintheta, dx, dy
         double precision b1,b2
         double precision factlog,binoco,tetraco
         complex(8) sum, eval0, eval1, eval2, eval3, eval4, eval5, eval, imaginaryi
         
         imaginaryi=(0.D0, 1.D0)
         sum=0.D0
         eval0=0.D0
         eval1=0.D0
         eval2=0.D0
         eval3=0.D0

         costheta=b1/sqrt(b1**2 + b2**2)
         sintheta=b2/sqrt(b1**2 + b2**2)
         cos2theta=(b1**2 - b2**2)/(b1**2 + b2**2)
         sin2theta=(2.D0*b1*b2)/(b1**2 + b2**2)
         lx=(b1*b2*dx)/sqrt(b1**2 + b2**2)
         ly=(b1*b2*dy)/sqrt(b1**2 + b2**2)
         l=sqrt(lx*lx+ly*ly)
         alpha=atan2(ly,lx)
         rmin=0
         rmax=n1 + n2 + (-m1 + abs(m1))/2 + (m2 + abs(m2))/2

         if(abs(l).lt.epsilon(l)) then
            if(m1.eq.m2) then
               vmin=0
               vmax=min(n1,(n1 + n2)/2)
               do v=vmin,vmax
                  r=abs(m1)+2*v
                  q=n1+n2-2*v
                  u=n1-v
                  p=q+r
                  if(abs(b1-b2).lt.epsilon(b1).and.q.eq.0) then
                     eval4=1.D0
                  else
                     eval4=cos2theta**q
                  end if
                  
                  eval0=((-1)**(q - u)*sin2theta**r*binoco(p,q)*binoco(q,u)*binoco(r,v)*eval4)/exp(factlog(p))
                  sum=sum+eval0
               end do
            end if
         else
            do r=rmin,rmax
               vmin=max(0,-n2 + r + (-m2 - abs(m2))/2)
               vmax=min(r,n1 + (-m1 + abs(m1))/2)
               eval3=0.D0
               do v=vmin,vmax
                  wmin=max(0,r - m1 - 2*v)
                  wmax=min(-r + rmax,n1 - v + (-m1 + abs(m1))/2)
                  eval2=0.D0
                  do w=wmin,wmax
                     
                     if(abs(b1-b2).lt.epsilon(b1)) then
                        zmin=max(m2 + 2*v - r,0,rmax - r - w)
                     else
                        zmin=max(m2 + 2*v - r,0)
                     end if
                     
                        zmax=min(-r + rmax - w,n2 - r + v + (m2 + abs(m2))/2)
                        eval1=0.D0
                        do z = zmin,zmax
                           s = m1 + 2*v + 2*w - r
                           t = -m2 - 2*v + 2*z + r
                           q = n1 + n2 + (abs(m1) - m1)/2 + (abs(m2) + m2)/2 - r - w - z
                           u = n1 + (abs(m1) - m1)/2 - v - w
                           p = q + r + s + t
                           if(abs(b1-b2).lt.epsilon(b1).and.q.eq.0) then
                              eval5=1.D0
                           else
                              eval5=cos2theta**q
                           end if

                           eval0 = ((-1)**(q + t - u)*sin2theta**r*costheta**t*sintheta**s*l**(s + t)*binoco(q,u)*binoco(r,v) &
                                &*binoco(s,w)*binoco(t,z)*eval5*tetraco(q,r,s,t))*exp(-imaginaryi*alpha*(s + t - 2*w - 2 *z))&
                                &/exp(factlog(p))

                           eval1=eval1+eval0
                        end do
                        eval2=eval2+eval1
                     end do
                     eval3=eval3+eval2
                  end do
                  sum=sum+eval3
               end do
            end if
            
            eval=sum*exp(-l**2/2)

            return
          end subroutine dolaps