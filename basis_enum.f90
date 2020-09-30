 !!!!!!!!!! this part is basis enumeration. !!!!!!!!!!!!!!!
!!!!!!!!!!!! the First part is the basis enumeration for the three particle Fock sector qqq !!!!!!!!!!!!!!!!!!!!
        subroutine mpstatesfill3p(nmax2,mj,kt,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat&
              &,nmpstate)
           implicit none
           integer nmax2,mj,kt,nmpstate
           integer,dimension(*) :: s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat
           integer i,k1,k2,k3,n1,m1,n2,m2,n3,m3,s1index,s2index,s3index,sum1,sum2,sum3,sumtot

           i=0

           do k2=0,kt-1
              do k1=0,kt-1-k2
                 do s3index=-1,1,2
                    do s2index=-1,1,2
                       do s1index=-1,1,2
                          do n3=0,floor(dble(nmax2)/2)
                             do n2=0,floor(dble(nmax2)/2)-n3
                                do n1=0,floor(dble(nmax2)/2)-n2-n3
                                   do m2=-floor(dble(nmax2-3-2*n1-2*n2-2*n3-(mj-s1index-s2index-s3index)/2)/2)&
                                        &,ceiling(dble(nmax2-3-2*n1-2*n2-2*n3+(dble(mj-s1index-s2index-s3index)/2))/2)
                                      do m1=-floor(dble(nmax2-3-2*n1-2*n2-2*n3-(dble(mj-s1index-s2index-s3index)/2-m2))/2)&
                                           &,ceiling(dble(nmax2-3-2*n1-2*n2-2*n3+(dble(mj-s1index-s2index-s3index)/2-m2))/2)


                                         k3=kt-k1-k2-1
                                         m3=nint(dble(mj-s1index-s2index-s3index)/2)-m1-m2
                                         sum1=2*n1+abs(m1)+1
                                         sum2=2*n2+abs(m2)+1
                                         sum3=2*n3+abs(m3)+1
                                         sumtot=sum1+sum2+sum3
                                         if(sumtot.le.nmax2) then
                                            i=i+1
                                            s1dat(i)=s1index
                                            s2dat(i)=s2index
                                            s3dat(i)=s3index
                                            n1dat(i)=n1
                                            m1dat(i)=m1
                                            n2dat(i)=n2
                                            m2dat(i)=m2
                                            n3dat(i)=n3
                                            m3dat(i)=m3
                                            k1dat(i)=k1
                                            k2dat(i)=k2
                                            k3dat(i)=k3
                                            ! Print*, k1dat(i),k2dat(i),k3dat(i),s1dat(i),s2dat(i),s3dat(i),n1dat(i),m1dat(i)&
                                            ! &,n2dat(i),m2dat(i),n3dat(i),m3dat(i)
                                         end if
                                      end do
                                   end do
                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
           
           nmpstate=i

           return
         end subroutine mpstatesfill3p
           
!!!!!!!!!!!! the second part is the basis enumeration for the four particle Fock sector qqqg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
        subroutine mpstatesfill4p(nmax2,mj,kt,s1dat,s2dat,s3dat,s4dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,n4dat,m4dat,&
          &k1dat,k2dat,k3dat,k4dat,nmpstate)
           implicit none
           integer nmax2,mj,kt,nmpstate
           integer,dimension(*) :: s1dat,s2dat,s3dat,s4dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,n4dat,m4dat,k1dat,&
           &k2dat,k3dat,k4dat
           integer i,k1,k2,k3,k4,n1,m1,n2,m2,n3,m3,n4,m4,s1index,s2index,s3index,s4index,sum1,sum2,sum3,sum4,sumtot
           integer ::m

           i=0

            do k2=0,kt-2
              do k1=0,kt-2-k2
                do k3=0,kt-2-k2-k1
                  k4=kt-k1-k2-k3-1
                  do s4index=-1,1,2
                    do s3index=-1,1,2
                      do s2index=-1,1,2
                        do s1index=-1,1,2
                          m=(mj-s1index-s2index-s3index)/2-s4index
                          do n4=0,floor(dble(nmax2)/2)
                            do n3=0,floor(dble(nmax2)/2)-n4
                              do n2=0,floor(dble(nmax2)/2)-n4-n3
                                do n1=0,floor(dble(nmax2)/2)-n4-n3-n2
                                  do m3=-floor(dble(nmax2-4-2*n1-2*n2-2*n3-2*n4-m)/2)&
                                        &,ceiling(dble(nmax2-4-2*n1-2*n2-2*n3-2*n4+m)/2)
                                    do m2=-floor(dble(nmax2-4-2*n1-2*n2-2*n3-2*n4-m+m3)/2)&
                                           &,ceiling(dble(nmax2-4-2*n1-2*n2-2*n3-2*n4+m-m3)/2)
                                      do m1=-floor(dble(nmax2-4-2*n1-2*n2-2*n3-2*n4-m+m2+m3)/2)&
                                           &,ceiling(dble(nmax2-3-2*n1-2*n2-2*n3-2*n4+m-m2-m3)/2)

                                         
                                         m4=m-m1-m2-m3
                                         sum1=2*n1+abs(m1)+1
                                         sum2=2*n2+abs(m2)+1
                                         sum3=2*n3+abs(m3)+1
                                         sum4=2*n4+abs(m4)+1
                                         sumtot=sum1+sum2+sum3+sum4
                                         if(sumtot.le.nmax2) then
                                            i=i+1
                                            s1dat(i)=s1index
                                            s2dat(i)=s2index
                                            s3dat(i)=s3index
                                            s4dat(i)=s4index
                                            n1dat(i)=n1
                                            m1dat(i)=m1
                                            n2dat(i)=n2
                                            m2dat(i)=m2
                                            n3dat(i)=n3
                                            m3dat(i)=m3
                                            n4dat(i)=n4
                                            m4dat(i)=m4
                                            k1dat(i)=k1
                                            k2dat(i)=k2
                                            k3dat(i)=k3
                                            k4dat(i)=k4
                                            ! Print*, k1dat(i),k2dat(i),k3dat(i),k4dat(i),s1dat(i),s2dat(i),s3dat(i),s4dat(i),&
                                            ! &n1dat(i),m1dat(i),n2dat(i),m2dat(i),n3dat(i),m3dat(i),n4dat(i),m4dat(i)
                                         end if
                                      end do
                                    end do
                                  end do
                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
           
           nmpstate=i

           return
         end subroutine mpstatesfill4p

!!!!!!!!!!!! the third one is the basis enumeration for five particle Fock sector qqqqqbar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
        subroutine mpstatesfill5p(nmax2,mj,kt,s1dat,s2dat,s3dat,s4dat,s5dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,n4dat,m4dat,&
          &n5dat,m5dat,k1dat,k2dat,k3dat,k4dat,k5dat,nmpstate)
           implicit none
           integer nmax2,mj,kt,nmpstate
           integer,dimension(*) :: s1dat,s2dat,s3dat,s4dat,s5dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,n4dat,m4dat,k1dat,&
           &k2dat,k3dat,k4dat,n5dat,m5dat,k5dat
           integer i,k1,k2,k3,k4,k5,n1,m1,n2,m2,n3,m3,n4,m4,n5,m5,s1index,s2index,s3index,s4index,s5index,sum1,sum2,sum3,&
           &sum4,sum5,sumtot
           integer :: m

           i=0

          do k2=0,kt-2
            do k1=0,kt-2-k2
              do k3=0,kt-2-k2-k1
                do k4=0,kt-2-k3-k2-k1
                  k5=kt-k1-k2-k3-k4-2
                  do s5index=-1,1,2
                    do s4index=-1,1,2
                      do s3index=-1,1,2
                        do s2index=-1,1,2
                          do s1index=-1,1,2
                            m=(mj-s1index-s2index-s3index-s4index-s5index)/2
                            do n5=0,floor(dble(nmax2-5)/2)
                              do n4=0,floor(dble(nmax2-5)/2)-n5
                                do n3=0,floor(dble(nmax2-5)/2)-n5-n4
                                  do n2=0,floor(dble(nmax2-5)/2)-n5-n4-n3
                                    do n1=0,floor(dble(nmax2-5)/2)-n5-n4-n3-n2
                                      do m5=-floor(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5-m)/2)&
                                        &,ceiling(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5+m)/2)
                                        do m4=-floor(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m5)/2)&
                                           &,ceiling(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m5)/2)
                                          do m3=-floor(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m5+m4)/2)&
                                           &,ceiling(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m5-m4)/2)
                                            do m2=-floor(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5-m+m5+m4+m3)/2)&
                                              &,ceiling(dble(nmax2-5-2*n1-2*n2-2*n3-2*n4-2*n5+m-m5-m4-m3)/2)

                                              m1=m-m2-m3-m4-m5
                                              sum1=2*n1+abs(m1)+1
                                              sum2=2*n2+abs(m2)+1
                                              sum3=2*n3+abs(m3)+1
                                              sum4=2*n4+abs(m4)+1
                                              sum5=2*n5+abs(m5)+1
                                              sumtot=sum1+sum2+sum3+sum4+sum5
                                              if(sumtot.le.nmax2) then
                                                i=i+1
                                                s1dat(i)=s1index
                                                s2dat(i)=s2index
                                                s3dat(i)=s3index
                                                s4dat(i)=s4index
                                                s5dat(i)=s5index
                                                n1dat(i)=n1
                                                m1dat(i)=m1
                                                n2dat(i)=n2
                                                m2dat(i)=m2
                                                n3dat(i)=n3
                                                m3dat(i)=m3
                                                n4dat(i)=n4
                                                m4dat(i)=m4
                                                n5dat(i)=n5
                                                m5dat(i)=m5
                                                k1dat(i)=k1
                                                k2dat(i)=k2
                                                k3dat(i)=k3
                                                k4dat(i)=k4
                                                k5dat(i)=k5
                                                ! Print*, k1dat(i),k2dat(i),k3dat(i),k4dat(i),k5dat(i),&
                                                ! &s1dat(i),s2dat(i),s3dat(i),s4dat(i),s5dat(i),n1dat(i),m1dat(i),&
                                                ! &n2dat(i),m2dat(i),n3dat(i),m3dat(i),n4dat(i),m4dat(i),n5dat(i),m5dat(i)
                                              end if
                                            end do
                                          end do
                                        end do
                                      end do
                                    end do
                                  end do
                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
           
           nmpstate=i

           return
        end subroutine mpstatesfill5p                       
           
                 
                             
