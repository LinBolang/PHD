subroutine mpstatesfillquark(nmax2,mj,kt,method,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat,nmpstate)
    implicit none
    integer nmax2,mj,kt,nmpstate,method
    integer,dimension(*) :: s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat
    integer i,k1,k2,k3,n1,m1,n2,m2,n3,m3,s1index,s2index,s3index,sum1,sum2,sum3,sumtot,color1,color2,color3

    i=0
    !!!!!!!! For the one particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    k1=Kt
    s1index=Mj 
    n1=0
    m1=(mj-s1index)/2
    color1=1
    sum1=2*n1+abs(m1)+1
    sumtot=sum1
    if(sumtot.le.nmax2) then
        i=i+1
        s1dat(i)=s1index
        s2dat(i)=0
        s3dat(i)=0
        n1dat(i)=n1
        m1dat(i)=m1
        n2dat(i)=0
        m2dat(i)=0
        n3dat(i)=0
        m3dat(i)=0
        k1dat(i)=k1
        k2dat(i)=0
        k3dat(i)=0
        ! Print*, k1dat(i),k2dat(i),k3dat(i),s1dat(i),s2dat(i),s3dat(i),n1dat(i),m1dat(i)&
        ! &,n2dat(i),m2dat(i),n3dat(i),m3dat(i)
    end if

    !!!!!!!!! For the two particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	If (method.eq.2) then
        do color2=1,2
		    do k1=0,kt-1
	            k2=Kt-k1
		        do s2index=-1,1,2
		        	do s1index=-1,1,2
		        		do n2=0,floor(dble(nmax2-2)/2)
		        			do n1=0,floor(dble(nmax2-2)/2)-n2
		        				do m2=-floor(dble(nmax2-2-2*n1-2*n2-((mj-s1index)/2-s2index))/2)&
		                            &,ceiling(dble(nmax2-2-2*n1-2*n2+(dble(mj-s1index)/2-s2index))/2)

		                            m1=nint(dble(mj-s1index)/2-s2index)-m2
		                            sum1=2*n1+abs(m1)+1
		                            sum2=2*n2+abs(m2)+1
		                            sumtot=sum1+sum2
		                            if(sumtot.le.nmax2) then
		                                i=i+1
		                                s1dat(i)=s1index
		                                s2dat(i)=s2index
		                                s3dat(i)=0
		                                n1dat(i)=n1
		                                m1dat(i)=m1
		                                n2dat(i)=n2
		                                m2dat(i)=m2
		                                n3dat(i)=0
		                                m3dat(i)=0
		                                k1dat(i)=k1
		                                k2dat(i)=k2
		                                k3dat(i)=0
		                                ! Print*, k1dat(i),k2dat(i),k3dat(i),s1dat(i),s2dat(i),s3dat(i),n1dat(i),m1dat(i)&
		                                ! &,n2dat(i),m2dat(i),n3dat(i),m3dat(i)
		                            end if
		                        enddo
		        			enddo
		        		enddo
		        	enddo
		        enddo
		    enddo
		enddo
	endif

    !!!!!!!!!!! For the three particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do color3=1,3
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
	                                    	enddo
	                                    end do
	                                end do
	                            end do
	                        end do
	                    end do
	                end do
	            end do
	        end do
	    end do
	enddo
           
           
    nmpstate=i

    return
end subroutine mpstatesfillquark

subroutine dimtotal1p(Nmax,Mj,Kt,dimtot)
    implicit none
    integer :: nmax,Mj,kt,dimtot

    dimtot=1

    return
end subroutine dimtotal1p

subroutine dimtotal2p(Nmax,Mj,Kt,dimtot)
    implicit none
    integer :: nmax,Mj,Kt,dimtot
    integer,dimension(5) :: SumSpin(5)

    call SumS2p(Nmax,Mj,SumSpin)

    dimtot=Kt*SumSpin(5)

    return
end subroutine dimtotal2p

subroutine SumS2p(Nmax,Mj,SumS)
    implicit none
    integer :: Nmax,Mj,nmax0
    integer,dimension(5) :: SumS

    SumS(1)=0
    nmax0=nmax-mod(nmax,2)

    SumS(2)=SumS(1)+((-2 + nmax0)*nmax0*(5 + nmax0))/24
    SumS(3)=SumS(2)+nmax0*(-4 + nmax0**2)/24
    SumS(4)=SumS(3)+nmax0*(2 + 3*nmax0 + nmax0**2)/24
    SumS(5)=SumS(4)+((-2 + nmax0)*nmax0*(2 + nmax0))/24

    return
end subroutine SumS2p

subroutine mpstatesfillPhoton(nmax2,mj,kt,method,s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat,nmpstate)
    implicit none
    integer nmax2,mj,kt,nmpstate,method
    integer,dimension(*) :: s1dat,s2dat,s3dat,n1dat,m1dat,n2dat,m2dat,n3dat,m3dat,k1dat,k2dat,k3dat
    integer i,k1,k2,k3,n1,m1,n2,m2,n3,m3,s1index,s2index,s3index,sum1,sum2,sum3,sumtot,color1,color2,color3

    i=0
    !!!!!!!! For the one particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    k1=Kt
	    s1index=Mj 
	    n1=0
	    m1=(mj-s1index)
	    color1=1
	    sum1=2*n1+abs(m1)+1
	    sumtot=sum1
	    if(sumtot.le.nmax2) then
	        i=i+1
	        s1dat(i)=s1index
	        s2dat(i)=0
	        s3dat(i)=0
	        n1dat(i)=n1
	        m1dat(i)=m1
	        n2dat(i)=0
	        m2dat(i)=0
	        n3dat(i)=0
	        m3dat(i)=0
	        k1dat(i)=k1
	        k2dat(i)=0
	        k3dat(i)=0
	        ! Print*, k1dat(i),k2dat(i),k3dat(i),s1dat(i),s2dat(i),s3dat(i),n1dat(i),m1dat(i)&
	        ! &,n2dat(i),m2dat(i),n3dat(i),m3dat(i)
	    end if

    !!!!!!!!! For the two particle Fock sector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !do color2=1,3
	    do k1=0,kt-1
            k2=Kt-1-k1
	        do s1index=-1,1,2
	        	do s2index=-1,1,2
	        		do n2=0,floor(dble(nmax2-2)/2)
	        			do n1=0,floor(dble(nmax2-2)/2)-n2
	        				do m1=-floor(dble(nmax2-2-2*n1-2*n2-((mj-s1index-s2index)/2))/2)&
	                            &,ceiling(dble(nmax2-2-2*n1-2*n2+(dble(mj-s1index-s2index)/2))/2)

	                            m2=nint(dble(2*mj-s1index-s2index)/2)-m1
	                            sum1=2*n1+abs(m1)+1
	                            sum2=2*n2+abs(m2)+1
	                            sumtot=sum1+sum2
	                            if(sumtot.le.nmax2) then
	                                i=i+1
	                                s1dat(i)=s1index
	                                s2dat(i)=s2index
	                                s3dat(i)=0
	                                n1dat(i)=n1
	                                m1dat(i)=m1
	                                n2dat(i)=n2
	                                m2dat(i)=m2
	                                n3dat(i)=0
	                                m3dat(i)=0
	                                k1dat(i)=k1
	                                k2dat(i)=k2
	                                k3dat(i)=0
	                                ! Print*, k1dat(i),k2dat(i),k3dat(i),s1dat(i),s2dat(i),s3dat(i),n1dat(i),m1dat(i)&
	                                ! &,n2dat(i),m2dat(i),n3dat(i),m3dat(i)
	                            end if
	                        enddo
	        			enddo
	        		enddo
	        	enddo
	        enddo
	    enddo
	!enddo           
           
    nmpstate=i

    return
end subroutine mpstatesfillPhoton

subroutine dimtotal2q(Nmax,Mj,Kt,dimtot)
    implicit none
    integer :: nmax,Mj,Kt,dimtot
    integer,dimension(5) :: SumSpin(5)

    call SumS2q(Nmax,Mj,SumSpin)

    dimtot=Kt*SumSpin(5)

    return
end subroutine dimtotal2q

    subroutine sums2q(Nmax,Mj,SumS)
        use basis_info
        implicit none
    integer :: Nmax, Mj
    integer, dimension(5) :: SumS ! the matrix length equal to 2^2+1
    integer, dimension(3)  :: dimNM !record different spin combination corresponding to the number of N&M combination
    integer :: s1index,s2index
    integer :: i,m

    i=1
    SumS(1)=0
    call Sumnm2q(Nmax,Mj,dimNM)

          do s2index=-1,1,2
            do s1index=-1,1,2
                i=i+1
                m=mj-(s1index+s2index)/2
                SumS(i)=SumS(i-1)+dimNM(m-(Mj-2))
            enddo
          enddo

        return
    end subroutine sums2q

    subroutine Sumnm2q(Nmax,Mj,dim)
        implicit none

    integer :: Nmax,Mj
    integer, dimension(3) :: dim ! the dimension come from 1-(-1)+1=3
    integer :: n1,m1,n2,m2,n3,m3,n4,m4
    integer :: i,j,m,sumN

    do j=1,3

        i=0
        m=(Mj-1)+j-1
        dim(j)=0

            do n2=0,floor(dble(nmax-2)/2.0D0)
              do n1=0,floor(dble(nmax-2)/2.0D0)-n2
                do m1=-floor(dble(nmax-2-2*n1-2*n2-m)/2.0D0)&
                   &,ceiling(dble(nmax-2-2*n1-2*n2+m)/2.0D0)

                    m2=m-m1
                    SumN=2*(n1+n2)+abs(m1)+abs(m2)+2

                    if(SumN.le.Nmax) then
                        i=i+1
                    endif
                enddo
              enddo
            enddo


      dim(j)=i

    enddo


        return
    end subroutine Sumnm2q