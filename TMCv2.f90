!
!
 !----------------------------------------------------------------------------
 !             Talmi-Moshinsky Transformation Coefficients                   ! 
 !                            for                                            !
 !                  2D Harmonic Oscillator Basis                             !
 !                                                                           !
 !    Yang Li <leeyoung@iastate.edu>, version 1, Nov. 18, 2013               !
 !                                    version 2, Feb. 22, 2015               !
 !    Contributors: Xingbo Zhao, Paul Wiecki, Guangyao Chen                  !
 !--                                                                         !
 !  Talmi-Moshinsy (TM) coefficients are the Clebsch-Gordan ceofficients for !
 !  the 2-dimensional harmonic oscillator functions, which can be used to    !
 !  separate the center-of-mass motion. See tmc.pdf for more details.        !
 !                                                                           !
 !- Definition of the TM transformation:                                     !
 !  Let psi(n,m, p) be the normalized 2D HO functions. Then TM transformation!
 !  relates the HO function in the single-particle coordinate with the       !
 !  relative coordinate:                                                     ! 
 !                                                                           !
 !  psi(n1,m1,p1)*psi(n2,m2,p2) =                                            !
 !      sum_{N,M,n,m} TMC(N,M,n,m;n1,m1,n2,m2;delta) * psi(N,M,P)*psi(n,m,p) !
 !  where P = cos(delta)*p1 + sin(delta)*p2 [vector],                        !
 !        p = sin(delta)*p1 - cos(delta)*p2 [vector].                        !
 !                                                                           !
 !  (By doing the TM transformation reccursively, one  can relate the single-!
 !  particle coordinate HO function to the Jacobi coordinate HO function.)   !
 !                                                                           !
 !- The HO function, expressed in the momentum space, reads,                 !
 !  b*psi(n,m,p) = sqrt(4*pi*n!/(n+|m|)!) * exp(i*m*phi-rho^2/2) * rho^|m| * !
 !                LaguerreL(n, |m|, rho^2)                                   !
 !   where LaguerreL(n, a, x) is the generalized Laguerre polynomial.        !
 !   rho = |p|/b, phi = arg p, b is the basis energy scale.                  !
 !   [recall for H = p^2/(2*m) + 1/2*m*w^2*r^2, b = sqrt(m*w), so            !
 !    in the transformation, B = sqrt((m1+m2)*w), b = sqrt(m1*m2/(m1+m2)*w)] !
 !                                                                           !
 !  psi(n,m,p) psisatisfy the orthonormality relations:                      !
 !                                                                           !
 !  int d^2 p/(2*pi)^2 psi(n,m,p) psi^*(n',m',p) = delta(n,n') delta(m,m')   !
 !  sum_{n,m} psi(n,m,p) psi(n,m,p') -> (2*pi)^2*delta^2(p-p')               !
 !                                                                           !
 !  where "^*" means the complex conjugate.                                  ! 
 !                                                                           !
 !- Properties of the TM coefficients (TMC):                                 !
 !                                                                           !
 !  * TMC is proportional to delta(2*n1+|m1|+2*n2+|m2|, 2*N+|M|+2*n+|m|);    !
 !  * TMC is proportional to delta(m1+m2, M+m);                              !
 !                                                                           !
 !--                                                                         !
 !  The new algorithm is based Guangyao Chen's method of identifying the     !
 !  terms in the generating functions.                                       ! 
 !                                                                           !
 !--                                                                         !
 ! This code is tested against the old algorithm up to Nmax=32. See          !
 ! new_tmc_timing.pdf new_tmc_normalization.pdf new_tmc_discrepancy_w_old.pdf!
 ! for details.                                                              !
 !                                                                           !
 !---------------------------------------------------------------------------!

    double precision function TMC(nn, mm, n, m, n1, m1, n2, m2, tandelta)
        implicit none
  !
  ! TMC(nn, mm, n, m; n1, m1, n2, m2; b2/b1)
  ! TMC(nn, mm, n, m; n1, m1, n2, m2; sqrt(x2/x1))
  !
  ! we actuall only need Binomial coefficient for n>0,0<=m<=n here
        integer, intent(in) :: nn, mm, n, m, n1, m1, n2, m2
        double precision, intent(in) :: tandelta
        integer :: a, b, EEp, EEm, Ep, Em, E1p, E1m
        double precision :: sindelta, cosdelta, s1, s2, t
        double precision, external :: Binomial, lognm
!
        if(2*nn+abs(mm)+2*n+abs(m)==2*n1+abs(m1)+2*n2+abs(m2).and.mm+m==m1+m2) then
!
            sindelta = tandelta/sqrt(1d0+tandelta**2)
            cosdelta = 1d0/sqrt(1d0+tandelta**2)
            t = -1d0/(tandelta*tandelta)

            EEp = nn + (abs(mm)+mm)/2
            EEm = nn + (abs(mm)-mm)/2
            Ep  = n  + (abs(m)+m)/2
            Em  = n  + (abs(m)-m)/2
            E1p = n1 + (abs(m1)+m1)/2
            E1m = n1 + (abs(m1)-m1)/2

            s1 = 0d0
            do a = max(0,E1m-Em), min(E1m,EEm)
                s1 = s1 + t**a * binomial(EEm,a) * binomial(Em,E1m-a) 
            enddo
            s2 = 0d0
            do b = max(0,E1p-Ep), min(E1p,EEp)
                s2 = s2 + t**b * binomial(EEp,b) * binomial(Ep,E1p-b) 
            enddo

            TMC = (1-2*mod(abs(nn+n+n1+n2+m+m1),2))*s1*s2*tandelta**(2*n1+abs(m1)) &
                * sindelta**(2*nn+abs(mm)) * cosdelta**(2*n+abs(m))                &
                * exp(lognm(n1,m1)+lognm(n2,m2)-lognm(nn,mm)-lognm(n,m))

        else 
            TMC = 0d0
        endif
! 
    end function
 
!==============================================================================
!   Auxiliary functions: factorials, binomials, multinomials etc              !
!==============================================================================
 
    Function LogNM(n, m)
    implicit none
! 	log(sqrt((n+abs(m))!n!))
    
    double precision :: logNM, logn
    integer :: n, m, i

    logn = 0D0;
    do i = 2, n
        logn = logn + log(dble(i));
    end do

    logNM = logn;
    do i = n+1, n+abs(m)
        logNM = logNM + log(dble(i));
    end do 
    
    logNM = 0.5D0 * (logNM + logn);
    
    End Function logNM

    double precision Function Binomial(n,m)
    implicit none
    ! binomial coefficients, n >= 0, m >= 0, and n >= m 
    ! Binomial(n, m) = n!/(m! * (n-m)!); 
    ! Binomial(n, m) = Binomial(n, n-m);
    integer, intent(in) :: n, m
    integer :: k, i
        
    k = min(m,n-m)
    Binomial = 1D0
    do i =  1, k
        Binomial = Binomial*dble(n-k+i)/dble(i)
    end do
               
    End Function

    !Function LogBinomial(n,m)
    !implicit none
    !! generalized binomial coefficients
    !! We follow the definition of Kronenburg ( arXiv:1105.3689v1 [math.CO], 2011)
    !!
    !! for n >=0, 
    !!	if 0 <= m <= n, Binomial(n, m) = n!/(m! * (n-m)!); 
    !!	otherwise Binomial(n, m) = 0;
    !! for n < 0, 
    !!	if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    !!	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
    !!	otherwise Binomial(n, m) = 0; 
    !!
    !! the symmetry Binomial(n, m) = Binomial(n, n-m) remains valid.
    !!
    !! the binomial theorem remains valid:
    !! ( x + y )^n = sum_{k=0}^{infty} Binomial(n, k) x^k y^{n-k};
    !!
    !! for n >=0, the summation terminates at k = n;
    !! for n < 0, the summation keeps going to infinity.
    !! of course, if n < 0, the above series only converges when |x| < |y|,
    !! when |x| > |y|, we should use:
    !! ( x + y )^n = sum_{k=0}^{infty} Binomial(n, n-k) x^{n-k} y^k 
    !! 		  = sum_{k'=-infty}^{n} Binomial(n,k') x^k' y^{n-k'}
    !! 
    !! This function evaluate log(|Binomial(n, m)|); namely, the sign of 
    !! binomial(n,m) is dropped;
    !! sign is given by a separate function BinomialSign(n, m);
    !! 
    !! tested via mathematica Binomial[] function, up to n, m = -100, 500.
    !!
    !! Young Li, Feb. 2, 2012, Groundhog Day. Phil saw his shadow, six more weeks
    !! winter.
            
        !double precision :: LogBinomial
        !integer :: n, m, k, ii, jj
!!        double precision, parameter :: iota = -1D12;
        
        !LogBinomial = 0D0;
        !if ( n >= 0 ) then
			!if( m >= 0 .and. m <= n) then
				
        		!do ii =  1, m
            		!LogBinomial = LogBinomial + 					&
            	!&		log(dble(n - m + ii)) - log(dble(ii));
        		!end do
				
			!else
				!return;
			!end if
        !else 
        
    !!	n < 0, if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    		!if ( m >= 0 ) then
    		
        		!do ii =  1, m
            		!LogBinomial = LogBinomial + 					&
            	!&		log(dble( ii - 1 - n )) - log(dble(ii));
        		!end do
   
    !!	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
        	!else if ( m <= n ) then
        		
        		!do ii = 1, -1 - n
        			!LogBinomial = LogBinomial + 					&
        		!&		log(dble( n - m + ii )) - log(dble(ii));
        		!end do 

    !!	otherwise Binomial(n, m) = 0;          		        		
        	!else
        		!return;
        	!end if
        
        !end if
        
    !End Function

	!Function BinomialSign(n, m)
	!implicit none
    !! sign of generalized binomial coefficients
    !! We follow the definition of Kronenburg ( arXiv:1105.3689v1 [math.CO], 2011)
    !!
    !! for n >=0, 
    !!	if 0 <= m <= n, Binomial(n, m) = n!/(m! * (n-m)!); 
    !!	otherwise Binomial(n, m) = 0;
    !! for n < 0, 
    !!	if m >=0, Binomial(n, m) = (-1)^m * Binomial(-n+m-1,m);
    !!	if m <=n, Binomial(n, m) = (-1)^(n-m) * Binomial(-m-1, n-m);
    !!	otherwise Binomial(n, m) = 0; 	
    !! see comments in LogBinomial() for more information.
    !integer :: n, m, BinomialSign;
    
    !if( n >= 0 ) then 
    	
    	!if( m >= 0 .and. m <= n ) then
    		
    		!BinomialSign = 1;
    		
    	!else
    		!BinomialSign = 0;
    	!end if
    !else ! n < 0
    	
    	!if ( m >= 0 ) then 
    		
    		!BinomialSign = (-1)**(mod(m,2));
    	
    	!else if ( m <= n) then
    		
    		!BinomialSign = (-1)**(mod(n-m,2));
    		
    	!else
    		
    		!BinomialSign = 0;
    	!end if
    	
    !end if
    
    !End Function            
    

