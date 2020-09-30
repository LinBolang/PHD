       module basis_info
         use numbers
         implicit none
         save

        !  integer :: nmax
         integer,dimension(:), allocatable :: s1dat,s2dat,s3dat
         integer,dimension(:), allocatable :: n1dat,m1dat,n2dat,m2dat,n3dat,m3dat
         integer,dimension(:), allocatable :: k1dat,k2dat,k3dat
         integer,dimension(:), allocatable :: s1dat1,s2dat1,s3dat1
         integer,dimension(:), allocatable :: n1dat1,m1dat1,n2dat1,m2dat1,n3dat1,m3dat1
         integer,dimension(:), allocatable :: k1dat1,k2dat1,k3dat1
         integer,dimension(:), allocatable :: s1dat2,s2dat2,s3dat2,s4dat2
         integer,dimension(:), allocatable :: n1dat2,m1dat2,n2dat2,m2dat2,n3dat2,m3dat2,n4dat2,m4dat2
         integer,dimension(:), allocatable :: k1dat2,k2dat2,k3dat2,k4dat2
         integer,dimension(:), allocatable :: s1dat3,s2dat3,s3dat3,s4dat3,s5dat3
         integer,dimension(:), allocatable :: n1dat3,m1dat3,n2dat3,m2dat3,n3dat3,m3dat3,n4dat3,m4dat3,n5dat3,m5dat3
         integer,dimension(:), allocatable :: k1dat3,k2dat3,k3dat3,k4dat3,k5dat3
         integer :: dimtot1,dimtot2,dimtot3,dimtot,seqeigen,select_distribution_function
         double precision :: ktot,bee
         double precision,dimension(:), allocatable :: masscountertermpool,zfactorpool
         double precision, dimension(:,:,:),allocatable :: integration,integration3t5
         real :: start_time,integrate_time,hamiltonian_time,time3
         double precision :: coupling3,coupling5,coupling3t5
         double precision,dimension(:),allocatable :: deltamq,deltamg

       end module basis_info
