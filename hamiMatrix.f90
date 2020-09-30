module hamiMatrix
    implicit none
    save

        !  integer :: nmax
!!!!!!!!!! For the first Fock sector !!!!!!!!!!!!!!!!!!!!!
    double precision,dimension(:),allocatable :: hamiltoniankineticvalue,hamiltonianinteractswvalue,hamiltonianinteractlongivalue 
    double precision,dimension(:),allocatable :: hamiltonianinteractogevalue,hamiltonianinteractinstvalue
    integer,dimension(:),allocatable :: i_nzk,j_nzk,i_nztsw,j_nztsw,i_nzlsw,j_nzlsw,i_nz_oge,j_nz_oge,i_nz_inst,j_nz_inst
  !!!!!!!!!! For the second Fock sector !!!!!!!!!!!!!!!!!!!!
    double precision,dimension(:),allocatable :: hamiltoniankineticvalue2,hamiltonianinteractvcvqtqg,hamiltonianinteractvcvgtqq
    integer,dimension(:),allocatable :: i_nzk2,j_nzk2,i_nz_vcvqtqg,j_nz_vcvqtqg,i_nz_vcvgtqq,j_nz_vcvgtqq
  !!!!!!!!!! For the third Fock Sector !!!!!!!!!!!!!!!!!!!!!
    double precision,dimension(:),allocatable :: hamiltoniankineticvalue3,hamiltonianinteractswvalue3,hamiltonianinteractlongivalue3
    double precision,dimension(:),allocatable :: hamiltonianinteractogevalue3,hamiltonianinteractogevalue3t5
    double precision,dimension(:),allocatable :: hamiltonianinteractinst3t5
    integer,dimension(:),allocatable :: i_nzk3,j_nzk3,i_nztsw3,j_nztsw3,i_nzlsw3,j_nzlsw3,i_nz_oge3,j_nz_oge3,i_nz_oge3t5 &
      & ,j_nz_oge3t5,i_nz_inst3t5,j_nz_inst3t5

end module hamiMatrix