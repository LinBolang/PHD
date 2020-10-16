module numbers
  implicit none
  save
  integer, parameter:: LI = SELECTED_INT_KIND(12)
  integer:: MAX_RENORM_ITER = 100
  double precision, parameter:: PI = acos(-1.D0), ARPACK_TOL = 1.D-10, MASS_RENORM_TOL = 1.D-10, NMCUTOFFWIDTH = 1.0D0
  double precision, parameter:: LAGRANGEMULTIPLIER = 10.1D0, B_DEFAULT = 0.600D0, binst_default = 100.6D0
  double precision, parameter:: Mu1_DEFAULT = 0.411D0, Md1_DEFAULT = 0.411D0, MS1_DEFAULT = 0.411D0, mg1_default = 0.05D0
  double precision, parameter:: Mu2_DEFAULT = 0.105D0, Md2_DEFAULT = 0.105D0, MS2_DEFAULT = 0.105D0, mg2_default = 0.05D0
  double precision, parameter:: Mu3_DEFAULT = 0.301D0, Md3_DEFAULT = 0.301D0, MS3_DEFAULT = 0.301D0, mg3_default = 0.05D0
  double precision, parameter:: Mu4_DEFAULT = 0.411D0, Md4_DEFAULT = 0.411D0, MS4_DEFAULT = 0.411D0, mg4_default = 0.05D0
  double precision, parameter:: M3t5_DEFAULT = 0.00D0
  double precision, parameter:: kappat1 = 1.00D0, kappal1 = 1.00D0, kappat3 = 0.0d0, kappal3 = 0.0d0
  double precision, parameter:: kappa4 = 1.0D0, coupling_default3 = 2.50D0, coupling_default5 = 1.70D0, &
    & coupling_default3t5 = 0.23D0
end module numbers
