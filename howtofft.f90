integer nt, ell, lmax, start_index, finish_index
parameter (nt=256, lmax = 255)
integer*8 invplan, fwdplan
real*8 arr(nt), nus(nt/2+1), total_time, pi, damp
parameter(pi = ACOS(-1.0D0), damp = 1e-6)
complex*16 coeffs(nt/2+1)

call dfftw_plan_dft_c2r_1d(invplan, nt, coeffs, arr, FFTW_ESTIMATE)

! TOTAL_TIME IN SECONDS
do i=1,nt/2+1
   nus(i) = (i-1.0)/total_time
enddo


do ell= 0,lmax
  coeffs(:) = 1./(-(2.*pi*nus)**2. + ell*(ell+1)*csun**2./Rsun**2. &
      + (-1)^0.5*nus*2.*pi*damp)
  
  !SUM OVER L 

   call dfftw_execute(invplan)
enddo


call dfftw_destroy_plan(invplan)
