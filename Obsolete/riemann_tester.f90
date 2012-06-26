program riemann_tester
use Streamer
implicit none
logical :: verbose = .false.
real(8), dimension(21) :: oneL, oneR, twoL, twoR, threeL, threeR,fourL, fourR, fiveL, fiveR
real(8), dimension(5) :: test
real(8), dimension(21) :: testervec
real(8), dimension(3,3) :: matrix

oneL = 0.D0 ; twoL = 0.D0 ; threeL = 0.D0 ; fourL = 0.D0 ; fiveL = 0.D0
oneR = 0.D0 ; twoR = 0.D0 ; threeR = 0.D0 ; fourR = 0.D0 ; fiveR = 0.D0

  oneL(2) = 1.0D0     ;   oneL(1) =    1.0D0   ;   oneL(3) = 0.0D0     ;   oneL(4) = 0.0D0     ;   oneL(5) =    0.0D0
  oneR(2) = 0.125D0   ;   oneR(1) =    0.1D0   ;   oneR(3) = 0.0D0     ;   oneR(4) = 0.0D0     ;   oneR(5) =    0.0D0
  twoL(2) = 1.0D0     ;   twoL(1) =    0.4D0   ;   twoL(3) =-2.0D0     ;   twoL(4) = 0.0D0     ;   twoL(5) =    0.0D0
  twoR(2) = 1.0D0     ;   twoR(1) =    0.4D0   ;   twoR(3) = 2.0D0     ;   twoR(4) = 0.0D0     ;   twoR(5) =    0.0D0
threeL(2) = 1.0D0     ; threeL(1) = 1000.0D0   ; threeL(3) = 0.0D0     ; threeL(4) = 0.0D0     ; threeL(5) =    0.0D0
threeR(2) = 1.0D0     ; threeR(1) =    0.01D0  ; threeR(3) = 0.0D0     ; threeR(4) = 0.0D0     ; threeR(5) =    0.0D0
 fourL(2) = 1.0D0     ;  fourL(1) =    0.01D0  ;  fourL(3) = 0.0D0     ;  fourL(4) = 0.0D0     ;  fourL(5) =    0.0D0
 fourR(2) = 1.0D0     ;  fourR(1) =  100.0D0   ;  fourR(3) = 0.0D0     ;  fourR(4) = 0.0D0     ;  fourR(5) =    0.0D0
 fiveL(2) = 5.99924D0 ;  fiveL(1) = 460.894D0  ;  fiveL(3) =19.5975D0  ;  fiveL(4) = 0.0D0     ;  fiveL(5) =    0.0D0
 fiveR(2) = 5.99242D0 ;  fiveR(1) =  46.0950D0 ;  fiveR(3) =-6.19633D0 ;  fiveR(4) = 0.0D0     ;  fiveR(5) =    0.0D0

test = riemann_solve(  oneL,   oneR, verbose_flag=verbose, t_out=.25d0 )
write(*,*) test
test = riemann_solve(  twoL,   twoR, verbose_flag=verbose, t_out=.15d0 )
write(*,*) test
test = riemann_solve(threeL, threeR, verbose_flag=verbose, t_out=.012d0)
write(*,*) test
test = riemann_solve( fourL,  fourR, verbose_flag=verbose, t_out=.035d0)
write(*,*) test
test = riemann_solve( fiveL,  fiveR, verbose_flag=verbose, t_out=.035d0)
write(*,*) test


testervec = dble([1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0])
call grid_coords(dble([1,1,0]),matrix(1,:),matrix(2,:),matrix(3,:))
call vels_transform(testervec,dble(reshape([1,0,0,0,1,0,0,0,1],[3,3])))
end program riemann_tester
