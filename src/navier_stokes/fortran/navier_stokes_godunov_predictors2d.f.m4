dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the minmod of two values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function minmod(a,b)
      implicit none
      REAL a,b
      minmod = 0.5d0*(sign(1.d0,a)+sign(1.d0,b))*min(abs(a),abs(b))
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the median of three values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function median(a,b,c)
      implicit none
      REAL a,b,c
      REAL minmod
      median = a + minmod(b-a,c-a)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the WENO5 interpolation of several values.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function WENO5_interp(Q)
      implicit none
      REAL Q(-2:2)
      REAL f(0:2)
      REAL IS(0:2)
      REAL omega_bar(0:2)
      REAL omega(0:2),omega_sum
      REAL alpha(0:2),alpha_sum
      INTEGER i
c     Compute the candidate interpolations.
      f(0) = (11.d0*Q( 0)-7.d0*Q(-1)+2.d0*Q(-2))/6.d0
      f(1) = ( 2.d0*Q(+1)+5.d0*Q( 0)-     Q(-1))/6.d0
      f(2) = (-1.d0*Q(+2)+5.d0*Q(+1)+2.d0*Q( 0))/6.d0
c     Compute the smoothness indicators.
      IS(0) = (13.d0/12.d0)*((Q( 0)-2.d0*Q(-1)+Q(-2))**2.d0) +
     &     0.25d0*((3.d0*Q( 0)-4.d0*Q(-1)+     Q(-2))**2.d0)
      IS(1) = (13.d0/12.d0)*((Q(+1)-2.d0*Q( 0)+Q(-1))**2.d0) +
     &     0.25d0*((     Q(+1)-                Q(-1))**2.d0)
      IS(2) = (13.d0/12.d0)*((Q(+2)-2.d0*Q(+1)+Q( 0))**2.d0) +
     &     0.25d0*((     Q(+2)-4.d0*Q(+1)+3.d0*Q( 0))**2.d0)
c     Compute the weights.
      omega_bar(0) = 0.1d0
      omega_bar(1) = 0.6d0
      omega_bar(2) = 0.3d0
      do i = 0,2
         alpha(i) = omega_bar(i)/(IS(i)+1.d-40)
      enddo
      alpha_sum = 0.d0
      do i = 0,2
         alpha_sum = alpha_sum + alpha(i)
      enddo
      do i = 0,2
         omega(i) = alpha(i)/alpha_sum
      enddo
c     Improve the accuracy of the weights (following the approach of
c     Henrick, Aslam, and Powers).
      do i = 0,2
         omega(i) = omega(i)*(omega_bar(i)+omega_bar(i)**2.d0
     &        -3.d0*omega_bar(i)*omega(i)+omega(i)**2.d0)/
     &        (omega_bar(i)**2.d0+omega(i)*(1.d0-2.d0*omega_bar(i)))
      enddo
      omega_sum = 0.d0
      do i = 0,2
         omega_sum = omega_sum + omega(i)
      enddo
      do i = 0,2
         omega(i) = omega(i)/omega_sum
      enddo
c     Compute the interpolant.
      WENO5_interp = 0.d0
      do i = 0,2
         WENO5_interp = WENO5_interp + omega_bar(i)*f(i)
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the sign of the input, returning zero if the absolute
c     value of x is less than a tolerance epsilon.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function sign_eps(x)
c
      implicit none
c
c     Constants.
c
      REAL EPSILON
      PARAMETER(EPSILON=1.0d-8)
c
c     Input.
c
      REAL x
c
c     Compute the sign of the input, returning zero if the absolute
c     value of x is less than a tolerance epsilon.
c
      if (dabs(x) .le. EPSILON) then
         sign_eps =  0.d0
      elseif (x  .ge. EPSILON) then
         sign_eps = +1.d0
      else
         sign_eps = -1.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     A Godunov predictor used to predict face and time centered values
c     from cell centered values using a Taylor expansion about each cell
c     center.
c
c     The predictor assumes that Q satisfies an equation of the form
c
c          dQ/dt + u * grad Q = F
c
c     i.e. Q satisfies an advection equation not in conservation form.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_godunov_predict2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,Q_L,Q_R,Q_6,Qscratch1,
     &     F,Fscratch1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0,qhalf1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx(0:NDIM-1),dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_6(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Qscratch1(ifirst1-nQgc1:ilast1+nQgc1,
     &               ifirst0-nQgc0:ilast0+nQgc0)

      REAL F(CELL2dVECG(ifirst,ilast,nFgc))
      REAL Fscratch1(ifirst1-nFgc1:ilast1+nFgc1,
     &               ifirst0-nFgc0:ilast0+nFgc0)

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     For ease of implementation, we make copies of Q and F with
c     permuted indices.
c
      do ic1 = ifirst1-nQgc1,ilast1+nQgc1
         do ic0 = ifirst0-nQgc0,ilast0+nQgc0
            Qscratch1(ic1,ic0) = Q(ic0,ic1)
         enddo
      enddo

      do ic1 = ifirst1-nFgc1,ilast1+nFgc1
         do ic0 = ifirst0-nFgc0,ilast0+nFgc0
            Fscratch1(ic1,ic0) = F(ic0,ic1)
         enddo
      enddo
c
c     Compute temporary predicted values on cell faces.
c
c     In this computation, transverse derivatives are not included.
c
      call navier_stokes_godunov_xsppm7_predict_normal2d( ! predict values on the x-faces
     &     dx(0),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,Q_L,Q_R,Q_6,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qtemp0)

      call navier_stokes_godunov_xsppm7_predict_normal2d( ! predict values on the y-faces
     &     dx(1),dt,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nQgc1,nQgc0,
     &     nFgc1,nFgc0,
     &     Qscratch1,dQ,Q_L,Q_R,Q_6,
     &     Fscratch1,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,
     &     qtemp1)
c
c     Compute final predicted values on cell faces.
c
c     This computation approximates transverse derivatives by centered
c     differences of the "temporary" predicted values.
c
      call navier_stokes_godunov_transverse_fix2d( ! update values on the x-faces
     &     dx(1),dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)

      call navier_stokes_godunov_transverse_fix2d( ! update values on the y-faces
     &     dx(0),dt,
     &     ifirst1,ilast1,ifirst0,ilast0,
     &     nugc1,nugc0,
     &     nqhalfgc1,nqhalfgc0,
     &     u1,u0,
     &     qtemp1,qtemp0,
     &     qhalf1)
c
      return
      end
      subroutine navier_stokes_godunov_plm_2nd_predict_normal2d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a PLM (piecewise linear method)
c     with unlimited second-order slopes.
c
      do ic1 = ifirst1-1,ilast1+1
         Qx = half*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         do ic0 = ifirst0-1,ilast0
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0  ,ic1)
            Qx = half*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0+1,ic1)
            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo
      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_plm_4th_predict_normal2d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qx,qL,qR
      REAL unorm
c
c     Predict face centered values using a PLM (piecewise linear method)
c     with unlimited fourth-order slopes.
c
      do ic1 = ifirst1-1,ilast1+1
         Qx = twothird*(Q(ifirst0-1+1,ic1)-Q(ifirst0-1-1,ic1))
     &        - sixth*half*(Q(ifirst0-1+2,ic1)-Q(ifirst0-1-2,ic1))
         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         do ic0 = ifirst0-1,ilast0
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0  ,ic1)
            Qx = twothird*(Q(ic0+1+1,ic1)-Q(ic0+1-1,ic1))
     &           - sixth*half*(Q(ic0+1+2,ic1)-Q(ic0+1-2,ic1))
            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0+1,ic1)
            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo
      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_plm_4th_limited_predict_normal2d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qx,qL,qR
      REAL dQQ_C,dQQ_L,dQQ_R,dQQ
      REAL unorm
c
c     Predict face centered values using a PLM (piecewise linear method)
c     with limited fourth-order slopes.
c
      do ic1 = ifirst1-1,ilast1+1
         do ic0 = ifirst0-1,ilast0+1
            dQQ_C = 0.5d0*(Q(ic0+1,ic1)-Q(ic0-1,ic1))
            dQQ_L =       (Q(ic0  ,ic1)-Q(ic0-1,ic1))
            dQQ_R =       (Q(ic0+1,ic1)-Q(ic0  ,ic1))
            if (dQQ_R*dQQ_L .gt. 0.d0) then
               dQQ = min(abs(dQQ_C),2.d0*abs(dQQ_L),2.d0*abs(dQQ_R))*
     c              sign(1.d0,dQQ_C)
            else
               dQQ = 0.d0
            endif
            dQ(ic0,ic1) = dQQ
         enddo

         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         do ic0 = ifirst0-1,ilast0
            Qx = dQ(ic0  ,ic1)
            qL = Q(ic0  ,ic1)
     &           + 0.5d0*(1.d0-unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0  ,ic1)
            Qx = dQ(ic0+1,ic1)
            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            qR = Q(ic0+1,ic1)
     &           - 0.5d0*(1.d0+unorm*dt/dx0)*Qx
     &           + 0.5d0*dt*F(ic0+1,ic1)
            qhalf0(ic0+1,ic1) =
     &           0.5d0*(qL+qR)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL-qR)
         enddo
      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_ppm_predict_normal2d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,Q_L,Q_R,Q_6,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_6(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL QQ,QQ_L,QQ_R
      REAL dQQ_C,dQQ_L,dQQ_R,dQQ
      REAL unorm,x,y
c
c     Predict face centered values using the standard PPM (piecewise
c     parabolic method).
c
      do ic1 = ifirst1-1,ilast1+1
         do ic0 = ifirst0-2,ilast0+2
            dQQ_C = 0.5d0*(Q(ic0+1,ic1)-Q(ic0-1,ic1))
            dQQ_L =       (Q(ic0  ,ic1)-Q(ic0-1,ic1))
            dQQ_R =       (Q(ic0+1,ic1)-Q(ic0  ,ic1))
            if (dQQ_R*dQQ_L .gt. 0.d0) then
               dQQ = min(abs(dQQ_C),2.d0*abs(dQQ_L),2.d0*abs(dQQ_R))*
     c              sign(1.d0,dQQ_C)
            else
               dQQ = 0.d0
            endif
            dQ(ic0,ic1) = dQQ
         enddo

         do ic0 = ifirst0-1,ilast0+1
            QQ = Q(ic0,ic1)
            QQ_L = 0.5d0*(Q(ic0-1,ic1)+Q(ic0  ,ic1)) -
     &           (1.d0/6.d0)*(dQ(ic0  ,ic1)-dQ(ic0-1,ic1))
            QQ_R = 0.5d0*(Q(ic0  ,ic1)+Q(ic0+1,ic1)) -
     &           (1.d0/6.d0)*(dQ(ic0+1,ic1)-dQ(ic0  ,ic1))
            if ((QQ_R-QQ)*(QQ-QQ_L) .le. 0.d0) then
               QQ_L = QQ
               QQ_R = QQ
            elseif ((QQ_R-QQ_L)*(QQ-0.5d0*(QQ_L+QQ_R)) .gt.
     &              +((QQ_R-QQ_L)**2.d0)/6.d0) then
               QQ_L = 3.d0*QQ-2.d0*QQ_R
            elseif ((QQ_R-QQ_L)*(QQ-0.5d0*(QQ_L+QQ_R)) .lt.
     &              -((QQ_R-QQ_L)**2.d0)/6.d0) then
               QQ_R = 3.d0*QQ-2.d0*QQ_L
            endif
            Q_L(ic0,ic1) = QQ_L
            Q_R(ic0,ic1) = QQ_R
         enddo

         do ic0 = ifirst0-1,ilast0+1
            QQ = Q(ic0,ic1)
            QQ_L = Q_L(ic0,ic1)
            QQ_R = Q_R(ic0,ic1)
            dQ(ic0,ic1) = QQ_R-QQ_L
            Q_6(ic0,ic1) = 6.d0*(QQ-0.5d0*(QQ_L+QQ_R))
         enddo

         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         do ic0 = ifirst0-1,ilast0
            y = +unorm*dt
            x = y/dx0
            QQ_L = Q_R(ic0  ,ic1) - 0.5d0*x*(
     &           dQ(ic0  ,ic1) - (1.d0-(2.d0/3.d0)*x)*Q_6(ic0  ,ic1))
     &           + 0.5d0*dt*F(ic0  ,ic1)
            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            y = -unorm*dt
            x = y/dx0
            QQ_R = Q_L(ic0+1,ic1) + 0.5d0*x*(
     &           dQ(ic0+1,ic1) + (1.d0-(2.d0/3.d0)*x)*Q_6(ic0+1,ic1))
     &           + 0.5d0*dt*F(ic0+1,ic1)
            qhalf0(ic0+1,ic1) =
     &           0.5d0*(QQ_L+QQ_R)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(QQ_L-QQ_R)
         enddo
      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_xsppm7_predict_normal2d(
     &     dx0,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     Q,dQ,Q_L,Q_R,Q_6,
     &     F,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL median,sign_eps,WENO5_interp
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nQgc0,nQgc1
      INTEGER nFgc0,nFgc1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx0,dt

      REAL Q(CELL2dVECG(ifirst,ilast,nQgc))
      REAL dQ(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_L(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_R(CELL2dVECG(ifirst,ilast,nQgc))
      REAL Q_6(CELL2dVECG(ifirst,ilast,nQgc))
      REAL F(CELL2dVECG(ifirst,ilast,nFgc))

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL QQ,QQ_L,QQ_R
      REAL QQ_star_L,QQ_star_R
      REAL QQ_WENO(-2:2)
      REAL QQ_WENO_L,QQ_WENO_R
      REAL QQ_4th_L,QQ_4th_R
      REAL dQQ_C,dQQ_L,dQQ_R,dQQ
      REAL unorm,x,y
      INTEGER i
      LOGICAL limit_value
c
c     Predict face centered values using the xsPPM7 scheme of Rider,
c     Greenough, and Kamm.
c
      do ic1 = ifirst1-1,ilast1+1
         do ic0 = ifirst0-2,ilast0+2
            dQQ_C = 0.5d0*(Q(ic0+1,ic1)-Q(ic0-1,ic1))
            dQQ_L =       (Q(ic0  ,ic1)-Q(ic0-1,ic1))
            dQQ_R =       (Q(ic0+1,ic1)-Q(ic0  ,ic1))
            if (dQQ_R*dQQ_L .gt. 0.d0) then
               dQQ = min(abs(dQQ_C),2.d0*abs(dQQ_L),2.d0*abs(dQQ_R))*
     c              sign(1.d0,dQQ_C)
            else
               dQQ = 0.d0
            endif
            dQ(ic0,ic1) = dQQ
         enddo

         do ic0 = ifirst0-1,ilast0+1
c
c     Compute a 7th order interpolation.
c
            QQ = Q(ic0,ic1)
            QQ_L = (1.0/420.d0)*(
     &           -   3.d0*Q(ic0+3,ic1)
     &           +  25.d0*Q(ic0+2,ic1)
     &           - 101.d0*Q(ic0+1,ic1)
     &           + 319.d0*Q(ic0  ,ic1)
     &           + 214.d0*Q(ic0-1,ic1)
     &           -  38.d0*Q(ic0-2,ic1)
     &           +   4.d0*Q(ic0-3,ic1))
            QQ_R = (1.0/420.d0)*(
     &           -   3.d0*Q(ic0-3,ic1)
     &           +  25.d0*Q(ic0-2,ic1)
     &           - 101.d0*Q(ic0-1,ic1)
     &           + 319.d0*Q(ic0  ,ic1)
     &           + 214.d0*Q(ic0+1,ic1)
     &           -  38.d0*Q(ic0+2,ic1)
     &           +   4.d0*Q(ic0+3,ic1))
            Q_L(ic0,ic1) = QQ_L
            Q_R(ic0,ic1) = QQ_R
c
c     Check for extrema or violations of monotonicity.
c
            QQ_L = median(QQ,QQ_L,Q(ic0-1,ic1))
            QQ_R = median(QQ,QQ_R,Q(ic0+1,ic1))

            QQ_star_L = median(QQ,QQ_L,3.d0*QQ-2.d0*QQ_R)
            QQ_star_R = median(QQ,QQ_R,3.d0*QQ-2.d0*QQ_L)

            limit_value =
     &           ((QQ_star_L-Q_L(ic0,ic1))**2.d0 .ge. 1.d-12) .or.
     &           ((QQ_star_R-Q_R(ic0,ic1))**2.d0 .ge. 1.d-12)

            if (limit_value) then
               QQ_L = Q_L(ic0,ic1) ! reset high-order values
               QQ_R = Q_R(ic0,ic1)

               do i = -2,2
                  QQ_WENO(i) = Q(ic0-i,ic1)
               enddo
               QQ_WENO_L = WENO5_interp(QQ_WENO)

               do i = -2,2
                  QQ_WENO(i) = Q(ic0+i,ic1)
               enddo
               QQ_WENO_R = WENO5_interp(QQ_WENO)

               QQ_4th_L = 0.5d0*(Q(ic0-1,ic1)+Q(ic0  ,ic1)) -
     &              (1.d0/6.d0)*(dQ(ic0  ,ic1)-dQ(ic0-1,ic1))

               QQ_4th_R = 0.5d0*(Q(ic0  ,ic1)+Q(ic0+1,ic1)) -
     &              (1.d0/6.d0)*(dQ(ic0+1,ic1)-dQ(ic0  ,ic1))

               if ( ((QQ_star_L-QQ)**2.d0 .le. 1.d-12) .and.
     &              ((QQ_star_R-QQ)**2.d0 .le. 1.d-12) ) then
c
c     Handle extrema.
c
                  QQ_WENO_L = median(QQ,QQ_WENO_L,QQ_L)
                  QQ_WENO_R = median(QQ,QQ_WENO_R,QQ_R)

                  QQ_star_L = median(QQ,QQ_WENO_L,Q(ic0-1,ic1))
                  QQ_star_R = median(QQ,QQ_WENO_R,Q(ic0+1,ic1))

                  QQ_star_L = median(
     &                 QQ,QQ_star_L,3.d0*QQ-2.d0*QQ_star_R)
                  QQ_star_R = median(
     &                 QQ,QQ_star_R,3.d0*QQ-2.d0*QQ_star_L)

                  Q_L(ic0,ic1) = median(QQ_WENO_L,QQ_star_L,QQ_L)
                  Q_R(ic0,ic1) = median(QQ_WENO_R,QQ_star_R,QQ_R)
               else
c
c     Handle steep slopes.
c
                  QQ_4th_L = median(QQ_4th_L,QQ_WENO_L,QQ_L)
                  QQ_4th_R = median(QQ_4th_R,QQ_WENO_R,QQ_R)

                  QQ_star_L = median(QQ,QQ_4TH_L,Q(ic0-1,ic1))
                  QQ_star_R = median(QQ,QQ_4TH_R,Q(ic0+1,ic1))

                  QQ_star_L = median(
     &                 QQ,QQ_star_L,3.d0*QQ-2.d0*QQ_star_R)
                  QQ_star_R = median(
     &                 QQ,QQ_star_R,3.d0*QQ-2.d0*QQ_star_L)

                  Q_L(ic0,ic1) = median(QQ_WENO_L,QQ_star_L,QQ_L)
                  Q_R(ic0,ic1) = median(QQ_WENO_R,QQ_star_R,QQ_R)
               endif
            endif
         enddo

         do ic0 = ifirst0-1,ilast0+1
            QQ = Q(ic0,ic1)
            QQ_L = Q_L(ic0,ic1)
            QQ_R = Q_R(ic0,ic1)
            dQ(ic0,ic1) = QQ_R-QQ_L
            Q_6(ic0,ic1) = 6.d0*(QQ-0.5d0*(QQ_L+QQ_R))
         enddo

         unorm = 0.5d0*(u0(ifirst0-1,ic1)+u0(ifirst0-1+1,ic1))
         do ic0 = ifirst0-1,ilast0
            y = +unorm*dt
            x = y/dx0
            QQ_L = Q_R(ic0  ,ic1) - 0.5d0*x*(
     &           dQ(ic0  ,ic1) - (1.d0-(2.d0/3.d0)*x)*Q_6(ic0  ,ic1))
     &           + 0.5d0*dt*F(ic0  ,ic1)
            unorm = 0.5d0*(u0(ic0+1,ic1)+u0(ic0+2,ic1))
            y = -unorm*dt
            x = y/dx0
            QQ_R = Q_L(ic0+1,ic1) + 0.5d0*x*(
     &           dQ(ic0+1,ic1) + (1.d0-(2.d0/3.d0)*x)*Q_6(ic0+1,ic1))
     &           + 0.5d0*dt*F(ic0+1,ic1)
            qhalf0(ic0+1,ic1) =
     &           0.5d0*(QQ_L+QQ_R)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(QQ_L-QQ_R)
         enddo
      enddo
c
      return
      end
c
      subroutine navier_stokes_godunov_transverse_fix2d(
     &     dx1,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     u0,u1,
     &     qtemp0,qtemp1,
     &     qhalf0)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Functions.
c
      REAL sign_eps
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1

      REAL dx1,dt

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qtemp0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qtemp1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Qy,qL_diff,qR_diff
      REAL vtan
c
c     Add transverse derivitives by taking centered differences of
c     temporary predicted values.
c
      do ic1 = ifirst1,ilast1
         vtan = 0.5d0*(u1(ic1,ifirst0-1)+u1(ic1+1,ifirst0-1))
         do ic0 = ifirst0-1,ilast0
            Qy = qtemp1(ic1+1,ic0)-qtemp1(ic1,ic0)
            qL_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1
            vtan = 0.5d0*(u1(ic1,ic0+1)+u1(ic1+1,ic0+1))
            Qy = qtemp1(ic1+1,ic0+1)-qtemp1(ic1,ic0+1)
            qR_diff =
     &           - 0.5d0*dt*vtan*Qy/dx1
            qhalf0(ic0+1,ic1) = qtemp0(ic0+1,ic1) +
     &           0.5d0*(qL_diff+qR_diff)+
     &           sign_eps(u0(ic0+1,ic1))*0.5d0*(qL_diff-qR_diff)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc