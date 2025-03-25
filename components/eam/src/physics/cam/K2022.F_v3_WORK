      module K2022
      implicit none
      private
      save
!lk+
!      integer, parameter :: r8 = selected_real_kind(12)
      public :: get_cirrus,get_cirrus_hom
      real w
      real sm,sp,sstar,ds,smid,xnstar,dummy
      real delta,xkappa
      real a1,b1,b2,xnwsat
      integer jinp,kinp
      integer inp_type(5)
      real xinp_tot(5),xinp_alpha(5)
      real xinp_b1(5),xinp_b2(5),xinp_delta(5),xinp_kappa(5)
      real b1ci,b2ci,deltaci,cnci,xkappaci
      real coolr,temp,rad,sig
!      real xwdown,xv,actg,actgint,phi,dphids,xlsub,xpsatw,xpsati
!      real xjhom,xshom,xfrz
!      external xjhom,xfrz,xv,actgint,FUNC
      common /frz/ coolr,temp,rad,sig
      common /upd/ w
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,dummy
      common /par/ delta,xkappa
      common /run/ a1,b1,b2,xnwsat
      common /inp/ jinp,inp_type,xinp_tot,xinp_alpha,
     >             xinp_delta,xinp_kappa,xinp_b1,xinp_b2
      common /cir/ b1ci,b2ci,cnci,deltaci,xkappaci
      contains
!lk+
!----------------xpsatw
      real function xpsatw(t)
!c
      real t,arg
      arg = 54.842763-6763.22/t-4.210*log(t)+0.000367*t
     >    + tanh(0.0415*(t-218.8))
     >    * (53.878-1331.22/t-9.44523*log(t)+0.014025*t)
      xpsatw = 0.01*exp(arg)
!c
      return
      end
!----------------phi
      real function phi(s)
      real sm,sp,sstar,ds,smid,xnstar,dummy
      integer kinp
      real s,z
!c
!c** evaluates ice activation functions
!c
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,dummy
!c
      if (kinp.eq.1) then
!c
       if (s.lt.sstar) then
        phi = 0.
       else
        phi = 1.
       endif
!c
      elseif (kinp.eq.2) then
!c
       if (s.le.sm) then
        phi = 0.
       elseif (s.ge.sp) then
        phi = 1.
       else
        phi = (s-sm)/(sp-sm)
       endif
!c
      elseif (kinp.eq.12) then
!c
       z = (s-smid)/ds
       phi = 0.5*(tanh(z)+1.)
!c
      endif
!c
      return
      end
!----------------actg
      real function actg(s)
      real voli,cpa,bk,rgas,tol,zero,pi,rstar,ctot,rliq
      real temp,press,wup,xnci,rci,alphaci,sigma
      real cnci,w,sp,sm,sstar,ds,smid,wdown
      real ssiout,xnhet,frachom,xnhom,xnhom0
      real gamma,gammad,xlsubk
      real theta, a1,eicgs,xnis,dvap,vth,xnwsat,flux
      real xkappaci
      integer i,k
      integer kinp,iin,ipar,iout,jinp
      real omega,v,ssihat,xnstar,dummy,delta,xkappa,b1,b2,coolr,tair
      real fluxci,b1ci,b2ci,deltaci,tauc,taug,taugci,xkqppaci,ssihom
      real s,ag1,ag2,ag3,arg,delf,delp1,deltas,fzm,fzp,rho,xks2,zm
      real zp
!c
!c** evaluates activation-growth integral
!c
      common /par/ delta,xkappa
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,dummy
!c
      actg = 1.e-8
      if ((s.le.0.).or.(kinp.gt.10.)) return
!c
      if (kinp.eq.1) then
!c
       if (s.le.sstar) return
!c
       rho = (1.+delta)*sqrt(1.+xkappa*(s**2-sstar**2)) - 1.
       actg = rho**2 / (rho+1.)
!c
      elseif (kinp.eq.2) then
!c
       if (s.le.sm) return
!c
       zm = sm / s
       zp = min((sp/s),1.)
       delp1 = 1.+ delta
       deltas = sp - sm
       xks2 = xkappa * s**2
!c
       arg = zm * sqrt(xks2/(1.+ xks2))
       fzm = asin(arg) / (s*sqrt(xkappa))
       arg = zp * sqrt(xks2/(1.+ xks2))
       fzp = asin(arg) / (s*sqrt(xkappa))
       delf = fzp - fzm
!c
       ag1 = (zp * sqrt(1.+xks2*(1.-zp**2))
     >     -  zm * sqrt(1.+xks2*(1.-zm**2)))
     >     + (1.+ xks2) * delf
       ag1 = ag1  * delp1/2. * s/deltas
       ag2 = delf / delp1    * s/deltas
       ag3 = - 2.*(phi(s)-phi(sm))
!c
       actg = ag1 + ag2 + ag3
!c
      endif

      return
      end

!----------------dphids
      real function dphids(s)
      real sm,sp,sstar,ds,smid,xnstar,dummy
      integer kinp
!c
!c** evaluates derivates of ice activation functions for numerical integration
!c
      real dels,s
      parameter (dels=0.005)
!c
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,dummy
      real phi,z
!c
      if (kinp.eq.12) then
!c
       z = (s-smid)/ds
       phi = 0.5*(tanh(z)+1.)
       dphids = 2.*phi*(1.-phi)/ds
!c
      endif
!c
      return
      end
!----------------actgint
      real function actgint(sigma)
      real delta,xkappa,sm,sp,sstar,ds,smid,xnstar
      integer kinp
      real s,sigma,rho
!c
      common /par/ delta,xkappa
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,s
!c
      actgint = 0.
      if (s.le.sigma) return
!c
      rho = (1.+delta)*sqrt(1.+xkappa*(s**2-sigma**2)) - 1.
      actgint = dphids(sigma) * rho**2/(rho+1.)
!c
      return
      end
!----------------xwdown
      real function xwdown(s)
!lk+
      real sm,sp,sstar,ds,smid,xnstar,dummy
      real delta,xkappa
      real a1,b1,b2,xnwsat
      integer jinp,kinp
      real b1ci,b2ci,deltaci,cnci,xkappaci
      real s,ssi,fac,rhoci,xwdci,xnci,rci
      real fourpi,slow,act,actgci,xwdinp
      integer k
!lk-
      parameter (fourpi=12.56637061)
      data slow /0./
!c
!c!!** the lower integration limit 'slow' may not always work when set close to 0
!c
      integer inp_type(5)
      real xinp_tot(5),xinp_alpha(5)
      real xinp_b1(5),xinp_b2(5),xinp_delta(5),xinp_kappa(5)
!c
      common /inp/ jinp,inp_type,xinp_tot,xinp_alpha,
     >             xinp_delta,xinp_kappa,xinp_b1,xinp_b2
      common /run/ a1,b1,b2,xnwsat
      common /par/ delta,xkappa
      common /cir/ b1ci,b2ci,xnci,deltaci,xkappaci,rci
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,ssi
!c
      !external actgint
!c
      ssi = s                         !for use in common block /sol/
!c
      xwdown = 0.
      if (s.le.0.) return
!c
      fac = (fourpi*s/xnwsat)/(a1*(s+1.))
!c
!c** evaluate total íce supersaturation loss term in terms of a quenching velocity
!c
      rhoci  = (1.+deltaci)*sqrt(1.+xkappaci*s**2) - 1.
      actgci = rhoci**2 / (rhoci+1.)
      xwdci  = (b1ci/b2ci**2)*xnci*actgci
      if ((rci.eq.0).and.(xnci.eq.0.)) then
       xwdci  = 0.
      end if
!c
      xwdinp = 0.
      do k = 1, jinp
!c
       xnstar = xinp_tot(k)
       if (xnstar.gt.0.) then
!c
        delta  = xinp_delta(k)
        xkappa = xinp_kappa(k)
        kinp   = inp_type(k)
!c
        if (kinp.le.10) then
!c
!c** analytical solutions
!c
         act = actg(s)
!c
        elseif (kinp.gt.10) then
!c
!c** numerical solutions
!c
         act = 0.
         if (s.gt.slow) call QSIMP(actgint,slow,s,act)
!c
        endif
        xwdinp = xwdinp +
     >           (xinp_b1(k)/xinp_b2(k)**2)*xinp_tot(k)*act
       endif
!c
      enddo
!c
      xwdown = fac*(xwdci+xwdinp)
!c
      return
      end
!----------------xv
      real function xv(s)
      real s,wup   !lk
!c
      common /upd/ wup

!c
      xv = wup - xwdown(s)
!c
      return
      end
!----------------xlsub
      real function xlsub(t)
!lk+
      real t
      real arg
!lk-
!c
      arg = -(t/123.75)**2
      xlsub = 46782.5+35.8925*t-0.07414*t**2+541.5*exp(arg)
!c
      return
      end
!----------------xpsati
      real function xpsati(t)
!c
      real t,arg
      arg = 9.550426-(5723.265/t)+3.53068*log(t)-0.00728332*t
      xpsati = 0.01*exp(arg)
!c
      return
      end
!----------------xjhom
      real function xjhom(aw,t)
!lk+
      real aw,t
      real t1,dmu,awi,daw,dum2
!lk-
!c
      t1    = 1./t
      dmu   = 210368.+131.438*t-3.32373e6*t1-41729.1*log(t)
      awi   = exp(dmu*t1/8.31441)
      daw   = aw-awi
      dum2  = -906.688+8502.28*daw-26924.4*daw**2+29179.6*daw**3
      xjhom = 10.**max(-20.,min(20.,dum2))
!c
      return
      end
!----------------xfrz
      real function xfrz(sw)
!c
!lk+
      real pi,sw
      real coolr,temp,rad,sig
      real cren,tau,zsig,vol,rate
      !real xjhom
!lk-
      parameter (pi=3.141592654)
!c
      common /frz/ coolr,temp,rad,sig
!c
      !external xjhom
!c
      cren = 304.4-2.*temp+4.e-3*temp**2
      tau  = 1./(cren*abs(coolr))
      zsig = (log(sig))**2
      vol  = (4.*pi/3.)*rad**3*exp(4.5*zsig)
      rate = vol*xjhom(sw,temp)
!c
      xfrz  = (rate*tau)-1.
!c
      return
      end

!----------------ZBRENT
      REAL FUNCTION ZBRENT(FUNC,X1,X2,TOL)
!lk+
      REAL X1,X2,TOL
      REAL EPS
      REAL A,B,FA,FB,FC,C,D,E
      REAL TOL1,XM,S,P,Q,R
      INTEGER ITER,ITMAX
      REAL FUNC
!lk-
      PARAMETER (ITMAX=100,EPS=3.E-8)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
!      IF(FB*FA.GT.0.) PAUSE 'Root not bracketed in ZBRENT.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
!      PAUSE 'Too many iterations (100) in ZBRENT.'
      ZBRENT=B
      RETURN
      END
!----------------TRAPZD
      SUBROUTINE TRAPZD(FUNC,A,B,S,N)
!lk+
      REAL A,B,S
      REAL TNM,DEL,X,SUM
      INTEGER J,IT,N
      REAL FUNC
!lk-
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
!----------------QSIMP
      SUBROUTINE QSIMP(FUNC,A,B,S)
!lk+
      REAL S,A,B
      REAL EPS,OST,OS,ST
      INTEGER J,JMAX
      REAL FUNC
!lk-
      PARAMETER (EPS=5.E-3,JMAX=20)
      EXTERNAL FUNC
      OST=-1.E30
      OS= -1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        S=(4.*ST-OST)/3.
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
        OS=S
        OST=ST
11    CONTINUE
!CCC   PAUSE 'Too many steps (20) in QSIMP.'
      END


!----------------xshom
      real function xshom(rad,sig,temp,wup)
!c
      !lk+
      real tol,simin,simax,gammad
      real sig,rad,temp,wup
      real coolr,tair,rliq,sigma
      real convert,swmin,swmax,swhom
      !real xfrz
      !lk-
      parameter (tol=1.e-3,simin=1.4,simax=1.7)
      parameter (gammad=9.761e-5)
!c
      common /frz/ coolr,tair,rliq,sigma
!c
      !external xfrz
!c
      tair  = temp
      sigma = sig
      rliq  = rad
      coolr = gammad*wup
!c
      convert = xpsati(temp)/xpsatw(temp)
      swmin = simin*convert
      swmax = simax*convert
!c
      swhom = ZBRENT(xfrz,swmin,swmax,tol)
      xshom = swhom/convert
!c
      return
      end
!--------------

!lk-

!lk-
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   MAIN SUBROUTINE  "get_cirrus"  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
      subroutine get_cirrus (zin,zpar,zout,iin,ipar,iout,print)
!c
!c** physical constants / parameters
!c
!lk+
      real rgasa,grav,voli,cpa,bk,rgas,tol,zero,pi,rstar,ctot,rliq
      real temp,press,wup,xnci,rci,alphaci,sigma
      real cnci,w,sp,sm,sstar,ds,smid,wdown
      real ssiout,xnhet,frachom,xnhom,xnhom0
      real gamma,gammad,xlsubk
      real theta, a1,eicgs,xnis,dvap,vth,xnwsat,flux
      real xkappaci
      integer i,k
      integer kinp,iin,ipar,iout,jinp
      real omega,v,ssihat,xnstar,dummy,delta,xkappa,b1,b2
      real fluxci,b1ci,b2ci,deltaci,tauc,taug,taugci,ssihom
      !real xv
!lk-
      parameter (rgasa=2.8704e6, grav=980.665,  voli=3.25e-23,
     >           cpa=1.00467e7,  bk=1.3807e-16, rgas=8.31447,
     >           tol=1.e-4,      zero=0.,       pi=3.141592654)
!c
!c** dry INP core radius set constant across all INP types
!c
      parameter (rstar=0.25e-4)
!c
!c** parameters used to evaluate homogeneous freezing solution
!c
      parameter (ctot=500.,rliq=0.25e-4,sigma=1.)
!c
      real zin(iin),zpar(ipar),zout(iout)
!c
      integer inp_type(5)
      real xinp_act(5),xinp_tot(5),xinp_alpha(5)
      real xinp_b1(5),xinp_b2(5),xinp_delta(5),xinp_kappa(5)
!c
      logical print
      !external xv
!c
      common /upd/ w
      common /sol/ sm,sp,sstar,ds,smid,xnstar,kinp,dummy
      common /par/ delta,xkappa
      common /run/ a1,b1,b2,xnwsat
      common /inp/ jinp,inp_type,xinp_tot,xinp_alpha,
     >             xinp_delta,xinp_kappa,xinp_b1,xinp_b2
      common /cir/ b1ci,b2ci,cnci,deltaci,xkappaci
!c
!c** retrieve input from in-array & INP-related parameters
!c
      temp    = zin(1)
      press   = zin(2)
      wup     = zin(3)
      xnci    = zin(4)
      rci     = zin(5)
      alphaci = zin(6)
      jinp    = zin(7)
      do k = 1, jinp
       i = 8+3*(k-1)
       inp_type(k) = int(zin(i))
       xinp_tot(k) = zin(i+1)
       xinp_alpha(k) = zin(i+2)
      enddo
!c
      cnci = xnci                   !for use in common block /cir/
      w    = wup                    !for use in common block /upd/
!c
      sp    = zpar(1)
      sm    = zpar(2)
      sstar = zpar(3)
      ds    = zpar(4)
      smid  = zpar(5)
!c
      if (wup.lt.0.01) then
       wdown   = 0.
       ssiout  = 0.
       xnhet   = 0.
       frachom = 0.
       xnhom   = 0.
       xnhom0  = 0.
       if (print) then
        write(6,*) 'no ice formation'
        write(6,*)
       endif
       goto 898
      endif
!c
      if ((jinp.eq.0).and.(xnci.eq.0.)) then
       wdown  = 0.
       xnhet  = 0.
       ssiout = xshom(rliq,sigma,temp,wup)-1.
       call get_cirrus_hom (wup,temp,press,ctot,rliq,sigma,xnhom0)
       xnhom   = xnhom0
       frachom = 100.
       if (print) then
        write(6,*) 'pure homogeneous freezing'
        write(6,*) 'shom        = ',ssiout
        write(6,*) 'n_hom [#/L] = ',xnhom*1.e3
        write(6,*)
       endif
       goto 898
      endif
!c
!c** constants & parameters
!c
      gamma  = cpa/rgasa
      gammad = grav/cpa
!c
      xlsubk = xlsub(temp)/rgas
      theta  = xlsubk/temp
      a1     = gammad*(theta-gamma)/temp
      eicgs  = 1000.*xpsati(temp)
      xnis   = eicgs/(bk*temp)
      dvap   = 4.e-3*temp**1.94/press
      vth    = 3.4234e3*sqrt(temp)
      xnwsat = voli*xnis
!c
!c** ice crystal growth parameters
!c
      do k = 1, jinp
       flux = 0.25*xinp_alpha(k)*vth
       xinp_b1(k) = flux*xnwsat
       xinp_b2(k) = flux/dvap
       xinp_delta(k) = xinp_b2(k)*rstar
      enddo
!c
      fluxci  = 0.25*alphaci*vth
      b1ci    = fluxci*xnwsat
      b2ci    = fluxci/dvap
      deltaci = b2ci*rci
!c
!c** growth regime parameter
!c
      tauc = 1./(a1*wup)
      do k = 1, jinp
       taug = (1.+xinp_delta(k))/(xinp_b1(k)/rstar)
       xinp_kappa(k) = (tauc/taug)*xinp_delta(k)/(1.+xinp_delta(k))
       if (print) then
        write(6,*) 'INP type',inp_type(k)
        if (xinp_kappa(k).gt.1.) then
         write(6,*) 'fast growth regime: kappa =',xinp_kappa(k)
        else
         write(6,*) 'slow growth regime: kappa =',xinp_kappa(k)
        endif
        write(6,*)
       endif
      enddo
!c
      taugci = (1.+deltaci)/(b1ci/rci)
      xkappaci = (tauc/taugci)*deltaci/(1.+deltaci)
!c
!c** conditions at homogeneous freezing-relaxation
!c
      ssihom = xshom(rliq,sigma,temp,wup)-1.
      wdown = xwdown(ssihom)
      do k  = 1, jinp
       kinp = inp_type(k)
       xinp_act(k) = xinp_tot(k)*phi(ssihom)
      enddo
!c
!c** regime identification
!c
      omega = wdown/wup
      if (print) then
       if (omega.lt.1.) then
        write(6,*) 'freezing regime:    omega =',omega
       else
        write(6,*) 'quenching regime:   omega =',omega
       endif
      endif
!c
      if (wdown.lt.wup) then
!c
       if (print) then
        write(6,*) '(homogeneous freezing activated)'
!c
        write(6,*)
        write(6,*) 'w_down (cm/s)          =',wdown
        write(6,*)
        if (xnci.ne.0.) then
         write(6,*) 'pre-existing ice (#/L) =',(xnci*1.e3)
         write(6,*)
        endif
        do k = 1, jinp
         write(6,*) 'INP type',inp_type(k)
         write(6,*) 'n_ice (het) (#/L)      =',(xinp_act(k)*1.e3)
         write(6,*)
        enddo
       endif
       xnhet = 0.
       do k = 1, jinp
        xnhet = xnhet + xinp_act(k)
       enddo
       if (print) then
        if (jinp.gt.1) then
         write(6,*)
         write(6,*)'n_ice (het,tot) (#/L)   =',(xnhet*1.e3)
        endif
       endif
!c
       ssiout = ssihom
       v      = wup-wdown
!c
       call get_cirrus_hom (wup,temp,press,ctot,rliq,sigma,xnhom0)
       call get_cirrus_hom (  v,temp,press,ctot,rliq,sigma,xnhom)
       frachom = 100.*(xnhom0-xnhom)/xnhom0
       if (print) write(6,*)'effective updraft speed v (cm/s) =',v
!c
       if (print) then
        write(6,*)
        write(6,*)'n_ice (hom) (#/L) =',(xnhom*1.e3)
        write(6,*)
        write(6,*)'n_ice (hom) (#/L) unperturbed    =',(xnhom0*1.e3)
        write(6,*)'reduction relat. to pure hom (%) =',frachom
       endif
!c
      else
!c
       if (print) write(6,*) '(homogeneous freezing suppressed)'
!c
       ssihat = ZBRENT(xv,zero,ssihom,tol)
       ssiout = ssihat
       do k  = 1, jinp
        kinp = inp_type(k)
        xinp_act(k) = xinp_tot(k)*phi(ssihat)
       enddo
!c
       if (print) then
        if (xnci.ne.0.) then
         write(6,*)
         write(6,*)'pre-existing ice (#/L)  =',(xnci*1.e3)
        endif
        write(6,*)
        do k = 1, jinp
         write(6,*) 'INP type',inp_type(k)
         write(6,*)'n_ice (het) (#/L)       =',(xinp_act(k)*1.e3)
        enddo
       endif
       xnhet = 0.
       do k = 1, jinp
        xnhet = xnhet + xinp_act(k)
       enddo
       if (print) then
        if (jinp.gt.1) then
         write(6,*)
         write(6,*)'n_ice (het,tot) (#/L)   =',(xnhet*1.e3)
        endif
       endif
!c
       xnhom = 0.
       xnhom0 = 0.
       frachom = 100.
!c
      endif
      if (print) then
        write(6,*)'_________________________'
        write(6,*)
        write(6,*)'*** peak ice supersaturation (%)   =',100.*ssiout
        write(6,*)'*** total ice crystal number (#/L) =',
     >             (xnhom+xnhet+xnci)*1.e3
        write(6,*)
      endif
!c
898   continue
!c
!c** set out-array and exit
!c
      zout(1) = wdown
      zout(2) = ssiout
      zout(3) = frachom
      zout(4) = xnhom
      zout(5) = xnhom0
      zout(6) = xnhet
      do k = 1, jinp
       zout(6+k) = xinp_act(k)
      enddo
!c
      return
      end
!C
!---------------------------------------------------------------
!c
!c
      subroutine get_cirrus_hom (w,temp,press,ctot,rliq,sigma,xnhom)
!c
!c** evaluated homogeneously nucleated ice crystal number for monodisperse aerosols
!c** based on Kärcher, B., DeMott, P. J., Jensen, E. J., & Harrington, J. Y.
!c** JGR 127, https://doi.org/10.1029/2021JD035805 (2022).
!c
!lk+
      real rgasa,grav,voli,cpa,bk,rgas,tol,zero,pi,rstar,ctot,rliq
      real temp,press,wup,xnci,rci,alphaci,sigma
      real w
      real xnhom
      real gamma,gammad,xlsubk
      real theta, a1,a2,eicgs,dvap,vth,flux,cisat
      real delta,b1,b2,d1,frim,rei,rim,shom,syk
      real tau,tdot,yk,zeta,awi,daw,dmu,dum2,t1
      real arg,ren
      real alpha
!lk-
      parameter (pi=3.141592654, voli=3.23e-23, grav=980.665,
     >           cpa=1.00467e7,  rgas=8.31447,  rgasa=2.8704e6,
     >           bk=1.3807e-16,  alpha=0.7)
      real phi
!c
      gamma  = cpa/rgasa
      gammad = grav/cpa
      tdot   = gammad*w
      zeta   = temp*(0.004*temp-2.)+304.4
      tau    = 1./(zeta*tdot)
      shom   = xshom(rliq,sigma,temp,w)
      xlsubk = xlsub(temp)/rgas
      theta  = xlsubk/temp
!c
      eicgs  = 1.e3*xpsati(temp)
      cisat  = eicgs/(bk*temp)
      dvap   = 4.0122e-3*temp**1.94/press
      vth    = 3.4234e3*sqrt(temp)
      flux   = 0.25*alpha*vth
!c
      a1     = gammad*(theta-gamma)/temp
      a2     = 1./cisat
      b1     = (voli*cisat)*flux*(shom-1.)
      b2     = flux/dvap
!c
!c rim is the monodisperse freezing-growth integral
!c
      delta  = b2*rliq
      d1     = 1.+delta
      yk     = 2.*b1*b2*tau/d1**2
      syk    = sqrt(yk)
      ren    = 3.*syk/(2.+sqrt(1.+9.*yk/pi))
      frim   = 4.*pi*b1/(b2**2*voli)
      rim    = (frim/d1)*(delta**2-1.+(1.+0.5*yk*d1**2)*(ren/syk))
      phi    = w*a1*shom/a2       
!c
      xnhom = 1./(1./(phi/rim)+1./ctot)     
!c
      return
      end
!c
!c
!C
!C+++++++++++++ Numerical Recipes
!C+++++++++++++ -----------------
!C
!C
!C
!C
      end module K2022
