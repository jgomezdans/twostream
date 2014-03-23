

SUBROUTINE TBarreUncoll_exact(tau,res)
   REAL(8), INTENT(IN) :: tau
   REAL(8), INTENT(OUT) :: res

   INTEGER :: j,ind
   REAL(8) :: iGammaloc
   INTEGER :: order

   iGammaloc = 0.0d0
   order=20

   DO j=0,order-1
      ind = order - j
      iGammaloc = ind / (1.0d0 + ind/(tau+iGammaloc))
   END DO
   iGammaloc=1.0d0/(tau + iGammaloc)

   res = exp(-tau)*(1.0d0 - tau + tau*tau*iGammaloc)
   
END SUBROUTINE TBarreUncoll_exact

SUBROUTINE gammas(rl,tl,mu0,gammaCoeff)
   REAL(8), INTENT(IN) :: rl,tl,mu0
   REAL(8), DIMENSION(4), INTENT(OUT) :: gammaCoeff

   REAL(8) :: wd,w0,w0half,wdsixth

   w0      = rl+tl
   wd      = rl-tl
   w0half  = w0*0.5_8
   wdsixth = wd/6.0_8

   gammaCoeff(1)=2.*(1. - w0half + wdsixth)
   gammaCoeff(2)=2.*(w0half + wdsixth)
   gammaCoeff(3)=2.*(w0half*0.5 + mu0*wdsixth)/w0
   gammaCoeff(4)=1.-gammaCoeff(3)

END SUBROUTINE gammas

SUBROUTINE dhrT1(rl,tl,gamma1,gamma2,gamma3,gamma4,tta,tau,AlbBS,Tdif,AbsVgt)
   REAL(8), INTENT(IN) :: rl,tl,gamma1,gamma2,gamma3,gamma4,tta,tau
   REAL(8), INTENT(OUT) :: AlbBS,Tdif,AbsVgt
!   LOGICAL :: dhrT1

   REAL(8) :: alpha1,alpha2,ksquare,k
   REAL(8) :: first_term,secnd_term1,secnd_term2,secnd_term3
   REAL(8) :: expktau,Tdir
   REAL(8) :: mu0,w0

   mu0=cos(tta)
   w0=rl+tl


   Tdir = exp(-tau/mu0)

   ! There is a difference between conservative and non-conservative scattering conditions */
   IF (w0 .ne. 1.0 .AND. w0 .ne. 0.0) THEN
      !NON_CONSERVATIVE SCATTERING

      ! some additional parameters 
      alpha1  = gamma1*gamma4 + gamma2*gamma3
      alpha2  = gamma1*gamma3 + gamma2*gamma4
      ksquare = gamma1*gamma1 - gamma2*gamma2
      k       = sqrt(ksquare)

      expktau = exp(k*tau)

      !Black Soil Albedo
      first_term  = ((1.0d0-ksquare*mu0*mu0)*((k+gamma1)*expktau + (k-gamma1)/expktau))
      IF (first_term .eq. 0.0) THEN
         !we will be dividing by zero : cannot continue.
!         dhrT1 = .false.
          RETURN
      ELSE
         first_term = 1.0d0/first_term
         secnd_term1 = (1.0d0 - k*mu0)*(alpha2 + k*gamma3)*expktau
         secnd_term2 = (1.0d0 + k*mu0)*(alpha2 - k*gamma3)/expktau
         secnd_term3 = 2.0d0 * k * (gamma3 - alpha2*mu0)*Tdir
         AlbBS = (w0 * first_term * (secnd_term1 - secnd_term2 - secnd_term3))

         !Transmission
         IF (ksquare .eq. 0.0) THEN
            first_term = 1.0d0
         ENDIF
         secnd_term1 = (1.0d0+k*mu0)*(alpha1+k*gamma4)*expktau
         secnd_term2 = (1.0d0-k*mu0)*(alpha1-k*gamma4)/expktau
         secnd_term3 = 2.0d0 * k * (gamma4 + alpha1*mu0)
         Tdif = - w0*first_term*(Tdir*(secnd_term1 - secnd_term2) - secnd_term3)

         !Absorption by vegetation
         AbsVgt = (1.0d0- (Tdif+Tdir) - AlbBS)
      ENDIF
   ELSE IF (w0 .eq. 0.) THEN
      !BLACK CANOPY
      AlbBS = 0.0d0
      Tdif  = 0.0d0
      AbsVgt = 1.0d0 - Tdir
   ELSE
      !CONSERATIVE SCATTERING
      AlbBS =  (1.0d0/(1.0d0 + gamma1*tau))*(gamma1*tau + (gamma3-gamma1*mu0)*(1.0d0-exp(-tau/mu0)));
      Tdif   = 1.0d0 - AlbBS - Tdir;
      AbsVgt = 0.0d0;
   ENDIF
END SUBROUTINE dhrT1


SUBROUTINE bhrT1(rl,tl,gamma1,gamma2,gamma3,gamma4,tau,AlbBS,Tdif,AbsVgt)
   REAL(8), INTENT(IN) :: rl,tl,gamma1,gamma2,gamma3,gamma4,tau
   REAL(8), INTENT(OUT) :: AlbBS,Tdif,AbsVgt
!   LOGICAL :: bhrT1

   REAL(8) :: a=0.705d0
   CALL dhrT1(rl,tl,gamma1,gamma2,gamma3,gamma4,acos(0.5d0/a),tau,AlbBS,Tdif,AbsVgt)
END SUBROUTINE bhrT1




SUBROUTINE twostream_solver(leaf_reflectance,leaf_transmittance,background_reflectance,&
   true_lai,structure_factor_zeta,structure_factor_zetaStar,sun_zenith_angle_degrees,&
   Collim_Alb_Tot,Collim_Tran_Tot,Collim_Abs_Tot,Isotrop_Alb_Tot,Isotrop_Tran_Tot,Isotrop_Abs_Tot)


   REAL(8), INTENT(IN) :: leaf_reflectance,leaf_transmittance,background_reflectance
   REAL(8), INTENT(IN) :: true_lai,structure_factor_zeta,structure_factor_zetaStar,sun_zenith_angle_degrees
   REAL(8), INTENT(OUT) :: Collim_Alb_Tot,Collim_Tran_Tot,Collim_Abs_Tot
   REAL(8), INTENT(OUT) :: Isotrop_Alb_Tot,Isotrop_Tran_Tot,Isotrop_Abs_Tot

   !internal and intermediate variables
!   LOGICAL :: ok
   REAL(8), DIMENSION(4) :: gammaCoeffs,gammaCoeffs_star
   REAL(8) :: tauprimetilde,tauprimestar,sun_zenith_angle_radians
   REAL(8) :: cosine_sun_angle

   !calculated fluxes
   REAL(8) :: Collim_Alb_BB, Collim_Tran_BB, Collim_Abs_BB
   REAL(8) :: Isotrop_Alb_BB,Isotrop_Tran_BB,Isotrop_Abs_BB
   REAL(8) :: Collim_Tran_BC, Isotrop_Tran_BC
   REAL(8) :: Collim_Tran_TotalOneWay,Isotrop_Tran_TotalOneWay
   REAL(8) :: Bellow_reinject_rad
   REAL(8), PARAMETER :: isotropic_cosine_constant=0.5/0.705

!   print*,leaf_reflectance,leaf_transmittance,background_reflectance
!   print*,true_lai,structure_factor_zeta,structure_factor_zetaStar,sun_zenith_angle_degrees
!   print*,"================="
   ! convert angular values
   sun_zenith_angle_radians = 3.14159265358979323846d0 * sun_zenith_angle_degrees / 180.00;
   cosine_sun_angle = cos(sun_zenith_angle_radians);

   ! calculate the 4 gamma coefficients both in isotropic and collimated illumination 
   call gammas(leaf_reflectance,leaf_transmittance,cosine_sun_angle,gammaCoeffs)

   call gammas(leaf_reflectance,leaf_transmittance,isotropic_cosine_constant,gammaCoeffs_star)
  
 
   ! estimate the effective value of the optical thickness 
   tauprimetilde = 0.5d0 * true_lai * structure_factor_zeta
   tauprimestar  = 0.5d0 * true_lai * structure_factor_zetaStar

   ! +++++++++++++ BLACK BACKGROUND ++++++++++++++++++++++ 
   ! * Apply the black-background 2stream solution 
   ! * These equations are written for the part of the incoming radiation 
   ! * that never hits the background but does interact with the vegetation 
   ! * canopy.
   ! * 
   ! * Note : the same routine dhrT1() is used for both the isotropic and
   ! * collimated illumination conditions but the calling arguments differ.
   ! * (especially the solar angle used).
   ! */

   !	/* 1) collimated source */
   CALL dhrT1(leaf_reflectance,leaf_transmittance,&
         gammaCoeffs(1),gammaCoeffs(2),gammaCoeffs(3),gammaCoeffs(4),&
         sun_zenith_angle_radians,tauprimetilde,&
         Collim_Alb_BB,Collim_Tran_BB,Collim_Abs_BB)
   ! 	/* 2) isotropic source */
!   print*,Collim_Alb_BB,Collim_Tran_BB,Collim_Abs_BB
   CALL dhrT1(leaf_reflectance,leaf_transmittance,&
         gammaCoeffs_star(1),gammaCoeffs_star(2),gammaCoeffs_star(3),gammaCoeffs_star(4),&
         acos(isotropic_cosine_constant),tauprimestar,&
          Isotrop_Alb_BB,Isotrop_Tran_BB,Isotrop_Abs_BB)
!   print*,Isotrop_Alb_BB,Isotrop_Tran_BB,Isotrop_Abs_BB
   ! +++++++++++++ BLACK CANOPY ++++++++++++++++++++++ 
   ! * Apply the black canopy solution.
   ! * These equations hold for the part of the incoming radiation 
   ! * that do not interact with the vegetation, travelling through 
   ! * its gaps.
   ! */

  ! 	/* 1) collimated source */
   Collim_Tran_BC  = exp( - tauprimetilde/cosine_sun_angle)
  ! 	/* 2) isotropic source */
   CALL TBarreUncoll_exact(tauprimestar,Isotrop_Tran_BC)

   ! /* Total one-way transmissions:
   ! * The vegetation canopy is crossed (one way) by the uncollided radiation
   ! * (black canopy) and the collided one (black background). */

   !	/* 1) collimated source */
   Collim_Tran_TotalOneWay  = Collim_Tran_BC  + Collim_Tran_BB 
   !	/* 2) isotropic source */
   Isotrop_Tran_TotalOneWay = Isotrop_Tran_BC + Isotrop_Tran_BB 
   
   ! * The Bellow_reinject_rad describes the process of reflecting toward the background
   ! * the upward travelling radiation (re-emitted from bellow the canopy). It appears in 
   ! * the coupling equations as the limit of the series: 
   ! *    1 + rg*rbv + (rg*rbv)^2 + (rg*rbv)^3 + ...
   ! *      where rg is the background_reflectance and rbv is Isotrop_Alb_BB
   ! *      (with Isotrop describing the Lambertian reflectance of the background).
   ! */
   Bellow_reinject_rad = 1.0d0 / (1.0d0 - background_reflectance*Isotrop_Alb_BB)

   !/* TOTAL ALBEDO */
   ! 	/* 1) collimated source */
   Collim_Alb_Tot = Collim_Alb_BB + &
      background_reflectance * Collim_Tran_TotalOneWay  * Isotrop_Tran_TotalOneWay * Bellow_reinject_rad 
   !	/* 2) isotropic source */
   Isotrop_Alb_Tot = Isotrop_Alb_BB + & 
      background_reflectance * Isotrop_Tran_TotalOneWay * Isotrop_Tran_TotalOneWay * Bellow_reinject_rad 
   
   !/* TOTAL TRANSMITION TO THE BACKGROUND LEVEL */
   !	/* 1) collimated source */
   Collim_Tran_Tot = Collim_Tran_TotalOneWay * Bellow_reinject_rad ;
   !	/* 2) isotropic source */
   Isotrop_Tran_Tot = Isotrop_Tran_TotalOneWay * Bellow_reinject_rad ;

   !/* TOTAL ABSORPTION BY THE VEGETATION LAYER */
   !	/* 1) collimated source */
   Collim_Abs_Tot = 1.0d0 - (Collim_Tran_Tot + Collim_Alb_Tot) + background_reflectance * Collim_Tran_Tot;
   !	/* 2) isotropic source */
   Isotrop_Abs_Tot = 1.0d0 - (Isotrop_Tran_Tot + Isotrop_Alb_Tot) + background_reflectance * Isotrop_Tran_Tot;

END SUBROUTINE twostream_solver

