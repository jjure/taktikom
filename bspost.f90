subroutine bspost(rlat,rlon,az,t,lgust,sog, sogc, twa,twac,twd, twso,twsb,tws30m,wg,uc,vc,tcs,rlsm,isc, iscb)

! procedura glede na lokacijo in cas izracuna
! twa, twac
!
  use ezspline_obj
  use ezspline
  
  use yomdata, only  : nwgtimeh, twsb_splineo, twsb30m_splineo
  use yomconf, only  : ltwsboost, lcurrent
  use yomconst, only : pi
  
  real,intent(in) :: rlat, rlon, az, t ! lat, lon, azimuth, cas
  logical,intent(in) :: lgust

  real,intent(out) :: sog, sogc, twa, twac, twd, twso, twsb, tws30m, wg, uc, vc, tcs
  real,intent(out) :: rlsm
  integer,intent(out) :: isc

  logical          :: lout
  real  :: ub, vb


  call wi(rlat,rlon,t,u,v, lout) ! intepolate  the wind value at lat, lon, t
  if (lcurrent) then
     ! call cui(rlat,rlon,uc,vc)
  else
     uc=0.
     vc=0.
  endif
  
  tcs=sqrt(uc**2+vc**2)*3600./1852.
  if (ti<=(nwgtimeh).and.lgust) then
     call wgi(rlat,rlon,t,wg)
  else
     wg=0.
  endif

  call check_rp(rlat,rlon,rlsm)

!  if (lout) rlsm = 1.

  tws=sqrt((u+uc)**2+(v+vc)**2)*3600./1852.      !speed in knots

 
  if(ltwsboost) then
     call ezspline_interp(twsb_splineo,tws,twsb,ier)
     call ezspline_error(ier)
     if (ier .ne. 0 ) then 
        write(*,*)'bspost twsb intepolations'
     endif
  else
     twsb=1.0
  endif
  twso=tws
  twsb=twsb*tws
  tws=twsb
 
  call ezspline_interp(twsb30m_splineo,tws,tws30m,ier)
  call ezspline_error(ier)
  if (ier .ne. 0 ) then 
     write(*,*)'bspost twsb30m intepolations'
     call exit(1)
  endif
  tws30m=tws30m*tws
  
  ub=-sin(az * pi /180)               !u,v componenets of the
  vb=-cos(az * pi /180)
  twa=wba(u,v,ub,vb)
  twd=180/pi*atan2(-u,-v)
  
 ! if ( rlsm > 0) then 
 !    tws = 0.
 ! endif
  if (lcurrent) then
     tcd=180/pi*atan2(-uc,-vc)
     tca=wba(ub,vb,uc,vc)
     dcs=tcs*cos(twa*pi/180.)
  else
     tca=0.
     dcs=0.
  endif
  
  
  if (twd<0) twd=360+twd
  
  call boatspeed_new (twa,tws,sog, sogc, twac,.true.,.false.,isc)

  sog=sog+dcs
  sogc=sogc + dcs 
!!!!!!!!!!!
!!!!!!!!!!
  tcs=dcs

end subroutine bspost
