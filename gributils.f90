! gributils.f90
! containing 
!   subroutine gribSetup
!   subroutine readGrib
!
subroutine gribSetup

  use yomconf, only : cpref, model, grib_uv, member
  use yomdata, only : u, v, lsm
  use yomexclzone, only   : lexclzone, nexclzone, cexclzone,texclzone, xy, npolypoints
  use yomgrib
  use EZspline_obj
  use EZspline

  implicit none
  

  integer             :: param_u, param_v, level, i, j, inpoly, mm, nn
  real                :: lat, lon
!  real(r8)            :: lat,lon, p_lat, p_lon, uu, vv, rtmp
  logical             :: lana

  namelist /namgrib/  grib_lsm
  
  interface
     subroutine readGrib(gribname,member2decode,param2decode,level2decode,lprint,lana,gribObj2d,gribObj3d,tsMin,tsMax)

       use yomgrib
       use EZspline_obj
       use EZspline

       character(256),intent(in)     :: gribname 
       integer,optional, intent(in)  :: member2decode
       integer, intent(in)           :: param2decode
       integer, intent(in)           :: level2decode
       logical, intent(in)           :: lprint
       logical, intent(in)           :: lana
       integer, intent(in), optional :: tsMin, tsMax
       type(decGribObj2d), intent(out),optional          :: gribObj2d
       type(decGribObj3d), intent(out),optional          :: gribObj3d
     end subroutine readGrib
  end interface

  interface
     subroutine gribInfo(gribObjDesc)
       use yomgrib
       use EZspline_obj
       use EZspline
       type(decGribDesc), intent(in)         :: gribObjDesc     
     end subroutine gribInfo
  end interface
  
 
  rewind(4)
  read (4,namgrib)

  write(*,*)'gribSetup'
 
!  grib_lsm=trim(cpref)//trim(grib_lsm)
!  grib_uv=trim(cpref)//trim(grib_uv)
  print*,'  grib_lsm     : ',trim(grib_lsm)
  print*,'  grib_uv      : ',trim(grib_uv)
  if (model.eq.'ecfc') then
!     grib_uv = '/home/jure/projects/gribRead/test_gribs/azores_fc_2008081112_0.5x0.5.grib'
     param_u=165
     param_v=166
     level=0
     member=0
     lana=.FALSE.
  elseif(model.eq.'gfsfc') then
!     grib_uv = '/home/jure/projects/gribRead/test_gribs/gfs-2010020500.grb'
     param_u=33
     param_v=34
     level=10
     member=0
     lana=.FALSE.
  elseif(model.eq.'ecana') then
!     gribname = '/home/jure/projects/gribRead/test_gribs/v7o_l1_an_1997_10.11_1.0x1.0.grib'
     param_u=165
     param_v=166
     level=0
     member=-1
     lana=.TRUE.
  elseif(model.eq.'eceps') then
!     gribname='/home/jure/SP/azores/gribs/ecmwf/eps/azores_eps_2009070200_0.5x0.5.grib'
     param_u=165
     param_v=166
     level=0
!     member=1
     lana=.FALSE.
  elseif(model.eq.'ecfcms') then
!     gribname='/home/jure/projects/gribRead/test_gribs/azores_fc_2008081600_0.5x0.5.grb'
     param_u=33
     param_v=34
     level=0
     member=0
     lana=.FALSE.
  else
     write(*,*)'Unknow model/grib type, no preset for this configuration!'
     call exit(1)
  end if

! Lans sea mask
!  grib_lsm = '/home/jure/projects/gribRead/test_gribs/rdr_lsm_0.25x0.25.grib'
  call readGrib(grib_lsm,0,172,0,.FALSE.,lana=.FALSE.,gribObj2d=lsm)
  write(*,*)'Grib object parameters for LSM'
  call gribInfo(lsm%desc)

  if (lexclzone) then
     nn=1
     do i=1, lsm%desc%nlon
        do j=1, lsm%desc%nlat
           lat=lsm%desc%rlata(j)
           lon=lsm%desc%rlona(i)
           call locpt(lon,lat,xy(1,1,1:npolypoints(nn)),xy(1,2,1:npolypoints(nn)),npolypoints(nn),inpoly,mm)
           write(45,'(I5,2F8.2)')inpoly,lon,lat
        end do
     end do
  end if
  

! Uncoment this part if you want fort.11 to be written out
!!$  write(11,*)"i j lat lon lsm"
!!$  do i=1, lsm%desc%nlon
!!$     do j=1, lsm%desc%nlat
!!$        lat=lsm%desc%rlata(j)
!!$        lon=lsm%desc%rlona(i)
!!$        call EZspline_interp(lsm%splineObj, lon, lat, uu, ier)
!!$        call EZspline_error(ier)
!!$        write(11,'(2I5,4F8.2)')i,j,lat, lon, uu
!!$     enddo
!!$  enddo

  print*,'member :', member
  call readGrib(grib_uv,member,param_u,level,lprint=.FALSE.,lana=lana,gribObj3d=u)
  print*,u%desc%nlon
  
  write(*,*)'Grib object parameters for U'
  call gribInfo(u%desc)
 
  call readGrib(grib_uv,member,param_v,level,lprint=.FALSE.,lana=lana,gribObj3d=v)
  write(*,*)'Grib object parameters for V'
  call gribInfo(v%desc)

end subroutine gribSetup


subroutine gribInfo(fieldObj)
  use grib_api
  use EZspline_obj
  use EZspline
 
  use yomgrib

  implicit none
  
  integer, external :: validTD

  type(decGribDesc), intent(in)          :: fieldObj

  write(*,100)'nlon   : ',fieldObj%nlon
  write(*,100)'nlat   : ',fieldObj%nlat
  write(*,101)'lat0   : ',fieldObj%lat0
  write(*,101)'lon0   : ',fieldObj%lon0
  write(*,101)'lat1   : ',fieldObj%lat1
  write(*,101)'lon1   : ',fieldObj%lon1
  write(*,101)'dlat   : ',fieldObj%dlat
  write(*,101)'dlon   : ',fieldObj%dlon
  write(*,100)'nstep  : ',fieldObj%nstep
  write(*,102)'ts min : ',fieldObj%ts_min,validTD(fieldObj%ts_min)
  write(*,102)'ts max : ',fieldObj%ts_max,validTD(fieldObj%ts_max)
100 format(A10,I4)
101 format(A10,f6.2)
102 format(A10,2I12)
end subroutine gribInfo

!---------------------------------------------------------------

subroutine readGrib(gribname,member2decode,param2decode,level2decode,lprint,lana,gribObj2d,gribObj3d,tsMin,tsMax)

! Subroutine returns decoded grib object (either 2d or 3d)
! If number2decode=0, then the forecast or analysis is decoded
! If number2decode > 0, eps forecast is decoded

  use grib_api
  use EZspline_obj
  use EZspline
 
  use yomgrib
  use yomexclzone, only   : lexclzone, nexclzone, cexclzone,texclzone, xy, npolypoints

  implicit none
 
  character(256),intent(in)     :: gribname  
  integer,optional, intent(in)  :: member2decode
  integer, intent(in)           :: param2decode
  integer, intent(in)           :: level2decode
  logical, intent(in)           :: lprint
  logical, intent(in)           :: lana

  type(decGribObj2d), intent(out),optional          :: gribObj2d
  type(decGribObj3d), intent(out),optional          :: gribObj3d  
  type(decGribObj3d)                                :: gribObj
  integer,optional, intent(in) :: tsMax, tsMin

  integer       :: idx,ifile,n, iret,i, j,  numberOfValues,date,m ,iigrib, numberOfPoints
  integer       :: paramIdSize, levelSize, stepSize, memberSize, dataDateSize, datatimeSize, nsteps
  
  integer       :: bcs1(2),bcs2(2),bcs3(2)
  integer       :: inpoly, mm, nn
  integer,dimension(:),allocatable    ::  eps, indicatorOfParameter, periodOfTime
  integer,dimension(:),allocatable    ::  members,levels
  real,dimension(:),allocatable       ::  rsteps
  integer,dimension(:),allocatable    ::  steps, times, dates, times_i, dates_i
  integer,dimension(:),allocatable    ::  paramId,fieldsToDecode
  logical,dimension(:),allocatable    ::  mask
  logical       :: ltwindow

  real, dimension(:),allocatable      :: data1d,lats,lons
  real, dimension(:,:),allocatable    :: data2d
  real, dimension(:,:,:), allocatable :: data3d

  real                                :: min_val, max_val, lat, lon, uu
  integer                             :: olevel, ostep, omember, oparamId, igrib
  integer                             :: iScansNegatively, jScansPositively,dataRepresentationType
 
  logical                             :: l3d
  integer                             :: jPointsAreConsecutive
  integer                             :: count,countFields,ii, jj
  integer                             :: ier, idummy1, idummy2
  integer                             :: year, month, day, hour
  character(255)                      :: errString
  integer                             :: vt, vta
  integer,external                    :: validTS
  integer                             :: k
  integer, external                   :: validTD
  logical                             :: ldecode
  
  print*,"---------------------------------------------------------------------"
  ltwindow = .false.
  if(present(tsMin).and.present(tsMax)) then
     ltwindow=.true.
     print*,tsMin,tsMax
  endif
  

  print*,'open grib file: ',trim(gribname)

  ! Decide which case we are: 2d or 3d or exit if we can't decide
  if (present(gribObj2d)) then
     write(*,*)'2d case!'
     l3d=.false.
  elseif(present(gribObj3d)) then
     write(*,*)'3d case!'
     l3d=.true.
  else
     write(*,*)'No output array given!'
     call exit(1)
  endif

!!! Create the index first depending on the configuration/type of the grib file
!!! Index is based on paramId, level, step and in the case of eps, number, 
!!! 

  if (member2decode>0) then 
     call grib_index_create(idx,gribname,'indicatorOfParameter,number,level,step')
     print*,'1: Eps Case  !!!'
  elseif(member2decode.eq.0) then
     call grib_index_create(idx,gribname,'indicatorOfParameter,level,step')
  elseif(member2decode.eq.-1) then
     call grib_index_create(idx,gribname,'indicatorOfParameter,level,dataDate,dataTime')
  endif

  gribObj%desc%member=member2decode


!!! Get the sizes of key arrays
!!! 
 
  call grib_index_get_size(idx,'indicatorOfParameter',paramIdSize)
  call grib_index_get_size(idx,'level',levelSize)
print*,levelSize
print*,paramIdSize
  if (member2decode>0) then
     call grib_index_get_size(idx,'number',memberSize)
     call grib_index_get_size(idx,'step',stepSize)
  elseif (member2decode>=0) then
     print*,member2decode
     print*,'3: Eps Case  !!!'
     call grib_index_get_size(idx,'step',stepSize)
  elseif(member2decode<0) then
     call grib_index_get_size(idx,'dataDate',dataDateSize)
     call grib_index_get_size(idx,'dataTime',dataTimeSize)
  endif

  print*,'CP 1'
!!! Allocate key's arrays
  print*,paramIdSize
  allocate(paramId(paramIdSize))
  allocate(levels(levelSize))
  if (member2decode>0) then
     allocate(members(memberSize))
     allocate(steps(stepSize))
  elseif (member2decode>=0) then
     allocate(steps(stepSize))
  elseif(member2decode<0) then
     allocate(times(dataTimeSize))
     allocate(dates(dataDateSize))
     allocate(times_i(dataTimeSize * dataDateSize))
     allocate(dates_i(dataTimeSize * dataDateSize))
  endif

  print*,'CP 2'
 
!!! Fill the key's arrays
  call grib_index_get(idx,'indicatorOfParameter',paramId)
  print*,'parameters : ', paramId

  print*,'CP 3'
  call grib_index_get(idx,'level',levels)
  print*,'levels     : ', levels
  if (member2decode>0) then
     call grib_index_get(idx,'number',members)
     print*,'members : ', members
     call grib_index_get(idx,'step',steps)
     print*,'steps   : ', steps
  elseif (member2decode>=0) then
     call grib_index_get(idx,'step',steps)
     print*,'steps   : ', steps
  elseif(member2decode<0) then
     call grib_index_get(idx,'dataDate',dates)
     print*,'dates   : ', dates
     call grib_index_get(idx,'dataTime',times)
     print*,'times   : ',times
  endif

  print*,'CP 4'

  if (lprint) then
     write(*,111)'Number of disnct paramters in grib file:',paramIdSize
     write(*,111)'levelSize:',levelSize
     if (member2decode>0) then
        write(*,111)'Number of distinct eps members :',memberSize 
     elseif (member2decode>=0) then
        write(*,111) 'stepSize=',stepSize
     elseif (member2decode<0) then
        write(*,111) 'dataDate=',dataDateSize
        write(*,111) 'dataTime=',dataTimeSize        
     endif
  endif

!!! Compute the number of time steps in the grib file for giving
!!! parameter

  if (member2decode>=0) then
     nsteps=stepSize 
  elseif (member2decode<0) then
     nsteps = dataDateSize * dataTimeSize
  endif

  
! We pass through member2decode also the information about grib type
! analysis(member2decode <0), forecast (member2decode=0) or ensamble (member2decode>0)  
  
  gribObj%desc%nstep=nsteps

! allocate the array to contain the list of distinct step times /analysis times
  if (.not.allocated(steps)) allocate(steps(nsteps))
  allocate(rsteps(nsteps))
  allocate(gribObj%desc%rstepa(nsteps))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LOOOOOOOOOOOOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  countFields=0
 
!!! Preselect indexes
  call grib_index_select(idx,'indicatorOfParameter',param2decode)
  call grib_index_select(idx,'level',level2decode)
 
  if (member2decode > 0) then
     call grib_index_select(idx,'number',member2decode)
  endif
  
  print*,'nsteps :',nsteps


!!! Prepare steps for analysis case
  if (member2decode < 0 ) then 
     k = 0
     do i = 1, dataDateSize
        do j = 1, dataTimeSize
           k = k + 1
           times_i(k)=j
           dates_i(k)=i
        enddo
     enddo
  endif

  do i=1,nsteps
     print*,'step : ',i, steps(i)

!!! Select the time index 
   
     if  (member2decode>=0) then
        call grib_index_select(idx,'step',steps(i))
     elseif (member2decode<0) then
        call grib_index_select(idx,'dataTime',times(times_i(i)),iret)
        call grib_index_select(idx,'dataDate',dates(dates_i(i)),iret)       
     endif
     
     call grib_new_from_index(idx,igrib,iret)   

!!! Decode date && time parameters
     if (member2decode >=0) then 
        call grib_get(igrib,'dataDate',idummy1,iret) ! analysis date
        call grib_get(igrib,'endStep',idummy2,iret)  ! forecast range
        call grib_get(igrib,'hour',hour)
        vta = validTS(idummy1,hour) 
        vt = vta + idummy2 * 60
     else
        call grib_get(igrib,'dataDate',idummy1,iret) ! analysis date
        call  grib_get(igrib,'dataTime',idummy2,iret) ! analysis time
        idummy2=idummy2/100
        vt=validTS(idummy1,idummy2)
     endif



     call grib_get(igrib,'typeOfGrid',gribObj%desc%typeOfGrid)
     if (trim(gribObj%desc%typeOfGrid).ne.'regular_ll') then 
        write(*,*) 'TypeOfGrid :', gribObj%desc%typeOfGrid
        write(*,*) 'Just regular lat lon grib implemented, exiting!'
        call exit(1)
     endif

     call grib_get(igrib,'numberOfPointsAlongAParallel',gribObj%desc%nlon) 
     call grib_get(igrib,'numberOfPointsAlongAMeridian',gribObj%desc%nlat) 
     call grib_get(igrib,'yFirst',gribObj%desc%lat0)
     call grib_get(igrib,'xFirst',gribObj%desc%lon0)
     call grib_get(igrib,'yLast',gribObj%desc%lat1)
     call grib_get(igrib,'xLast',gribObj%desc%lon1)     
     call grib_get(igrib,'jDirectionIncrementInDegrees',gribObj%desc%dlat)
     call grib_get(igrib,'iDirectionIncrementInDegrees',gribObj%desc%dlon)
     call grib_get(igrib,'iScansNegatively',iScansNegatively)
     call grib_get(igrib,'jScansPositively',jScansPositively)
     call grib_get(igrib,'jPointsAreConsecutive',jPointsAreConsecutive)
     call grib_get(igrib,'dataRepresentationType',dataRepresentationType)
     call grib_get(igrib,'month',month)
     call grib_get(igrib,'year',year)
     call grib_get(igrib,'day',day)
     call grib_get(igrib,'hour',hour)
    
    ! temporary fix: to work on west longitude 
!    print*,gribObj%desc%lon0,gribObj%desc%lon1
    if (gribObj%desc%lon0 > 0 .and. gribObj%desc%lon1 > 180) then 
            gribObj%desc%lon0=gribObj%desc%lon0-360.
            gribObj%desc%lon1=gribObj%desc%lon1-360.
    endif
!    print*,gribObj%desc%lon0,gribObj%desc%lon1
    !    call exit(0)
     if (lprint) then 
        write(*,111)'dataRepresentationType  :',dataRepresentationType
        write(*,111)'iScansNegatively        :',iScansNegatively
        write(*,111)'jScansPositively        :',jScansPositively
        write(*,111)'jPointsAreConsecutive   :',jPointsAreConsecutive
        write(*,111)'day                     :',day
        write(*,111)'month                   :',month
        write(*,111)'year                    :',year
        write(*,111)'hour                    :',hour
        write(*,111)'timestamp               :',vt
        write(*,111)'decoded TS              :',validTD(vt)
     endif
     ! Note: i direction is defined as west to east along a parallel of latitude, or left to right along an x axis.
     !j direction is defined as south to north along a meridian of longitude, or bottom to top along a y axis.
     !  	0       Adjacent points in i direction are consecutive
     !		(FORTRAN: (I,J))
     !	1 	Adjacent points in j direction are consecutive
     !		(FORTRAN: (J,I)) 

     if (iScansNegatively.ne.0.and.jScansPositively.ne.0) then
        write(*,*)'Unimpleneted way of Scanning either in i either in j direction!'
        write(*,*)'Should be relatively easy to implement!!'
        call exit(1)
     endif

     if (dataRepresentationType.ne.0) then
        write(*,*)'Unimpleneted dataRepresentationType (not a lat/lon grib)!'
        call exit(1)
     endif

     !     get the size of the values array
        call grib_get_size(igrib,'values',numberOfPoints)
        
        if(lprint) write(*,111) 'numberOfPoints          :',numberOfPoints    
        
        if (.not.allocated(data1d)) then
           allocate(data1d(numberOfPoints))
           allocate(lats(numberOfPoints))
           allocate(lons(numberOfPoints)) 
           allocate(data2d(gribObj%desc%nlon,gribObj%desc%nlat))
           allocate(data3d(gribObj%desc%nlon,gribObj%desc%nlat,nsteps))
        endif

       
! Definicijsko obmocje     
        if (.not.allocated(gribObj%desc%rlata)) then
           allocate(gribObj%desc%rlata(gribObj%desc%nlat))
           allocate(gribObj%desc%rlona(gribObj%desc%nlon)) 
           do ii=1,gribObj%desc%nlon
              gribObj%desc%rlona(ii)=gribObj%desc%lon0+((ii-1)*gribObj%desc%dlon)
           enddo
! Pozor: normalna smer v medridionalni smeri je N-S, zato zacnemo z lat1 in ne lat0
           do ii=1,gribObj%desc%nlat
              gribObj%desc%rlata(ii)=gribObj%desc%lat1+((ii-1)*gribObj%desc%dlat)
           enddo
        end if
     
       
        call grib_get(igrib,'paramId',oparamId)
        if(member2decode>0) then
           call grib_get(igrib,'number',omember)
        elseif(member2decode<0) then
           call grib_get(igrib,'step',ostep)
        endif
        call grib_get(igrib,'level',olevel)
        call grib_get_data(igrib,lats,lons,data1d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Rotate data in a way, that data(1,1) coresponds to SW corner
     ! data(i,j)  : i goes along parallel
     !            : j goes along meridian =>  data(l
     
        count=0
        do ii=gribObj%desc%nlat,1,-1
           do jj=1,gribObj%desc%nlon
              count=count+1
              data2d(jj,ii)=data1d(count)
           enddo
        enddo
       
        ldecode=.false.

        if (ltwindow) then
           if (vt.ge.tsMin.and.vt.le.tsMax) then
              ldecode=.true.
              print*,'=============>'
           endif
        else
           ldecode=.true.
        endif
        
        if(ldecode) then
           countFields=countFields+1
           data3d(:,:,countFields)=data2d
           
           call grib_get(igrib,'min',min_val)
           call grib_get(igrib,'max',max_val)
           write(*,'(A,i3,A,i2,A,i4,A,i3,A,f10.2,A,f10.2)') ' paramId=',oparamId,&
                '   number=',gribObj%desc%member,&
                '   level=' ,olevel, &
                '   step='  ,ostep, &
                '   min_val=',min_val, &
                '   max_val= ',max_val
 
           print*,i,countFields
           gribObj%desc%rstepa(countFields)=vt
        endif
     enddo

  print*,'Number of fields decoded :', countFields
  call grib_index_release(idx)

!!! Get the time interval
  gribObj%desc%ts_min=minval(gribObj%desc%rstepa(1:countFields))
  gribObj%desc%ts_max=maxval(gribObj%desc%rstepa(1:countFields))
  
! Spline object setup -  Ezspline setup 
  bcs1=(/ 0, 0/)
  bcs2=(/ 0, 0/)
  bcs3=(/ 0, 0/)
  
  if (countFields==1) then
     gribObj2d%desc=gribObj%desc     
     write(*,*)'EZspline init gribObj',paramId
     call EZspline_init(gribObj2d%splineObj, gribObj2d%desc%nlon, gribObj2d%desc%nlat, &
          & bcs1, bcs2, ier)
     call EZspline_error(ier)
  
     gribObj2d%splineObj%x1=gribObj%desc%rlona
     gribObj2d%splineObj%x2=gribObj%desc%rlata
     nn=1
     if (lexclzone) then 
        do i=1,gribObj2d%desc%nlon
           do j=1,gribObj2d%desc%nlat
              lat=gribObj%desc%rlata(j)
              lon=gribObj%desc%rlona(i)
              call locpt(lon,lat,xy(1,1,1:npolypoints(nn)),xy(1,2,1:npolypoints(nn)),npolypoints(nn),inpoly,mm)
              if (inpoly==-1) then
                 data2d(i,j)=1
              endif
!              write(46,'(I4,","I4,",",F8.3,",",F8.3,",",F8.3)'),i,j,lon,lat,data2d(i,j)
           enddo
        enddo
     endif
     call EZspline_setup(gribObj2d%splineObj, data2d, ier)
     call EZspline_error(ier)
     
     call EZspline_isGridRegular(gribObj2d%splineObj,ier)
     call EZspline_error(ier)
  elseif(nsteps.gt.0) then
     gribObj3d=gribObj
     gribObj3d%desc%nstep=countFields

!     call EZspline_init(gribObj3d%splineObj, gribObj3d%desc%nlon, gribObj3d%desc%nlat,&
!          & gribObj3d%desc%nstep,bcs1, bcs2, bcs3,ier)
     call EZlinear_init(gribObj3d%splineObj, gribObj3d%desc%nlon, gribObj3d%desc%nlat,&
          & gribObj3d%desc%nstep,ier)
     call EZspline_error(ier) 

     
     
     gribObj3d%splineObj%x1=gribObj%desc%rlona
     gribObj3d%splineObj%x2=gribObj%desc%rlata
     gribObj3d%splineObj%x3=gribObj%desc%rstepa(1:countFields)   
  
     call EZspline_setup(gribObj3d%splineObj, data3d(:,:,1:countFields), ier)
     call EZspline_error(ier)
 
     call EZspline_isGridRegular(gribObj3d%splineObj,ier)
     call EZspline_error(ier) 
  else
     write(*,*)'Hmmm, internal error(1)'
     call exit(1)
  endif
  
  deallocate(paramId)
  deallocate(levels)
  deallocate(data1d)
  deallocate(data2d)
  deallocate(data3d)
  deallocate(lats)
  deallocate(lons)
  if (allocated(times_i)) deallocate(times_i)
  if (allocated(dates_i)) deallocate(dates_i)
  if (allocated(times)) deallocate(times)
  if (allocated(members)) deallocate(members)
  
111 format(A50,I10)
end subroutine readGrib




