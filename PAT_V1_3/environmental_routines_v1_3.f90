MODULE environmental_routines_v1_3
 USE general_routines
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,ALLOCATABLE::table5(:,:,:)
 REAL(iwp),ALLOCATABLE::table1(:,:),table2(:,:),table3(:,:),table4(:,:)
 REAL(iwp),PARAMETER::pi=3.141592653589793238_iwp,zero=0.0_iwp
 REAL(iwp)::t_m1,t_p1,delay_1,t_m2,t_p2,delay_2,delay_3
 REAL(iwp)::wl,ts,ao,tauo,ce,e1,e2,e3,e4,e5
 REAL(iwp)::aa,bb,Io
 PRIVATE::iwp
 PRIVATE::t_m1,t_p1,delay_1,t_m2,t_p2,delay_2
 PRIVATE::wl,ts,ao,tauo,ce,e1,e2,e3,e4,e5
 PRIVATE::aa,bb,Io
 PUBLIC::pi,zero
 PUBLIC::table1,table2,table3,table4,table5
CONTAINS
!
!-------------------------------------------------------------------------
! Counting of time
!-------------------------------------------------------------------------
SUBROUTINE day_of_year(jd,d)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::jd
 REAL(iwp),INTENT(out)::d
 INTEGER::year,month,jdate,jhour
 REAL(iwp)::jd_i,day
 CALL civil_date(jd,year,month,day) 
 jdate=(year-1)*10000+1200+31
 jhour=2400
 CALL julian_day(jdate,jhour,jd_i) 
 d=jd-jd_i
END SUBROUTINE day_of_year
!
SUBROUTINE civil_date(jd,year,month,day)
 IMPLICIT NONE
 REAL(iwp),INTENT(IN)::jd
 INTEGER,INTENT(OUT)::year,month
 REAL(iwp),INTENT(OUT)::day
 INTEGER::z,w,a,b,c,d,e,f
 REAL(iwp)::q
 q=jd+0.5
 z=INT(q)
 w=INT((jd-1867216.25)/36524.25)
 a=z+1+w-INT(w/4)
 b=a+1524
 c=INT((b-122.1)/365.25)
 d=INT(c*365.25)
 e=INT((b-d)/30.6001)
 f=INT(e*30.6001)
 day=b-d-f+(q-z)
 IF(e<14) THEN
   month=e-1
 ELSE
   month=e-13
 END IF
 IF(month>2) THEN
   year=c-4716
 ELSE
   year=c-4715
 END IF
!Corretion for leap years
 IF(month==2.and.day>=31)THEN
   day=day-2
   year=year-1
 END IF

END SUBROUTINE civil_date
!
SUBROUTINE julian_day(jdate,jhour,jd)
 IMPLICIT NONE
 INTEGER::jdate,jhour,y,m,d,a,b,c,e,f
 REAL(iwp)::jd,hour,minute
 y=INT(jdate/10000)
 m=INT((jdate-y*10000)/100)
 d=INT(jdate-y*10000-m*100)
 hour=FLOAT(jhour/100)
 minute=jhour-hour*100
 IF(m<3) THEN
   y=y-1
   m=m+12
 END IF
 a=INT(y/100)
 b=INT(a/4)
 c=2-a+b
 e=INT(365.25*(y+4716))
 f=INT(30.6001*(m+1))
 jd=FLOAT(c+d+e+f)-1524.5+hour/24.+minute/1440.
RETURN
END SUBROUTINE julian_day
!
FUNCTION date_hour(jd) RESULT(char_num)
 IMPLICIT NONE
 INTEGER::jdate,jhour,yyyy,mm
 REAL(iwp),INTENT(in)::jd
 REAL(iwp)::ut,minutes,seconds,dd
 CHARACTER(len=14)::char_num
 CHARACTER(len=8)::char_idate
 CHARACTER(len=6)::char_ihour
 ut=jd+0.5_iwp-AINT(jd+0.5_iwp)+1e-6
 CALL civil_date(jd,yyyy,mm,dd)
 jdate=yyyy*10000+mm*100+INT(dd)
 minutes=(ut*24_iwp-INT(ut*24_iwp))*60_iwp
 seconds=(minutes-INT(minutes))*60_iwp
 jhour=INT(ut*24_iwp)*10000+INT(minutes)*100+INT(seconds)
 WRITE(char_idate,'(i8)')jdate
 WRITE(char_ihour,'(i6)')jhour
 char_ihour=ADJUSTL(char_ihour)
 IF(jhour>=10000.and.jhour<100000)char_ihour='0'//TRIM(char_ihour)
 IF(jhour>=1000.and.jhour<10000)char_ihour='00'//TRIM(char_ihour)
 IF(jhour>=100.and.jhour<1000)char_ihour='000'//TRIM(char_ihour)
 IF(jhour>=10.and.jhour<100.)char_ihour='0000'//TRIM(char_ihour)
 IF(jhour<10)char_ihour='00000'//TRIM(char_ihour)
 char_num=char_idate//char_ihour
END FUNCTION date_hour
!=========================================================================
! Solar time
!-------------------------------------------------------------------------
SUBROUTINE solar_time(d,jd,long,timezone,saving,d2)
! This routine performs the conversion between clock time and solar time
! It requires the knowledge of the location given by the longitude (long)
! in degrees (positive to the east of the Prime Meridiam) and 
! the local standards to which local clocks are set given by the time zone
! (timezone) in hours from UTC and, when applicable, the daylight saving 
! time (saving)
! 
 IMPLICIT NONE
 INTEGER,INTENT(in)::saving
 INTEGER::year,month,day_of_week,jdate,jhour
 REAL(iwp),INTENT(in)::d,jd,long,timezone
 REAL(iwp),INTENT(out)::d2
 REAL(iwp)::tz,sintz,costz,sin2tz,cos2tz,sin3tz,cos3tz,eqr,eqm
 REAL(iwp)::day,jd1,jd2
 !---> IF saving==1, define sumer time
 ! starts on the last Sunday in March (jd1)
 ! ends on the las Sunay in October (jd2)
 IF(saving==1)THEN
   CALL civil_date(jd,year,month,day)
   jdate=year*10000+331
   jhour=1200
   CALL julian_day(jdate,jhour,jd1)
   day_of_week=INT(jd1-INT(jd1/7.0)*7)+1
   IF(day_of_week/=7)jd1=jd1-day_of_week
   jdate=year*10000+1031
   jhour=1200
   CALL julian_day(jdate,jhour,jd2)
   day_of_week=INT(jd2-INT(jd2/7.0)*7)+1
   IF(day_of_week/=7)jd2=jd2-day_of_week
 END IF
 !
 tz=2.0_iwp*pi*d/365
 sintz=SIN(tz)
 costz=COS(tz)
 sin2tz=2.*sintz*costz
 cos2tz=costz*costz-sintz*sintz
 sin3tz=sintz*cos2tz+costz*sin2tz
 cos3tz=costz*cos2tz-sintz*sin2tz
!---> Equation of time  in radians
 eqr=0.000075+0.001868*costz-0.032077*sintz-0.014615*cos2tz-0.040849*sin2tz
!---> convert equation of time to minutes:
 eqm=229.18*eqr
!---> calculate local hour (in days):
 IF(saving==1)THEN
   IF(jd>=jd1.and.jd<jd2)THEN
     d2=d+(eqm+4.0_iwp*long-60.0_iwp*(timezone+1))/1440.0_iwp
     IF(d2<0.0_iwp)d2=365.0_iwp+d2
   ELSE
     d2=d+(eqm+4.0_iwp*long-60.0_iwp*timezone)/1440.0_iwp
     IF(d2<0.0_iwp)d2=365.0_iwp+d2
   END IF
 ELSE
   d2=d+(eqm+4.0_iwp*long-60.0_iwp*timezone)/1440.0_iwp
   IF(d2<0.0_iwp)d2=365.0_iwp+d2
 END IF
RETURN
END SUBROUTINE solar_time

!
!-------------------------------------------------------------------------
! Reservoir temperature boundary conditions
!
!------------------Glossary of input variables----------------------------
!
! Scalar reals:
! wl               reservoir level
! ts               yearly mean water temperature at the surface
! ao               amplitude of annual variation at the surface
! tauo             day of maximumwater temperature at the surface
! tb               yearly mean water temperature at the bottom
! c1 to c5         Bofang's formula parameters
!-------------------------------------------------------------------------
SUBROUTINE reservoir_parameters
  IMPLICIT NONE
  READ(10,*)wl,ts,ao,tauo,ce,e1,e2,e3,e4,e5
  WRITE(11,'(/a,E11.4)')"Reservoir level                                = ",wl
  WRITE(11,'(a,E11.4)')"Yearly mean water temperature at the surface   = ",ts
  WRITE(11,'(a,E11.4)')"Amplitude of annual variation at the surface   = ",ao
  WRITE(11,'(a,E11.4)')"Day of maximumwater air temperature            = ",tauo
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter c                   = ",ce
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter e1                  = ",e1
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter e2                  = ",e2
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter e3                  = ",e3
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter e4                  = ",e4
  WRITE(11,'(a,E11.4)')"Bofang's formula parameter e5                  = ",e5
END SUBROUTINE reservoir_parameters
!
FUNCTION reservoir_temperature(z,d)RESULT(t)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::z,d
 REAL(iwp)::t,y,tm,a,t0
 y=wl-z
 IF(y<zero) y=zero
 tm=ce+(ts-ce)*EXP(-e1*y)
 a=-ao*EXP(-e2*y)
 t0=(e3-e4*EXP(-e5*y))*365.0/12.0+tauo
 t=tm+a*COS(2*pi/365.0*(d-t0))
END FUNCTION reservoir_temperature
!
!-------------------------------------------------------------------------
! Other temperature boundary conditions
!
!------------------Glossary of input variables----------------------------
!
! Scalar integer:
! jdate            date in yyyymmdd format
! jhour            time in hhmm format
! 
! Scalar real:
! temp             temperature
!-------------------------------------------------------------------------
SUBROUTINE temp_table(ifile)
  IMPLICIT NONE
  INTEGER::i,icont,io,jdate,jhour
  REAL(iwp)::temp,jd
  CHARACTER(len=100),INTENT(in)::ifile
  OPEN(20,FILE=ifile,STATUS='OLD')
  icont=0
  DO
    READ(20,*,IOSTAT=io)jdate,jhour,temp
    IF(io<0)THEN
      EXIT
    ELSE
      icont=icont+1
    END IF
  END DO
  ALLOCATE(table4(icont,2))
  REWIND 20
  DO i=1,icont
    READ(20,*,IOSTAT=io)jdate,jhour,temp
    CALL Julian_Day(jdate,jhour,jd)
    table4(i,1)=jd
    table4(i,2)=temp
  END DO
END SUBROUTINE temp_table
!
!-------------------------------------------------------------------------
! Solar radiation boundary conditions
!-------------------------------------------------------------------------
SUBROUTINE radiation_parameters(utemp)
 CHARACTER(len=1),INTENT(in)::utemp
 WRITE(11,'(/a)')"Exponential function"
 READ(10,*)aa,bb
 WRITE(11,'(a,E11.4,a,E11.4,3a)')                                        &
 "Ib/cos Z = Io exp(",aa,"+",bb," cos Z) in J/(",utemp," m2)"
 SELECT CASE(utemp)
 CASE("d")
   Io=1367*86400.0_iwp
 CASE("h")
   Io=1367*3600.0_iwp
 CASE("s")
   Io=1367
 END SELECT
END SUBROUTINE radiation_parameters
!
FUNCTION exponential_function(cos_z)RESULT(ih)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::cos_z
 REAL(iwp)::ih
 ih=Io*EXP(aa+bb*cos_z)
END FUNCTION exponential_function
!
SUBROUTINE rad_table(ifile)
  IMPLICIT NONE
  INTEGER::i,icont,io,jdate,jhour
  REAL(iwp)::ibeam,idiffuse,jd
  CHARACTER(len=100),INTENT(in)::ifile
  OPEN(20,FILE=ifile,STATUS='OLD')
  icont=0
  DO
    READ(20,*,IOSTAT=io)jdate,jhour,ibeam,idiffuse
    IF(io<0)THEN
      EXIT
    ELSE
      icont=icont+1
    END IF
  END DO
  ALLOCATE(table3(icont,3))
  REWIND 20
  DO i=1,icont
    READ(20,*,IOSTAT=io)jdate,jhour,ibeam,idiffuse
    CALL julian_day(jdate,jhour,jd)
    table3(i,1)=jd
    table3(i,2)=ibeam
    table3(i,3)=idiffuse
  END DO
END SUBROUTINE rad_table
!
SUBROUTINE solar_position(d,fi,delta,omega,cos_z)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::d,fi
 REAL(iwp)::dr,tz,cos_z,ut
 REAL(iwp),INTENT(out)::delta,omega
 REAL(iwp),PARAMETER::pi=3.141592653589793238_iwp,zero=0.0_iwp
 ut=d-AINT(d)
 dr=pi/180.0_iwp
!----day angle
 tz=2.0_iwp*pi*d/365.0_iwp
!----solar declination
 delta=rdecl(tz)
!----solar hour angle
 omega=15.0*(ut*24.0-12.)*dr
!---cosine of solar zenith angle
 cos_z=SIN(fi)*SIN(delta)+COS(fi)*COS(delta)*COS(omega)
 IF(cos_z<zero)cos_z=0.0
END SUBROUTINE solar_position
!
SUBROUTINE ray_parametric_equation(phi,azimuth,delta,omega,cos_z,dir)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::phi,azimuth,delta,omega,cos_z
 REAL(iwp)::alpha_r,azimuth_s,cos_a,sin_a,sin_z,solar_zenith_angle
 REAL(iwp),INTENT(out)::dir(:)
 REAL(iwp),PARAMETER::pi=3.141592653589793238_iwp,zero=0.0_iwp
   solar_zenith_angle=ACOS(cos_z) 
   sin_z=SIN(solar_zenith_angle)
   azimuth_s=(COS(omega)*COS(delta)*SIN(phi)-SIN(delta)*COS(phi))/sin_z
   azimuth_s=ACOS(MAX(-1.0_iwp,MIN(1.0_iwp,azimuth_s)))
   IF(omega<zero)azimuth_s=-azimuth_s
   alpha_r=azimuth-azimuth_s
   sin_a=SIN(alpha_r)
   cos_a=COS(alpha_r)
   dir(1)=-sin_a*sin_z
   dir(2)=cos_a*sin_z
   dir(3)=cos_z
END SUBROUTINE ray_parametric_equation
!
FUNCTION incident_angle(delta,omega,phi,yy,beta)RESULT(cos_alpha)
 IMPLICIT NONE
 REAL(iwp),INTENT(IN)::delta,omega,phi,yy,beta
 REAL(iwp)::aa,bb,cc,cos_alpha
 REAL(iwp),PARAMETER::pi=3.141592653589793238_iwp,zero=0.0_iwp 
 aa=COS(yy)*SIN(phi)-SIN(yy)*COS(phi)*COS(beta)
 bb=COS(yy)*COS(phi)+SIN(yy)*SIN(phi)*COS(beta)
 cc=SIN(yy)*SIN(beta)
 cos_alpha=(aa*SIN(delta)+bb*COS(omega)*COS(delta)+cc*SIN(omega)*        &
   COS(delta)) 
 IF(cos_alpha<zero) cos_alpha=zero 
END FUNCTION incident_angle
!
SUBROUTINE Kumar_model(d,cos_z,ibeam,idiffuse)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::d,cos_z
 REAL(iwp)::epsilon,G0,m,Tb,Td
 REAL(iwp),INTENT(out)::ibeam,idiffuse
!----Sun-Earth distance
!----correction factor
 epsilon=1.0+0.0344*COS(2.0*pi/365.25*d-0.0488869)
!----extraterrestrial irradiance
 G0=epsilon*1367.0  
!---optical air mass
 m=SQRT(1229.0+(614.0*cos_z)*(614.0*cos_z))-614.0*cos_z  
!---atmosferic atenuation
 Tb=0.56*(EXP(-0.65*m)+EXP(-0.095*m))  
!
 IF(cos_z < zero) THEN
   ibeam=zero
 ELSE
   ibeam=G0*Tb
 END IF
 Td=0.271-0.294*Tb
 idiffuse=1367.0*Td*cos_z
 RETURN
END SUBROUTINE Kumar_model
!
FUNCTION rdecl(tz)RESULT(decl)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::tz
 REAL(iwp)::decl,sintz,costz,sin2tz,cos2tz,sin3tz,cos3tz
  sintz = SIN(tz)
  costz = COS(tz)
  sin2tz = 2.*sintz*costz
  cos2tz = costz*costz-sintz*sintz
  sin3tz = sintz*cos2tz + costz*sin2tz
  cos3tz = costz*cos2tz - sintz*sin2tz
!---> declination in radians
  decl = 0.006918 - 0.399912*costz  + 0.070257*sintz   &
                  - 0.006758*cos2tz + 0.000907*sin2tz  &
                  - 0.002697*cos3tz + 0.001480*sin3tz  
END FUNCTION rdecl
!
SUBROUTINE form_shadow_table(filename,azimuth,phi,g_num,g_coord,iflux,   &
 inci_s,points_s)
 IMPLICIT NONE
 INTEGER,INTENT(in)::g_num(:,:),iflux(:,:),inci_s(:,:)
 REAL(iwp),INTENT(in)::azimuth,phi,g_coord(:,:),points_s(:,:)
 CHARACTER(len=100),INTENT(in)::filename
 INTEGER::hfbc,i,ibound,icont,iel,iflag,iside,j,jdate,jhour,ndim,nip_s,  &
   nod_s,hfbc2
 REAL(iwp)::beta,cos_alpha,cos_z,d1,d2,delta,omega,det,n1,n2,n3,yy
 INTEGER,ALLOCATABLE::num_s(:),iflux2(:,:)
 REAL(iwp),ALLOCATABLE::coord_s(:,:),dir(:),fun_s(:),gc(:),jac_s(:,:),   &
    der_s(:,:)
 CHARACTER(len=100)::ifile
 LOGICAL::flag
!-------------------------------------------------------------------------
 INTEGER,ALLOCATABLE::incidencia2(:,:)
 INTEGER::ivec1(4),ivec2(4),ivec3(4),ivec4(4),ivec5(4),ivec6(4),cell_type2
 INTEGER::num(20),nn
!-------------------------------------------------------------------------
 nod_s=UBOUND(inci_s,2)  
 nip_s=UBOUND(points_s,1)
 hfbc=UBOUND(iflux,1)
 nn=UBOUND(g_coord,2)
 ndim=3
 ALLOCATE(num_s(nod_s),coord_s(nod_s,ndim),dir(ndim),fun_s(nod_s),gc(ndim),&
   jac_s(ndim-1,ndim),der_s(ndim-1,nod_s))
 ALLOCATE(table5(8760,hfbc,nip_s))
 READ(10,*)hfbc2
 ALLOCATE(iflux2(hfbc2,2))
 IF(hfbc2>0)THEN
   WRITE(11,'(/a,/a)')"Extra elements for shadow","   Element     iside"
   READ(10,*)(iflux2(i,:),i=1,hfbc2)
   WRITE(11,'(2i10)')(iflux2(i,:),i=1,hfbc2)
   OPEN(22,file='shadow.vtk',status='replace',action='write')  
   ALLOCATE(incidencia2(hfbc+hfbc2,4))
   ivec1=(/13,1,3,15/)
   ivec2=(/7,19,17,5/)
   ivec3=(/1,7,5,3/)
   ivec4=(/19,13,15,17/)
   ivec5=(/13,19,7,1/)
   ivec6=(/3,5,17,15/)
   boundary0: DO ibound=1,hfbc
     iel=iflux(ibound,1) 
     iside=iflux(ibound,2) 
     num=g_num(:,iel)
     SELECT CASE ( iside)
       CASE(1)
         incidencia2(ibound,:)=num(ivec1)
       CASE(2)
         incidencia2(ibound,:)=num(ivec2)
       CASE(3)
         incidencia2(ibound,:)=num(ivec3)
       CASE(4)
         incidencia2(ibound,:)=num(ivec4)
       CASE(5)
         incidencia2(ibound,:)=num(ivec5)
       CASE(6)
         incidencia2(ibound,:)=num(ivec6)
     END SELECT 
   END DO boundary0
   boundary1: DO ibound=1,hfbc2
     iel=iflux2(ibound,1) 
     iside=iflux2(ibound,2) 
     num=g_num(:,iel)
     SELECT CASE ( iside)
       CASE(1)
         incidencia2(hfbc+ibound,:)=num(ivec1)
       CASE(2)
         incidencia2(hfbc+ibound,:)=num(ivec2)
       CASE(3)
         incidencia2(hfbc+ibound,:)=num(ivec3)
       CASE(4)
         incidencia2(hfbc+ibound,:)=num(ivec4)
       CASE(5)
         incidencia2(hfbc+ibound,:)=num(ivec5)
       CASE(6)
         incidencia2(hfbc+ibound,:)=num(ivec6)
     END SELECT 
   END DO boundary1
   cell_type2=9
   WRITE(22,'(a)') "# vtk DataFile Version 3.0"
   WRITE(22,'(a)') "radiation"
   WRITE(22,'(a)') "ASCII"
   WRITE(22,'(a)') "DATASET UNSTRUCTURED_GRID"
   WRITE(22,'(a,i7,a)')"POINTS",nn," float" 
   DO i=1,nn
    WRITE(22,'(3f12.4)')g_coord(:,i)
   END DO
   WRITE(22,'(a,2i10)') "CELLS",hfbc+hfbc2,5*(hfbc+hfbc2)
   DO i=1,hfbc+hfbc2
     WRITE(22,'(a,4i6)') "  4",incidencia2(i,1:4)-1 
   END DO
   WRITE(22,'(a,i5)') "CELL_TYPES",hfbc+hfbc2
   DO i=1,hfbc+hfbc2
     WRITE(22,'(i5)') cell_type2
   END DO
   WRITE(22,'(a,i5)')"CELL_DATA",hfbc+hfbc2
   WRITE(22,'(a)')"SCALARS radiation float 1"
   WRITE(22,'(a)')"LOOKUP_TABLE default"
   DO ibound=1,hfbc
     WRITE(22,'(a)')" 1"
   END DO
   DO ibound=1,hfbc2
     WRITE(22,'(a)')" 2"
   END DO
   CLOSE(22)
 END IF

 icont=0
 DO j=1,365
   DO jhour=1,24
     icont=icont+1
     d1=FLOAT(j)+FLOAT(jhour-1)/24.0_iwp  
     CALL solar_position(d1,phi,delta,omega,cos_z)
     IF(cos_z>zero)THEN
       CALL ray_parametric_equation(phi,azimuth,delta,omega,cos_z,dir)
       boundary: DO ibound=1,hfbc
         iel=iflux(ibound,1) ; iside=iflux(ibound,2) 
         num_s=g_num(inci_s(iside,:),iel) 
         coord_s=TRANSPOSE(g_coord(:,num_s))
         gauss_pts: DO i=1,nip_s
           iflag=1
           CALL shape_fun(fun_s,points_s,i)
           CALL shape_der(der_s,points_s,i)
           jac_s=MATMUL(der_s,coord_s)
           CALL det_s(jac_s,n1,n2,n3,det)
           beta=azimuth+(pi*0.5_iwp-ATAN2(n2,n1))
           yy=ACOS(n3)
           cos_alpha=incident_angle(delta,omega,phi,yy,beta)
           IF(cos_alpha>zero) THEN
             gc=MATMUL(fun_s,coord_s) 
             CALL ray_face_intersection(g_num,g_coord,iflux,iflux2,inci_s,ibound,&
               gc,dir,flag)  
             IF(flag)iflag=0 
             table5(icont,ibound,i)=iflag
           ELSE
             table5(icont,ibound,i)=0
           END IF
         END DO gauss_pts
       END DO boundary
     ELSE
       iflag=0
       table5(icont,:,:)=iflag
     END IF
   END DO
 END DO
 ifile=TRIM(filename)//'_shadow.dat'
 OPEN(50,file=ifile,form='unformatted',status='replace',action='write')
 WRITE(50)table5
 CLOSE(50)
END SUBROUTINE form_shadow_table
!  
SUBROUTINE shadow_table(ifile,hfbc,nip_s)
  IMPLICIT NONE
  INTEGER,INTENT(in)::hfbc,nip_s
  CHARACTER(len=100),INTENT(in)::ifile
  ALLOCATE(table5(8760,hfbc,nip_s))
  OPEN(50,FILE=ifile,form='unformatted',status='old')
  READ(50)table5
  CLOSE(50)
END SUBROUTINE shadow_table
!
!-------------------------------------------------------------------------
! Convection boundary conditions
!
!------------------Glossary of input variables----------------------------
!
! Scalar integer:
! jdate            date in yyyymmdd format
! jhour            time in hhmm format
! 
! Scalar real:
! t_m1             yearly mean temperature
! t_p1             amplitude of annual temperature variation
! delay_1          yearly temperature phase difference
! t_m2             yearly mean amplitude 
! t_p2             amplitude of annual variation of the daily amplitude
! delay_2          yearly phase difference of the daily amplitude
! delay_3          daily temperature phase difference
! temp             temperature
! htco             heat convection coefficient
!-------------------------------------------------------------------------
SUBROUTINE air_parameters
  WRITE(11,'(/a)')"Harmonic function"
  READ(10,*)t_m1,t_p1,delay_1,t_m2,t_p2,delay_2,delay_3
  WRITE(11,'(a,E11.4)')"Yearly mean temperature                        = ",t_m1
  WRITE(11,'(a,E11.4)')"Amplitude of annual temperature variation      = ",t_p1
  WRITE(11,'(a,E11.4)')"Yearly temperature phase difference            = ",delay_1
  WRITE(11,'(a,E11.4)')"Yearly mean amplitude                          = ",t_m2
  WRITE(11,'(a,E11.4)')"Amplitude of the daily amplitude variation     = ",t_p2
  WRITE(11,'(a,E11.4)')"Yearly phase difference of the daily amplitude = ",delay_2
  WRITE(11,'(a,E11.4)')"Daily temperature phase differencee            = ",delay_3
END SUBROUTINE air_parameters
!
SUBROUTINE air_table1(ifile,airin)
  IMPLICIT NONE
  INTEGER,INTENT(in)::airin
  INTEGER::i,icont,io,jdate,jhour
  REAL(iwp)::temp,jd,htco
  CHARACTER(len=100),INTENT(in)::ifile
  OPEN(20,FILE=ifile,STATUS='OLD')
  icont=0
  DO
    READ(20,*,IOSTAT=io)jdate,jhour,temp
    IF(io<0)THEN
      EXIT
    ELSE
      icont=icont+1
    END IF
  END DO
  SELECT CASE(airin)
  CASE(2)
    ALLOCATE(table1(icont,3))
  CASE(3)
    ALLOCATE(table1(icont,2)) 
  END SELECT
  REWIND 20
  DO i=1,icont
    SELECT CASE(airin)  
    CASE(2)
      READ(20,*,IOSTAT=io)jdate,jhour,temp,htco 
      CALL Julian_Day(jdate,jhour,jd)
      table1(i,1)=jd
      table1(i,2)=temp
      table1(i,3)=htco
    CASE(3)
      READ(20,*,IOSTAT=io)jdate,jhour,temp 
      CALL Julian_Day(jdate,jhour,jd)
      table1(i,1)=jd
      table1(i,2)=temp
    END SELECT
  END DO
END SUBROUTINE air_table1
!
SUBROUTINE air_table2(ifile)
  IMPLICIT NONE
  INTEGER::i,icont,io,jdate,jhour
  REAL(iwp)::temp,jd,htco
  CHARACTER(len=100),INTENT(in)::ifile
  OPEN(20,FILE=ifile,STATUS='OLD')
  icont=0
  DO
    READ(20,*,IOSTAT=io)jdate,jhour,temp
    IF(io<0)THEN
      EXIT
    ELSE
      icont=icont+1
    END IF
  END DO
  ALLOCATE(table2(icont,3))
  REWIND 20
  DO i=1,icont
    READ(20,*,IOSTAT=io)jdate,jhour,temp,htco 
    CALL Julian_Day(jdate,jhour,jd)
    table2(i,1)=jd
    table2(i,2)=temp
    table2(i,3)=htco
  END DO
END SUBROUTINE air_table2
!
FUNCTION harmonic_function(d)RESULT(t_a)
 IMPLICIT NONE
 REAL(iwp),INTENT(in)::d
 REAL(iwp)::t_a,omega,t_0,t_b,t_e,t_p
 t_0=d
 omega=2.0_iwp*pi/(365.0_iwp)
 t_b=t_m1-t_p1*(COS(omega*(t_0-delay_1)))
!
 t_e=t_m2-t_p2*(COS(omega*(t_0-delay_2)))
!
 t_p=t_e/2.0_iwp
 omega=2.0_iwp*pi
 t_a=t_b-t_p*COS(omega*(t_0-delay_3))
END FUNCTION harmonic_function
!
!-------------------------------------------------------------------------
! Extract values from tables
!-------------------------------------------------------------------------
FUNCTION find_in_table(table,jd)RESULT(dd)
 IMPLICIT NONE
 INTEGER::i,n
 REAL(iwp),INTENT(in)::jd,table(:,:)
 REAL(iwp)::dd
 n=UBOUND(table,1)
 DO i=1,n
   IF(jd>table(i+1,1))cycle
   dd=table(i,2)+(table(i+1,2)-table(i,2))*(jd-table(i,1))/              &
      (table(i+1,1)-table(i,1))
   EXIT   
 END DO
END FUNCTION find_in_table
!
FUNCTION find_in_table2(table,jd)RESULT(dd)
 IMPLICIT NONE
 INTEGER::i,n
 REAL(iwp),INTENT(in)::jd,table(:,:)
 REAL(iwp)::dd
 n=UBOUND(table,1)
 DO i=1,n
   IF(jd>table(i+1,1))cycle
   dd=table(i,3)+(table(i+1,3)-table(i,3))*(jd-table(i,1))/              &
      (table(i+1,1)-table(i,1))
   EXIT   
 END DO
END FUNCTION find_in_table2
!
FUNCTION find_in_table5(table,d2,ibound,i)RESULT(dd)
 IMPLICIT NONE
 INTEGER,INTENT(in)::ibound,i
 INTEGER::irow
 REAL(iwp),INTENT(in)::d2
 INTEGER,INTENT(in)::table(:,:,:)
 REAL(iwp)::dd,dt,one=1.0_iwp
 irow=INT(d2*24.0_iwp+1) !; write(*,*)"d2",d2,"irow",irow , UBOUND(table,1)
 dt=(d2-(FLOAT(irow)-one)/24.0_iwp)*24 !; write(*,*)d2, (FLOAT(irow)-one)/24.0_iwp+one,dt
 IF(irow<8760)THEN
   dd=table(irow,ibound,i)+(table(irow+1,ibound,i)-table(irow,ibound,i))*dt
 ELSE
   dd=table(8760,ibound,i)+(table(1,ibound,i)-table(8760,ibound,i))*dt 
 END IF  
END FUNCTION find_in_table5

!
SUBROUTINE ray_triangle_intersect(orig,dir,v0,v1,v2,flag)
! From: 
!https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(in)::orig(:),dir(:),v0(:),v1(:),V2(:)
 LOGICAL,INTENT(out)::flag
 REAL(iwp)::d,Epsilon=1e-10_iwp,NdotRayDirection,t,zero=0.0_iwp
 REAL(iwp)::c(3),edge(3),n(3),p(3),v0v1(3),v0v2(3),vp(3)
 ! compute the plane's normal
   v0v1=v1-v0
   v0v2=v2-v0
 ! no need to normalize
   n = Vector_Product(v0v1,v0v2)
 ! Step 1: finding P
 ! Check if the ray and plane are parallel.
   NdotRayDirection=DOT_PRODUCT(n,dir)
   IF(ABS(NdotRayDirection)<Epsilon)THEN ! almost 0
     flag=.FALSE.  !they are parallel, so they don't intersect! 
     RETURN
   END IF
 ! compute d parameter 
   d=-DOT_PRODUCT(n,v0)
 ! compute t 
   t=-(DOT_PRODUCT(n,orig) + d)/NdotRayDirection
 ! check if the triangle is behind the ray
   IF(t < zero)THEN
     flag=.FALSE.  ! the triangle is behind 
     RETURN
   END IF
 ! compute the intersection point 
    p=orig+t*dir
 ! Step 2: inside-outside test
 ! edge 0
   edge=v1-v0
   vp=p-v0
   c=Vector_Product(edge,vp)
   IF(DOT_PRODUCT(n,c)<zero)THEN
     flag=.FALSE.  ! P is on the right side   
     RETURN
   END IF
 ! edge 1
   edge=v2-v1
   vp=p-v1;
   c=Vector_Product(edge,vp)
   IF(DOT_PRODUCT(n,c)<zero) THEN
     flag=.FALSE. ! P is on the right side   
     RETURN
   END IF
 ! edge 2
   edge=v0-v2
   vp=p-v2
   c=Vector_Product(edge,vp)
   IF(DOT_PRODUCT(n,c)<zero)THEN
     flag=.FALSE. !P is on the right side
     RETURN
   END IF
   flag=.TRUE.
END SUBROUTINE ray_triangle_intersect
!
SUBROUTINE ray_face_intersection(g_num,g_coord,iflux,iflux2,inci_s,ibound,orig,dir,flag)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(in)::ibound,g_num(:,:),iflux(:,:),iflux2(:,:),inci_s(:,:)
 REAL(iwp),INTENT(in)::g_coord(:,:),orig(:),dir(:)
 LOGICAL,INTENT(out)::flag
 INTEGER::jbound,hfbc,hfbc2,iel,iside,ndim,nod_s
 REAL(iwp)::zero=0.0_iwp
 REAL(iwp),ALLOCATABLE::num(:),coord(:,:)
 nod_s=UBOUND(inci_s,2)  
 hfbc=UBOUND(iflux,1)
 hfbc2=UBOUND(iflux2,1)
 ndim=3
 ALLOCATE(num(nod_s),coord(nod_s,ndim))
   flag=.FALSE.
   jbound=0
   IF(ibound>1)jbound=1
   DO WHILE(.not.flag)
     jbound=jbound+1 
     IF(jbound==ibound)jbound=jbound+1
     IF(jbound>hfbc)EXIT
     iel=iflux(jbound,1) ; iside=iflux(jbound,2)
     num=g_num(inci_s(iside,:),iel); coord=TRANSPOSE(g_coord(:,num))
     !triangle 1
     CALL ray_triangle_intersect(orig,dir,coord(1,:),coord(5,:),coord(3,:),flag)
     IF(flag)RETURN
     !triangle 2
     CALL ray_triangle_intersect(orig,dir,coord(1,:),coord(7,:),coord(5,:),flag)
     IF(flag)RETURN
   END DO
   IF(hfbc2>0)THEN
     flag=.FALSE.
     jbound=0
     DO WHILE(.not.flag)
       jbound=jbound+1 
       IF(jbound>hfbc2)EXIT
       iel=iflux2(jbound,1) ; iside=iflux2(jbound,2)
       num=g_num(inci_s(iside,:),iel); coord=TRANSPOSE(g_coord(:,num))
       !triangle 1
       CALL ray_triangle_intersect(orig,dir,coord(1,:),coord(5,:),coord(3,:),flag)
       IF(flag)RETURN
       !triangle 2
       CALL ray_triangle_intersect(orig,dir,coord(1,:),coord(7,:),coord(5,:),flag)
       IF(flag)RETURN
     END DO

   END IF    
 END SUBROUTINE ray_face_intersection 
!
END MODULE environmental_routines_v1_3