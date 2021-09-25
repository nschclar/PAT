MODULE environmental_routines
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
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
 PUBLIC::table1,table2,table3,table4
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
 ut=d-AINT(d)
 dr=pi/180.0_iwp
!----day angle
 tz=2.0_iwp*pi*d/365.0_iwp
!----solar declination
 delta=rdecl(tz)
!----solar hour angle
 omega=-15.0*(ut*24.0-12.)*dr
!---cosine of solar zenith angle
 cos_z=SIN(fi)*SIN(delta)+COS(fi)*COS(delta)*COS(omega)
 IF(cos_z<zero)cos_z=0.0
END SUBROUTINE solar_position
!
FUNCTION incident_angle(delta,omega,phi,yy,beta)RESULT(cos_alpha)
 IMPLICIT NONE
 REAL(iwp),INTENT(IN)::delta,omega,phi,yy,beta
 REAL(iwp)::aa,bb,cc,cos_alpha
 aa=COS(yy)*SIN(phi)-SIN(yy)*COS(phi)*COS(beta)
 bb=COS(yy)*COS(phi)+SIN(yy)*SIN(phi)*COS(beta)
 cc=SIN(yy)*SIN(beta)
 cos_alpha=(aa*SIN(delta)+bb*COS(omega)*COS(delta)-cc*SIN(omega)*        &
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
  WRITE(11,'(a,E11.4)')"Amplitude of the daily amplitude variation     = ",tp_2
  WRITE(11,'(a,E11.4)')"Yearly phase difference of the daily amplitude = ",delay_2
  WRITE(11,'(a,E11.4)')"Daily temperature phase differencee            = ",delay_3
END SUBROUTINE air_parameters
!
SUBROUTINE air_table1(ifile)
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
  ALLOCATE(table1(icont,3))
  REWIND 20
  DO i=1,icont
    READ(20,*,IOSTAT=io)jdate,jhour,temp,htco 
    CALL Julian_Day(jdate,jhour,jd)
    table1(i,1)=jd
    table1(i,2)=temp
    table1(i,3)=htco
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
END MODULE environmental_routines