PROGRAM PAT_v1_3
!-------------------------------------------------------------------------
! THEMAL ANALISYS OF STRUCTURES SUBJECTED TO ENVIRONMETAL ACTION
!-------------------------------------------------------------------------
!
! Created by
!
! Noemi Alejandra Schclar Leitao
! Eloisa Maria Castilho dos Santos
!
! Laboratorio Nacional de Engenharia Civil
! Portugal
!------------------Version 1_3--------------------------------------------
! Modify by Noemi Alejandra Schclar Leitao
!    
! 06/03/2022 - Local to solar time convertion 
! 26/05/2023 - Shadow (3D only)
! 30/05/2023 - Discrete air temperature with constant total thermal 
!              coefficient (airin=3)
!------------------Glossary of variable names-----------------------------
!
! Scalar integers:
! airin            type of input air temperatures
! fixed_freedoms_1 number of fixed temperatures (reservoir temperatures)
! fixed_freedoms_2 number of fixed temperatures (other cases)
! fluxin           type of input flux data
! hfbc             number of boundaries with prescibed flux (solar 
!                  radiation)
! htbc             number of boundaries with convection heat transfer
! i,ibound,iel     simple counters
! igroup           group identifier for convection data
! indic            if indic=0 all nodes have the same initial temperature,
!                  if indic>0 variavel initial temperature
! iout             unit identifier of the unit connected to the file
! ishadow          shadow index    
! iside            number of the loaded side
! iwp              SELECTED_REAL_KIND(15)
! j                simple counter
! jdate            date in format yyyymmdd
! jhour            hour in format hhmm
! ndim             number of dimensions
! nels             number of elements
! neq              number of degree of freedom in the mesh
! nfiles           number of convection boundaries 
! nip              number of integrating points per element
! nip_s            number of integrating points of the associated line or
!                  surface element
! nn               number of nodes in the mesh
! nod              number of nodes per element
! nod_s            number of nodes of the associated element
! npri             output printed every npri time steps
! nprops           number of material properties
! np_types         number of different property types
! nsensors         number of points with output for every time steps
! nside            number of element faces
! nstep            number of time steps
! saving           if saving=1, summer time is considered
! shadowin         if shadowin=1, compute table,
!                  if shadowin=2, read values from a file
! tempin           type of input fixed temperatures (other temperatures)
!
! Scalar reals
! ab               absorption coefficient
! azimuth          azimuth of the reference axes
! beta             azimuth of the surface normal
! cos_z            cosine of the solar zenith angle
! cos_alpha        cosine of the incidence angle
! cp               specific heat
! d1,d2            decimal day of the year
! delta            declination angle
! det              determinant of the Jacobian matrix
! dtim             calculation time step
! htco,htco1,htco2 convection coefficient
! ibeam            beam radiation component
! idiffuse         diffuse radiation component
! jd,jd_i,jd_f     Julian days
! long             Earth's longitude in decimal degrees
! n1,n2,n3         unit normal components
! omega            solar hour angle
! penalty          set to 1.0E+20
! phi              Earth's latitude in decimal degrees
! rho              density
! t_a,t_a1,t_a2    air temperature
! theta            time-integration weighting parameter
! timezone         local standards to which local clocks are set given by
!                  the time zone (timezone) in hours from UTC
! t_m,t_p,delay    harmonic function parameters
! ucte             time conversion factor 
! val0, val1       initial temperatures
! val2             output point temperature
! yy               tilt angle
! 
! Scalar characters
! element          element type
! element_s        associated line or surface element
! filename         file name
! ifile            file name and extension
! char_num         keep an integer as a string varaivel
! utemp            time units
!
! Scalar logical
! shadow           set to .TRUE. if the ray-tracing algorithm is applied    
!
! Dynamic integer arrays
! etype            element property types
! g_num            node number for all elements
! ielsen           output point element numbers
! iflux            element and element side with prescribed flux
! inci_s           node order for the associated line or surface element
! itrans           element and element side with convection heat transfer
! kdiag            diagonal term locations
! node_1           nodes with fixed temperatures (reservoir temperatures)
! node_2           nodes with fixed temperatures (other temperatures)
! node_3           nodes with variavel intial temperature
! num              element node numbers
! num_s            element node numbers for the associated line or surface
!                  element
!
! Dynamic real arrays
! bk               global heat stiffness matrix
! bkcopy           stores time indepent global conductivity matrix
! bp               global capacitance matrix
! bpcopy           stores initial global capacitance matrix
! coord            element nodal coordinates
! coord_s          associated element nodal coordinates
! der              shape function derivatives with respect to local
!                  coordinates
! der_s            shape function derivatives with respect to local
!                  coordinates for associated elements
! deriv            shape function derivatives with respet to global
!                  coordinates
! dir              direction of the solar ray
! elt              element nodal temperatures
! elt_s            associated element nodal temperatures
! fun              shape functions
! fun_s            shape functions of the associated element
! gc               integrating point coordinates
! g_coord          nodal coordinates for all elements
! jac              Jacobian matrix
! jac_s            Jacobian matrix of the associated element
! kay              thermal conductivity matrix
! kp               element heat stiffness matrix
! kp_s             associated element heat stiffness matrix
! loads            temperature values
! loadsh           convective load vector
! loadsq           radiation load vector
! newload          load vector for the present time step
! ntn              cross product of shape functions
! ntn_s            cross product of shape functions of the associated
!                  element
! oldlo            load vector for the previous time step
! pm               element capacitance matrix
! points           integrating point local coordinates
! points_s         integrating point local coordinates for the associated
!                  element
! prop             elemt properties
! sensors          output points coordinates
! storbp_1         stores augmented diagonal terms (reservoir temperatures)
! storbp_2         stores augmented diagonal terms (other temperatures)
! value_1          fixed boundary values (reservoir temperatures)
! value_2          fixed boundary values (other temperatures)
! value_3          initial temperature values
! weights          weighting coefficients
! weights_s        weighting coefficients for the associated element
!-------------------------------------------------------------------------
 USE general_routines  ;  USE environmental_routines_v1_3
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::airin,fixed_freedoms_1,fixed_freedoms_2,fluxin,hfbc,htbc,i,    &
   ibound,iel,igroup,indic,iout,iside,ishadow,j,jdate,jhour,ndim,nels,   &
   neq,nfiles,nip,nip_s,nn,nod,nod_s,npri,nprops,np_types,nsensors,nside,&
   nstep,saving,shadowin,tempin
 INTEGER::iflag,icont
 REAL(iwp)::ab,azimuth,beta,cos_alpha,cos_z,cp,d1,d2,delta,det,dtim,htco,&
   htco1,htco2,ibeam,idiffuse,jd,jd_i,jd_f,long,n1,n2,n3,omega,          &
   penalty=1.0e20_iwp,phi,rho,t_a,t_a1,t_a2,theta,timezone,ucte,val0,    &
   val1,val2,yy
 CHARACTER(len=15)::element,element_s
 CHARACTER(len=100)::filename,ifile
 CHARACTER(len=1)::utemp
 CHARACTER(len=14)::char_num
 LOGICAL::shadow
!-------------------dynamic arrays----------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g_num(:,:),ielsen(:),iflux(:,:),          &
   inci_s(:,:),itrans(:,:),kdiag(:),node_1(:),node_2(:),node_3(:),num(:),&
   num_s(:)
 REAL(iwp),ALLOCATABLE::bk(:),bkcopy(:),bp(:),bpcopy(:),coord(:,:),      &
   coord_s(:,:),der(:,:),der_s(:,:),deriv(:,:),dir(:),elt(:),elt_s(:),   &
   fun(:),fun_s(:),gc(:),g_coord(:,:),jac(:,:),jac_s(:,:),kay(:,:),      &
   kp(:,:),kp_s(:,:),loads(:),loadsh(:),loadsq(:),newlo(:),ntn(:,:),     &
   ntn_s(:,:),oldlo(:),pm(:,:),points(:,:),points_s(:,:),prop(:,:),      &
   sensors(:,:),storbp_1(:),storbp_2(:),value_1(:),value_2(:),value_3(:),&
   weights(:),weights_s(:)
!-------------------extra variables for vtk files-------------------------
 INTEGER,ALLOCATABLE::incidencia2(:,:),jvec1(:)
 INTEGER::ivec(20),ivec1(4),ivec2(4),ivec3(4),ivec4(4),ivec5(4),ivec6(4),&
   cell_type,cell_type2,nod_vtk
 REAL(iwp)::rr
!------------------input and initialisation-------------------------------
 WRITE(*,*)'Input the filename:'
 READ(*,*)filename
 ifile=TRIM(filename)//'.dat'
 OPEN(10,FILE=ifile,STATUS='OLD',ACTION='READ')
 ifile=TRIM(filename)//'.res' 
 OPEN(11,FILE=ifile,STATUS='REPLACE',ACTION='WRITE')
!
 READ(10,*)element,nod,nels,nn,nip,ndim,dtim,theta
 WRITE(11,'(a,a,/a,i15,/a,i15,/a,i15,/a,i15,/a,i15,/a,E15.4,/a,E15.4)')  &
 "Element type................................  ",ADJUSTR(element),      &
 "Number of nodes per element.................  ",nod,                   &
 "Number of elements..........................  ",nels,                  &
 "Number of nodes.............................  ",nn,                    &
 "Number of integrating points................  ",nip,                   &
 "Number of dimensions........................  ",ndim,                  &
 "Calculation time step.......................  ",dtim,                  &
 "Time-integration weighting parameter........  ",theta
 READ(10,*)jdate,jhour ; CALL julian_day(jdate,jhour,jd_i)
 READ(10,*)jdate,jhour ; CALL julian_day(jdate,jhour,jd_f)
 char_num=date_hour(jd_i)
 WRITE(11,'(6a,1x,3a)')                                                  &
 "Initial time................................ ",char_num(1:4),'-',      &
 char_num(5:6),'-',char_num(7:8),char_num(9:10),':',char_num(11:12)
 char_num=date_hour(jd_f)
 WRITE(11,'(6a,1x,3a)')                                                  &
 "Final time.................................. ",char_num(1:4),'-',      &
 char_num(5:6),'-',char_num(7:8),char_num(9:10),':',char_num(11:12)
 READ(10,*)utemp
 WRITE(11,'(a,a15)')                                                     &
 "Time units..................................  ",utemp
 SELECT CASE(utemp)
 CASE("d")
   ucte=1.0_iwp
 CASE("h")
   ucte=24.0_iwp
 CASE("s")
   ucte=86400.0_iwp
 CASE DEFAULT
   WRITE(*,*)"the units of time must be day, hour or second"
   STOP
 END SELECT
 nstep=NINT((jd_f-jd_i)*ucte/dtim)
 dtim=(jd_f-jd_i)*ucte/nstep
 READ(10,*)phi,azimuth
 WRITE(11,'(a,E15.4/a,E15.4)')                                           &
 "Earth's latitude............................  ",phi,                   &
 "Azimuth of the reference axes...............  ",azimuth
 phi=phi*pi/180.0_iwp
 azimuth=azimuth*pi/180.0_iwp
 nprops=ndim+4
 neq=nn
 ALLOCATE(coord(nod,ndim),der(ndim,nod),deriv(ndim,nod),dir(ndim),       &
  elt(nod),etype(nels),fun(nod),gc(ndim),g_coord(ndim,nn),               &
  g_num(nod,nels),jac(ndim,ndim),kay(ndim,ndim),kdiag(neq),kp(nod,nod),  &
  num(nod),ntn(nod,nod),pm(nod,nod),points(nip,ndim),weights(nip))
!
 CALL surface(element,nod,element_s,nod_s,nip_s,nside)
 ALLOCATE(coord_s(nod_s,ndim),der_s(ndim-1,nod_s),fun_s(nod_s),          &
  inci_s(nside,nod_s),jac_s(ndim-1,ndim),kp_s(nod_s,nod_s),              &
  ntn_s(nod_s,nod_s),num_s(nod_s),points_s(nip_s,ndim-1),elt_s(nod_s),   &
  weights_s(nip_s))
 CALL num_surface(element,nod,inci_s)
!
 READ(10,*)np_types
 WRITE(11,'(a,i15)')                                                     &
 "Number of different property types..........  ",np_types
 ALLOCATE(prop(nprops,np_types)) ; READ (10,*) prop
 WRITE(11,'(//a)')"Thermal properties"
 IF(ndim==2)WRITE(11,'(a,a)')                                            &
 "Material         kx             ky            rho             cp     ",&
 "       htco            ab"
 IF(ndim==3)WRITE(11,'(/a,a)')                                           &
 "Material         kx             ky             kz            rho     ",&
 "        cp            htco            ab"
 DO i=1,np_types
   WRITE(11,'(i8,7E15.4)')i,prop(:,i)
 END DO
 etype=1
 IF(np_types>1) READ(10,*) etype
 READ(10,*)g_coord 
 READ(10,*)g_num 
 WRITE(11,'(//a)')"Node coordinates"
 IF(ndim==2)THEN
   WRITE(11,'(a)')"    Node          x              y"
   WRITE(11,'(i8,2E15.4)')(i,g_coord(:,i),i=1,nn)
 ELSE
   WRITE(11,'(a)')"    Node          x              y              z"
   WRITE(11,'(i8,3E15.4)')(i,g_coord(:,i),i=1,nn)
 END IF
 WRITE(11,'(//a)')"Element connectivities"
 WRITE(11,'(a)')" Element      Nodes"
 DO i=1,nels
   WRITE(11,'(i8,3x,20i8)')i,g_num(:,i)
 END DO
!------------------integrating point local coordinates--------------------
 CALL sample(element,points,weights)
 CALL sample(element_s,points_s,weights_s)
!------------------vtk boundary conditions data---------------------------
 ivec=(/1,7,19,13,3,5,17,15,8,12,20,9,4,11,16,10,2,6,18,14/)
 ALLOCATE(jvec1(nn))
!------------------reservoir temperature boundary condition---------------
 READ(10,*)fixed_freedoms_1
 IF(fixed_freedoms_1/=0)THEN
   WRITE(11,'(//a)')"Reservoir temperature boundary conditions"
   ALLOCATE(node_1(fixed_freedoms_1),value_1(fixed_freedoms_1),          &
      storbp_1(fixed_freedoms_1))
   CALL reservoir_parameters
   WRITE(11,'(/a,/a)')"It is applied to","    Node"
   READ(10,*)(node_1(i),i=1,fixed_freedoms_1)
   WRITE(11,'(i8)')(node_1(i),i=1,fixed_freedoms_1)
   value_1=g_coord(ndim,node_1)
 !----------------- vtk file----------------------------------------------
   jvec1=0
   jvec1(node_1)=1
   OPEN(22,file='Fixed_freedoms_1.vtk')
   nod_vtk=8
   cell_type=12
   WRITE(22,'(a)') "# vtk DataFile Version 3.0"
   WRITE(22,'(a)') "fixed_freedoms_1"
   WRITE(22,'(a)') "ASCII"
   WRITE(22,'(a)') "DATASET UNSTRUCTURED_GRID"
   WRITE(22,'(a,i7,a)') "POINTS",nn," float"
   DO i=1,nn ; WRITE(22,'(3f12.4)')g_coord(:,i) ; END DO
   WRITE(22,'(a,2i10)') "CELLS",nels,(nod_vtk+1)*(nels)
   DO i=1,nels 
     WRITE(22,'(22i6)') nod_vtk,g_num(ivec(1:8),i)-1 
   END DO
   WRITE(22,'(a,i5)') "CELL_TYPES",nels
   DO i=1,nels ; WRITE(22,'(i5)') cell_type ; END DO
   WRITE(22,'(a,i10)') "POINT_DATA",nn
   WRITE(22,'(a)')"VECTORS vectors float"
   DO i=1,nn
     IF(jvec1(i)==0)THEN
       WRITE(22,'(a)')"0 0 0"
     ELSE
       WRITE(22,'(a)')"1 1 1"   
     END IF
   END DO
   CLOSE(22) 
 ELSE
   WRITE(11,'(//a)')"There is no reservoir temperatures boundary conditions"
 END IF
!------------------other temperature boundary condition-------------------
 READ(10,*)fixed_freedoms_2
 IF(fixed_freedoms_2/=0)THEN
   WRITE(11,'(//a)')"Other temperatures boundary conditions"
   ALLOCATE(node_2(fixed_freedoms_2),value_2(fixed_freedoms_2),           &
     storbp_2(fixed_freedoms_2))
   READ(10,*)tempin
   SELECT CASE(tempin)
   CASE(1)
     WRITE(11,'(a)')"Fixed temperature"
     WRITE(11,'(/a,/a)')"It is applied to","    Node      Value"
     READ(10,*)(node_2(i),value_2(i),i=1,fixed_freedoms_2)
     WRITE(11,'(i8,E11.4)')(node_2(i),value_2(i),i=1,fixed_freedoms_2)
   CASE(2)
     READ(10,*)ifile
     WRITE(11,'(2a)')"Discrete values from ",TRIM(ifile)
     CALL temp_table(ifile)
     WRITE(11,'(/a,/a)')"It is applied to","    Node"
     READ(10,*)(node_2(i),i=1,fixed_freedoms_2)
     WRITE(11,'(i8)')(node_2(i),i=1,fixed_freedoms_2)
   END SELECT
 !----------------- vtk file----------------------------------------------
   jvec1=0
   jvec1(node_2)=1
   OPEN(22,file='Fixed_freedoms_2.vtk')
   nod_vtk=8
   cell_type=12
   WRITE(22,'(a)') "# vtk DataFile Version 3.0"
   WRITE(22,'(a)') "fixed_freedoms_1"
   WRITE(22,'(a)') "ASCII"
   WRITE(22,'(a)') "DATASET UNSTRUCTURED_GRID"
   WRITE(22,'(a,i7,a)') "POINTS",nn," float"
   DO i=1,nn ; WRITE(22,'(3f12.4)')g_coord(:,i) ; END DO
   WRITE(22,'(a,2i10)') "CELLS",nels,(nod_vtk+1)*(nels)
   DO i=1,nels 
     WRITE(22,'(22i6)') nod_vtk,g_num(ivec(1:8),i)-1 
   END DO
   WRITE(22,'(a,i5)') "CELL_TYPES",nels
   DO i=1,nels ; WRITE(22,'(i5)') cell_type ; END DO
   WRITE(22,'(a,i10)') "POINT_DATA",nn
   WRITE(22,'(a)')"VECTORS vectors float"
   DO i=1,nn
     IF(jvec1(i)==0)THEN
       WRITE(22,'(a)')"0 0 0"
     ELSE
       WRITE(22,'(a)')"0 0 1"   
     END IF
   END DO
   CLOSE(22) 
 ELSE
   WRITE(11,'(//a)')"There is no other temperatures boundary conditions"
 END IF
!------------------radiation boundary condition---------------------------
 shadow=.FALSE.
 READ(10,*)hfbc
 IF(hfbc/=0)THEN
   WRITE(11,'(//a)')"Radiation boundary conditions"
   READ(10,*)fluxin,long,timezone,saving,shadow
   IF(ndim==2)shadow=.FALSE.
   WRITE(11,'(/a,E11.4)')                                                &
     "Earth's longitude in decimal degrees           = ",long
   WRITE(11,'(a,E11.4)')                                                 &
     "Time zone                                      = ",timezone
   WRITE(11,'(a,I11)')                                                   &
     "Saving, if saving=1, summer time is considered = ",saving         
   WRITE(11,'(a,L11)')                                                   &
     "Shadow                                         = ",shadow
   SELECT CASE(fluxin)
   CASE(1)
     CALL radiation_parameters(utemp)
   CASE(3)
     READ(10,*)ifile
     WRITE(11,'(2a)')"Discrete values from ",TRIM(ifile)
     CALL rad_table(ifile) 
   END SELECT 
   ALLOCATE(iflux(hfbc,2))
   WRITE(11,'(/a,/a)')"It is applied to","   Element     iside"
   READ(10,*)(iflux(i,:),i=1,hfbc)
   WRITE(11,'(2i10)')(iflux(i,:),i=1,hfbc)
!------------------vtk file----------------------------------------------
   OPEN(22,file='radiation.vtk',status='replace',action='write')  
   ALLOCATE(incidencia2(hfbc,4))
   ivec1=(/13,1,3,15/)
   ivec2=(/7,19,17,5/)
   ivec3=(/1,7,5,3/)
   ivec4=(/19,13,15,17/)
   ivec5=(/13,19,7,1/)
   ivec6=(/3,5,17,15/)
   boundary: DO ibound=1,hfbc
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
   END DO boundary
   cell_type2=9
   WRITE(22,'(a)') "# vtk DataFile Version 3.0"
   WRITE(22,'(a)') "radiation"
   WRITE(22,'(a)') "ASCII"
   WRITE(22,'(a)') "DATASET UNSTRUCTURED_GRID"
   WRITE(22,'(a,i7,a)')"POINTS",nn," float" 
   DO i=1,nn
    WRITE(22,'(3f12.4)')g_coord(:,i)
   END DO
   WRITE(22,'(a,2i10)') "CELLS",hfbc,5*hfbc
   DO i=1,hfbc 
     WRITE(22,'(a,4i6)') "  4",incidencia2(i,1:4)-1 
   END DO
   WRITE(22,'(a,i5)') "CELL_TYPES",hfbc
   DO i=1,hfbc
     WRITE(22,'(i5)') cell_type2
   END DO
   WRITE(22,'(a,i5)')"CELL_DATA",hfbc
   WRITE(22,'(a)')"SCALARS radiation float 1"
   WRITE(22,'(a)')"LOOKUP_TABLE default"
   DO ibound=1,hfbc
     WRITE(22,'(a)')" 1"
   END DO
   CLOSE(22)
 ELSE
   WRITE(11,'(//a)')"There is no radiation boundary conditions"
 END IF
!-------------------------------------------------------------------------
 IF(shadow)THEN
   READ(10,*)shadowin
   IF(shadowin==1)THEN
     CALL form_shadow_table(filename,azimuth,phi,g_num,g_coord,iflux,   &
       inci_s,points_s)
   ELSE
     READ(10,*)ifile
     WRITE(11,'(2a)')"Shadow table ",TRIM(ifile)
     CALL shadow_table(ifile,hfbc,nip_s)
   END IF
 END IF
! corte
!------------------convection boundary condition--------------------------
 airin=0
 idiffuse=zero
 READ(10,*)htbc
 IF(htbc/=0)THEN
   WRITE(11,'(//a)')"Convection boundary conditions"
   READ(10,*)airin
   SELECT CASE(airin)
   CASE(1)
     CALL air_parameters
     nfiles=1
   CASE(2)
     READ(10,*)nfiles
     READ(10,*)ifile
     IF(nfiles==1)THEN
       WRITE(11,'(2a)')"Discrete values from ",TRIM(ifile)
     ELSE
       WRITE(11,'(2a)')"Discrete values from (1) ",TRIM(ifile)
     END IF
     CALL air_table1(ifile,airin)
     IF(nfiles==2)THEN
       READ(10,*)ifile
       WRITE(11,'(2a)')"Discrete values from (2) ",TRIM(ifile)
       CALL air_table2(ifile)
     ELSE IF(nfiles>2)THEN
       WRITE(11,'(2a)')"The program only considers 2 different groups",  &
         " of convection boundary conditions"
       STOP
     END IF
!-------------------------------------------------------------------------
! For variavel heat transfer coefficient it works with implicit scheme
! theta is set to 1 
!-------------------------------------------------------------------------
     IF(theta<1.0_iwp)THEN
       WRITE(11,'(2a)')"The time-integration weighting parameter was ",  &
         "set to 1 --> Implicit scheme" 
       theta=1.0_iwp
     END IF
   CASE(3)
     READ(10,*)ifile
     nfiles=1
     WRITE(11,'(2a)')"Discrete values from ",TRIM(ifile)
     CALL air_table1(ifile,airin) 
   END SELECT
   SELECT CASE(nfiles)
   CASE(1)
     ALLOCATE(itrans(htbc,2))
     WRITE(11,'(/a,/a)')"It is applied to","   Element     iside"
     READ(10,*)(itrans(i,:),i=1,htbc)
     WRITE(11,'(2i10)')(itrans(i,:),i=1,htbc)
!------------------vtk file-----------------------------------------------
     OPEN(22,file='convection.vtk',status='replace',action='write') 
     IF(hfbc/=0) DEALLOCATE (incidencia2)
     ALLOCATE(incidencia2(htbc,4))
     boundary_a: DO ibound=1,htbc
       iel=itrans(ibound,1) 
       iside=itrans(ibound,2) 
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
     END DO boundary_a
     cell_type2=9
     WRITE(22,'(a)') "# vtk DataFile Version 3.0"
     WRITE(22,'(a)') "convection"
     WRITE(22,'(a)') "ASCII"
     WRITE(22,'(a)') "DATASET UNSTRUCTURED_GRID"
     WRITE(22,'(a,i7,a)')"POINTS",nn," float" 
     DO i=1,nn
      WRITE(22,'(3f12.4)')g_coord(:,i)
     END DO
     WRITE(22,'(a,2i10)') "CELLS",htbc,5*htbc
     DO i=1,htbc 
       WRITE(22,'(a,4i6)') "  4",incidencia2(i,1:4)-1 
     END DO
     WRITE(22,'(a,i5)') "CELL_TYPES",htbc
     DO i=1,htbc
       WRITE(22,'(i5)') cell_type2
     END DO
     WRITE(22,'(a,i5)')"CELL_DATA",htbc
     WRITE(22,'(a)')"SCALARS radiation float 1"
     WRITE(22,'(a)')"LOOKUP_TABLE default"
     DO ibound=1,htbc
       WRITE(22,'(a)')" 1"
     END DO
     CLOSE(22)
   CASE(2)
     ALLOCATE(itrans(htbc,3))
     WRITE(11,'(/a,/a)')"It is applied to","   Element     iside    igroup"
     READ(10,*)(itrans(i,:),i=1,htbc)
     WRITE(11,'(3i10)')(itrans(i,:),i=1,htbc)
   END SELECT
 ELSE
   WRITE(11,'(//a)')"There is no convection boundary conditions"
 END IF
!------------------initial temperature------------------------------------
 WRITE(11,'(//a)')"Initial temperatures"
 READ(10,*)indic
 IF (indic==0) THEN
   READ(10,*)val0
   WRITE(11,'(a,E11.4)')"Constant value = ",val0
 ELSE
   ALLOCATE(node_3(indic), value_3(indic))
   WRITE(11,'(/a)')"    Node      Value"
   READ(10,*)(node_3(i),value_3(i),i=1,indic)
   WRITE(11,'(i8,E11.4)')(node_3(i),value_3(i),i=1,indic)
 END IF
!------------------output data--------------------------------------------
 READ(10,*)npri
 WRITE(11,'(//a,i5,a)')"The results are printed every ",npri," time steps"
 READ(10,*)nsensors
 WRITE(11,'(//a,i5)')"The number of required point results is =",nsensors
 IF(nsensors/=0)THEN
   ALLOCATE(ielsen(nsensors),sensors(nsensors,ndim))
   DO j=1,nsensors
     READ(10,*)ifile,sensors(j,:)
     CALL find_points(j,g_num,g_coord,iel,sensors)
     ielsen(j)=iel
     iout=1000+j
     ifile=TRIM(filename)//'_'//TRIM(ifile)//'.res'
     OPEN(iout,FILE=ifile)
   END DO
 END IF
!------------------loop the elements to find global arrays sizes----------
 kdiag=zero
 elements_1: DO iel=1,nels
   num=g_num(:,iel)
   CALL fkdiag(kdiag,num)
 END DO elements_1
 DO i=2,neq
   kdiag(i)=kdiag(i)+kdiag(i-1)
 END DO
 WRITE(11,'(//2(a,i10))')                                                &
 "There are",neq,"  equations and the skyline storage is ",kdiag(neq)
 ALLOCATE(bp(kdiag(neq)),bk(kdiag(neq)),loads(0:neq),newlo(0:neq),       &
   loadsq(0:neq),oldlo(0:neq),loadsh(0:neq))
 bp=zero ; bk=zero ;  loadsh=zero
!------------------load vector initialisation-----------------------------
 loads=zero
 IF(indic==0)THEN
   loads=val0
 ELSE
   DO i=1,indic ; loads(node_3(i))=value_3(i) ; END DO
 END IF
!------------------element heat stiffness and capacitance matrix assembly-
 elements_2: DO iel=1,nels 
   kay=zero ; DO i=1,ndim ; kay(i,i)=prop(i,etype(iel)); END DO 
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num))
   j=etype(iel)
   rho=prop(ndim+1,j) ; cp=prop(ndim+2,j)
   kp=zero ; pm=zero 
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i) ; CALL shape_fun(fun,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac) ; CALL invert(jac) ; deriv=MATMUL(jac,der)
     kp=kp+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
     CALL cross_product(fun,fun,ntn) ; pm=pm+ntn*det*weights(i)*rho*cp
   END DO gauss_pts_1
   CALL fsparv(bk,kp,num,kdiag) ; CALL fsparv(bp,pm,num,kdiag)
 END DO elements_2
!
 IF(airin/=2)THEN
!-------------------------------------------------------------------------
!for constant convective heat transfer coefficient only
!-------------------------------------------------------------------------
!------------------boundary convection matrix and load vector assembly----
   IF(htbc /= 0)THEN 
     boundary_1: DO ibound=1,htbc
       iel=itrans(ibound,1) ; iside=itrans(ibound,2)
       num_s=g_num(inci_s(iside,:),iel)
       coord_s=TRANSPOSE(g_coord(:,num_s))
       elt_s=zero ; kp_s=zero
       htco=prop(ndim+3,etype(iel))
       gauss_pts_2: DO i=1,nip_s
         CALL shape_der(der_s,points_s,i)
         CALL shape_fun(fun_s,points_s,i)
         jac_s=MATMUL(der_s,coord_s) 
         CALL det_s(jac_s,n1,n2,n3,det)
         elt_s=elt_s+fun_s*det*weights_s(i)*htco  
         CALL cross_product(fun_s,fun_s,ntn_s)
         kp_s=kp_s+ntn_s*det*weights_s(i)*htco
       END DO gauss_pts_2
       loadsh(num_s)=loadsh(num_s)+elt_s  
       CALL fsparv(bk,kp_s,num_s,kdiag)
     END DO boundary_1
   END IF
!------------------apply Finite Diference Method--------------------------
   bk=bk*theta*dtim ; bp=bp+bk ; bk=bp-bk/theta
!------------------store augmented diagonal terms-------------------------
   IF(fixed_freedoms_1/=0)THEN
     bp(kdiag(node_1))=bp(kdiag(node_1))+penalty
     storbp_1=bp(kdiag(node_1))
   END IF
   IF(fixed_freedoms_2/=0)THEN
     bp(kdiag(node_2))=bp(kdiag(node_2))+penalty
     storbp_2=bp(kdiag(node_2))
   END IF
!------------------factorise equations------------------------------------
   CALL sparin(bp,kdiag)
 ELSE
!-------------------------------------------------------------------------
!for varaiavel convective heat transfer coefficient the boundary 
!convection matrix and the load vector are computed in each time step
!-------------------------------------------------------------------------
   ALLOCATE(bpcopy(kdiag(neq)),bkcopy(kdiag(neq)))
   bkcopy=bk
   bpcopy=bp
  END IF
!------------------step 0 initialisation----------------------------------
 CALL day_of_year(jd_i,d1) 
!------------------boundary radiation load vector assembly----------------
 loadsq=zero
 IF(hfbc/=0)THEN
   CALL solar_time(d1,jd_i,long,timezone,saving,d2) 
   CALL solar_position(d2,phi,delta,omega,cos_z) 
   SELECT CASE(fluxin)
   CASE(1)
     IF(cos_z<zero) THEN
       ibeam=zero
     ELSE
       ibeam=exponential_function(cos_z) 
     END IF
     idiffuse=zero
   CASE(2)
     CALL Kumar_model(d2,cos_z,ibeam,idiffuse)
   CASE(3)
     ibeam=find_in_table(table3,jd_i)
     IF(cos_z>0.26_iwp)THEN
       ibeam=ibeam/cos_z
     ELSE
       ibeam=zero
     END IF
     idiffuse=find_in_table2(table3,jd_i)
   END SELECT 
   boundary_2: DO ibound=1,hfbc
     iel=iflux(ibound,1) ; iside=iflux(ibound,2) ! ;WRITE(*,*)"ibound",ibound,"iel",iel,"iside",iside
     num_s=g_num(inci_s(iside,:),iel)
     coord_s=TRANSPOSE(g_coord(:,num_s))
     ab=prop(ndim+4,etype(iel))
     elt_s=zero
     gauss_pts_3: DO i=1,nip_s
       CALL shape_der(der_s,points_s,i) ! ; write(*,*)"i",i
       CALL shape_fun(fun_s,points_s,i)
       jac_s=MATMUL(der_s,coord_s)
       CALL det_s(jac_s,n1,n2,n3,det)
       IF(ndim==2)THEN
         IF(n1>zero)THEN
           beta=azimuth+3.0_iwp*pi*0.5_iwp
         ELSE
           beta=azimuth+pi*0.5_iwp
         END IF
         yy=ACOS(n2)
       ELSE
         beta=azimuth+(pi*0.5_iwp-ATAN2(n2,n1))
         yy=ACOS(n3)
       END IF
       cos_alpha=incident_angle(delta,omega,phi,yy,beta)
       IF(shadow.and.cos_alpha>zero)THEN
          ishadow=find_in_table5(table5,d2,ibound,i) 
          cos_alpha=cos_alpha*ishadow
       END IF
       elt_s=elt_s-fun_s*det*weights_s(i)*ab*(ibeam*cos_alpha+idiffuse*  &
         (1.0_iwp+COS(yy))/2.0_iwp) 
     END DO gauss_pts_3
     loadsq(num_s)=loadsq(num_s)+elt_s
   END DO boundary_2 
 END IF
!------------------air temperature for step 0-----------------------------
 SELECT CASE(airin)
 CASE(0)
   t_a=zero
 CASE(1)
   t_a=harmonic_function(d1)
 CASE(2)
!-------------------------------------------------------------------------
! For variavel heat transfer coefficient it works with implicit scheme 
!-------------------------------------------------------------------------
   t_a=zero
   loadsh=zero
   loadsq=zero
 CASE(3)
   t_a=find_in_table(table1,jd_i)  
 END SELECT
 oldlo=loadsh*t_a+loadsq
!------------------time stepping loop-------------------------------------
 timesteps: DO j=1,nstep
   jd=jd_i+j*dtim/ucte
   CALL day_of_year(jd,d1) 
!------------------boundary convection load vector------------------------
   SELECT CASE(airin)
   CASE(0)
     t_a=zero
   CASE(1)
     t_a=harmonic_function(d1)
   CASE(2)  !for varaiavel convective heat transfer coefficient only
     t_a1=find_in_table(table1,jd)
     htco1=find_in_table2(table1,jd)
     IF(nfiles==2)THEN
       t_a2=find_in_table(table2,jd)
       htco2=find_in_table2(table2,jd)
     END IF
     bk=bkcopy
     bp=bpcopy
     loadsh=zero
     boundary_4: DO ibound=1,htbc
       iel=itrans(ibound,1)
       iside=itrans(ibound,2)
       IF(nfiles==2)THEN
         igroup=itrans(ibound,3)
       ELSE
         igroup=1
       END IF
       IF(igroup==1)THEN
         t_a=t_a1
         htco=htco1
       ELSE
         t_a=t_a2
         htco=htco2
       END IF
       num_s=g_num(inci_s(iside,:),iel)
       coord_s=TRANSPOSE(g_coord(:,num_s))
       elt_s=zero ; kp_s=zero
       gauss_pts_5: DO i=1,nip_s
         CALL shape_der(der_s,points_s,i)
         CALL shape_fun(fun_s,points_s,i)
         jac_s=MATMUL(der_s,coord_s)
         CALL det_s(jac_s,n1,n2,n3,det)
         elt_s=elt_s+fun_s*det*weights_s(i)*htco
         CALL cross_product(fun_s,fun_s,ntn_s)
         kp_s=kp_s+ntn_s*det*weights_s(i)*htco
       END DO gauss_pts_5
       loadsh(num_s)=loadsh(num_s)+elt_s*t_a
       CALL fsparv(bk,kp_s,num_s,kdiag)
     END DO boundary_4
     bk=bk*theta*dtim ; bp=bp+bk ; bk=bp-bk/theta 
     IF(fixed_freedoms_1/=0)THEN
       bp(kdiag(node_1))=bp(kdiag(node_1))+penalty
       storbp_1=bp(kdiag(node_1))
     END IF
     IF(fixed_freedoms_2/=0)THEN
       bp(kdiag(node_2))=bp(kdiag(node_2))+penalty
       storbp_2=bp(kdiag(node_2))
     END IF
     CALL sparin(bp,kdiag)
     t_a=1.0_iwp
   CASE(3)  
     t_a=find_in_table(table1,jd)
   END SELECT
!------------------matrix-by-vector multiplication------------------------
   CALL linmul_sky(bk,loads,newlo,kdiag)
!------------------boundary radiation load vector assembly----------------
   loadsq=zero
   IF(hfbc/=0)THEN
     CALL solar_time(d1,jd,long,timezone,saving,d2)
     CALL solar_position(d2,phi,delta,omega,cos_z)
     SELECT CASE(fluxin)
     CASE(1)
       IF(cos_z<=zero) THEN
         ibeam=zero
       ELSE
         ibeam=exponential_function(cos_z) 
       END IF
       idiffuse=zero
     CASE(2)
       CALL Kumar_model(d2,cos_z,ibeam,idiffuse)
     CASE(3)
       ibeam=find_in_table(table3,jd)
       IF(cos_z>0.26_iwp)THEN
         ibeam=ibeam/cos_z
       ELSE
         ibeam=zero
       END IF
       idiffuse=find_in_table2(table3,jd)
     END SELECT 
     boundary_5: DO ibound=1,hfbc
       iel=iflux(ibound,1) ; iside=iflux(ibound,2) 
       num_s=g_num(inci_s(iside,:),iel)
       coord_s=TRANSPOSE(g_coord(:,num_s))
       ab=prop(ndim+4,etype(iel))
       elt_s=zero
       gauss_pts_6: DO i=1,nip_s
         CALL shape_der(der_s,points_s,i)
         CALL shape_fun(fun_s,points_s,i)
         jac_s=MATMUL(der_s,coord_s)
         CALL det_s(jac_s,n1,n2,n3,det)
         IF(ndim==2)THEN
           IF(n1>zero)THEN
             beta=azimuth+3.0_iwp*pi*0.5_iwp
           ELSE
             beta=azimuth+pi*0.5_iwp
           END IF
           yy=ACOS(n2)
         ELSE
           beta=azimuth+(pi*0.5_iwp-ATAN2(n2,n1))
           yy=ACOS(n3)
         END IF
         cos_alpha=incident_angle(delta,omega,phi,yy,beta)
         IF(shadow.and.cos_alpha>zero)THEN
           ishadow=find_in_table5(table5,d2,ibound,i)
           cos_alpha=cos_alpha*ishadow  ;              ! write(102,*)"d1",d1,"d2",d2,"ishadow",ishadow
         END IF
         elt_s=elt_s-fun_s*det*weights_s(i)*ab*(ibeam*cos_alpha+idiffuse &
          *(1+COS(yy))/2.) 
       END DO gauss_pts_6
       loadsq(num_s)=loadsq(num_s)-elt_s
     END DO boundary_5 
   END IF
!------------------load vector--------------------------------------------
   newlo=newlo+dtim*(theta*(loadsh*t_a+loadsq)+(1-theta)*oldlo)
   oldlo=loadsh*t_a+loadsq
!------------------fixed boundary conditions------------------------------
   IF(fixed_freedoms_1/=0)THEN
     DO i=1,fixed_freedoms_1
       newlo(node_1(i))=storbp_1(i)*reservoir_temperature(value_1(i),d1)
     END DO
   END IF
   IF(fixed_freedoms_2/=0)THEN
     IF(tempin==2)THEN
       val1=find_in_table(table4,jd)
       value_2=val1
     END IF
       newlo(node_2)=storbp_2*value_2
   END IF
!------------------equations solution-------------------------------------
   CALL spabac(bp,newlo,kdiag) ; loads=newlo
!------------------write sensors temperatures-----------------------------
   char_num=date_hour(jd)
   IF(nsensors/=0)THEN
     elements_6: DO  i=1,nsensors
       iel=ielsen(i)
       num=g_num(:,iel)
       elt=loads(num)
       CALL shape_fun(fun,sensors,i)
       val2=DOT_PRODUCT(fun,elt)
       iout=1000+i 
       WRITE(iout,'(5a,3x,5a,E15.4,2x,E15.4)')char_num(1:4),'-',char_num(5:6),'-',&  !,2x,3E15.4
         char_num(7:8),char_num(9:10),':',char_num(11:12),':',           &
         char_num(13:14),val2,t_a           ! ,d1,reservoir_temperature(120.00_iwp,d1),t_a
     END DO elements_6
   END IF 
!------------------write node temperaures---------------------------------
   IF(j/npri*npri==j)THEN
     ifile=TRIM(filename)//'_'//char_num//'.res'
     OPEN(50,file=ifile,status='replace',action='write')
     WRITE(50,'(5a,3x,5a)')char_num(1:4),'-',char_num(5:6),'-',          &
       char_num(7:8),char_num(9:10),':',char_num(11:12),':',             &
       char_num(13:14)
     WRITE(50,'(i6,e15.4)')(i,loads(i),i=1,nn)
     CLOSE(50)
   END IF
 END DO timesteps
END PROGRAM PAT_v1_3

