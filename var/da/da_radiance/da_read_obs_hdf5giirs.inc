subroutine da_read_obs_hdf5giirs (iv, infile_tb, infile_clw)

   !--------------------------------------------------------
   !  Purpose:
   !     Read FY-4B GIIRS Level-1 HDF5 data and form
   !     the WRFDA radiance innovation structure.
   !
   !     Current simplified version:
   !       1. Read Only ES_RealLW
   !       2. Read Only WN_LW
   !       3. Convert spectral radiance to brightness temperature
   !          using Planck's law
   !       4. No cloud information is read
   !
   !--------------------------------------------------------

   implicit none

   character(len=*), intent(in)    :: infile_tb, infile_clw
   type(iv_type),    intent(inout) :: iv

#if defined(HDF5)
! fixed parameter values
   integer,parameter::nwn_lw=725       ! Number of WN_LW channels
   integer,parameter::nobs_giirs=128   ! Maximum allowed NumberOfScans
   integer,parameter::time_dims=6      ! Time dimension
   integer,parameter::nsubfile=27      ! low  resolution pixel width
   integer,parameter::nfile_max = 972  ! 27 hdf files contains ~5min of data
                                       ! at most 972 files for a 3-h time window
! interface variable
   integer iret                        ! return status
   integer(HID_T) fhnd1, fhnd2         ! file handle
   integer(HID_T) ahnd1, ahnd2         ! attribute handle
   integer(HID_T) dhnd1, dhnd2         ! dataset handle
   integer(HID_T) shnd1, shnd2         ! dataspace handle
   integer(HSIZE_T) sz1(1)             ! array size 1
   integer(HSIZE_T) sz2(2)             ! array size 2
   integer(HID_T) :: atype
   integer(SIZE_T) :: str_len
   integer(HSIZE_T) :: dims(1)

!  array data
   real(4) :: es_real_lw(nwn_lw,nobs_giirs)
   real(4) :: wn_lw(nwn_lw)
   real(4) :: clw(nwn_lw,nobs_giirs)     
   real(4) :: latitude_lw(nobs_giirs)
   real(4) :: longitude_lw(nobs_giirs)

   real(4) :: sat_zenith_lw(nobs_giirs)
   real(4) :: sat_azimuth_lw(nobs_giirs)
   real(4) :: sun_zenith_lw(nobs_giirs)
   real(4) :: sun_azimuth_lw(nobs_giirs)

   real(r_kind) :: tbb_lw(nwn_lw,nobs_giirs)
   real(r_kind),parameter  :: tbmin  = 50._r_kind
   real(r_kind),parameter  :: tbmax  = 550._r_kind

!  Planck constants
!
!  L = c1 * nu^3 / [ exp(c2*nu/T) - 1 ]
!
!  nu : cm-1
!  L  : mW/(m2 sr cm-1)
!  T  : K

   real(r_kind), parameter :: c1 = 1.191042e-5_r_kind
   real(r_kind), parameter :: c2 = 1.4387752_r_kind

   integer(i_kind)   :: idate5(6)
   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id
   real(r_kind)      :: tb, crit
   real(r_kind)      :: radiance_lw

   real(kind=8)                   :: obs_time
   type(datalink_type), pointer   :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse

   integer           :: i,j,k,l,m,n, ifile, isubfile, landsea_mask
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth
   integer(i_kind)  :: num_giirs_local, num_giirs_global, num_giirs_used, num_giirs_thinned
   integer(i_kind)  :: num_giirs_used_tmp, num_giirs_file
   integer(i_kind)  :: num_giirs_local_local, num_giirs_global_local, num_giirs_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1, filename2
   character(10)    :: begin_date
   character(12)    :: begin_time
   integer          :: nchan,ilongitude,ilatitude,ichannels
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max)
   character(80)    :: fname_clw(nfile_max)
   logical          :: fexist, got_clw_file

! Allocatable arrays
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)

   if (trace_use) call da_trace_entry("da_read_obs_hdf5giirs")

   
!  0.0  Initialize variables
!-----------------------------------
   platform_id  = 52  ! Table-2 Col 1 corresponding to 'fy4'
   satellite_id = 2   ! Table-2 Col 3
   sensor_id    = 98  ! Table-3 Col 2 corresponding to 'giirs'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_giirs_file    = 0
   num_giirs_local   = 0
   num_giirs_global  = 0
   num_giirs_used    = 0
   num_giirs_thinned = 0
   
   sz1(1)=nobs_giirs
   sz2(1) = nwn_lw
   sz2(2) = nobs_giirs

   do i = 1, rtminit_nsensor
      if (platform_id == rtminit_platform(i) .and. &
          satellite_id == rtminit_satid(i) .and. &
          sensor_id == rtminit_sensor(i)) then
         inst = i
         exit
      end if
   end do

   if (inst == 0) then
      call da_warning(__FILE__,__LINE__, &
         (/"The combination of Satellite_Id and Sensor_Id for GIIRS is not found"/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")
      return
   end if

   ! Initialize HDF5 library and Fortran90 interface
   call H5open_f(iret)
   if (iret.lt.0) then
      call da_warning(__FILE__,__LINE__,(/"Problems in Initializing HDF5 library. Can not read GIIRS HDF5 data."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")
      return
   end if

   nchan = iv%instid(inst)%nchan
   write(unit=stdout,fmt=*) 'GIIRS WRFDA nchan = ', nchan
   allocate(data_all(1:nchan))

! 1.0 Assign file names and prepare to read amsr2 files
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   ! first check if L1GIIRS.HDF5 is available
   filename1 = trim(infile_tb)//'.HDF5'
   filename2 = trim(infile_clw)//'.HDF5'
   inquire (file=filename1, exist=fexist)
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
      fname_clw(nfile) = filename2
   else
      ! check if L1GIIRS.HDF5 is available for multiple input files
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max
         write(filename1,fmt='(A,A,I2.2,A)') trim(infile_tb),'-',i,'.HDF5'
         write(filename2,fmt='(A,A,I2.2,A)') trim(infile_clw),'-',i,'.HDF5'
         inquire (file=filename1, exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_tb(nfile)  = filename1
            fname_clw(nfile) = filename2
         write(unit=stdout,fmt=*)'fname_tb=',nfile,filename1
         else
            exit
         end if
      end do
   end if

   if ( nfile == 0 ) then
      call da_warning(__FILE__,__LINE__, &
         (/"No valid GMI L1B.HDF5 or L1GIIRS.HDF5 file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")
      return
   end if


   infile_loop:  do ifile = 1, nfile
      num_giirs_file_local    = 0
      num_giirs_local_local   = 0
      num_giirs_global_local  = 0

   ! open HDF5 file for read
      call H5Fopen_f(fname_tb(ifile),H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"Cannot open GIIRS HDF5 file: "//trim(infile_tb)/))
         cycle infile_loop
      end if
      got_clw_file = .false.
      call H5Fopen_f(fname_clw(ifile),H5F_ACC_RDONLY_F,fhnd2,iret,H5P_DEFAULT_F)
      if ( iret == 0 ) then
         got_clw_file = .true.
      endif

      ! Read Observing Beginning Date
      call h5aopen_f(fhnd1, "Observing Beginning Date", ahnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__, (/"HDF5 attribute open error: Beginning Date"/))
         call h5fclose_f(fhnd1, iret)
         return
      end if
      call h5aget_type_f(ahnd1, atype, iret)
      dims(1) = 1 
      call h5aread_f(ahnd1, atype, begin_date, dims, iret)
      call h5aclose_f(ahnd1, iret)

      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__, (/"HDF5 attribute read error: Beginning Date"/))
         call h5fclose_f(fhnd1, iret)
         return
      end if

      write(unit=stdout,fmt=*) 'begin_date: ', begin_date 
      read(begin_date(1:4), '(I4)') idate5(1)
      read(begin_date(6:7), '(I2)') idate5(2)
      read(begin_date(9:10), '(I2)') idate5(3)

      ! Read Observing Beginning Time
      call h5aopen_f(fhnd1, "Observing Beginning Time", ahnd2, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__, (/"HDF5 attribute open error: Beginning Time"/))
         call h5fclose_f(fhnd1, iret)
         return
      end if
      call h5aget_type_f(ahnd2, atype, iret)
      dims(1) = 1
      call h5aread_f(ahnd2, atype, begin_time, dims, iret)
      call h5aclose_f(ahnd2, iret)

      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__, (/"HDF5 attribute read error: Beginning Time"/))
         call h5fclose_f(fhnd1, iret)
         return
      end if

      write(unit=stdout,fmt=*) 'begin_time: ', begin_time
      read(begin_time(1:2), '(I2)') idate5(4)
      read(begin_time(4:5), '(I2)') idate5(5)
      read(begin_time(7:8), '(I2)') idate5(6)
      
   !  Read LW wavenumber
      call H5Dopen_f(fhnd1, '/Data/WN_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Data/WN_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, wn_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Data/WN_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read LW real radiance
      call H5Dopen_f(fhnd1, '/Data/ES_RealLW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Data/ES_RealLW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, es_real_lw, sz2, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Data/ES_RealLW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read latitude
      call H5Dopen_f(fhnd1, '/Geolocation/Latitude_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Latitude_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, latitude_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Latitude_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read longitude
      call H5Dopen_f(fhnd1, '/Geolocation/Longitude_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Longitude_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, longitude_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Longitude_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read sensor zenith angle
      call H5Dopen_f(fhnd1, '/Geolocation/Sensor_Zenith_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Sensor_Zenith_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, sat_zenith_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Sensor_Zenith_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read sensor azimuth angle
      call H5Dopen_f(fhnd1, '/Geolocation/Sensor_Azimuth_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Sensor_Azimuth_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, sat_azimuth_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Sensor_Azimuth_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read solar zenith angle
      call H5Dopen_f(fhnd1, '/Geolocation/Solar_Zenith_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Solar_Zenith_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, sun_zenith_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Solar_Zenith_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   !  Read solar azimuth angle
      call H5Dopen_f(fhnd1, '/Geolocation/Solar_Azimuth_LW', dhnd1, iret)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 open error for: /Geolocation/Solar_Azimuth_LW"/))
         call H5Fclose_f(fhnd1, iret)
         return
      end if

      call H5Dread_f(dhnd1, H5T_NATIVE_REAL, sun_azimuth_lw, sz1, iret, &
                     H5S_ALL_F, H5S_ALL_F)
      if (iret.lt.0) then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: /Geolocation/Solar_Azimuth_LW"/))
      endif
      call H5Dclose_f(dhnd1, iret)

   ! close infile_tb and HDF5
      call H5Fclose_f(fhnd1,iret)

      if ( got_clw_file ) then
      ! read CLW from infile_clw:
         call H5Dopen_f(fhnd2,'/S1/cloudWaterPath',dhnd2,iret)
         call H5Dread_f(dhnd2,H5T_NATIVE_REAL,clw,sz2,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning(__FILE__,__LINE__, &
               (/"HDF5 read error for: CLW data"/))
         endif
         call H5Dclose_f(dhnd2,iret)
      ! sample display
         write(unit=stdout,fmt=*)'clw(pixel=1,scan=1): ',clw(100,1),clw(2,100),clw(100,2)
    
      ! read SurfaceType from infile_clw:
         call H5Dopen_f(fhnd2,'/S1/surfaceTypeIndex',dhnd2,iret)
         call H5Dread_f(dhnd2,H5T_NATIVE_INTEGER,landsea_mask,sz1,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning(__FILE__,__LINE__, &
               (/"HDF5 read error for: CLW data"/))
         endif
         call H5Dclose_f(dhnd2,iret)

         ! close infile_clw and HDF5
         call H5Fclose_f(fhnd2,iret)         
      end if

   !  Print test values
      write(unit=stdout,fmt=*) 'GIIRS LW test:'
      write(unit=stdout,fmt=*) 'WN_LW(1)       = ', wn_lw(1)
      write(unit=stdout,fmt=*) 'ES_RealLW(1,1) = ', es_real_lw(1,1)


! 2.0 Loop to read hdf file and assign information to a sequential structure
!-------------------------------------------------------------------------

   !  Allocate sequential structure
      if (.not. head_allocated) then
         allocate(head)
         nullify(head%next)
         p => head
         head_allocated = .true.
      end if

      ! start scan_loop
      scan_loop:     do isubfile=1, nsubfile	 

         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)

         write(stdout,*) 'DEBUG idate5 = ', idate5
         write(stdout,*) 'DEBUG obs_time = ', obs_time
         write(stdout,*) 'DEBUG time_slots = ', time_slots(0), time_slots(num_fgat_time)

         if ( obs_time < time_slots(0) .or. obs_time >= time_slots(num_fgat_time) ) cycle scan_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and. obs_time  < time_slots(ifgat) ) exit
         end do

      !  start fov_loop: longitude
         fov_loop:      do j = 1, nobs_giirs

            if ( sat_zenith_lw(j) > 65.0 ) cycle fov_loop
            
            ! Assign geolocation before calling da_llxy.
            info%lat = real(latitude_lw(j), r_kind)
            info%lon = real(longitude_lw(j), r_kind)

            call da_llxy(info, loc, outside, outside_all)
            if (outside_all) cycle fov_loop

            num_giirs_global       = num_giirs_global + 1
            num_giirs_global_local = num_giirs_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE

            num_giirs_local       = num_giirs_local + 1
            num_giirs_local_local = num_giirs_local_local + 1
            write(unit=info%date_char, &
            fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
               idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
               ':', idate5(5), ':', idate5(6)
            info%elv = 0.0

! 3.0  Make Thinning
! Map obs to thinning grid
!-------------------------------------------------------------------
            if (thinning) then
               dlat_earth = info%lat !degree
               dlon_earth = info%lon
               if (dlon_earth<zero)  dlon_earth = dlon_earth+r360
               if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
               dlat_earth = dlat_earth*deg2rad !radian
               dlon_earth = dlon_earth*deg2rad
               crit = 1.
               call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
               if (.not. iuse) then
                  num_giirs_thinned = num_giirs_thinned+1
                  cycle fov_loop
               end if
            end if

            num_giirs_used = num_giirs_used + 1
            data_all = missing_r

         !  Calculate LW brightness temperature
            do i = 1, nwn_lw
               radiance_lw = real(es_real_lw(i,j),r_kind)
               if (radiance_lw > 0.0_r_kind) then
                  tbb_lw(i,j) = c2 * real(wn_lw(i),r_kind) / &
                     log(1.0_r_kind + &
                     c1 * real(wn_lw(i),r_kind)**3 / radiance_lw)
               else
                  tbb_lw(i,j) = missing_r
               end if
            end do
            write(unit=stdout,fmt=*) 'TBB_LW(1,1)    = ', tbb_lw(1,1)

   
            if (i <= nwn_lw) then
               if (tbb_lw(i,j) >= 100.0_r_kind .and. tbb_lw(i,j) <= 450.0_r_kind) then
                     data_all(i) = tbb_lw(i,j)
               endif
            end if
  

! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate(p%tb_inv(nchan))
            p%info = info
            p%loc  = loc      
            p%landsea_mask = 1
            p%scanpos = 1.0_r_kind
            p%satzen = real(sat_zenith_lw(j), r_kind)
            p%satazi = real(sat_azimuth_lw(j), r_kind)
            p%solzen = real(sun_zenith_lw(j), r_kind)
            p%solazi = real(sun_azimuth_lw(j), r_kind)
            p%sensor_index = inst
            p%ifgat = 1
            p%tb_inv(1:nchan) = data_all(1:nchan)

            allocate(p%next)
            p => p%next
            nullify(p%next)

         end do fov_loop
      end do scan_loop

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_file    : ',num_giirs_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_global  : ',num_giirs_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_local   : ',num_giirs_local_local

   end do infile_loop

   call H5close_f(iret)

   deallocate(data_all) ! Deallocate data arrays

   if (thinning .and. num_giirs_global > 0 ) then
#ifdef DM_PARALLEL
   ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         j = j + thinning_grid(inst,ifgat)%itxmax
      end do 

      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            in(j) = thinning_grid(inst,ifgat)%score_crit(i)
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)
	  
      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            if ( ABS(out(j)-thinning_grid(inst,ifgat)%score_crit(i)) > 1.0E-10 ) &
            thinning_grid(inst,ifgat)%ibest_obs(i) = 0
         end do
      end do

      deallocate( in  )
      deallocate( out )

#endif

   ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_giirs_used_tmp = num_giirs_used
      do j = 1, num_giirs_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

      ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_giirs_thinned = num_giirs_thinned + 1
            num_giirs_used = num_giirs_used - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

!  Allocate innovation structure
   iv%total_rad_pixel   = iv%total_rad_pixel + num_giirs_used
   iv%total_rad_channel = iv%total_rad_channel + num_giirs_used*nchan
   
   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_giirs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_giirs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'GIIRS data counts: '
   write(stdout,fmt='(a,i7)') ' In file: ',num_giirs_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_giirs_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_giirs_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_giirs_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_giirs_thinned


!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------
   if (num_giirs_used > 0) then

      iv%instid(inst)%num_rad = num_giirs_used
      iv%instid(inst)%info%nlocal = num_giirs_used

      call da_allocate_rad_iv(inst, nchan, iv)
   end if

!  6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------
   p => head

   do i = 1, num_giirs_used

      call da_initialize_rad_iv(p%sensor_index, i, iv, p)

      current => p
      p => p%next

      deallocate(current%tb_inv)
      deallocate(current)

   end do

   deallocate(p)
   deallocate (ptotal)

   if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")
#else
   call da_error(__FILE__,__LINE__,(/"Needs to be compiled with HDF5 library"/))
#endif

end subroutine da_read_obs_hdf5giirs