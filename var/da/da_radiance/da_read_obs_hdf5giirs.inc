subroutine da_read_obs_hdf5giirs (iv, infile_tb)

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

   character(len=*), intent(in)    :: infile_tb
   type(iv_type),    intent(inout) :: iv

#if defined(HDF5)

!-----------------------------------------------------------------------
!  HDF5 variables
!-----------------------------------------------------------------------

   integer :: iret

   integer(HID_T) :: fhnd
   integer(HID_T) :: dhnd

   integer(HSIZE_T) :: sz_lw(2)
   integer(HSIZE_T) :: sz_wn_lw(1)
   integer(HSIZE_T) :: sz_geo(1)

!-----------------------------------------------------------------------
!  GIIRS dimensions
!-----------------------------------------------------------------------

   integer, parameter :: nobs_giirs = 128
   integer, parameter :: nwn_lw     = 725

!-----------------------------------------------------------------------
!  HDF5 data
!-----------------------------------------------------------------------

   real(4) :: es_real_lw(nwn_lw,nobs_giirs)
   real(4) :: wn_lw(nwn_lw)

!-----------------------------------------------------------------------
!  Geolocation and observation angles
!-----------------------------------------------------------------------

   real(4) :: latitude_lw(nobs_giirs)
   real(4) :: longitude_lw(nobs_giirs)

   real(4) :: sensor_zenith_lw(nobs_giirs)
   real(4) :: sensor_azimuth_lw(nobs_giirs)

   real(4) :: solar_zenith_lw(nobs_giirs)
   real(4) :: solar_azimuth_lw(nobs_giirs)

!-----------------------------------------------------------------------
!  Brightness temperature
!-----------------------------------------------------------------------

   real(r_kind) :: tbb_lw(nwn_lw,nobs_giirs)

!-----------------------------------------------------------------------
!  Planck constants
!
!  L = c1 * nu^3 / [ exp(c2*nu/T) - 1 ]
!
!  nu : cm-1
!  L  : mW/(m2 sr cm-1)
!  T  : K
!-----------------------------------------------------------------------

   real(r_kind), parameter :: c1 = 1.191042e-5_r_kind
   real(r_kind), parameter :: c2 = 1.4387752_r_kind

!-----------------------------------------------------------------------
!  WRFDA variables
!-----------------------------------------------------------------------

   integer(i_kind) :: platform_id
   integer(i_kind) :: satellite_id
   integer(i_kind) :: sensor_id
   integer(i_kind) :: inst

   integer :: i,j
   integer :: nchan

   real(r_kind) :: radiance_lw

   logical :: fexist

!-----------------------------------------------------------------------
!  Sequential structure
!-----------------------------------------------------------------------

   type(datalink_type), pointer :: head, p, current
   type(info_type) :: info
   type(model_loc_type) :: loc

   logical :: head_allocated
   logical :: outside, outside_all

!-----------------------------------------------------------------------
!  Statistics
!-----------------------------------------------------------------------

   integer(i_kind) :: num_giirs_local
   integer(i_kind) :: num_giirs_global
   integer(i_kind) :: num_giirs_used

!-----------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_read_obs_hdf5giirs")

   head_allocated = .false.

   num_giirs_local  = 0
   num_giirs_global = 0
   num_giirs_used   = 0

!-----------------------------------------------------------------------
!  Initialize geolocation and observation angles
!-----------------------------------------------------------------------

   latitude_lw       = 0.0
   longitude_lw      = 0.0
   sensor_zenith_lw  = 0.0
   sensor_azimuth_lw = 0.0
   solar_zenith_lw   = 0.0
   solar_azimuth_lw  = 0.0

!-----------------------------------------------------------------------
!  0. GIIRS instrument identification
!-----------------------------------------------------------------------

   platform_id  = 52
   satellite_id = 2
   sensor_id    = 98

   inst = 0

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

!-----------------------------------------------------------------------
!  Number of channels defined by WRFDA
!-----------------------------------------------------------------------

   nchan = iv%instid(inst)%nchan

   write(unit=stdout,fmt=*) 'GIIRS WRFDA nchan = ', nchan

!-----------------------------------------------------------------------
!  Check input file
!-----------------------------------------------------------------------

   inquire(file=trim(infile_tb),exist=fexist)

   if (.not. fexist) then

      call da_warning(__FILE__,__LINE__, &
         (/"GIIRS HDF5 file does not exist: "//trim(infile_tb)/))

      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")

      return

   end if

!-----------------------------------------------------------------------
!  1. Initialize HDF5
!-----------------------------------------------------------------------

   call H5open_f(iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"Problems in initializing HDF5 library."/))

      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")

      return

   end if

!-----------------------------------------------------------------------
!  2. Open HDF5 file
!-----------------------------------------------------------------------

   call H5Fopen_f(trim(infile_tb), H5F_ACC_RDONLY_F, fhnd, iret, &
                  H5P_DEFAULT_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"Cannot open GIIRS HDF5 file: "//trim(infile_tb)/))

      if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")

      return

   end if

!-----------------------------------------------------------------------
!  3. Read LW wavenumber
!-----------------------------------------------------------------------

   sz_wn_lw(1) = nwn_lw

   call H5Dopen_f(fhnd, '/Data/WN_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Data/WN_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, wn_lw, sz_wn_lw, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Data/WN_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  4. Read LW real radiance
!-----------------------------------------------------------------------

   sz_lw(1) = nwn_lw
   sz_lw(2) = nobs_giirs

   call H5Dopen_f(fhnd, '/Data/ES_RealLW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Data/ES_RealLW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, es_real_lw, sz_lw, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Data/ES_RealLW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5. Read geolocation and observation angles
!-----------------------------------------------------------------------

   sz_geo(1) = nobs_giirs

!-----------------------------------------------------------------------
!  5.1 Read latitude
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Latitude_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Latitude_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, latitude_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Latitude_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5.2 Read longitude
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Longitude_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Longitude_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, longitude_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Longitude_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5.3 Read sensor zenith angle
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Sensor_Zenith_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Sensor_Zenith_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, sensor_zenith_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Sensor_Zenith_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5.4 Read sensor azimuth angle
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Sensor_Azimuth_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Sensor_Azimuth_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, sensor_azimuth_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Sensor_Azimuth_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5.5 Read solar zenith angle
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Solar_Zenith_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Solar_Zenith_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, solar_zenith_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Solar_Zenith_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  5.6 Read solar azimuth angle
!-----------------------------------------------------------------------

   call H5Dopen_f(fhnd, '/Geolocation/Solar_Azimuth_LW', dhnd, iret)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 open error for: /Geolocation/Solar_Azimuth_LW"/))

      call H5Fclose_f(fhnd, iret)

      return

   end if

   call H5Dread_f(dhnd, H5T_NATIVE_REAL, solar_azimuth_lw, sz_geo, iret, &
                  H5S_ALL_F, H5S_ALL_F)

   if (iret < 0) then

      call da_warning(__FILE__,__LINE__, &
         (/"HDF5 read error for: /Geolocation/Solar_Azimuth_LW"/))

   endif

   call H5Dclose_f(dhnd, iret)

!-----------------------------------------------------------------------
!  6. Close HDF5 file
!-----------------------------------------------------------------------

   call H5Fclose_f(fhnd, iret)

!-----------------------------------------------------------------------
!  7. Calculate LW brightness temperature
!-----------------------------------------------------------------------

   do j = 1, nobs_giirs

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

   end do

!-----------------------------------------------------------------------
!  8. Print test values
!-----------------------------------------------------------------------

   write(unit=stdout,fmt=*) 'GIIRS LW test:'
   write(unit=stdout,fmt=*) 'WN_LW(1)       = ', wn_lw(1)
   write(unit=stdout,fmt=*) 'ES_RealLW(1,1) = ', es_real_lw(1,1)
   write(unit=stdout,fmt=*) 'TBB_LW(1,1)    = ', tbb_lw(1,1)

!-----------------------------------------------------------------------
!  9. Allocate sequential structure
!-----------------------------------------------------------------------

   if (.not. head_allocated) then

      allocate(head)

      nullify(head%next)

      p => head

      head_allocated = .true.

   end if

!-----------------------------------------------------------------------
!  10. Loop over GIIRS observations
!-----------------------------------------------------------------------

   do j = 1, nobs_giirs

      ! Assign geolocation before calling da_llxy.
      info%lat = real(latitude_lw(j), r_kind)
      info%lon = real(longitude_lw(j), r_kind)

      call da_llxy(info, loc, outside, outside_all)

      if (outside_all) cycle

      num_giirs_global = num_giirs_global + 1

      if (outside) cycle

      num_giirs_local = num_giirs_local + 1

      num_giirs_used = num_giirs_used + 1

      allocate(p%tb_inv(nchan))
      p%info = info
      p%loc  = loc      
      p%landsea_mask = 1
      p%scanpos = 1.0_r_kind
      p%satzen = real(sensor_zenith_lw(j), r_kind)
      p%satazi = real(sensor_azimuth_lw(j), r_kind)
      p%solzen = real(solar_zenith_lw(j), r_kind)
      p%solazi = real(solar_azimuth_lw(j), r_kind)
      p%sensor_index = inst
      p%ifgat = 1
      p%tb_inv(1:nchan) = tbb_lw(1:nchan, j)

      allocate(p%next)
      p => p%next
      nullify(p%next)

   end do

!-----------------------------------------------------------------------
!  11. Allocate innovation structure
!-----------------------------------------------------------------------

   iv%total_rad_pixel   = iv%total_rad_pixel + num_giirs_used
   iv%total_rad_channel = iv%total_rad_channel + num_giirs_used*nchan
   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_giirs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_giirs_global

   if (num_giirs_used > 0) then

      iv%instid(inst)%num_rad = num_giirs_used
      iv%instid(inst)%info%nlocal = num_giirs_used

      call da_allocate_rad_iv(inst, nchan, iv)

   end if

!-----------------------------------------------------------------------
!  12. Transfer sequential structure to innovation structure
!-----------------------------------------------------------------------

   p => head

   do i = 1, num_giirs_used

      call da_initialize_rad_iv(p%sensor_index, i, iv, p)

      current => p
      p => p%next

      deallocate(current%tb_inv)
      deallocate(current)

   end do

   deallocate(p)

!-----------------------------------------------------------------------
!  13. Cleanup
!-----------------------------------------------------------------------

   if (trace_use) call da_trace_exit("da_read_obs_hdf5giirs")

#endif

end subroutine da_read_obs_hdf5giirs