subroutine da_qc_giirs (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for GIIRS HDF5 data.
   !
   ! Simplified HDF5 GIIRS QC:
   !   1. Mixed surface check
   !   2. Land/sea check
   !   3. iuse check
   !   4. Absolute innovation check
   !   5. Relative innovation check
   !
   ! Cloud-related QC is NOT performed because the current
   ! FY-4B GIIRS HDF5 product does not provide the required
   ! cloud information.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it
   integer, intent(in)             :: i
   integer, intent(in)             :: nchan
   type (y_type),  intent(in)      :: ob
   type (iv_type), intent(inout)   :: iv

   ! local variables
   logical :: lmix
   integer :: n, k, isflg
   integer :: fgat_rad_unit, ios

   integer :: ngood(nchan)
   integer :: nrej(nchan)
   integer :: nrej_omb_abs(nchan)
   integer :: nrej_omb_std(nchan)

   integer :: num_proc_domain
   integer :: nrej_mixsurface
   integer :: nrej_land

   real :: inv_grosscheck

   character(len=30) :: filename

   if (trace_use) call da_trace_entry("da_qc_giirs")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0

   nrej_mixsurface = 0
   nrej_land       = 0
   num_proc_domain = 0

   write(*,*) "GIIRS HDF5 Radiance Data quality control"

   do n = iv%instid(i)%info%n1, iv%instid(i)%info%n2

      if (iv%instid(i)%info%proc_domain(1,n)) &
            num_proc_domain = num_proc_domain + 1

      !---------------------------------------------------------------
      ! 0. Initialize QC flags assuming all observations are good
      !---------------------------------------------------------------

      iv%instid(i)%tb_qc(:,n) = qc_good

      !---------------------------------------------------------------
      ! 1. Reject all channels over mixed surface type
      !---------------------------------------------------------------

      isflg = iv%instid(i)%isflg(n)

      lmix = (isflg == 4) .or. &
             (isflg == 5) .or. &
             (isflg == 6) .or. &
             (isflg == 7)

      if (lmix) then

         iv%instid(i)%tb_qc(:,n) = qc_bad

         if (iv%instid(i)%info%proc_domain(1,n)) &
            nrej_mixsurface = nrej_mixsurface + 1

      end if

      !---------------------------------------------------------------
      ! 2. Land/sea check
      !
      !    Keep channels 2, 3 and 4 over land.
      !    Other channels are rejected when only_sea_rad is enabled.
      !---------------------------------------------------------------

      if (isflg > 0) then

         do k = 1, nchan

            if (k /= 2 .and. k /= 3 .and. k /= 4) then

               if (only_sea_rad) then

                  iv%instid(i)%tb_qc(k,n) = qc_bad

                  if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_land = nrej_land + 1

               end if

            end if

         end do

      end if

      !---------------------------------------------------------------
      ! 3. Check iuse
      !---------------------------------------------------------------

      do k = 1, nchan

         if (satinfo(i)%iuse(k) == -1) then
            iv%instid(i)%tb_qc(k,n) = qc_bad
         end if

      end do

      !---------------------------------------------------------------
      ! 4. Absolute innovation check
      !
      !    Reject if |O-B| > 15 K
      !---------------------------------------------------------------

      do k = 1, nchan

         inv_grosscheck = 15.0

         if (use_satcv(2)) inv_grosscheck = 100.0

         if (abs(iv%instid(i)%tb_inv(k,n)) > inv_grosscheck) then

            iv%instid(i)%tb_qc(k,n) = qc_bad

            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_omb_abs(k) = nrej_omb_abs(k) + 1

         end if

      end do

      !---------------------------------------------------------------
      ! 5. Relative innovation check
      !
      !    Reject if |O-B| > 3 * observation error
      !---------------------------------------------------------------

      do k = 1, nchan

         if (abs(iv%instid(i)%tb_inv(k,n)) > &
             3.0 * iv%instid(i)%tb_error(k,n)) then

            iv%instid(i)%tb_qc(k,n) = qc_bad

            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_omb_std(k) = nrej_omb_std(k) + 1

         end if

         !------------------------------------------------------------
         ! Final QC decision
         !------------------------------------------------------------

         if (iv%instid(i)%tb_qc(k,n) == qc_bad) then

            iv%instid(i)%tb_error(k,n) = 500.0

            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej(k) = nrej(k) + 1

         else

            if (iv%instid(i)%info%proc_domain(1,n)) &
               ngood(k) = ngood(k) + 1

         end if

      end do

   end do

   !-----------------------------------------------------------------
   ! MPI statistics
   !-----------------------------------------------------------------

   call da_proc_sum_int  (num_proc_domain)
   call da_proc_sum_int  (nrej_land)
   call da_proc_sum_int  (nrej_mixsurface)

   call da_proc_sum_ints (nrej_omb_abs)
   call da_proc_sum_ints (nrej_omb_std)
   call da_proc_sum_ints (nrej)
   call da_proc_sum_ints (ngood)

   !-----------------------------------------------------------------
   ! Write QC statistics
   !-----------------------------------------------------------------

   if (rootproc) then

      if (num_fgat_time > 1) then

         write(filename,'(i2.2,a,i2.2)') &
            it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_', &
            iv%time

      else

         write(filename,'(i2.2,a)') &
            it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)

      end if

      call da_get_unit(fgat_rad_unit)

      open(fgat_rad_unit, file=trim(filename), &
           form='formatted', iostat=ios)

      if (ios /= 0) then

         write(unit=message(1),fmt='(A,A)') &
            'error opening the output file ', filename

         call da_error(__FILE__,__LINE__,message(1:1))

      end if

      write(fgat_rad_unit,fmt='(/a/)') &
         ' Quality Control Statistics for '// &
         iv%instid(i)%rttovid_string

      if (num_proc_domain > 0) &
         write(fgat_rad_unit,'(a20,i7)') &
         ' num_proc_domain  = ', num_proc_domain

      if (nrej_land > 0) &
         write(fgat_rad_unit,'(a20,i7)') &
         ' nrej_land        = ', nrej_land

      if (nrej_mixsurface > 0) &
         write(fgat_rad_unit,'(a20,i7)') &
         ' nrej_mixsurface  = ', nrej_mixsurface

      write(fgat_rad_unit,'(a20)') ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)') nrej_omb_abs(:)

      write(fgat_rad_unit,'(a20)') ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)') nrej_omb_std(:)

      write(fgat_rad_unit,'(a20)') ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)') nrej(:)

      write(fgat_rad_unit,'(a20)') ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)') ngood(:)

      close(fgat_rad_unit)

      call da_free_unit(fgat_rad_unit)

   end if

   if (trace_use) call da_trace_exit("da_qc_giirs")

end subroutine da_qc_giirs