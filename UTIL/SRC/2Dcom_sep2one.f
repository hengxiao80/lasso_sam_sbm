      program main

      implicit none

      character(120) filename
      character(80) long_name
      character(10) units
      character(8)  name
      character(1)  blank/" "/
      character(10) c_min, c_max
      character (12)c_dx, c_dy, c_time
      character(4) rankchar

      integer(2), allocatable :: byte(:)
      real dx,dy,time
      integer nsubs,nsubsx,nsubsy,nx,ny,nz,nfields,nstep
      integer i,k,n,nx_gl,ny_gl,ifields

      ! External functions:

      integer, external :: iargc
      real fldmin, fldmax

      ! Read the file-name from the comman line:

      n=COMMAND_ARGUMENT_COUNT()
      if(n.eq.0) then
        print*,'no input-file name is specified.'
        print*,'Format: 2Dcom_sep2one input.2Dcom (without trailing _*)'
        stop
      end if
      call getarg(1,filename)

      ! Read files; merge data from different subdomains;
      ! save as a netcdf file.

      open(1,file=trim(filename)//"_0",status='old',form='unformatted')
      open(3,file=trim(filename),status='new',form='unformatted')

      read(1,end=3333,err=3333) nstep
      read(1) long_name(1:32)
      write(3) nstep
      write(3) long_name(1:32)
      print*,long_name(1:32)
      read(long_name,'(8i4)') nx,ny,nz,nsubs,nsubsx,nsubsy,nfields
      read(1) c_dx
      read(1) c_dy
      read(1) c_time
      print*, c_time,c_dx,c_dy
      write(3) c_dx
      write(3) c_dy
      write(3) c_time

      allocate(byte(nx*ny))

      do ifields=0,nfields-1
      
        do n=0,nsubs-1

          if(n.ne.0) then 
            if (ifields .eq. 0) then
              write(rankchar,'(i4)') n
              open(8+n,file=trim(filename)//"_"//trim(adjustl(rankchar)),status='old',form='unformatted')
              read(8+n) ! nstep
              read(8+n) ! long_name(1:32)
              read(8+n) ! c_dy
              read(8+n) ! c_dy
              read(8+n) ! c_time
            endif
            read(8+n) (byte(k),k=1,nx*ny)
          else
            read(1) name,blank,long_name,blank,units,blank,c_max,c_min
            write(3) name,blank,long_name,blank,units,blank,c_max,c_min
            print*,ifields+1,trim(long_name),' ',trim(units), '  ',c_min,c_max
            read(1) (byte(k),k=1,nx*ny)
          end if

          write(3) (byte(k),k=1,nx*ny)

        end do ! n

      end do ! ifields

 3333 continue

      end



