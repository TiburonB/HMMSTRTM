module gamma_io
    use drct2profile
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_gamma(xgamma, modelfile, n, nres, code, chain)
          ! write gammas to a drct file with t entries of size n + 6
          ! code = 4
          ! chain = 1
          ! time = 1
          !  nx1 gamma array at time t 
          implicit none
          real(8),dimension(:,:), intent(in) :: xgamma
          character(len=250),intent(in) :: modelfile
          character(len=250) :: gamma_file
          integer :: titr, ios
          integer, intent(in) :: nres
          integer, intent(in)::n
          character(len=11) :: format_string
          character(len=3) :: nnodes_string
          character(len=4),intent(in) :: code
          character(len=1),intent(in) :: chain
          integer :: vrecsize
          character(len=100) :: model_path, model_name
          character(len=3) :: mitr, mitrm1
          integer :: model_iteration, dind, slind, gdl
          logical :: dir_e
          slind = index(modelfile, "/", .true.)
          dind = index(modelfile, ".", .true.)
          call get_model_str(modelfile, model_path, model_name, model_iteration, slind, dind)
          !write(*,*) modelfile
          write(mitrm1,'(I3)') model_iteration
          !write(*,*) mitrm1
          ! make sure directory exists
          if ( model_iteration == 1 ) then
               gamma_file = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitrm1//"GAMMA/"
               gdl = dind + 9
          else
               gamma_file = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitrm1//"GAMMA/"
               gdl = dind + 5
          endif
          inquire(file=trim(gamma_file(:gdl)), exist = dir_e)

          if ( dir_e .eqv. .false. ) then
          !      write(*,*) "no trim = " //gamma_file(:gdl)
                write(*,*) "mkdir '"//trim(gamma_file(:gdl))//"'"
                call system("mkdir '"//trim(gamma_file(:gdl))//"'")
          endif
          gamma_file = gamma_file(:gdl)//this_profile%code
          gamma_file = trim(gamma_file)
          !write(*,*) gamma_file

          !write(*,*) "gamma file = " // trim(gamma_file)
          inquire(file=trim(gamma_file(:gdl+4)), exist = dir_e)

          if ( model_iteration == 1 ) then
               gamma_file = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitrm1//"GAMMA/"//this_profile%code
               gdl = dind + 9
          else
               gamma_file = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitrm1//"GAMMA/"//this_profile%code
               gdl = dind + 5
          endif

          vrecsize = n*8 + 9
          ! write(*,*) ("WRITING GAMMAS to " // trim(gamma_file))
          ! write(*,*) "VRECSIZE = ", vrecsize

           open(3, file=trim(gamma_file),status='replace', form = 'unformatted', &
                   access = 'direct', recl= vrecsize, iostat = ios)
           if (ios /= 0) then
                write(*,*) ios
                write(*,*) "Error 4. Improperly formatted gamma output file."
                stop
           endif
           do titr = 1, nres
!              write(*,*) xgamma(titr, :)
              write(3, rec = titr, iostat = ios) code, chain, titr, xgamma(titr, :)
              if ( ios /= 0) then
                      write(*,*) "ERROR writing to gamma file."
              endif
           enddo
           !write(*,*) "FINISHED WRITING TO GAMMA FILE " // trim(gamma_file)
           close(3)

    end subroutine write_gamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_gamma(xgamma, gamma_file, nres, n)
          implicit none
          real(8),dimension(:,:), intent(inout) :: xgamma
          character(len=*),intent(in) :: gamma_file
          integer :: ttitr, titr, ios
          integer, intent(in) :: nres
          integer, intent(in)::n
          character(len=11) :: format_string
          character(len=3) :: nnodes_string
          character(len=4) :: code
          character(len=1) :: chain
          integer :: vrecsize
          vrecsize = n*8 + 9
          write(*,*) vrecsize
           !write(*,*) ("WRITING GAMMAS to " // gamma_file)
           open(3, file=trim(gamma_file),status='OLD', form = 'unformatted', &
                   access = 'direct', recl= vrecsize, iostat = ios)
           if (ios /= 0) then
                write(*,*) ios
                write(*,*) "Error 4. Improperly formatted gamma output file."
                stop
           endif
           do titr = 1, nres
              read(3, rec = titr, iostat = ios) code, chain, ttitr, xgamma(titr, :)
              if ( ios /= 0) then
                      write(*,*) "ERROR reading from gamma file."
              endif
           enddo

           close(3)

    end subroutine read_gamma

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine loop_gamma(xgamma, gamma_file, n, nres, record, chain_id, retcode)
          implicit none
          real(8),dimension(:,:),allocatable,intent(inout) :: xgamma
          character(len=*),intent(in) :: gamma_file
          integer :: ttitr, titr, ios
          integer, intent(in) :: n
          integer, intent(inout) :: record, nres
          character(len=11) :: format_string
          character(len=3) :: nnodes_string
          character(len=4) :: code, old_code
          character(len=1) :: chain
          character(len=1),intent(inout) :: chain_id
          character(len=4),intent(inout) :: retcode
          integer :: vrecsize, end_record
          real(8),dimension(n) :: temp_gamma
          
          vrecsize = n*8 + 9
          !write(*,*) vrecsize
           !write(*,*) ("READ GAMMAS from " // trim(gamma_file))
           open(3, file=trim(gamma_file),status='OLD', form = 'unformatted', &
                   access = 'direct', recl= vrecsize, iostat = ios)
           if (ios /= 0) then
                write(*,*) ios
                write(*,*) "Error 4. Improperly formatted gamma output file."
                stop
           endif
           titr = record

           do 
              !write(*,*) "SIZEOF READ = " , ( sizeof(code) + sizeof(chain) + sizeof(titr) + sizeof(temp_gamma))
               ! write(3, rec = titr, iostat = ios) code, chain, titr, xgamma(titr, :)
              read(3, rec=titr, iostat = ios) code, chain, ttitr, temp_gamma(:)

              !write(*,*) code, chain, titr, "RECORD =",record
              if ( ios /= 0 ) then
              !        write(*,*) "IOS /=0", ios
                      end_record = titr
                      nres = end_record - record
                      record = -1 ! reset record value, signal eof
                      exit
              endif
              if ( titr == record ) then
                      retcode = code
                      chain_id = chain
                      old_code = code
              endif
              if ( code /= old_code ) then
                      end_record = titr
                      nres = end_record - record
                      record = end_record ! reset record value for next read
                      exit
              endif
              titr = titr + 1
           enddo
           
           !write(*,*) "ALLOCATING XGAMMA, populating array." 
           !write(*,*) "NRES = " , nres
           if (allocated(xgamma)) deallocate(xgamma)
           allocate(xgamma(nres, n))
           !write(*,*) "ALLOCATED XGAMMA, populating array." 
           do titr = record, end_record-1
              !write(*,*) sizeof(xgamma(titr-record+1,:)) + sizeof(ttitr) + sizeof(chain) + sizeof(code)
              read(3, rec = titr, iostat = ios) code, chain, ttitr, xgamma(1+(titr-record), :)
              if (ios /= 0 ) then
                    record = -1 ! end condition
                    close(3)
                    return
              endif
           enddo

           close(3)

    end subroutine loop_gamma



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------
  subroutine get_model_str(model_file, model_path, model_name, model_iteration, slind, dind)
          implicit none
          character(len=*),intent(in) :: model_file
          integer,intent(in) :: slind, dind
          character(len=slind),intent(inout) :: model_path
          character(len=dind-slind-1),intent(inout):: model_name
          integer, intent(inout) :: model_iteration
          integer :: str_len, ind
          integer :: ios

          str_len = len_trim(model_file)
          model_path(:slind) = model_file(:slind)
          model_name(:dind-slind-1) = model_file(slind+1:dind)
          read(model_file(dind-3:dind-1),*,iostat=ios) model_iteration
          if ( ios/= 0) then
            model_iteration = 1
          endif
          !write(*,*) "MODEL_FILE = " , model_file(:len_trim(model_file))
          !write(*,*) "MODEL_PATH = " , model_path(:len_trim(model_path))
          !write(*,*) "MODEL_NAME = " , model_name
          !write(*,*) "MODEL_ITERATION = " , model_iteration
  end subroutine get_model_str

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

subroutine get_gamma_profile(database,last_record, tgamma, n, tcode)
      ! -----------------
      implicit none
      ! ----------------- Params
      character(len=80),intent(in) :: database
      integer,intent(in) :: n
      real(8),dimension(:,:),intent(inout) :: tgamma
      character(len=4),intent(in) :: tcode
      real(8),dimension(n) :: ttgamma
      integer,intent(inout) :: last_record
      ! ----------------- Return is this_profile
      integer :: ios
      integer :: ritr, pos
      character(len=4) :: oldcode
      character(len=4) :: code
      character(len=1) :: chain
      integer :: t, titr
      integer :: vrecsize
      ! -----------------

      vrecsize = n*8 + 9
      tgamma = 0
      oldcode = "NONE"
      open(4, file = trim(database), status = 'OLD', form='unformatted', &
              access = 'direct', recl = vrecsize, iostat = ios)
      if (ios /= 0) stop 'error opening file'
      ritr = last_record + 1
      t = 0
      do
          read(4, rec=ritr, iostat= ios) code, chain, titr, ttgamma(:)
          if ( ios /= 0) then ! last seq in file
             last_record = -1
             close (4)
             exit
          endif
          if ( code /= tcode ) then
               if (oldcode /= "NONE") then
                     last_record = ritr-1
                     close(4)
                     exit
               endif
               ritr = ritr + 1
               cycle
          else
               tgamma(titr,:) = ttgamma
               ritr = ritr + 1
               oldcode = tcode
          endif
      enddo
      close (4)
end subroutine get_gamma_profile


!------------------------------------------------------------------------------------------------------------

end module gamma_io
