! rama2drct.f95
! read in a drct file, write out a drct file with a 'ramachar' added to each entry... new entry size = 137
program rama2drct
    implicit none
    ! this program steps through all available drct files in a directory

    integer,parameter :: MAXLINE = 500
    integer,parameter :: MAXSEQLENGTH = 500
    character(len=4) :: ocode
    character(len=1) :: chainID, ochain, sschar, aa, mrchar
    character(len=1) :: rotamer, domainID, context
    integer(2) :: iseq, resseq, i
    integer :: nrec, prec, J, K, ios
    integer :: vrecsize = 136
    integer,dimension(50) :: A, B
    integer(2),dimension(3) :: clusterID
    real(4) :: csum, x, gap_freq, ins_freq
    real(4),dimension(3) :: calpha, dih
    type Profile ! type which stores profile, dih, aa for each !position and seq code
      real(4),dimension(20, MAXSEQLENGTH) :: profiles
      real(4),dimension(3, MAXSEQLENGTH) :: dihs
      character(len=1),dimension(MAXSEQLENGTH) :: seq
      character(len=1),dimension(MAXSEQLENGTH) :: ss_seq
      character(len=1),dimension(MAXSEQLENGTH) :: mr_seq
      integer,dimension(MAXSEQLENGTH) :: seq_code
      character(len=4) :: code
      integer :: nres
      integer :: last_record
    end type
   integer :: jarg
   type(Profile) :: this_profile
   character(len=100) :: infile, outfile
   
   jarg = command_argument_count() ! number of commands from command line

   call get_command_argument(1, infile )
   call get_command_argument(2, outfile )

   ! (main)
   this_profile = getramaseq(infile, outfile, 0)
   stop
   contains

type(Profile) function getramaseq(database, new_database, last_record) result(this_profile)
      
      implicit none
      ! ----------------- Params
      character(len=*),intent(in) :: database, new_database      
      integer,intent(in) :: last_record
      ! ----------------- Return is this_profile
      ! -----------------
      integer :: ios
      integer :: ritr, pos
      ! ----------------- values contained at each record 
      integer(2) :: iseq, nseq 
      character(len=1) :: chainID, sschar, aa, rotamer
      real(4) :: calpha(3), dih(3)
      real(4) :: csum, gap_freq, ins_freq
      character(len=1) :: context, domainID
      real(4),dimension(20) :: t_profile
      integer(2) :: clusterID(3)
      character(len=4) :: oldcode, code
      character(len=1) :: rama
      character(len=1) :: this_rama
      ! -----------------

      this_profile%profiles = 0
      this_profile%dihs = 0
      this_profile%seq = ''
      this_profile%seq_code = 0 
      this_profile%code = ''
      this_profile%nres = 0
      this_profile%last_record = -1
      code = "    "
      this_rama = " " 
      write(*,*) "OPENING ", database
      open(2, file = trim(database), status = 'OLD', form='unformatted', &
              access = 'direct', recl = vrecsize, iostat = ios)
      if (ios /= 0) stop 'error opening input drct file'
      open(3, file=trim(new_database),status='replace',form='unformatted', access='direct',recl=vrecsize, iostat=ios)
      if (ios /= 0) stop 'error opening output drct file'
      
      ritr = last_record + 1
      do
          read(2, rec=ritr, iostat= ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
                  calpha, dih, csum, t_profile , gap_freq, ins_freq, domainID, mrchar, clusterID
          if (ios /= 0) stop 'error reading input drct file'
          rama = getrama(dih) 
          !write(*,*) rama
          this_rama = rama
          write(3, rec = ritr, iostat = ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
                  calpha, dih, csum, t_profile, gap_freq, ins_freq, this_rama, mrchar, clusterID
          if (ios /= 0) stop 'error writing to drct file'
          !write(*,*) dih
          !if (oldcode /= code .and. code /= "    ") then
          !        stop
          !endif
          code = oldcode
          ritr = ritr + 1
      enddo

      close (2)
      close (3)
  end function getramaseq
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------
  character(len=1) function getrama(dih) result(rama)
    implicit none
    character(len=10) :: rama_string
    real(4),intent(in) :: dih(3)
    real(4),dimension(2,10) :: phipsi
    !character(len=1) :: rama
    integer :: i 
    real :: min_distance 
    integer :: min_index 
    real :: delx, dely, skipx, skipy, this_distance
    
    ! parallel lists... phipsi(i) corresponds to rama_string(i)
    !phipsi = ( (-61.91, -45.20), (-109.78, 20.88), (-70.58, 147.22), (-132.89, 142.43),
    !           (-135.03, 77.26), (-85.03, 72.26), (-165.00, 175.00), (55.88, 38.62),
    !           (85.82, -0.03), (80.00, -170.00))
    
    ! H = (-61.91, -45.20)
    ! G = (-109.78, 20.88)
    ! B = (-70.58, 147.22)
    ! E = (-132.89, 142.43)
    ! d = (-135.03, 77.26)
    ! b = (-85.03, 72.26)
    ! e = (-165.00, 175.00)
    ! L = (55.88, 38.62)
    ! l = (85.82, -0.03)
    ! x = (80.00, -170.00))

    phipsi = reshape((/-61.91, -45.20, -109.78, 20.88, -70.58, 147.22, -132.89, 142.43, -135.03, 77.26, -85.03, 72.26,&
                                 -165.00, 175.00, 55.88, 38.62, 85.82, -0.03, 80.00, -170.00 /), shape(phipsi))
    
    rama_string = "HGBEdbeLlx"
    min_distance = 1000
    min_index = -1

    ! if |dih[2]|< 90 then return 'c'
    if ((dih(1) > 998 .and. dih(2) > 998 .and. dih(3) > 998) .or. &
            ( ABS(dih(1)) < 0.1 .and. ABS(dih(2)) < 0.1 .and. ABS(dih(3)) < 0.1)) then
        rama = '?'
    else if (ABS(dih(1)) > 0 .and. ABS(dih(2)) < 0.1 .and. ABS(dih(3)) < 0.1) then
        ! C terminus... set to undefined
        rama = '?'
    else if (ABS(dih(1)) < 0.1 .and. ABS(dih(2)) > 0 .and. ABS(dih(3)) > 0) then
        ! N terminus... set to undefined
        rama = '?'
    else if ( ABS(dih(3)) < 90) then
        ! cis 
        rama = 'c'
    else
       i = 1
       do while (i <= size(phipsi(2,:)))

            dely = phipsi(2,i) - dih(2)
            if (dely < -180) dely = dely + 360
            if (dely > 180) dely = dely - 360

            delx = phipsi(1,i) - dih(1)
            if (delx < -180) delx = delx + 360
            if (delx > 180) delx = delx - 360

            this_distance = distance(delx,dely) 
            if (this_distance < min_distance) then
                min_distance = this_distance
                min_index = i 
            endif
            i = i + 1
        end do
        if (min_index == -1) then
                write(*,*) "MIN_INDEX ERROR"
                stop
        endif

        rama = rama_string(min_index:min_index)
    endif

  end function getrama
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------
  function distance(delx,dely) result(hypo)
    implicit none
    real,intent(in) :: delx, dely
    real :: hypo
    hypo = (delx ** 2 + dely ** 2) ** 0.5
    !write(*,*) hypo
  end function distance
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------

end program rama2drct
