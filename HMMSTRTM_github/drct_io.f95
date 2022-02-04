program drct_io
    implicit none
    ! this program steps through all available drct files in a directory

    integer,parameter :: MAXLINE = 5000
    integer,parameter :: MAXSEQLENGTH = 5000
    character(len=4) :: ocode
    character(len=1) :: chainID, ochain, sschar, aa, mrchar, ramachar
    character(len=1) :: rotamer, domainID, context
    integer(2) :: iseq, resseq
    integer(4) ::  i,itr, record
    integer :: nrec, prec, J, K, ios
    integer :: vrecsize = 136
    integer :: trecsize
    integer,dimension(50) :: A, B
    integer(2),dimension(3) :: clusterID
    real(4) :: csum, x, gap_freq, ins_freq
    real(4),dimension(3) :: calpha, dih
    character(len=30) :: drct_fname
    type Profile ! type which stores profile, dih, aa for each !position and seq code
      real(4),dimension(20, MAXSEQLENGTH) :: profiles
      real(4),dimension(3, MAXSEQLENGTH) :: dihs
      character(len=1),dimension(MAXSEQLENGTH) :: seq
      character(len=1),dimension(MAXSEQLENGTH) :: ss_seq
      character(len=1),dimension(MAXSEQLENGTH) :: mr_seq
      character(len=1),dimension(MAXSEQLENGTH) :: rama_seq
      integer,dimension(MAXSEQLENGTH) :: seq_code
      character(len=4) :: code
      character(len=1) :: chain
      integer :: nres
      integer :: last_record
    end type
   integer :: jarg
   type(Profile) :: this_profile
   character(len=4):: en_string, code
   logical :: nan_profile
   integer :: entry_number = 0
   character(len=7) :: format_string
   character(len=4) :: szstr
   integer :: count_match
   character(len=100) :: out_f
   jarg = command_argument_count() ! number of commands from command line
   if ( jarg >= 3 ) then
        call get_command_argument(1, drct_fname )
        call get_command_argument(2, code )
        call get_command_argument(2, out_f )
   else if (jarg >= 2) then
        call get_command_argument(1, drct_fname )
        call get_command_argument(2, code )
        out_f = './dio.txt'
   else if (jarg >= 1) then
        call get_command_argument(1, drct_fname )
        code = 'null'
        out_f = './dio.txt'
   else
        write(*,*) 'Usage: read_drct drctd'
        ! write(*,*) 'Example: master2drct ./dssp/ ./membrane_relations/ ./profiles/ ./drcts/' 
        stop 'read_drct.f95'
   endif
   
   open(3,file=out_f,iostat = ios)
   if ( ios /= 0 ) then
          stop
   endif

   ! database loop script !!!!!!!!!!!!!!!!!!!!!!
   i = 0
   record = 0
   count_match = 0 
   this_profile = get_profile(drct_fname, record)
   record = this_profile%last_record
   do while (this_profile%last_record /=  -1 .and. record /= 0)
           ! write if we've reached desired code, or if code n/a 
           if ( code == 'null' .or. code == this_profile%code) then
           ! write if ! TM PROTEINS
           !if ( any(this_profile%mr_seq(:this_profile%nres) /= '1' )) then 
                   count_match = count_match + 1
                   ! output script !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   write(szstr,'(I4)') size(this_profile%seq)
                   format_string = '(' // szstr // 'A)'
                   write(3,format_string) this_profile%seq(:this_profile%nres)      ! AA SEQ
                   write(3,format_string) this_profile%mr_seq(:this_profile%nres)   ! TM SEQ
                   write(3,format_string) this_profile%ss_seq(:this_profile%nres)   ! SS SEQ
                   write(3,format_string) this_profile%rama_seq(:this_profile%nres) ! RAMA SEQ
                   write(3,'(A6)') this_profile%code // ":" // this_profile%chain   ! CODE:CHAIN
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           endif
           !write(*,*) this_profile%code
           ! cycle script !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           i = i + 1 
           record = this_profile%last_record
           this_profile = get_profile(drct_fname, record)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end do
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,'(a6,i5,a10)') "FOUND ", count_match , " ENTRIES."
   close(3)
   stop
   contains
  
           

type(Profile) function get_profile(database, last_record) result(this_profile)

      implicit none
      ! ----------------- Params
      character(len=*),intent(in) :: database
      integer(4),intent(in) :: last_record
      ! ----------------- Return is this_profile
      ! -----------------
      integer :: ios
      integer :: ritr, pos
      ! ----------------- values contained at each record 
      integer(2) :: iseq, nseq
      character(len=1) :: chainID, sschar, aa, rotamer, ramachar, old_chainID
      real(4) :: calpha(3), dih(3)
      real(4) :: csum, gap_freq, ins_freq
      real(4),dimension(20) :: t_profile
      integer(2) :: clusterID(3)
      character(len=4) :: oldcode
      integer(2) :: hrow
      character(len=1) :: last_pdbtm
      ! -----------------

      this_profile%profiles = 0
      this_profile%dihs = 0
      this_profile%seq = ''
      this_profile%seq_code = 0
      this_profile%code = ''
      this_profile%nres = 0
      this_profile%last_record = -1
      last_pdbtm = ""
      !write(*,*) "OPENING ", database
      open(2, file = trim(database), status = 'OLD', form='unformatted', &
              access = 'direct', recl = vrecsize, iostat = ios)

      if (ios /= 0) stop 'error opening file'

      ritr = last_record + 1
      hrow = 0
      do
!                      write(2, rec = nrec, iostat = reason) iseq, iseq, code, chainID, aa, sschar, rotamer, &
!                        calpha, dih, csum, profile, gap_freq, ins_freq, ramachar, mrchar, clusterID
!  # record size: 136 bytes
!  # (0)seq_list (1)resSeq (2)structure_name (3)chainID (4)seq_resid (5)ss (6)rotamer (7-9)xyzca (10-12)ang (13)csum (14-33)vall_pro (34)gap_freq (35)ins_freq (36)domainID (37)ctx (38-40)clusterID
!  # $template = "s s A4 A A A A f f f f f f f f f f f f f f f f f f f f f f f f f f f f f A A s3";
          t_profile = 0
          read(2, rec=ritr, iostat= ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
                  calpha, dih, csum, t_profile , gap_freq, ins_freq, ramachar, mrchar, clusterID

          !read(2, rec=ritr, iostat= ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
          !        calpha, dih, csum, t_profile , gap_freq, ins_freq, ramachar, mrchar, clusterID
          !write(*,*) ins_freq
          ! GET C PERSISTENCE LENGTH?
          if ( .false. ) then
                  if (last_pdbtm == "C") then
                          if (last_pdbtm == mrchar) then
                                hrow = hrow + 1
                          else
                                write(*,*) hrow
                                hrow = 0
                         endif
                  endif
          endif

          if (ritr - last_record > MAXSEQLENGTH .or. ios /= 0) then ! last seq in file
             !write(*,*) "GOT HERE"
             this_profile%nres = ritr - last_record -1
             this_profile%profiles = this_profile%profiles(:, :this_profile%nres)
             this_profile%dihs = this_profile%dihs(:, :this_profile%nres)
             this_profile%seq = this_profile%seq(:this_profile%nres)
             this_profile%ss_seq = this_profile%ss_seq(:this_profile%nres)    ! Addded 4/7/21 
             this_profile%mr_seq = this_profile%mr_seq(:this_profile%nres) !  TLB
             this_profile%rama_seq = this_profile%rama_seq(:this_profile%nres)
             this_profile%seq_code = get_seq_code(this_profile%seq, this_profile%nres)
             this_profile%chain = old_chainID
             exit
          endif
          if (ritr == last_record + 1) then ! if first time through loop, set this_profile%code = code
              this_profile%code = oldcode
          endif
          if (oldcode /= this_profile%code) then
             this_profile%nres = ritr - last_record - 1
             this_profile%profiles = this_profile%profiles(:,:this_profile%nres)
             this_profile%dihs = this_profile%dihs(:,:this_profile%nres)
             this_profile%seq = this_profile%seq(:this_profile%nres)
             this_profile%ss_seq = this_profile%ss_seq(:this_profile%nres)    ! Addded 4/7/21 
             this_profile%mr_seq = this_profile%mr_seq(:this_profile%nres) !  TLB
             this_profile%rama_seq = this_profile%rama_seq(:this_profile%nres)
             this_profile%seq_code = get_seq_code(this_profile%seq, this_profile%nres)
             this_profile%last_record = ritr
             this_profile%chain = old_chainID
             exit
          endif
          if (oldcode == this_profile%code) then
             pos = ritr - last_record
             this_profile%profiles(:, pos) = t_profile
             this_profile%dihs(:, pos) = dih
             this_profile%seq(pos) = aa
             this_profile%ss_seq(pos) = sschar  ! Addded 4/7/21 
             this_profile%mr_seq(pos) = mrchar !  TLB
             !this_profile%rama_seq(pos) = ramatype(dih(1),dih(2),dih(3)) ! code to handle old format
             this_profile%rama_seq(pos) = ramachar ! new format. 
             this_profile%chain = chainID
             last_pdbtm = mrchar
          endif
          ritr = ritr + 1
          old_chainID = chainID
      enddo

      close (2)
  end function get_profile
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------
  function get_seq_code(seq,nres) result(seq_code)
    implicit none
    character(len=20),parameter :: aa1="ACDEFGHIKLMNPQRSTVWY"
    character,dimension(:) :: seq
    integer,dimension(nres) :: seq_code
    integer :: itr, nres
    do itr=1, nres
        seq_code(itr) = index(aa1,seq(itr))
    enddo
  end function get_seq_code
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------

end program drct_io
