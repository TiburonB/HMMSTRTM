module profile_io
      implicit none
      ! methods list 
      ! ------------------------------------------------------------
      public :: get_profile
      ! ------------------------------------------------------------
      public :: get_seq_code
      private :: ramatype! private data  ! TLB 5/29/2021 Make this public. HMMSTR can access Profile directly
        integer, parameter :: MAXSEQLENGTH = 5000 ! current largest 
        integer, parameter :: vrecsize = 136 ! 136 bytes... machine dependent, may also be 34 bytes
      ! ------------------------------------------------------------
      public
      type Profile ! type which stores profile, dih, aa for each !position and seq code
        real(4),dimension(20, MAXSEQLENGTH) :: profiles
        real(4),dimension(3, MAXSEQLENGTH) :: dihs
        character(len=1),dimension(MAXSEQLENGTH) :: seq
        character(len=1),dimension(MAXSEQLENGTH) :: tm_seq
        character(len=1),dimension(MAXSEQLENGTH) :: ss_seq
        character(len=1),dimension(MAXSEQLENGTH) :: rama_seq
        integer,dimension(:),allocatable :: seq_code
        character(len=4) :: code
        character(len=1) :: chain
        integer :: nres
        integer :: last_record
      end type
      ! ------------------------------------------------------------ 
      !-------------------------------------------------------------
        type(Profile),public :: this_profile
      !-------------------------------------------------------------

contains

! read from a formatted amino acid profile file
type(Profile) function read_profile(profile_file) result(this_profile)
     implicit none
     character(len=100) :: profile_file
     character(len=300) :: aline
     character(len=30),dimension(30) :: args
     integer :: titr, bitr, ios
     integer :: nres

     open(1, file=profile_file,iostat=ios)
     if (ios/=0) then
        write(*,*) "Error 2. Improperly formatted profile. Error on open."
        stop
     endif
     titr = 1
     do 
        read(1,'(a)') aline
        args = get_substrings(aline)
        select case (trim(args(1)))
        case ('nres')
             read(args(2),*) nres
             this_profile%nres = nres
        case ('code')
             this_profile%code = trim(args(2))
        case ('chain')  
             this_profile%chain = trim(args(2))
        case ('eof')
             allocate(this_profile%seq_code(nres))
             this_profile%seq(:nres) = profile2seq(this_profile%profiles, this_profile%nres)
             this_profile%seq_code(:nres) = get_seq_code(this_profile%seq, nres)
             exit
        case default
            do bitr = 1,20
                read(args(bitr), *) this_profile%profiles(bitr,titr) 
            enddo
            titr = titr + 1
        end select
     enddo
     close(1)
end function read_profile
   
function profile2seq(profile, nres) result(seq)
        real(4),dimension(:,:),intent(in) :: profile
        integer,intent(in) :: nres
        character(len=1),dimension(:),allocatable :: seq
        integer :: titr, i 
        character(len=20) :: alpha

        alpha = "ACDEFGHIKLMNPQRSTVWY"
        allocate(seq(nres))
        do titr = 1,nres
           i =  maxloc(profile(:,titr),1)
           seq(titr) = alpha(i:i)
        enddo
        
end function profile2seq

! write the current profile into a formatted, keyworded file
subroutine write_profile(this_profile, profile_file)

     implicit none
     integer :: titr, ios
     character(len=*) :: profile_file
     type(Profile) :: this_profile
     open(1, file=profile_file,iostat=ios)
     if (ios/=0) then
        write(*,*) "Error 2. Improperly formatted profile. Error on open."
        stop
     endif

     write(1,'(a5,i4)') 'nres ', this_profile%nres
     write(1,'(a5,a4)') 'code ', this_profile%code
     write(1,'(a6,a1)') 'chain ', this_profile%chain
     do titr = 1, this_profile%nres
        write(1, '(20F8.5)') this_profile%profiles(:,titr)
     enddo
     write(1,'(a3)') 'eof'

     close(1)
end subroutine write_profile

! given the name of a binary file containing profile info and the position of the last observed record,
! return the next sequence's profile
! passing 0 as last_record will reutrn 1st sequence's profile.
! if the profile was the last sequence in the database, then last_record == -1 in return type
type(Profile) function get_profile(database, last_record) result(this_profile)
      
      ! -----------------
      implicit none
      ! ----------------- Params
      character(len=100) :: database      
      integer :: last_record
      ! ----------------- Return is this_profile
      ! -----------------
      integer :: ios,i,j 
      integer :: ritr, pos
      ! ----------------- values contained at each record 
      integer(2) :: iseq, nseq 
      character(len=1) :: chainID, sschar, aa, rotamer
      real(4) :: calpha(3), dih(3)
      real(4) :: csum, gap_freq, ins_freq
      character(len=1) :: ramachar, tmchar
      real(4),dimension(20) :: t_profile
      integer(2) :: clusterID(3)
      character(len=4) :: oldcode
      logical :: first
      ! -----------------

      this_profile%profiles = 0
      this_profile%dihs = 0
      this_profile%seq = ''
      this_profile%rama_seq = ''
      this_profile%tm_seq = ''
      this_profile%ss_seq = '' 
      !this_profile%seq_code = 0 
      this_profile%code = ''
      this_profile%nres = 0
      this_profile%last_record = -1

      open(2, file = trim(database), status = 'OLD', form='unformatted', &
              access = 'direct', recl = vrecsize, iostat = ios)

      if (ios /= 0) stop 'error opening file'
      
      ritr = last_record + 1
      first = .true. 
      do
!          read(2, rec=ritr, iostat= ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
!               calpha, dih, csum, t_profile , gap_freq, ins_freq, ramachar, mrchar, clusterID
          t_profile = 0
          read(2, rec=ritr, iostat= ios) iseq, nseq, oldcode, chainID, aa, sschar, rotamer, &
                  calpha, dih, csum, t_profile , gap_freq, ins_freq, ramachar, tmchar, clusterID
          if ( ios/=0 .and. first) then
                  this_profile%nres = -1
                  exit
          else if (ritr - last_record > MAXSEQLENGTH .or. ios/= 0  ) then ! last seq in file
             !write(*,*) "GOT HERE"
             this_profile%nres =  ritr - last_record
             this_profile%profiles(:,this_profile%nres+1:) = 0
             !this_profile%profiles = this_profile%profiles(:20, :this_profile%nres)
             !this_profile%dihs = this_profile%dihs(:, :this_profile%nres)
             !this_profile%seq = this_profile%seq(:this_profile%nres)
             !this_profile%ss_seq = this_profile%ss_seq(:this_profile%nres)
             !this_profile%tm_seq = this_profile%tm_seq(:this_profile%nres)
             !this_profile%rama_seq = this_profile%rama_seq(:this_profile%nres)
             allocate(this_profile%seq_code(this_profile%nres))
             this_profile%seq_code(:) = get_seq_code(this_profile%seq, this_profile%nres)
             !call write_profile()
      !       write(*,*) "GOT PROFILE"
             exit
          endif
          if (ritr == last_record + 1) then ! if first time through loop, set this_profile%code = code
              this_profile%code = oldcode
          endif
          if (oldcode /= this_profile%code) then
             this_profile%nres = ritr - last_record
             !this_profile%profiles = this_profile%profiles(:20,:this_profile%nres)
             !this_profile%dihs = this_profile%dihs(:,:this_profile%nres)
             !this_profile%seq = this_profile%seq(:this_profile%nres)
             !this_profile%ss_seq = this_profile%ss_seq(:this_profile%nres)
             !this_profile%tm_seq = this_profile%tm_seq(:this_profile%nres)
             !this_profile%rama_seq = this_profile%rama_seq(:this_profile%nres)
             allocate(this_profile%seq_code(this_profile%nres))
             this_profile%seq_code(:) = get_seq_code(this_profile%seq, this_profile%nres)
             this_profile%last_record = ritr - 1
             !call print_profile()
       !      write(*,*) "GOT PROFILE"
             exit
          endif
          if (oldcode == this_profile%code) then
             pos = ritr - last_record
             this_profile%chain = chainID
        !     write(*,*) pos
             this_profile%profiles(:, pos) = t_profile(:)
             this_profile%dihs(:, pos) = dih(:)
             this_profile%seq(pos) = aa
             this_profile%tm_seq(pos) = tmchar
             if ( database == "vall16L.drct" ) then 
                     this_profile%rama_seq(pos) = ramatype(dih(1),dih(2),dih(3)) ! code to handle old format
             else
                 this_profile%rama_seq(pos) = ramachar ! new format. 
             endif
             !this_profile%rama_seq(pos) = ramachar
             this_profile%ss_seq(pos) = sschar
          endif
          ritr = ritr + 1
          first = .false.
      enddo

      close (2)
end function get_profile
! ---------------------------------------------------------------
subroutine print_profile()
    implicit none
    integer :: i 
    !write(*,*) this_profile%nres, this_profile%code
    !write(*,*) this_profile%seq(:this_profile%nres)
    !write(*,*) this_profile%ss_seq(:this_profile%nres)
    !write(*,*) this_profile%tm_seq(:this_profile%nres)
    !write(*,*) this_profile%rama_seq(:this_profile%nres)
    do i=1,this_profile%nres
        write(*,*) this_profile%profiles(19,i)
        !write(*,*) this_profile%dihs(:,i)
    enddo

end subroutine print_profile
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
    !write(*,*) "GOT SEQ_CODE"
end function get_seq_code
! ---------------------------------------------------------------
character(len=1) function ramatype(ph,ps,om)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   real,intent(in)::ph
   real,intent(in)::ps
   real,intent(in)::om
   !--------------------------------------------
   integer::i
   integer::icen
   real::d
   real::dmin
   real::ds
   real::df
   real,dimension(11,2)::ramacen
   character(len=11)::rama_string="HGBEdbeLlxc"

   !--------------------------------------------
   ramacen(1,1)=-61.91   ;ramacen(1,2)=-45.20
   ramacen(2,1)=-109.78  ;ramacen(2,2)=20.88
   ramacen(3,1)=-70.58   ;ramacen(3,2)=147.22
   ramacen(4,1)=-132.89  ;ramacen(4,2)=142.43
   ramacen(5,1)=-135.03  ;ramacen(5,2)=77.26
   ramacen(6,1)=-85.03   ;ramacen(6,2)=72.26
   ramacen(7,1)=-165.00  ;ramacen(7,2)=175.00
   ramacen(8,1)=55.88    ;ramacen(8,2)=38.62
   ramacen(9,1)=85.82    ;ramacen(9,2)=-.03
   ramacen(10,1)=80.     ;ramacen(10,2)=-170.00
   ramacen(11,1)=-70.0   ;ramacen(11,2)=150.00
   if(ph>=998.)then; ramatype="?"; return; endif
   if(ps>=998.)then; ramatype="?"; return; endif
   if(om<90.and.om>-90)then; ramatype="c"; return; endif
   dmin=999
   do i=1,11
     df=abs(ph-ramacen(i,1)); if(df>180.) df= 360.-df
     ds=abs(ps-ramacen(i,2)); if(ds>180.) ds= 360.-ds
     d=sqrt(ds*ds+df*df)
     if(d<dmin) then
       icen=i; dmin=d
     endif
   enddo
   ! rama_string="HGBEdbeLlxc"
   ramatype=rama_string(icen:icen)
  end function ramatype
function get_substrings(aline) result(args)

    character(len=300),intent(in) :: aline
    character(len=30),dimension(30) :: args
    integer :: b1,b2
    integer :: itr, ind
    character(len=30) :: t_arg
    b1 = 1
    b2 = 1
    ind = 1
    do itr = 1, len(aline)
        b2 = index(aline(b1:), ' ') + b1
        if (b1 == b2 - 1) then
           b1 = b2
           cycle
        endif
        t_arg = aline(b1:b2-1)
        if (ind > 30) exit
        args(ind) = t_arg
        ind = ind + 1
        b1 = b2
    enddo
endfunction get_substrings
! ---------------------------------------------------------------
end module profile_io
