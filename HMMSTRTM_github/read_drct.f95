
program read_drct
    implicit none
    ! this program steps through all available drct files in a directory

    integer,parameter :: MAXLINE = 2000
    integer,parameter :: MAXSEQLENGTH = 2000
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
   real :: ss_H_count 
   real :: ss_B_count 
   integer :: helix_count 
   jarg = command_argument_count() ! number of commands from command line
   if (jarg >= 2) then
        call get_command_argument(1, drct_fname )
        call get_command_argument(2, code )
   else if (jarg >= 1) then
        call get_command_argument(1, drct_fname )
   else
        write(*,*) 'Usage: read_drct drctd'
        stop 'read_drct.f95'
   endif
   open(3,file="t.txt",iostat = ios)
   if ( ios /= 0 ) then
          stop
   endif 
   i = 0
   record = 0
   this_profile = get_profile(drct_fname, record)
   record = this_profile%last_record
   do while (this_profile%last_record /=  -1 .and. record /= 0)
           
           if (this_profile%nres == 0) then
                   record = this_profile%last_record
                   this_profile = get_profile(drct_fname, record)
                   i = i + 1
                   exit
           endif
           nan_profile = .false. 
           helix_count = 0
           write(*,*) this_profile%code, this_profile%chain
           
           ss_H_count = 0
           ss_B_count = 0
           do itr=1,this_profile%nres
              if ( this_profile%ss_seq(itr) == "H" .or. this_profile%ss_seq(itr) == "G") then
                      ss_H_count = ss_H_count + 1.0
              endif
              if ( this_profile%ss_seq(itr) == "E") then 
                      ss_B_count = ss_B_count + 1.0
              endif
           enddo
           ss_H_count = ss_H_count / this_profile%nres
           ss_B_count = ss_B_count / this_profile%nres
           if ( any(this_profile%mr_seq(:this_profile%nres) == 'L' )) then ! TM PROTEINS
           !if ( ss_H_count > 0 .and. ss_B_count == 0 ) then ! ALPHA PROTEINS
           !if ( ss_H_count == 0 .and. ss_B_count > 0 ) then ! BETA PROTEINS
           !if ( ss_H_count > 0 .and. ss_B_count > 0 ) then ! ALPHA/BETA PROTEINS
           !if ( this_profile%code == code ) then       
           ! size(this_profile%seq) == size(this_profile%mr_seq) .and. &
                  ! ((any(this_profile%mr_seq(:this_profile%nres) == ' ')) .eqv. .false.)  ) then
                   !write(*,*) this_profile%seq(:this_profile%nres)
                   !write(*,*) "len(SEQ) = ", size(this_profile%seq)
                   !write(*,*) "SS_SEQ   = ", this_profile%ss_seq(:this_profile%nres)
                   !write(*,*) "len(SS_SEQ) = ", size(this_profile%ss_seq)
                   !write(*,*) 
                   write(*,'(a5)') this_profile%code  // this_profile%chain
                   helix_count = 0
                   if ( .false. ) then
                           do itr=1, this_profile%nres
                               if (this_profile%mr_seq(itr) == "H" .and. this_profile%mr_seq(itr-1) /= "H") then
                                    helix_count = helix_count + 1
                                    write(*,*) "HELIX ", helix_count
                               endif     
                               if ( this_profile%mr_seq(itr) == "H") then
                                    write(*,'(20F8.5)') this_profile%profiles(:,itr)
                               endif
                           enddo
                   endif
                   write(szstr,'(I4)') size(this_profile%seq)
                   format_string = '(' // szstr // 'A)'
                   !write(*,format_string,iostat=ios) this_profile%seq(:this_profile%nres)
                   write(*,format_string,iostat=ios) this_profile%mr_seq(:this_profile%nres)
                   !write(*,format_string,iostat=ios) this_profile%ss_seq(:this_profile%nres)
                   !write(*,format_string,iostat=ios) this_profile%rama_seq(:this_profile%nres)
                   if ( ios /= 0 ) then
                           write(*,*) this_profile%code
                   endif
               !    write(*,format_string,iostat=ios) this_profile%mr_seq(:this_profile%nres)
                   if ( ios /= 0 ) then
                           write(*,*) this_profile%code
                   endif
                   !write(*,*) "RAMA_SEQ = ", this_profile%rama_seq(:this_profile%nres)
                   !write(*,*) "len(MR_SEQ) = ", size(this_profile%mr_seq(:this_profile%nres))
                   !!!write(*,*) "SEQ_CODE = ", this_profile%seq_code(:this_profile%nres)
                   !write(*,*) "len(SEQ_CODE) = ", size(this_profile%seq_code)
                   !write(*,*) "CHAIN     = ", this_profile%chain
                   !write(*,*) " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
                   do itr=1,this_profile%nres
                      do j =1,20
                         if (this_profile%profiles(j,itr) /= this_profile%profiles(j,itr)) then
                                 nan_profile = .true.
                         endif
                      enddo 
                      if (nan_profile) then 
                          !write(*,*) this_profile%profiles(:,itr)
               !           write(*,'(a4)') this_profile%code
                      endif 
                  enddo
           endif
           if (nan_profile) then
               !write(*,"(A5)") this_profile%code//this_profile%chain
           endif
                  !write(*,*) "Dih = ", this_profile%dihs(:,:this_profile%nres)
           i = i + 1 
           record = this_profile%last_record
           this_profile = get_profile(drct_fname, record)
  end do
    write(*,*) i  
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
             this_profile%last_record = ritr 
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
  subroutine set_flags(ramaseq,nres,use_str,angles,iflag)
   !----------------------------------------
   implicit none
   !----------------------------------------
   character(len=1),dimension(:),intent(out)::ramaseq
   integer,intent(in)::nres
   integer,dimension(:),intent(out)::iflag
   logical,intent(in)::use_str
   real,dimension(:,:),intent(in)::angles
   !----------------------------------------
   integer::i
   !----------------------------------------
   ramaseq="?"
   iflag=1
   if(use_str) then
     do i=1,nres
        if(all(angles(:,i)/=999.) .and. all(angles(:,i)/=0.00000000)) then
           iflag(i)=3
        endif
        if(iflag(i)==3) then
           ramaseq(i)=ramatype(angles(1,i),angles(2,i),angles(3,i))
        endif
     enddo
   endif
  end subroutine set_flags
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


end program read_drct
