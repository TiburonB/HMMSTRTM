module hmm_io ! handles input / output for .hmm file format. 
     implicit none
     ! methods list
     ! --------------------------------------------------
     public :: write_hmm
     ! --------------------------------------------------
     public :: read_hmm
     ! --------------------------------------------------
     ! --------------------------------------------------
     public 
     type HMM_
        character(len=250) :: model_file
        real,allocatable :: priors (:)
        character(1),dimension(:),allocatable :: AA_seq, SS_seq, TM_seq, RAMA_seq
        real::AA_EPS  ! learning rate parameter. 
        real::SS_EPS  ! learning rate parameter. 
        real::TM_EPS  ! learning rate parameter. 
        real::RAMA_EPS! learning rate parameter. 
        real,dimension(:,:),allocatable :: AA
        real,dimension(:,:),allocatable :: SS
        real,dimension(:,:),allocatable :: TM
        real,dimension(:,:),allocatable :: RAMA
        real,dimension(:,:),allocatable :: dihs
        integer :: n,k,u,nta,ntb,ntab
        character(len=4),dimension(20) :: emissions
        character(len=30),dimension(20) :: emission_alphabet
        real,dimension(:),allocatable :: AA_ground, SS_ground, TM_ground, RAMA_ground !
        real,allocatable :: a(:,:) ! TRANSITIONS
        real,dimension(:),allocatable :: RAMA_reweight, SS_reweight, TM_reweight
        integer,dimension(:,:),allocatable :: aties,bties,abties ! array is nt x 2
        integer,dimension(:,:),allocatable :: btie_lookup ! array is n x 400
        integer,dimension(:),allocatable :: atie_lookup,abtie_lookup ! array is n
        !real,dimension(20) :: bground 
        !real,dimension(6)  :: ssground 
        !real,dimension(10) :: tmground 
        !real,dimension(11) :: ramaground 
    endtype
     private :: get_substrings
!    character(len=250) :: hmmfile
!    type(HMM) :: HMMSTR
!    ! MAIN 
!    iarg = command_argument_count()
!    if (iarg < 1) then
!        write(*,*) 'ERROR 1. Usage: read_hmm hmm_file.hmm'
!        stop
!    else if ( iarg == 1) then 
!         ! LOAD FILE, CREATE COPY @ ./HMM_OUT.hmm for testing purposes.
!        call get_command_argument(1,hmmfile)
!
!        ! READ IN FILE TO HMMSTR
!        call read_hmm('./HMMS/'//hmmfile, HMMSTR)
!        ! WRITE OUT HMMSTR TO FILE
!        call write_hmm('./HMMS/'//hmmfile,'./HMM_WAREHOUSE/HMM_OUT.hmm', HMMSTR)
!
!    end if

contains 

integer function write_hmm(hmmfile_in, hmmfile_out, HMMSTR) result(ios)
   implicit none
   character(len=*),intent(in) :: hmmfile_in, hmmfile_out
   type(HMM_),intent(inout) :: HMMSTR
   character(len=300) :: aline
   character(len=30),dimension(30) :: args
   integer :: itr, nitr, nitr_j
   ! write all comment lines from _in to _out
   write(*,*) ("WRITING HMM to " // hmmfile_out)
   open(1, file=hmmfile_in, iostat = ios)
   if (ios /= 0) then
        write(*,*) "Error 2. Improperly formatted input file."
        stop
   endif
   open(2, file=hmmfile_out, iostat = ios)
   if (ios /= 0) then
        write(*,*) "Error 2. Improperly formatted output file."
        stop
   endif
   
   do ! copy the comment lines from input file to output file. 
      read(1, '(a)', iostat = ios) aline
      if (ios/=0) then
          write(*,*) "Error 2. Improperly formatted hmm."
          stop
      endif
      args = get_substrings(aline)
      select case(trim(args(1)))
      case('COMMENT', 'comment')
          write(2,*) adjustl(trim(aline))
      case default
          exit
      end select 
   enddo    
   close(1)

   write(2, "(A7,F20.15)") "AA_EPS ", HMMSTR%AA_EPS
   write(2, "(A7,F20.15)") "SS_EPS ", HMMSTR%SS_EPS
   write(2, "(A7,F20.15)") "TM_EPS ", HMMSTR%TM_EPS
   write(2, "(A9,F20.15)") "RAMA_EPS ", HMMSTR%RAMA_EPS
   ! n_a_ties
   write(2,"(A9,I7)") "n_a_ties ", HMMSTR%nta
   ! n_b_ties
   write(2,"(A9,I7)") "n_b_ties ", HMMSTR%ntb
   ! n_ab_ties
   write(2,"(A10,I7)") "n_ab_ties ", HMMSTR%ntab
   ! n_nodes
   write(2,"(A8,I4)") "n_nodes ", HMMSTR%n 
   ! k_nodes
   write(2,"(A8,I2)") "k_nodes ", HMMSTR%k!-HMMSTR%u
   ! u_nodes
   write(2,"(A8,I2)") "u_nodes ", HMMSTR%u 
   ! emissions
   write(2,"(A10,4A4)") "emissions ", HMMSTR%emissions(:4)
   do itr =1,4
       write(2,*) adjustl(trim(HMMSTR%emission_alphabet(itr))) 
   enddo
   write(2,"(A10,20F8.5)") "AA_ground ", HMMSTR%AA_ground(:)
   write(2,"(A10,6F8.5)") "SS_ground ", HMMSTR%SS_ground(:)
   write(2,"(A12,11F8.5)") "RAMA_ground ", HMMSTR%RAMA_ground(:)
   write(2,"(A10,10F8.5)") "TM_ground ", HMMSTR%TM_ground(:)
   write(2,"(A12,3F12.5)") "SS_reweight ", HMMSTR%SS_reweight(:)
   write(2,"(A14,11F12.5)") "RAMA_reweight ", HMMSTR%RAMA_reweight(:)
   write(2,"(A14,8F12.5)") "TM_reweight ", HMMSTR%TM_reweight(:)
   call update_seqs(HMMSTR)
   do nitr=0,HMMSTR%n-1 ! write node data
      write(2,"(A5,I4,F8.5)") "node ", nitr, HMMSTR%priors(nitr+1) ! Prior line 1
      write(2,"(A5,I4,A4,A1,20F8.5)")"node ", nitr, " AA ", HMMSTR%AA_seq(nitr+1), HMMSTR%AA(nitr+1,:) ! AA line 2
      write(2,"(A5,I4,A4,A1,6F8.5)") "node ", nitr, " SS ", HMMSTR%SS_seq(nitr+1), HMMSTR%SS(nitr+1,:) ! SS line 3
      write(2,"(A5,I4,A4,A1,10F8.5)") "node ", nitr, " TM ", HMMSTR%TM_seq(nitr+1),HMMSTR%TM(nitr+1,:) ! TM line 4
      write(2,"(A5,I4,A6,A1,11F8.5)") "node ", nitr, " RAMA ", HMMSTR%RAMA_seq(nitr+1),HMMSTR%RAMA(nitr+1,:) ! RAMA line 5
      write(2,"(A5,I4,A5,3F11.5)") "node ",nitr," DIH ",HMMSTR%dihs(nitr+1,:) ! DIH line 6
      write(2,"(A7)") "endnode"
   enddo
   ! transition probs.
   write(2,"(A14)") "START_OUTTRANS"
   do nitr=1,HMMSTR%n
      do itr=1,HMMSTR%n
         if (HMMSTR%a(nitr,itr) > 0.0) then
             write(2,"(I5,I5,F11.5)") nitr-1,itr-1,HMMSTR%a(nitr,itr)
         endif
      enddo
   enddo
   write(2,"(A12)") "END_OUTTRANS"
   write(2,"(A10)") "START_TIES"
   write(2,"(A6)")  "A_TIES"
   !if ( HMMSTR%nta > 0 ) then
           do itr=1,HMMSTR%nta 
              nitr = HMMSTR%aties(itr,1)-1
              nitr_j = HMMSTR%aties(itr,2)-1
              write(2,"(2I5)") nitr, nitr_j
           enddo
   !endif
   write(2,"(A6)") "B_TIES"
   !if ( HMMSTR%ntb > 0 ) then
           do itr=1,HMMSTR%ntb
              nitr = HMMSTR%bties(itr,1)-1
              nitr_j = HMMSTR%bties(itr,2)-1
   !           write(*,*) "BTIE OUTPUT ", nitr, nitr_j
              write(2,"(2I5)",iostat=ios) nitr, nitr_j
              if ( ios /= 0) then
                      write(*,*) ios
                      stop
              endif
           enddo 
   !endif
   write(2,"(A7)") "AB_TIES"
   !if ( HMMSTR%ntab > 0 ) then
           do itr=1,HMMSTR%ntab
              nitr = HMMSTR%abties(itr,1)-1
              nitr_j = HMMSTR%abties(itr,2)-1
              write(2,"(2I5)",iostat=ios) nitr, nitr_j
              if ( ios /= 0) then
                      write(*,*) ios
                      stop
              endif
           enddo
   !endif
   write(2,"(A8)") "END_TIES"
   write(2,"(A3)") "EOF"
   close(1)
   close(2)
   ios = 1
end function write_hmm
subroutine update_seqs(HMMSTR)
    implicit none
    Type(HMM_),intent(inout) :: HMMSTR
    integer  :: i , j, arg_max
    write(*,*) "UPDATING"
    do i=1, HMMSTR%n
      do j = 1, 4
        if ( j == 1 ) then
                arg_max = maxloc(HMMSTR%AA(i,:), 1)
                HMMSTR%AA_seq(i) = HMMSTR%emission_alphabet(j)(arg_max: arg_max)
        else if  ( j == 2 ) then
                arg_max = maxloc(HMMSTR%SS(i,:), 1)
                HMMSTR%SS_seq(i) = HMMSTR%emission_alphabet(j)(arg_max: arg_max)
        else if  ( j == 3 ) then
                arg_max = maxloc(HMMSTR%TM(i,:), 1)
                HMMSTR%TM_seq(i) = HMMSTR%emission_alphabet(j)(arg_max: arg_max)
        else if  ( j == 4 ) then
                arg_max = maxloc(HMMSTR%RAMA(i,:), 1)
                HMMSTR%RAMA_seq(i) = HMMSTR%emission_alphabet(j)(arg_max: arg_max)
                call update_dih(HMMSTR, i)
        endif
      enddo
    enddo
end subroutine update_seqs

subroutine update_dih(HMMSTR, i)
   implicit none
   Type(HMM_), intent(inout) :: HMMSTR
   integer, intent(in) :: i
   real,dimension(11,2)::ramacen
   character(len=11),parameter :: rama_string="HGBEdbeLlxc"
   integer :: loc
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
 
   loc = index(rama_string, HMMSTR%RAMA_seq(i))
   HMMSTR%dihs(i,:2) = ramacen(loc,:)
   if ( HMMSTR%RAMA_seq(i) /= 'c' ) then
           HMMSTR%dihs(i,3) = 180.0
   else
           HMMSTR%dihs(i,3) = 0.00
   endif
end subroutine update_dih

type(HMM_) function read_hmm(hmmfile, QUIET) result(HMMSTR)
    implicit none
    character(len=250),intent(in) :: hmmfile
    character(len=300) :: aline
    logical,intent(in), optional :: QUIET
    integer :: ios, iarg, i, j, temp
    character(len=20) :: arg
    !type(HMM_),intent(in) :: HMMSTR
    character(len=30),dimension(30) :: args
    integer :: lb, itr, nitr,titr ! loop break
    real,dimension(20) :: sumAA = 0 
    real,dimension(6)  :: sumSS = 0 
    real,dimension(10) :: sumTM = 0
    real,dimension(11) :: sumRAMA = 0
    logical,dimension(3) :: tie_type
    tie_type = [ .false., .false., .false.]
  ! READ MODEL
    HMMSTR%model_file = hmmfile
    nitr = 1
    open(1, file=hmmfile,iostat=ios)
    if (ios/=0) then
        write(*,*) "Error 2. Improperly formatted hmm. Error on open."
        stop
    endif
    do
      read(1,'(a)',iostat=ios) aline
      if (ios/=0) then
          write(*,*) "Error 2. Improperly formatted hmm. Error on read."
          write(*,*) ios
          stop
      endif
      args = get_substrings(aline)
      !do itr=1,size(args,dim=1)
      !write(*,*) args(itr)
      !enddo
      
      select case (trim(args(1)))
      case ('COMMENT', 'comment')
          continue
      case ('AA_EPS')
          read(args(2),*) HMMSTR%AA_EPS
      case ('SS_EPS')
          read(args(2),*) HMMSTR%SS_EPS
      case ('TM_EPS')
          read(args(2),*) HMMSTR%TM_EPS
      case ('RAMA_EPS')
          read(args(2),*) HMMSTR%RAMA_EPS
      case ('n_a_ties')
          read(args(2),*) HMMSTR%nta
          allocate(HMMSTR%aties(HMMSTR%nta,2)); HMMSTR%aties = 0 
      case ('n_b_ties')
          read(args(2),*) HMMSTR%ntb
          allocate(HMMSTR%bties(HMMSTR%ntb,2)); HMMSTR%bties = 0 
      case ('n_ab_ties', 'n_ties')
          read(args(2),*) HMMSTR%ntab
          allocate(HMMSTR%abties(HMMSTR%ntab,2)); HMMSTR%abties = 0 
      case ('n_nodes')
          read(args(2),*)HMMSTR%n
          allocate(HMMSTR%AA(HMMSTR%n,20))   ; HMMSTR%AA = 0
          allocate(HMMSTR%SS(HMMSTR%n,6))    ; HMMSTR%SS = 0
          allocate(HMMSTR%TM(HMMSTR%n,10))   ; HMMSTR%TM = 0
          allocate(HMMSTR%RAMA(HMMSTR%n,11)) ; HMMSTR%RAMA = 0
          allocate(HMMSTR%dihs(HMMSTR%n,3))  ; HMMSTR%dihs = 0
          allocate(HMMSTR%priors(HMMSTR%n))  ; HMMSTR%priors = 0
          allocate(HMMSTR%AA_seq(HMMSTR%n))  
          allocate(HMMSTR%SS_seq(HMMSTR%n))  
          allocate(HMMSTR%TM_seq(HMMSTR%n))  
          allocate(HMMSTR%RAMA_seq(HMMSTR%n))
          allocate(HMMSTR%AA_ground(20))  ; HMMSTR%AA_ground = 0
          allocate(HMMSTR%SS_ground(6))   ; HMMSTR%SS_ground = 0 
          allocate(HMMSTR%SS_reweight(3))   ; HMMSTR%SS_reweight = 0 
          allocate(HMMSTR%TM_ground(10))  ; HMMSTR%TM_ground = 0 
          allocate(HMMSTR%RAMA_ground(11)); HMMSTR%RAMA_ground = 0 
          allocate(HMMSTR%RAMA_reweight(6)); HMMSTR%RAMA_reweight = 0 
          allocate(HMMSTR%TM_reweight(8)); HMMSTR%TM_reweight = 0 
          allocate(HMMSTR%a(HMMSTR%n,HMMSTR%n)); HMMSTR%a = 0 
          allocate(HMMSTR%atie_lookup(HMMSTR%n)); HMMSTR%atie_lookup = 0 
          allocate(HMMSTR%abtie_lookup(HMMSTR%n)); HMMSTR%abtie_lookup = 0 
          allocate(HMMSTR%btie_lookup(HMMSTR%n,400)); HMMSTR%btie_lookup = 0 
          do itr = 1, HMMSTR%n
             HMMSTR%atie_lookup(itr) = 0
             HMMSTR%abtie_lookup(itr) = 0 
             HMMSTR%btie_lookup(itr,1) = 1
          enddo
  !        write(*,*) "LOADED N = ",HMMSTR%n
      case ('k_nodes')
          read(args(2),*)HMMSTR%k
  !        write(*,*) "LOADED K = ", HMMSTR%k
      case ('u_nodes')
          read(args(2),*)HMMSTR%u
          !HMMSTR%k = HMMSTR%k + HMMSTR%u
  !        write(*,*) "LOADED U = ", HMMSTR%u
      case ('AA_ground')
          do itr=2,21
                read(args(itr),*)HMMSTR%AA_ground(itr-1)
          enddo
  !        write(*,*) "LOADED backgroud = ", HMMSTR%bground
      case ('SS_ground')
          do itr=2,7
                read(args(itr),*)HMMSTR%SS_ground(itr-1)
          enddo
      case ('RAMA_ground')
          do itr=2,12
                read(args(itr),*)HMMSTR%RAMA_ground(itr-1)
          enddo
      case ('RAMA_reweight')
          do itr = 2, 7
                read(args(itr),*) HMMSTR%RAMA_reweight(itr-1)
          enddo
      case ('SS_reweight')
          do itr = 2, 4
                read(args(itr),*) HMMSTR%SS_reweight(itr-1)
          enddo
      case ('TM_reweight')
          do itr = 2, 9
                read(args(itr),*) HMMSTR%TM_reweight(itr-1)
          enddo
      case ('TM_ground')
          do itr=2,11
                read(args(itr),*)HMMSTR%TM_ground(itr-1)
          enddo
      case ('EMISSIONS', 'emissions') ! LOADING EMISSIONS
          itr = 1
          do
            if (itr > 4) exit
            read(args(1+itr),*)HMMSTR%emissions(itr)
            itr = itr + 1
          enddo
  !        write(*,*) "LOADED EMISSION TYPES: ", HMMSTR%emissions
          itr = 1
          do 
            if (itr > 4) exit
            read(1,'(a)',iostat=ios) aline
            if (ios/=0) then
                write(*,*) "Error 2. Improperly formatted hmm."
                stop
            endif
            args = get_substrings(aline)
            read(args(1),*)HMMSTR%emission_alphabet(itr)
            itr = itr + 1
          enddo 
  !        write(*,*) "LOADED EMISSION abcdetc: ", HMMSTR%emission_alphabet
      case ('NODE', 'node') ! LOADING NODES
          read(args(3),*)HMMSTR%priors(nitr)
          do ! get all data from this node.
              read(1,'(a)',iostat=ios) aline
              if (ios/=0) then
                  write(*,*) aline
                  write(*,*) "Error 2. Improperly formatted hmm."
                  stop
              endif
              args = get_substrings(aline)
              if (trim(args(1)) == 'endnode' .or. trim(args(1)) == 'ENDNODE') then 
                      nitr = nitr + 1 
                      exit ! eol
              endif
              select case (args(3))
              case('AA','aa')
                      read(args(4),*)HMMSTR%AA_seq(nitr)
                      do itr = 5,24
                           read(args(itr),*)HMMSTR%AA(nitr,itr-4)
                      enddo
            !          sumAA = sumAA + HMMSTR%priors(nitr) * HMMSTR%AA(nitr,:)
   !                  write(*,*) "Loaded node ", nitr, " AA ", HMMSTR%nodes(nitr)%AA(:)
              case('SS','ss')
                      read(args(4),*)HMMSTR%SS_seq(nitr)
                      do itr=5,10
                          read(args(itr),*)HMMSTR%SS(nitr,itr-4)
                      enddo            
           !           sumSS = sumSS + HMMSTR%priors(nitr) * HMMSTR%SS(nitr,:)   
   !                   write(*,*) "Loaded node ", nitr, " SS ", HMMSTR%nodes(nitr)%SS(:)
              case('TM','tm')
                      read(args(4),*)HMMSTR%TM_seq(nitr)
                      do itr=5,14
                          read(args(itr),*)HMMSTR%TM(nitr,itr-4)
                      enddo
           !           sumTM = sumTM + HMMSTR%priors(nitr) * HMMSTR%TM(nitr,:)   
   !                   write(*,*) "Loaded node ", nitr, " TM ", HMMSTR%nodes(nitr)%TM(:)
              case('RAMA','rama')
                      read(args(4),*)HMMSTR%RAMA_seq(nitr)
                      do itr=5,15
                           read(args(itr),*)HMMSTR%RAMA(nitr,itr-4)
                      enddo
           !           sumRAMA = sumRAMA + HMMSTR%priors(nitr) * HMMSTR%RAMA(nitr, :)
   !                   write(*,*) "Loaded node ", nitr, " RAMA ", HMMSTR%nodes(nitr)%RAMA(:)
              case('DIH','dih')
                      do itr=4,6
                           read(args(itr),*)HMMSTR%dihs(nitr,itr-3)
                      enddo
              case('label','LABEL')
                      continue
   !                   write(*,*) "Loaded node ", nitr, " DIH ", HMMSTR%nodes(nitr)%dihs(:)
              case default
                  write(*,*) "Error 3. Unhandled Keyword (while parsing a node).", args(3)
                  stop
              end select
              !write(*,*) "Loaded node ", nitr 
          enddo
      case ('START_OUTTRANS','start_outtrans')
           HMMSTR%a = 0
           do
              read(1,'(a)',iostat=ios) aline
              if (ios/=0) then
                exit
              endif
              args = get_substrings(aline)
              if (trim(args(1)) == 'END_OUTTRANS' .or. trim(args(1)) == 'end_outtrans') exit
              read(args(1),*) i
              i = i + 1
              read(args(2),*) j
              j = j + 1
              read(args(3),*,iostat=ios)HMMSTR%a(i,j)
              !write(*,*) i,j, HMMSTR%transitions(i,j)
           enddo
     !     write(*,*)HMMSTR%transitions
       case ('START_TIES', 'start_ties')
           do itr=1,HMMSTR%nta+HMMSTR%ntb+HMMSTR%ntab+4
               !write(*,*) itr
               titr = titr + 1
               read(1,'(a)',iostat=ios) aline
               args = get_substrings(aline)
               if (trim(args(1)) == "A_TIES") then
                 titr = 0
                 tie_type = [ .true., .false., .false. ]
                 cycle
               else if (trim(args(1)) == "B_TIES") then
                 titr = 0
                 tie_type = [ .false., .true., .false. ]
                 cycle
               else if (trim(args(1)) == "AB_TIES") then
                 titr = 0
                 tie_type = [ .false., .false., .true. ]
                 cycle
               endif
               if (trim(args(1)) == "END_TIES" .or. trim(args(1)) == "end_ties") then
                       exit
               endif
               !write(*,*) args
               read(args(1),*) i
               read(args(2),*) j
               i = i + 1
               j = j + 1
               if ( tie_type(1) ) then
                       HMMSTR%aties(titr,1) =  i
                       HMMSTR%aties(titr,2) =  j
                       HMMSTR%atie_lookup(i) = titr
               else if ( tie_type(2) ) then
               !        write(*,*) "BTIE", i, j
                       HMMSTR%bties(titr,1) =  i
                       HMMSTR%bties(titr,2) =  j
                       HMMSTR%btie_lookup(j,1) = HMMSTR%btie_lookup(j,1) + 1
                       HMMSTR%btie_lookup(j,HMMSTR%btie_lookup(j,1)) = i
               else if ( tie_type(3) ) then
                       HMMSTR%abties(titr,1) =  i
                       HMMSTR%abties(titr,2) =  j
                       HMMSTR%abtie_lookup(i) = titr
               endif
            enddo
       case ('EOF')
          if (QUIET .eqv. .false.) then
                write(*,*) "Successfully parsed " // trim(hmmfile) // " in hmm_io.f95."
          endif
          exit ! <-- normal termination of do loop
       case default
          write(*,*) args
          write(*,*) "Error 3. Unhandled Keyword."
          stop 
      end select 
    enddo
    !HMMSTR%AA_ground = sumAA / sum(sumAA)
    !HMMSTR%SS_ground = sumSS / sum(sumSS)
    !HMMSTR%TM_ground = sumTM / sum(sumTM)
    !HMMSTR%RAMA_ground = sumRAMA / sum(sumRAMA)
    if (QUIET .eqv. .false.) then
           ! write(*,*) "default bg freqs = ", HMMSTR%bground
            write(*,*) "AA_ground = ", HMMSTR%AA_ground
            write(*,*) "SS_ground = ", HMMSTR%SS_ground
            write(*,*) "TM_ground = ", HMMSTR%TM_ground
            write(*,*) "RAMA_ground = ", HMMSTR%RAMA_ground
    endif
    !write(*,*) HMMSTR%a
    close(1)
end function read_hmm

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

!-------------------------------------------------------------------------------------------------------------
type(HMM_) function release_HMM(HMMSTR)
   implicit none
   type(HMM_),intent(inout)::HMMSTR
   if(allocated(HMMSTR%AA))         deallocate(HMMSTR%AA)
   if(allocated(HMMSTR%SS))         deallocate(HMMSTR%SS)
   if(allocated(HMMSTR%TM))         deallocate(HMMSTR%TM)
   if(allocated(HMMSTR%RAMA))       deallocate(HMMSTR%RAMA)
   if(allocated(HMMSTR%dihs))       deallocate(HMMSTR%dihs)
   if(allocated(HMMSTR%priors))     deallocate(HMMSTR%priors)
   if(allocated(HMMSTR%AA_seq))     deallocate(HMMSTR%AA_seq)
   if(allocated(HMMSTR%SS_seq))     deallocate(HMMSTR%SS_seq)
   if(allocated(HMMSTR%TM_seq))     deallocate(HMMSTR%TM_seq)
   if(allocated(HMMSTR%RAMA_seq))   deallocate(HMMSTR%RAMA_seq)
   if(allocated(HMMSTR%AA_ground))  deallocate(HMMSTR%AA_ground)
   if(allocated(HMMSTR%SS_ground))  deallocate(HMMSTR%SS_ground)
   if(allocated(HMMSTR%TM_ground))  deallocate(HMMSTR%TM_ground)
   if(allocated(HMMSTR%RAMA_ground))deallocate(HMMSTR%RAMA_ground)
   if(allocated(HMMSTR%a))          deallocate(HMMSTR%a)
  end function release_HMM
!------------------------------------------------------------------------------------------------------------
end module hmm_io
