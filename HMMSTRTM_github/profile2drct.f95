program profile2drct
    implicit none
    integer,parameter :: MAXLINE = 500
    character(len=4) :: code, ocode
    character(len=1) :: chainID, ochain, sschar, aa, mrchar, ramachar
    character(len=1) :: rotamer, domain ! future use
    integer(2) :: iseq, resseq, i
    integer :: nrec, prec, vrecsize, J, K, ios, trecsize
    integer,dimension(50) :: A, B
    integer(2),dimension(3) :: clusterID
    real :: csum, x
    real,dimension(3) :: calpha, dih
    real,dimension(20) :: profile
    character(len=100) filename ! profile file
    character (len=MAXLINE) :: aline
    character(len=23) :: res1
    logical :: gapchar, putchar
    integer :: jarg
    integer*4 :: N
    character(len=100) :: profilefile
    character(len=5) :: code_chain
    integer :: c, offset
    ! AA ground frequencies.
    real ,dimension(20) :: Fi=(/ 0.08279,0.01937,0.05855,0.05992, &
        0.04014,0.08089,0.02275,0.05552,0.05959,0.08020,0.02070,0.04729, &
        0.04599,0.03728,0.04640,0.06246,0.05888,0.06866,0.01507,0.03756/) 

    character :: caps
    res1 = 'ACDEFGHIKLMNPQRSTVWYBXZ'
    
    vrecsize = 136 ! 136 bytes... machine dependent .. also try 34 bytes
    
    jarg = command_argument_count() ! number of commands from command line
    
    rotamer = ' '
    clusterID = 0
    offset = 0
    if (jarg < 3) then
        write(*,*) 'Usage: xprofile2drct profilefile codechain outfile'
        stop
    else
        call getarg(1, profilefile)
        call getarg(2, code_chain)
        call getarg(3, filename)
    endif

    code = code_chain(1:4)
    ochain = code_chain(5:5)

    write(*,*) "Output file is ", trim(filename)
    open(3, file=profilefile)
    open(2, file=filename,status='replace',form='unformatted',&
            access='direct',recl=vrecsize)
    
    chainID = ochain
    aa = 'X'
    nrec = 0
    prec = 0
    csum = 1.0
    sschar = 'U'
    rotamer = ' '
    calpha = 0.0
    dih = 0.0
    iseq = 0
    aline = ' '

    do
        read (3, '(a)', iostat = ios) aline
        if (ios/=0) exit
        if (aline(1:1) == "#") cycle ! comment line
        if (aline(1:1) == ">") cycle ! comment line
        N = getfields(50, aline, A, B)
        !! write(*,*) aline, N
        !! first field : residue number
        read(aline(A(2):B(2)),*,iostat=ios) iseq
        if (ios/=0) stop 'ERROR reading first field: sequence number.'
        !! write(*,*) 'iseq=', iseq
        !! second field : amino acid character
        aa = aline(A(3):A(3))
        !! write(*,*) 'aa=', aa
        !! field 3 ~ 22 : profile
        !! diagnostic
        !! write(*, '(a)') aline(A(3):)
        read(aline(A(4):),*,iostat=ios) profile(1:20)
        if (ios /=0) stop 'ERROR reading profile, check format.'
        x = 0. ! sum of profile
        do i = 1,20 ! loop over val in profile
            x = x + profile(i)
        enddo
        profile = profile/x ! normalize values in profile
        nrec = nrec + 1
        sschar = ' ' 
        ramachar = ' '
        mrchar = ' '
        trecsize = sizeof(iseq)*2 + sizeof(code) + sizeof(chainID) + sizeof(aa) + sizeof(sschar) + sizeof(rotamer) + &
                sizeof(calpha) + sizeof(dih) + sizeof(csum) + sizeof(profile) + sizeof(0.) + sizeof(0.) + sizeof(ramachar) + &
                sizeof(mrchar) + sizeof(clusterID)
        if (trecsize /= vrecsize) then
                write(*,*) "Attempting to write record of size ", trecsize
                write(*,*) "sizeof(iseq)", sizeof(iseq)
                write(*,*) "sizeof(iseq)", sizeof(iseq)
                write(*,*) "sizeof(code)", sizeof(code)
                write(*,*) "sizeof(chainID)", sizeof(chainID)
                write(*,*) "sizeof(aa)", sizeof(aa)
                write(*,*) "sizeof(sschar)", sizeof(sschar)
                write(*,*) "sizeof(rotamer)", sizeof(rotamer)
                write(*,*) "sizeof(calpha)", sizeof(calpha)
                write(*,*) "sizeof(dih)", sizeof(dih)
                write(*,*) "sizeof(csum)", sizeof(csum)
                write(*,*) "sizeof(profile)", sizeof(profile)
                write(*,*) "sizeof(0.)", sizeof(0.)
                write(*,*) "sizeof(0.)", sizeof(0.)
                write(*,*) "sizeof(ramachar)", sizeof(ramachar)
                write(*,*) "sizeof(mrchar)", sizeof(mrchar)
                write(*,*) "sizeof(clusterID)", sizeof(clusterID)
                stop
        endif
!  # record size: 136 bytes
!  # (0)seq_list (1)resSeq (2)structure_name (3)chainID (4)seq_resid (5)ss (6)rotamer (7-9)xyzca (10-12)ang (13)csum (14-33)vall_pro (34)gap_freq (35)ins_freq (36)domainID (37)ctx (38-40)clusterID
!  # $template = "s s A4 A A A A f f f f f f f f f f f f f f f f f f f f f f f f f f f f f A A s3";
        write(2, rec = nrec) iseq, iseq, code, chainID, aa, ' ', rotamer, &
                   calpha, dih, csum, profile, 0., 0., ' ', '1', clusterID
   enddo

   write(*, '('' Sequence Length: '', i10)') iseq
   write(*, '('' Last Record: '', i10)') nrec
   prec = 4 * nrec
   write (*, '('' Total Size: '' , i10, '' bytes. '')') prec
   close(2)
contains
! ---------------------------------------------------------------------------------
! given character,dimension(:)::card, return card as uppercase
function to_ucase(card) result(ucase)
    implicit none
    character,dimension(:),intent(in) :: card
    integer :: I, icap
    character,dimension(len(card)) :: ucase
    icap = ichar('A') - ichar('a')
    do I = 1, len(card)
       if (ichar('a') <= ichar(card(I)) .and. ichar(card(I)) <= ichar('z')) then
           ucase(I:I) = char(ichar(card(I:I))+icap)
       else
           ucase(I:I) = card(I:I)
       endif
    enddo
end function to_ucase
! ---------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------
function to_lcase(card) result(lcase)
    implicit none
    character,dimension(:),intent(in) :: card
    integer :: I, icap
    character,dimension(len(card)) :: lcase
    icap = ichar('A') - ichar('a')
    do I = 1, len(card)
       if (ichar('A') <= ichar(card(I)) .and. ichar(card(I)) <= ichar('Z')) then
           lcase(I:I) = char(ichar(card(I:I))-icap)
       else
           lcase(I:I) = card(I:I)
       endif
    enddo
end function to_lcase
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! find the limits of up to N fields of non-blank characters in 'aline'
! The start of field I is A(I), the end is B(I)
! e.g. A(1) is starting position of field 1, B(1) is ending position of field 1
! return value for getfields is the number of fields...
integer function getfields(N, aline, A, B) 
    implicit none
    integer,dimension(50) :: A, B
    character(len=*),intent(in) :: aline
    integer :: ret
    integer :: N
    integer :: I, FF, LL
    
    I = 1 ! position in aline
    FF = 0 ! Field # 
    LL = len(aline) ! number of characters in aline
    do while ( FF < N .and. I <= LL ) 
        if ( I > LL ) exit ! early termination of loop 
        do while ( aline(I:I) == ' ' .and. I <= LL ) 
           I = I + 1
        enddo
        if (I <= LL) then
           FF = FF + 1
           A(FF) = I
           do while ( aline(I:I) /= ' ' .and. I <= LL )
               I = I + 1
           enddo
           B(FF) = I - 1
        end if
    enddo
    getfields = FF
end function getfields
! --------------------------------------------------------------------------------
end program profile2drct
