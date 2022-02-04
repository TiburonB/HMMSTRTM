program hmm_inject
! -------------------------------------------------------------------
    use profile_io
    implicit none
! -------------------------------------------------------------------
    character(len=5000) :: seq
    character,dimension(:),allocatable :: aa_seq
    integer,dimension(:),allocatable :: seq_code
    real(8),dimension(:,:),allocatable :: aa_profile
    type(Profile) :: seq_profile
    character(len=4) :: code
    character(len=30) :: HMMFILE
    character(len=1) :: chain
    integer :: nres, titr
    integer :: iarg
! -------------------------------------------------------------------

    code = 'null'
    chain = 'X'

    iarg = command_argument_count()
    if (iarg < 2) then
        write(*,*) 'Usage: xinject HMMFILE sequence <code> <chain>'
        stop  
    else     
        call getarg(1,HMMFILE)
        call getarg(2, seq)
    end if
    if (iarg >= 3) then
            call getarg(3, code)
    endif
    if (iarg >= 4) then
            call getarg(4, chain)
    endif
    seq = trim(seq)
    !write(*,*) trim(seq)
    nres = len(trim(seq))
    !write(*,*) nres
    allocate(aa_seq(nres))
    do titr=1,nres
        aa_seq(titr) = seq(titr:titr)
    enddo
    !write(*,*) aa_seq
    allocate(seq_code(nres))
    seq_code(:) = get_seq_code(aa_seq, nres)
    allocate(aa_profile(20,nres))
    aa_profile = 0
    allocate(seq_profile%seq_code(nres))
    do titr =1, nres
         aa_profile(seq_code(titr),titr) = 1
         seq_profile%profiles(:,titr) = aa_profile(:,titr)
         seq_profile%seq_code(titr) = seq_code(titr)
    enddo
    seq_profile%nres = nres
    seq_profile%code = code
    seq_profile%chain = chain


    write(*,*) "Writing to file"
    call write_profile(seq_profile, './tmp/'//code//chain//'.profile')
    deallocate(seq_code)
    deallocate(aa_profile)

    write(*,*) './xscratch '//HMMFILE//' ./tmp/'//code//chain//'.profile'
    call execute_command_line('./xscratch '//HMMFILE//' ./tmp/'//code//chain//'.profile')

end program hmm_inject
