program hmm_run
! -------------------------------------------------------------------
    use hmm
    implicit none
! -------------------------------------------------------------------
    character(len=250) :: modelfile
    character(len=250) :: drctfile
    character(len=4) :: code
    integer :: iarg
! -------------------------------------------------------------------

    code = 'null'
    iarg = command_argument_count()
    if (iarg < 3) then
        !write(*,*) 'To run using database: xscratch hmmfile drctfile code'''
        !stop       ! ./a.out ./HMMSTR/HMMSTR.hmm ./final.drct 7bpo
    end if         ! ./a.out ./HMMSTR/HMMSTR.hmm ./tmp/7bpoA.profile
    call getarg(1, modelfile)
    call getarg(2, drctfile)
    if (iarg >= 3) then
            call getarg(3, code)
    endif

    write(*,*) "SCRATCHING CODE :", code
    call scratch_hmm(modelfile, drctfile, code)

end program hmm_run
