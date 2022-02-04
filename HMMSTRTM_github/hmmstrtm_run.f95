program hmmstrtm_run
! -------------------------------------------------------------------
    use hmmstrtm
    implicit none
! -------------------------------------------------------------------
    character(len=250) :: modelfile
    character(len=250) :: drctfile
    character(len=4) :: code
    integer :: iarg
! -------------------------------------------------------------------

    iarg = command_argument_count()
    if (iarg < 3) then
        write(*,*) 'Usage: xhmmstr_train model_R.hmm drctfile epochs '
        stop       ! ./a.out ./HMMS/HMMSTR.hmm ./final.drct 2
    end if
    call getarg(1, modelfile)
    call getarg(2, drctfile)
    call getarg(3, code)

    write(*,*) "SCRATCHING CODE :", code
    call scratch_hmmstr(modelfile, drctfile, code)

end program hmmstrtm_run
