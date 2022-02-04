program hmm_train
! -------------------------------------------------------------------
    use hmm
    implicit none
! -------------------------------------------------------------------
    character(len=250) :: modelfile
    character(len=250) :: drctfile
    character(len=4) :: ec
    integer :: epochs ! number of epochs to train for
    integer :: iarg
! -------------------------------------------------------------------

    iarg = command_argument_count()
    if (iarg < 3) then
        write(*,*) 'Usage: xhmmstr_train model_R.hmm drctfile epochs '
        stop       ! ./a.out ./HMMS/HMMSTR.hmm ./final.drct 2
    end if
    call getarg(1, modelfile)
    call getarg(2, drctfile)
    call getarg(3, ec)

    write(*,*) "EC ", ec
    read(ec, '(i2)') epochs
    write(*,*) "TRAINING FOR ", epochs, " EPOCHS."
    call train_hmm(modelfile, drctfile, epochs)

end program hmm_train
