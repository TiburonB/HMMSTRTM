program gamma_heat
   use gamma_io
   use hmm_io
   implicit none
    
   character(len=250) :: modelfile
   character(len=250) :: gammadrctfile
   character(len=4) :: ec
   integer :: epochs ! number of epochs to train for
   integer :: iarg
   type(HMM_) :: HMMSTR
   integer :: record, nitr, nres, titr, itr
   real(8),dimension(:,:),allocatable :: xgamma
   real(8),dimension(:),allocatable :: sum_gamma
   character(len=1) :: chain_id
   character(len=4) :: code
    iarg = command_argument_count()
    if (iarg < 2) then
        write(*,*) 'Usage: xhmmstr_train model_R.hmm drctfile '
        stop    
    end if
    call getarg(1, modelfile)
    call getarg(2, gammadrctfile)

    HMMSTR = read_hmm(modelfile, .false.) ! READ THE MODEL
    allocate(sum_gamma(HMMSTR%n))
    sum_gamma = 0 
    record = 1
    itr = 1
    do while (record /= -1 )
        call loop_gamma(xgamma, gammadrctfile, HMMSTR%n, nres, record, chain_id, code)
        do nitr = 1, HMMSTR%n
           sum_gamma(nitr) = sum_gamma(nitr) + sum(xgamma(:,nitr))
        enddo
        itr = itr + 1
    enddo
    do nitr = 1, HMMSTR%n
            write(*,*) nitr, sum_gamma(nitr)
    enddo
    deallocate(sum_gamma)
end program gamma_heat
