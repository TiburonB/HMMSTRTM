program max_gamma
   use gamma_io
   use hmm_io
   implicit none
    
   character(len=250) :: modelfile
   character(len=250) :: gammadrctfile
   character(len=250) :: outfile
   character(len=4) :: ec
   integer :: ios
   integer :: epochs ! number of epochs to train for
   integer :: iarg
   type(HMM_) :: HMMSTR
   integer :: record, nitr, nres, titr, itr
   real(8),dimension(:,:),allocatable :: xgamma
   real(8) :: sum_gamma
   real(8),dimension(:),allocatable :: max_gamma_v
   character(len=5),dimension(:),allocatable :: max_gamma_code_chain
   integer,dimension(:),allocatable :: max_gamma_ind
   character(len=1) :: chain_id
   character(len=4) :: code

    iarg = command_argument_count()
    if (iarg < 2) then
        write(*,*) 'Usage: xhmmstr_train model_R.hmm drctfile '
        stop    
    end if
    call getarg(1, modelfile)
    call getarg(2, gammadrctfile)
    call getarg(3, outfile)
    
    !write(*,*) "GAMMADRCTFILE = ", gammadrctfile
    HMMSTR = read_hmm(modelfile, .false.) ! READ THE MODEL
    allocate(max_gamma_v(HMMSTR%n))
    allocate(max_gamma_code_chain(HMMSTR%n))
    allocate(max_gamma_ind(HMMSTR%n))
    do nitr =1 , HMMSTR%n
        max_gamma_code_chain(nitr) = "NONEA"
        max_gamma_ind(nitr) = -1
    enddo
    max_gamma_v = 0 
    sum_gamma = 0 
    record = 1
    itr = 1
    do while (record /= -1 )
        !write(*,*) itr
        call loop_gamma(xgamma, gammadrctfile, HMMSTR%n, nres, record, chain_id, code)
        do titr = 1, nres
              sum_gamma = sum(xgamma(titr, :))
              do nitr = 1, HMMSTR%n
                 if ( xgamma(titr,nitr) / sum_gamma > max_gamma_v(nitr) ) then
                      max_gamma_v(nitr) = xgamma(titr, nitr) / sum_gamma
                      max_gamma_code_chain(nitr) = code // chain_id
                      max_gamma_ind(nitr) = titr
                 endif
                 !sum_gamma(nitr) = sum_gamma(nitr) + sum(xgamma(:,nitr)) ! gamma_heat
              enddo
        enddo
        itr = itr + 1
    enddo
    open(4, file=outfile, iostat = ios)
    if ( ios /= 0) then
            write(*,*) ios
            write(*,*) 'unable to open output of max_gamma file'
    endif
    do nitr = 1, HMMSTR%n
            write(4,*) nitr, max_gamma_code_chain(nitr), max_gamma_ind(nitr)
            !write(*,*) nitr, sum_gamma(nitr) ! gamma heat
    enddo
    close(4)
    deallocate(max_gamma_v)
    deallocate(max_gamma_code_chain)
    deallocate(max_gamma_ind)
end program max_gamma
