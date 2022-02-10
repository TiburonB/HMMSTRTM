program gamma_read
       ! this program reads a specific time point of a specific gamma file and prints it 
! -------------------------------------------------------------------
    use hmm_io
    use gamma_io
    implicit none
! -------------------------------------------------------------------
    character(len=250) :: modelfile
    character(len=250) :: gammafile
    character(len=4) :: titr_str, szstr
    real(8),dimension(:, :),allocatable :: tgamma ! (nres, HMMSTR%n)
    character(len=20) :: format_string
    integer :: titr, nres
    integer :: iarg
    type(HMM_) :: HMMSTR
! -------------------------------------------------------------------

    iarg = command_argument_count()
    if (iarg < 3) then
        !write(*,*) 'To run using database: xscratch hmmfile drctfile code'''
        !stop       ! ./a.out ./HMMSTR/HMMSTR.hmm ./final.drct 7bpo
    end if         ! ./a.out ./HMMSTR/HMMSTR.hmm ./tmp/7bpoA.profile
    call getarg(1, modelfile)
    call getarg(2, gammafile)
    call getarg(3, titr_str)
    read(titr_str,*) titr
    call getarg(4, titr_str)
    read(titr_str,*) nres
    HMMSTR = read_hmm(modelfile, .true.) 
    allocate(tgamma(nres,HMMSTR%n))
    tgamma = 0
    call read_gamma(tgamma, gammafile, nres, HMMSTR%n)
    write(szstr,'(I4)') HMMSTR%n
    format_string = '('//szstr//'F8.5)'
    write(*, format_string) tgamma(titr,:)
    deallocate(tgamma)

end program gamma_read
