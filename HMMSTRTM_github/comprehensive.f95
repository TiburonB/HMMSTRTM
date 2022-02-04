! comprehensive.f95
! TLB 10/12/21 
program comprehensive
! -------------------------------------------------------------------
    use profile_io
   ! use drct2gamma ! DIAGNOSTICSArchived
    use hmm_io
    use hmm
    implicit none
! -------------------------------------------------------------------
    character(len=250) :: gamma_dir
    character(len=250) :: modelfile
    character(len=250) :: drctfile
    character(len=4) :: chain_id
    integer :: iarg
    integer :: t, titr, nitr, nres, itr, bitr
    character(len=1),dimension(:),allocatable :: seq1, tmseq1, ssseq1, ramaseq1
    type(HMM_) :: HMMSTR1
    integer :: sitr    
    real,dimension(:,:),allocatable::aaprofile
    real,dimension(:,:),allocatable::angles
  !character(len=1),dimension(:),allocatable::aaseq
  character(len=1),dimension(:),allocatable::ramaseq
  character(len=1),dimension(:),allocatable::ssseq
  character(len=1),dimension(:),allocatable::tmseq
  ! mda stuff ! --------------------------------------- 9/8/21
  integer,dimension(4) :: training_flags
  real(8),dimension(:),allocatable::ct
  real(8),dimension(:,:),allocatable::alpha
  real(8),dimension(:,:),allocatable::beta
  !real(8),dimension(:, :, :),allocatable :: zeta ! (nres, HMMSTR%n, HMMSTR%n)
  real(8),dimension(:, :),allocatable :: tgamma ! (nres, HMMSTR%n)
  logical :: viterbi_voting
  ! TM_ROC
  real :: TMness
  integer :: TM_correct, TM_bool
  real(8) :: sum_TMgamma
  integer :: tm_reference_ind
  character(len=mtm),parameter   :: reference_tm_string="12BHCILFU_"
  character(len=1) :: tm_c
  !!!!!!!!!!!!!!!!!!!!!!VOTING METHOD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,dimension(4,4) :: Q3_confusion
    integer,dimension(6,6) :: RAMA_confusion
    integer,dimension(8,8) :: TM_confusion
    character(len=3),dimension(4) :: Q3_cats
    character(len=5),dimension(6) :: RAMA_cats
    character(len=3),dimension(8) :: TM_cats
     real(8),dimension(20) :: aa_prof
     real(8),dimension(6) :: ss_prof
     real(8),dimension(11) :: rama_prof
     real(8),dimension(10) :: tm_prof
    real(8),dimension(:,:),allocatable :: tm_prof_c, rama_prof_c, ss_prof_c, aa_prof_c
    real(8) :: acc
  character(len=1),dimension(:,:),allocatable :: paradigm, ground_truth_profiles
  real,dimension(2) :: mda_wa ! mda_weight/accuracy
  character(len=1) :: run_gamma
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! model file organizer ..! 1 10.12.21
  character(len=100) :: model_path, model_name, gdir, gfile
  character(len=3) :: mitr, mitrm1
  integer :: model_iteration, dind, slind, gdl
  logical :: dir_e, gf_e ! gdir exists, gfile exists
  logical :: tm
  integer :: last_g_record 
  real(8),dimension(:,:),allocatable :: zeta
  integer,dimension(:,:),allocatable :: zeta_dict
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !character(len=1),dimension(:),allocatable :: mda
  integer,dimension(:),allocatable :: mda

    !integer,parameter::m=20
    !integer,parameter::mdssp=6
    !integer,parameter::mtm=10
    !integer,parameter::mrama=11
    !character(len=m),parameter     :: aa="ACDEFGHIKLMNPQRSTVWY"
    !character(len=mrama),parameter :: rama_string="HGBEdbeLlxc"
    !character(len=mdssp),parameter :: ss_string="HEGST_"
! -------------------------------------------------------------------
    Q3_cats = [ "S/T", "G/H", " E " , " _ " ]
    RAMA_cats = [ " G/H ", " L/l ", " B/b ", "E/e/d", "  x  ", "  c  " ]
    TM_cats = [ "1/2", " B ", " H ", " C ", " I ", " L ", " F ", " U " ]
    Q3_confusion = 0
    RAMA_confusion = 0
    TM_confusion = 0

    training_flags = [1, 0, 0, 0]
    iarg = command_argument_count()
    if (iarg == 3) then
        !write(*,*) 'Usage: modelfile drctfile '
        call getarg(1, modelfile)
        call getarg(2, drctfile)
       call getarg(3, run_gamma)
        chain_id = "NONE"
    else if ( iarg == 4) then
        !write(*,*) 'Usage: modelfile drctfile '
     call getarg(1, modelfile)
     call getarg(2, drctfile)
     call getarg(3, chain_id)
     call getarg(4, run_gamma)
   end if
  slind = index(modelfile, "/", .true.)
  dind = index(modelfile, ".", .true.)
  call model_str(modelfile, model_path, model_name, model_iteration, slind, dind)
    HMMSTR1 = read_hmm(modelfile, .true.) ! READ THE MODEL
    call get_intrans(intrans,HMMSTR1)   ! get indeces of nodes transferring into each node
    call get_outtrans(outtrans,HMMSTR1) ! get indeces of nodes transferred to from each node
    
    if ( index(model_name, 'TM') /= 0) then
        tm = .true.
    else
        tm = .false.
    endif
    last_g_record = 0
     write(mitrm1,'(I3)') model_iteration 

    !write(*,*) chain_id(:4)
    this_profile = get_profile(drctfile, 0)
    sitr = 1
    
    do while (this_profile%last_record /= -1 ) 
           
            ! we skip if chain_id is specified and this code != chain_id 
            if ( chain_id .ne. "NONE" .and. this_profile%code .ne. chain_id ) then
                sitr = sitr + 1
                this_profile = get_profile(drctfile, this_profile%last_record)
                cycle
           ! else if ( chain_id .eq. "NONE" .and. mod(sitr,100) == 0 ) then
           !        write(*,*) "COMPLETED ", sitr, " CHAINS" 
            endif
            
            t = this_profile%nres - 1
            nres = t
            allocate(aaprofile(m,t)) ; aaprofile = 0
            allocate(angles(3,t)) ; angles=0.      !  (Phi, Psi, Omega)
            allocate(seq1(t))
            allocate(tmseq1(t))
            allocate(ssseq1(t))
            allocate(ramaseq1(t))
            allocate(tmseq(t))
            allocate(ssseq(t))
            allocate(ramaseq(t))
            allocate(aa_prof_c(t, m))
            allocate(ss_prof_c(t, mdssp))
            allocate(tm_prof_c(t, mtm))
            allocate(rama_prof_c(t, mrama))

            ! set profile for 1 .. nres
            do titr=1,t
                if (all(this_profile%profiles(:,titr) == 0)) then
                    if (this_profile%seq_code(titr) == 0 ) then  ! if there is an empty profile, empty seq_code... set to b_ground
                        aaprofile(:,titr) = HMMSTR%AA_ground
                    else
                        aaprofile(this_profile%seq_code(titr),titr) = 1.00
                    endif
                else
                    aaprofile(:,titr) = this_profile%profiles(:,titr)
                endif
            enddo
   
            angles = this_profile%dihs(:,:t)
            call set_flags(ramaseq,t,angles)
            do itr=1,t
                    ssseq(itr) = this_profile%ss_seq(itr)
                    tmseq(itr) = this_profile%tm_seq(itr)
            enddo
            
          !!!!!!!!!!!!GAMMA VOTING PROCEDURE, UPDATED 11.8.21 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          allocate(ct(nres)); ct = 0
          allocate(alpha(HMMSTR1%n,nres)); alpha = 0
          allocate(beta(HMMSTR1%n,nres));  beta = 0
          allocate(tgamma(nres,HMMSTR1%n)); tgamma = 0
          call zeta_init(zeta,zeta_dict,HMMSTR1,nres)
          !allocate(zeta(nres,HMMSTR1%n,HMMSTR1%n)); zeta = 0
          
          if ( model_iteration == 1 ) then
              gfile = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitrm1//"GAMMA.drct"
              gdl = dind + 9
          else
               gfile = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitrm1//"GAMMA.drct"
              gdl = dind + 5
          endif
          inquire(file=trim(gfile(:gdl+4)), exist = gf_e)
          gf_e = .false. ! this line lets us evaluate on AA + RAMA input
          if ( gf_e .eqv. .false. ) then
                  ! TLB 5/28/21 checking getalpha/getbeta...
                  ! calculate alpha, beta, gamma
                  call getalpha(ct,alpha,nres,HMMSTR1,aaprofile,intrans,outtrans,ramaseq,ssseq,tmseq,training_flags,1,nres) ! 6/22 Added training flags
                  call getbeta(ct,beta,nres,HMMSTR1,aaprofile,intrans,outtrans,ramaseq,ssseq,tmseq,training_flags,1,nres)
                  call get_gamma_zeta(ct, beta, alpha, nres, HMMSTR1, aaprofile, ramaseq, ssseq, tmseq, &
                                 outtrans, training_flags, tgamma, zeta, zeta_dict) ! 6/22 Added flags
          else
                  call get_gamma_profile(gfile,last_g_record,tgamma,HMMSTR1%n,this_profile%code)
          endif 
          ! TLB 10/12 moved gamma_voting logic to main loop...
          do titr = 1, t
             aa_prof = 0
             ss_prof = 0
             tm_prof = 0
             rama_prof = 0
             do nitr = HMMSTR1%k+1, HMMSTR1%n
                  aa_prof    = aa_prof  + tgamma(titr,nitr)*HMMSTR1%AA(nitr,:)
                  ss_prof    = ss_prof  + tgamma(titr,nitr)*HMMSTR1%SS(nitr,:)
                  tm_prof    = tm_prof  + tgamma(titr,nitr)*HMMSTR1%TM(nitr,:)
                  rama_prof  = rama_prof  + tgamma(titr,nitr)*HMMSTR1%RAMA(nitr,:)
             enddo
             aa_prof_c(titr,:) = aa_prof(:)
             ss_prof_c(titr,:) = ss_prof(:)
             tm_prof_c(titr,:) = tm_prof(:)
             rama_prof_c(titr,:) = rama_prof(:)
             call voting(HMMSTR1, seq1, tmseq1, ssseq1, ramaseq1, aa_prof, tm_prof, ss_prof, rama_prof, titr)
         enddo
     !     call viterbi(nres, HMMSTR1, aaprofile, intrans, outtrans, seq1, ramaseq, ssseq, tmseq, training_flags, paradigm) 
     !     stop
          !!! TM ROC LOGIC !!!!!!!!!!!!!!!!!!!
          do titr = 1, t-1
                  sum_TMgamma = 0
                  do nitr = HMMSTR1%k+1, HMMSTR1%n
                      tm_reference_ind = maxloc(HMMSTR1%TM(nitr,:),1)
                      tm_c = reference_tm_string(tm_reference_ind:tm_reference_ind)
          !write(*,*) "HERE", tm_reference_ind, tm_c
                      if (any(["H", "C", "I", "L", "F", "B"] == tm_c)) then
                          sum_TMgamma = sum_TMgamma + tgamma(titr, nitr)
                      endif
                  enddo
                  TMness = (sum_TMgamma / (sum(tgamma(titr,:)))) 
                  if ( any(["H", "C", "I", "L", "F", "B"] == tmseq(titr) ) ) then 
                      TM_bool = 1
                  else
                      TM_bool = 0
                  endif 
                  if ((TMness > .5 .and. TM_bool==1) .or. (TMness <= 0.5 .and. TM_bool==0)) then
                       TM_correct = 1           
                  else 
                       TM_correct = 0
                  endif
                  ! PRINTOUT STATEMENT !!!!!!!!!!!!!!!
                  !write(*,*) TMness, TM_correct, TM_bool, tmseq(titr)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          enddo
          ! get ocnfusion matrices
          call get_Q3_confusion(this_profile%ss_seq, ssseq1, t, Q3_confusion)
          call get_RAMA_confusion(this_profile%rama_seq, ramaseq1, t, RAMA_confusion)
          call get_TM_confusion(this_profile%tm_seq, tmseq1, t, TM_confusion)
          ! Deallocate ct, alpha, beta
          if(allocated(ct))               deallocate(ct)
          if(allocated(alpha))            deallocate(alpha)
          if(allocated(zeta))             deallocate(zeta)
          if(allocated(zeta_dict))             deallocate(zeta_dict)
          if(allocated(beta))             deallocate(beta)
          ! Deallocate seq, profile, angles, ramaseq, iflag
          if(allocated(aaprofile))          deallocate(aaprofile)
          if(allocated(angles))           deallocate(angles)
          !if(allocated(aaseq))          deallocate(aaseq)

        
        allocate(ground_truth_profiles(4,t))
        ground_truth_profiles(1,:) = this_profile%seq(:t)
        ground_truth_profiles(2,:) = ssseq
        ground_truth_profiles(3,:) = ramaseq
        ground_truth_profiles(4,:) = tmseq
        allocate(paradigm(4,t))
        paradigm(1,:) = seq1
        paradigm(2,:) = ssseq1
        paradigm(3,:) = ramaseq1
        paradigm(4,:) = tmseq1
        mda_wa = 0
        !write(*,*) " GOT HERE"
        !mda(profile, paradigm, mda, t, t_mda)

        call get_mda(ground_truth_profiles, paradigm, mda_wa, t, mda)
        write(*,*) " GOT MDA = " , mda_wa
        
         if ( run_gamma .eq. 'T') then
                write(*,*) "START"
                write(*,*) this_profile%code
                write(*,*) nres
                write(*,*) "1" ! profile weight ( to be implemented ) 
                write(*,*) HMMSTR1%n
                 do itr = 1, t-1
                     write(*,*) tgamma(itr, : )
                 enddo
         else ! normal comprehensive output
                ! if we only ran one chain then save gamma matrix to gamma_dir ( './gamma/' ) 
                !if ( chain_id .ne. "NONE") then
                !     call write_gamma(tgamma, './gamma/'//this_profile%code, HMMSTR1%n, nres)
                !endif
                sitr = sitr + 1
                !write(*,*) sitr
                if ( mod(sitr, 10) == 0 ) then ! do someting every 10 chains..
                endif
                !!!!!!!!!!!!!!!!!!!!!!!!
                write(*,*) "START"
                write(*,*) this_profile%code
                write(*,*) nres
                write(*,*) mda_wa(2)
                write(*,*) mda
                !write(*,*) Q3_cats
                do itr = 1, 4
                    !write(*,*) Q3_cats(itr),Q3_confusion(itr,:)
                    write(*,*) Q3_confusion(itr,:)
                enddo
                !write(*,*) RAMA_cats
                do itr = 1, 6
                    !write(*,*) RAMA_cats(itr),RAMA_confusion(itr,:)
                    write(*,*) RAMA_confusion(itr,:)
                enddo
                !write(*,*) TM_cats
                do itr = 1, 8
                    !write(*,*) TM_cats(itr),TM_confusion(itr,:)
                    write(*,*) TM_confusion(itr,:)
                enddo
                ! GAMMA-weighted RAMA PROFILE
                do titr = 1, nres
                    write(*,*) rama_prof_c(titr,:) 
                enddo
                ! GAMMA-weighted TM PROFILE
                do titr = 1, nres
                    write(*,*) tm_prof_c(titr,:) 
                enddo
                write(*,*) this_profile%seq
                write(*,*) ramaseq
                write(*,*) ramaseq1
                write(*,*) tmseq
                write(*,*) tmseq1
        endif ! end printout section 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LOOP ITR!!!!!!!!!!!!!!!!!!!!!!        
        if(allocated(aa_prof_c))  deallocate(aa_prof_c)
        if(allocated(ss_prof_c))  deallocate(ss_prof_c)
        if(allocated(tm_prof_c))  deallocate(tm_prof_c)
        if(allocated(rama_prof_c))  deallocate(rama_prof_c)
        if(allocated(ramaseq))          deallocate(ramaseq)
        if(allocated(tmseq))          deallocate(tmseq)
        if(allocated(ramaseq1))          deallocate(ramaseq1)
        if(allocated(tmseq1))          deallocate(tmseq1)
        if(allocated(ssseq))          deallocate(ssseq)
        if(allocated(seq1))          deallocate(seq1)
        if(allocated(ssseq1))          deallocate(ssseq1)
        if(allocated(tgamma))             deallocate(tgamma)
        if(allocated(paradigm))   deallocate(paradigm)
        if(allocated(ground_truth_profiles))          deallocate(ground_truth_profiles)
        Q3_confusion = 0
        RAMA_confusion = 0
        TM_confusion = 0
        this_profile = get_profile(drctfile, this_profile%last_record)
        !stop
    enddo

    contains

subroutine get_TM_confusion(gt, prediction, t, TM_confusion) 
    implicit none
    character(len=*),dimension(:) :: gt
    character(len=*),dimension(:) :: prediction
    integer,intent(in) :: t
    integer,dimension(8,8),intent(inout) :: TM_confusion
    integer :: i , j
    !Q3_confusion = 0 

    do titr = 1, t
        if ( any(["1", "2"] == prediction(titr) ) ) then
                j = 1
        else if ( "B" == prediction(titr) )then
                j = 2
        else if ( "H" == prediction(titr) )then
                j = 3
        else if ( "C" == prediction(titr) )then
                j = 4
        else if ( "I" == prediction(titr) )then
                j = 5
        else if ( "L" == prediction(titr) )then
                j = 6
        else if ( "F" == prediction(titr) )then
                j = 7
        else
                j = 8 ! U
        endif
        if ( any(["1", "2"] == gt(titr) ) ) then
                i = 1
        else if ( "B" == gt(titr) ) then
                i = 2
        else if ( "H" == gt(titr) ) then
                i = 3
        else if ( "C" == gt(titr) ) then
                i = 4
        else if ( "I" == gt(titr) ) then
                i = 5
        else if ( "L" == gt(titr) ) then
                i = 6
        else if ( "F" == gt(titr) ) then
                i = 7
        else
                i = 8 ! U
        endif
        TM_confusion(i,j) = TM_confusion(i,j) + 1
    enddo
end subroutine get_TM_confusion


subroutine get_Q3_confusion(gt, prediction, t, Q3_confusion) 
    implicit none
    character(len=*),dimension(:) :: gt
    character(len=*),dimension(:) :: prediction
    integer,intent(in) ::  t
    integer,dimension(4,4),intent(inout) :: Q3_confusion
    integer :: i , j
    !Q3_confusion = 0

    do titr = 1, t
        if ( any(["G","H"] == prediction(titr) ) )then
                j = 1
        else if ( any(["S", "T"] == prediction(titr) ) ) then
                j = 2
        else if ( prediction(titr) == "E" ) then
                j = 3
        else if ( prediction(titr) == "-") then
                j = 2
        else 
                cycle
        endif
        if ( any(["G","H"] == gt(titr) )) then
                i = 1
        else if ( any(["S", "T"] == gt(titr) ) ) then
                i = 2
        else if ( gt(titr) == "E" ) then
                i = 3
        else if ( gt(titr) == "-" ) then
                i = 2
        else 
                cycle
        endif
        Q3_confusion(i,j) = Q3_confusion(i,j) + 1
    enddo
end subroutine get_Q3_confusion

subroutine get_RAMA_confusion(gt, prediction, t, RAMA_confusion) 
    implicit none
    character(len=*),dimension(:) :: gt
    character(len=*),dimension(:) :: prediction
    integer,intent(in) :: t
    integer,dimension(6,6) :: RAMA_confusion
    integer :: i , j
    !RAMA_confusion = 0

    do titr = 1, t
       if ( any(["G","H"]==prediction(titr)) ) then
               j = 1
       else if ( any(["B","b"]==prediction(titr))) then
               j = 2
       else if ( any(["E","e","d"]==prediction(titr))) then
               j = 3
       else if ( any(["L","l"]==prediction(titr))) then
               j = 4
       else if (prediction(titr) == "c") then
               j = 6
       else if ( prediction(titr) == "x") then
               j = 5
       else
               cycle
       endif
       if ( any(["G","H"]==gt(titr)) ) then
               i = 1
       else if ( any(["B","b"]==gt(titr))) then
               i = 2
       else if ( any(["E","e","d"]==gt(titr))) then
               i = 3
       else if ( any(["L","l"]==gt(titr))) then
               i = 4
       else if (gt(titr) == "c") then
               i = 6
       else if ( gt(titr) == "x") then
               i = 5
       else
               cycle
       endif
        RAMA_confusion(i,j) = RAMA_confusion(i,j) + 1
    enddo
end subroutine get_RAMA_confusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_substrings(aline) result(args)

    character(len=50000),intent(in) :: aline
    character(len=30),dimension(3000) :: args
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
end function get_substrings

end program comprehensive
