module hmm
        use hmm_io!,  only : read_hmm, write_hmm
        use profile_io!, only : get_profile
        use gamma_io
        implicit none
        type(HMM_) :: HMMSTR
!        type(Profile) :: this_profile
        !real,dimension(:,:),allocatable :: profile
        !real,parameter::AA_EPS  = 0.0001 ! learning rate parameter. 
        !real,parameter::RAMA_EPS= 0.000001 ! learning rate parameter. 
        !real,parameter::TM_EPS  = 0.000000001 ! learning rate parameter. 
        real,parameter :: ncount = 0.0001
        type trans
          integer,dimension(:),allocatable::c
        end type
        integer::c=0
        integer,parameter::m=20
        integer,parameter::mdssp=6
        integer,parameter::mtm=10
        integer,parameter::mrama=11
        character(len=m),parameter     :: aa="ACDEFGHIKLMNPQRSTVWY"
        character(len=mrama),parameter :: rama_string="HGBEdbeLlxc"
        character(len=mdssp),parameter :: ss_string="HEGST-"
        character(len=mtm),parameter   :: tm_string="12BHCILFU_"
        type(trans),dimension(:),allocatable::intrans,outtrans
        real(8),dimension(:),allocatable :: memobt ! memoization objective function values per same t across nodes 1..n
! NEW HMM_model object:
!     type HMM
!        character(len=250) :: model_file
!        real,allocatable :: priors (:)
!        character(1),dimension(:),allocatable :: AA_seq, SS_seq, TM_seq, RAMA_seq
!        real,dimension(:,:),allocatable :: AA
!        real,dimension(:,:),allocatable :: SS
!        real,dimension(:,:),allocatable :: TM
!        real,dimension(:,:),allocatable :: RAMA
!        real,dimension(:,:),allocatable :: dihs
!        integer :: n,k,u,nt
!        character(len=4),dimension(20) :: emissions
!        character(len=30),dimension(20) :: emission_alphabet
!        real,dimension(:),allocatable :: AA_ground, SS_ground, TM_ground, RAMA_ground !
!        real,allocatable :: a(:,:) ! TRANSITIONS
!        real,allocatable :: ties(:,:) ! TRANSITIONS
!    endtype

     !hmmfile = "./HMMSTR.hmm"
     !HMMSTR = read_hmm('./HMMS/'//hmmfile)

contains

! write .seq output which holds predicted sequences from viterbi and gamma-weight 3rd column is CI --> B factor .pdb
subroutine write_seq(HMMSTR, tgamma, nres, code, ground_truth_profiles, paradigm)
  implicit none!write_seq(HMMSTR, tgamma, nres, code, ground_truth_variables, paradigm)
  real(8),dimension(:,:), intent(in) :: tgamma
  character(len=1),dimension(:,:),allocatable,intent(in) :: paradigm, ground_truth_profiles
  type(HMM_),intent(in) :: HMMSTR
  integer,intent(in) :: nres
  character(len=4),intent(in) :: code
  integer :: titr, nitr
  character(len=20) format_string
  character(len=4) szstr
  character(len=100) :: model_path, model_name
  character(len=250) :: modelfile, seqfile
  character(len=3) :: mitr, mitrm1
  integer :: model_iteration, dind, slind, gdl, ios
  !  gamma-voted sequences
  character(len=1),dimension(:),allocatable :: seq1
  character(len=1),dimension(:),allocatable :: ssseq1
  character(len=1),dimension(:),allocatable :: ramaseq1
  character(len=1),dimension(:),allocatable :: tmseq1
  real(8),dimension(2) :: mda
  integer,dimension(:),allocatable :: tmda
  
  modelfile = HMMSTR%model_file
  slind = index(modelfile, "/", .true.)
  dind = index(modelfile, ".", .true.)
  call model_str(modelfile, model_path, model_name, model_iteration, slind, dind)
  write(mitrm1,'(I3)') model_iteration
  seqfile = model_path(1:slind)//model_name(:dind-slind-1)//code//".seq"

  gdl = dind + 7 
  write(*,*) ".seq file = ", seqfile(:gdl)
  open(3, file=seqfile(:gdl), iostat = ios)
     if (ios /= 0) then
        write(*,*) "Error 2. Improperly formatted output file."
        stop
   endif

   write(3, '(A4)') code
   write(szstr,'(I4)') nres
   format_string = '(' // szstr // 'A)'
   write(3,'(A12)') "GROUND TRUTH"
   write(3,format_string,iostat=ios) ground_truth_profiles(1,:nres)
   write(3,format_string,iostat=ios) ground_truth_profiles(2,:nres)
   write(3,format_string,iostat=ios) ground_truth_profiles(3,:nres)
   write(3,format_string,iostat=ios) ground_truth_profiles(4,:nres)
   write(3,'(A16)') "VITERBI PARADIGM"
   write(3,format_string,iostat=ios) paradigm(1,:)
   write(3,format_string,iostat=ios) paradigm(2,:)
   write(3,format_string,iostat=ios) paradigm(3,:)
   write(3,format_string,iostat=ios) paradigm(4,:)
   allocate(seq1(nres))
   allocate(ssseq1(nres))
   allocate(tmseq1(nres))
   allocate(ramaseq1(nres))
   call get_mda(ground_truth_profiles, paradigm, mda, nres, tmda)
   call gamma_voting(nres, HMMSTR, tgamma, seq1, ssseq1, ramaseq1, tmseq1)

   write(3,'(A11)') "GAMMA VOTED"
   do titr=1, this_profile%nres-1
      format_string = '(I4,4A4,I3,3F11.5)' 
      write(3, format_string,iostat=ios) titr, seq1(titr), ssseq1(titr), ramaseq1(titr), tmseq1(titr), &
              tmda(titr), this_profile%dihs(:,titr)
   enddo
   deallocate(seq1)
   deallocate(ssseq1)
   deallocate(ramaseq1)
   deallocate(tmseq1)
   deallocate(tmda)
  close(3)
end subroutine write_seq

subroutine scratch_hmm(modelfile, drctfile, code)
   implicit none
  character(len=*),intent(inout) :: modelfile
  character(len=*) :: drctfile
  character(len=4) :: code
  integer,dimension(4) :: training_flags
  real,dimension(6) :: eval_flags
  real(8),dimension(:, :),allocatable :: tgamma ! (nres, HMMSTR%n)
  real(8),dimension(:,:),allocatable :: zeta
  integer,dimension(:,:),allocatable :: zeta_dict
  integer :: sitr
  training_flags = [1, 0, 0, 0]
                   !AA RAMA SS TM
  ! [ mda, Viterbi, pt, gamma, seq]
  eval_flags = [ 0., 0., 0., 0., 1., 0.]
  HMMSTR = read_hmm(modelfile, .true.) ! READ THE MODEL
    call get_intrans(intrans,HMMSTR)   ! get indeces of nodes transferring into each node
    call get_outtrans(outtrans,HMMSTR) ! get indeces of nodes transferred to from each node
  if ( index(drctfile, '.profile') /= 0 ) then ! .profile file as input

          this_profile = read_profile(drctfile)
          ! ---------------------------------------------------------------------------
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call run_hmm(HMMSTR, tgamma, zeta, zeta_dict, training_flags, eval_flags)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ---------------------------------------------------------------------------
           
          if(allocated(zeta))             deallocate(zeta)
          if(allocated(tgamma))             deallocate(tgamma)

  else if ( index(drctfile, '.drct') /= 0 ) then
          this_profile = get_profile(drctfile, 0)
          do while(this_profile%nres /= -1)
                  if ( this_profile%code /= code  ) then  ! debug loop cycler
                          sitr = sitr + 1
                          this_profile = get_profile(drctfile, this_profile%last_record)
                          cycle
                  endif

                  ! ---------------------------------------------------------------------------
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  call run_hmm(HMMSTR, tgamma, zeta, zeta_dict, training_flags, eval_flags)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! ---------------------------------------------------------------------------
                  
                  if(allocated(zeta))             deallocate(zeta)
                  if(allocated(tgamma))             deallocate(tgamma)
                  
                  sitr = sitr + 1
                  this_profile = get_profile(drctfile, this_profile%last_record)
          enddo
   endif

end subroutine scratch_hmm

! TLB 1/9/22 exhaust_hmm runs forward backward algorithm at all possible start and end timepoints. 
!  This is meant to 'fix' possible off-by-X errors for prediction purposes.
subroutine exhaust_hmm(modelfile, drctfile, code)
   implicit none
  character(len=*),intent(inout) :: modelfile
  character(len=*) :: drctfile
  character(len=4) :: code
  integer,dimension(4) :: training_flags
  real,dimension(6) :: eval_flags
  real(8),dimension(:, :),allocatable :: tgamma ! (nres, HMMSTR%n)
  real(8),dimension(:,:),allocatable :: zeta
  integer,dimension(:,:),allocatable :: zeta_dict
  integer :: sitr
  training_flags = [1, 0, 0, 0]
                   !AA RAMA SS TM
  ! [ MDA, Viterbi, pt, gamma, .seq, exhaust]
  eval_flags = [ 0., 0., 0., 0., 1., 1.]
  HMMSTR = read_hmm(modelfile, .false.) ! READ THE MODEL
    call get_intrans(intrans,HMMSTR)   ! get indeces of nodes transferring into each node
    call get_outtrans(outtrans,HMMSTR) ! get indeces of nodes transferred to from each node

  this_profile = get_profile(drctfile, 0)
  do while(this_profile%nres /= -1)
          if ( this_profile%code /= code  ) then  ! debug loop cycler
                  sitr = sitr + 1
                  this_profile = get_profile(drctfile, this_profile%last_record)
                  cycle
          endif

          ! ---------------------------------------------------------------------------
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call run_hmm(HMMSTR, tgamma, zeta, zeta_dict, training_flags, eval_flags)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ---------------------------------------------------------------------------

          if(allocated(zeta))             deallocate(zeta)
          if(allocated(tgamma))             deallocate(tgamma)

          sitr = sitr + 1
          this_profile = get_profile(drctfile, this_profile%last_record)
  enddo
end subroutine exhaust_hmm

! TLB 12/14/21 moved bulk of train_hmmstr's inner-loop into run_hmm
subroutine run_hmm(HMMSTR, tgamma, zeta, zeta_dict, alg_flags, eval_flags)

  ! tgamma = probability state j at time i
  ! zeta = probability state zeta_dict(i) goes to state zeta_dict(j) at time t
  !      = zeta(zeta_dict(i,j),t)
  ! mp = model probability = - log sum ct
  real(8),dimension(:, :),allocatable,intent(inout) :: tgamma ! (nres, HMMSTR%n)
  real(8),dimension(:,:),allocatable,intent(inout) :: zeta
  integer,dimension(:,:),allocatable,intent(inout) :: zeta_dict
  character(len=1),dimension(:,:),allocatable :: ground_truth_profiles
  real :: mp ! model probability is the output
  type(HMM_),intent(in) :: HMMSTR
  integer :: nres
  ! evaluation / objective function criteria
  ! [ AA, RAMA, SS, TM ] 
  integer,dimension(4),intent(in) :: alg_flags
  ! [ MDA, Viterbi, pt, gamma, seq, exhaust]
  !  if pt == 1, set eval_flags[3] = mp 
  real,dimension(6),intent(inout) :: eval_flags
  ! mda stuff ! --------------------------------------- 9/8/21
  character(len=1),dimension(:,:),allocatable :: paradigm
  real(8),dimension(2) :: mda_wa
  logical :: mda_acc
  ! profile sequences 
  integer,dimension(:),allocatable :: seq
  real,dimension(:,:),allocatable::profile
  real,dimension(:,:),allocatable::angles
  character(len=1),dimension(:),allocatable::aaseq
  character(len=1),dimension(:),allocatable::ramaseq
  character(len=1),dimension(:),allocatable::ssseq
  character(len=1),dimension(:),allocatable::tmseq
  character(len=1),dimension(:),allocatable::tmseq_unlabel
  ! Gamma-voting prediction arrays
  character(len=1),dimension(:),allocatable :: seq1
  character(len=1),dimension(:),allocatable :: ssseq1
  character(len=1),dimension(:),allocatable :: ramaseq1
  character(len=1),dimension(:),allocatable :: tmseq1
  ! Forward-Backward pass variables
  real(8),dimension(:),allocatable::ct
  real(8),dimension(:,:),allocatable::alpha
  real(8),dimension(:,:),allocatable::beta
  ! Exhaust Forward-Backward
  real(8),dimension(:,:),allocatable::alpha_sum
  real(8),dimension(:,:),allocatable::beta_sum
  ! counters, temporaries, misc
  integer :: titr, titr2, titr3, nitr
  real :: treal
  character(len=100) :: model_path, model_name
  character(len=250) :: modelfile, seqfile
  character(len=3) :: mitr, mitrm1
  integer :: model_iteration, dind, slind, gdl, ios
  integer,dimension(:),allocatable :: tmda
  logical :: unlabel_bool

  unlabel_bool = .false. ! TMHMM unlableing scheme 1/3/22
  modelfile = HMMSTR%model_file 
  slind = index(modelfile, "/", .true.)
  dind = index(modelfile, ".", .true.)
  call model_str(modelfile, model_path, model_name, model_iteration, slind, dind)
  !write(*,*) "RUNNING ", modelfile
  !write(*,*) this_profile%code
  ! set sequence (numbered aa)
  if ( allocated(profile) )  deallocate(profile)
  if ( allocated(angles) )  deallocate(angles)
  if ( allocated(seq) )  deallocate(seq)
  if ( allocated(aaseq) )  deallocate(aaseq)
  if ( allocated(ramaseq) )  deallocate(ramaseq)
  if ( allocated(tmseq) )  deallocate(tmseq)
  if ( allocated(seq1) )  deallocate(seq1)
  if ( allocated(ssseq1) )  deallocate(ssseq1)
  if ( allocated(ramaseq1) )  deallocate(ramaseq1)
  if ( allocated(tmseq1) )  deallocate(tmseq1)
  if ( allocated(tmseq_unlabel) )  deallocate(tmseq_unlabel)
  if ( allocated(ssseq) )  deallocate(ssseq)
  nres = this_profile%nres-1
  !write(*,*) nres
  allocate(seq(nres)) ; seq = 0
  seq = this_profile%seq_code(:nres)
  ! allocate seq, profile, angles, ramaseq, iflag
  allocate(profile(m,nres)) ; profile = 0
  allocate(angles(3,nres)) ; angles=0.      !  (Phi, Psi, Omega)
  allocate(aaseq(nres))
  allocate(ramaseq(nres))
  allocate(tmseq(nres))
  allocate(tmseq_unlabel(nres))
  allocate(ssseq(nres))
  allocate(seq1(nres))
  allocate(ssseq1(nres))
  allocate(ramaseq1(nres))
  allocate(tmseq1(nres)) 
  !write(*,*) "DONE ALLOCATING."
  ! --------------------------------------------------------------------------
  ! gather information from this_profile
  ! set profile for 1 .. nres
  do titr=1,nres
      if (all(this_profile%profiles(:,titr) == 0)) then
          if (seq(titr) == 0 ) then  ! if there is an empty profile, empty seq_code... set to b_ground
              profile(:,titr) = HMMSTR%AA_ground
          else
              profile(seq(titr),titr) = 1.00
          endif
      else
          profile(:,titr) = this_profile%profiles(:,titr)
      endif
  enddo

  ! set angles
  angles = this_profile%dihs(:,:nres)

  ! using angles, set ramaseq
  call set_flags(ramaseq,nres,angles)
  aaseq(:nres) = this_profile%seq(:nres)
  ssseq(:nres) = this_profile%ss_seq(:nres)
  tmseq(:nres) = this_profile%tm_seq(:nres)
  !write(*,*) tmseq(:nres)
  !call unlabel(tmseq, tmseq_unlabel, nres)
  !write(*,*) tmseq_unlabel(:nres)
  allocate(ground_truth_profiles(4,nres))
  ground_truth_profiles(1,:) = aaseq(:nres)
  ground_truth_profiles(2,:) = ssseq(:nres)
  ground_truth_profiles(3,:) = ramaseq(:nres)
  ground_truth_profiles(4,:) = tmseq(:nres)
  ! ---------------------------------------------------------------------------
  ! allocate ct, alpha, beta, gamma
  if ( allocated(ct) )  deallocate(ct)
  allocate(ct(nres)); ct = 0
  if ( allocated(alpha) )  deallocate(alpha)
  allocate(alpha(HMMSTR%n,nres)); alpha = 0
  if ( allocated(beta) )  deallocate(beta)
  allocate(beta(HMMSTR%n,nres));  beta = 0
  if ( allocated(alpha_sum) )  deallocate(alpha_sum)
  allocate(alpha_sum(HMMSTR%n,nres)); alpha_sum = 0
  if ( allocated(beta_sum) )  deallocate(beta_sum)
  allocate(beta_sum(HMMSTR%n,nres));  beta_sum = 0
  if ( allocated(tgamma) )  deallocate(tgamma)
  allocate(tgamma(nres,HMMSTR%n)); tgamma = 0
  call zeta_init(zeta,zeta_dict,HMMSTR,nres)
  ! TLB 5/28/21 checking getalpha/getbeta...
  ! 6/22 Added training flags
  !write(*,*) "Calculating Rabiner Variables ..."
  if ( unlabel_bool ) then ! 1/3/22 TMHMM TM unlabeling scheme, experimental
          call getalpha(ct,alpha,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq_unlabel,alg_flags,1,nres) 
          call getbeta(ct,beta,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq_unlabel,alg_flags,1,nres)
          call get_gamma_zeta(ct, beta, alpha, nres, HMMSTR, profile, ramaseq, ssseq, tmseq_unlabel, & 
                         outtrans, alg_flags, tgamma, zeta, zeta_dict) ! 6/22 Added flags ! 11/18 zeta_dict update
  else
         if (eval_flags(6) == 1. ) then ! exhaustive evaluation, experimental 
              do titr = 1 , nres - 1
                 do titr2 = titr+1, nres
                  call getalpha(ct,alpha,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq,alg_flags,titr,titr2) 
                  call getbeta(ct,beta,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq,alg_flags,titr,titr2)
                  do titr3 = titr, titr2
                      do nitr = 1, HMMSTR%n
                          alpha_sum(nitr, titr3) = alpha_sum(nitr, titr3) + alpha(nitr, titr3)
                          beta_sum(nitr, titr3) = beta_sum(nitr, titr3) + beta(nitr, titr3)
                      enddo
                  enddo
                 enddo
              enddo
              call getalpha(ct,alpha,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq,alg_flags,1,nres) 
              call get_gamma_zeta(ct, beta_sum, alpha_sum, nres, HMMSTR, profile, ramaseq, ssseq, tmseq, & 
                 outtrans, alg_flags, tgamma, zeta, zeta_dict) ! 6/22 Added flags ! 11/18 zeta_dict update
         else ! MAIN BRANCH ... Rabiner variables calculated here:
              call getalpha(ct,alpha,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq,alg_flags,1,nres) 
              call getbeta(ct,beta,nres,HMMSTR,profile,intrans,outtrans,ramaseq,ssseq,tmseq,alg_flags,1,nres)
              call get_gamma_zeta(ct, beta, alpha, nres, HMMSTR, profile, ramaseq, ssseq, tmseq, & 
                 outtrans, alg_flags, tgamma, zeta, zeta_dict) ! 6/22 Added flags ! 11/18 zeta_dict update
         endif

  endif
  ! ---------------------------------------------------------------------------
  if ( eval_flags(1) == 1. ) then ! Run Viterbi, get mda metric for given observation
           call viterbi(nres, HMMSTR, profile, intrans, outtrans, aaseq, ramaseq, ssseq, tmseq, alg_flags, paradigm, .false. )
          call get_mda(ground_truth_profiles, paradigm, mda_wa, nres, tmda)
          write(*,*) "VITERBI MDA ACCURACY = ", mda_wa(2)
   endif
   if ( eval_flags(2) == 1. ) then ! Run Viterbi
          call viterbi(nres, HMMSTR, profile, intrans, outtrans, aaseq, ramaseq, ssseq, tmseq, alg_flags, paradigm, .false.)
   endif
   if ( eval_flags(3) == 1. ) then ! gather model probability of the given observation
       treal = 1
       do titr=1, nres
          treal = treal + log(ct(titr))
       enddo
       treal = -treal 
       mp = treal
   endif
   if ( eval_flags(4) == 1. ) then ! Write gamma matrix to output drct file
       call write_gamma(tgamma, modelfile, HMMSTR%n, nres, this_profile%code, this_profile%chain)
   endif
   if ( eval_flags(5) == 1. ) then ! Gather and print Viterbi, Gamma-voted prediction sequences, and get mda 
           call viterbi(nres, HMMSTR, profile, intrans, outtrans, aaseq, ramaseq, ssseq, tmseq, alg_flags, paradigm, .false. )
          call get_mda(ground_truth_profiles, paradigm, mda_wa, nres, tmda)
          write(*,*) "VITERBI MDA ACCURACY = ", mda_wa(2)
           call gamma_voting(nres, HMMSTR, tgamma, seq1, ssseq1, ramaseq1, tmseq1)
           write(*,*) "GAMMA VOTED"
           write(*,*) seq1(:)
           write(*,*) ssseq1(:)
           write(*,*) ramaseq1(:)
           write(*,*) tmseq1(:)
             paradigm(1,:) = seq1(:nres) ! change paradigm to fit gamma-voted profiles.
  paradigm(2,:) = ssseq1(:nres)
  paradigm(3,:) = ramaseq1(:nres)
  paradigm(4,:) = tmseq1(:nres)
          call write_seq(HMMSTR, tgamma, nres, this_profile%code, ground_truth_profiles, paradigm)
          call get_mda(ground_truth_profiles, paradigm, mda_wa, nres, tmda)
          write(*,*) "GAMMA VOTED MDA ACCURACY = ", mda_wa(2)
   endif 
  !write(*,*) " Resolving Ties." ! TLB 12/7
  call untie(HMMSTR, tgamma, zeta, zeta_dict, outtrans, nres) 
                                              ! if tie states (i,j) exist in B or AB ties: 
                                              ! for (tiitr, ties i --> j) , update j gamma += sum (i gamma for i --> j )
                                              ! for i in ( i --> j)
                                              !      i gamma = j gamma
                                                  ! if tie states(i,j) exist in A or AB ties:
                                              ! for i,i2 in Atie(tiitr):
                                              !    for j in i --> j
                                              !       if j in Atie(tiitr[1,:])
                                              !         j2 = Atie(tiitr[2,Atie(tiitr[1,:]).index(j)])
                                              !         for titr in 1,T 
                                              !              t = sum(zeta(i,j,titr),zeta(i2,j2,titr))
                                              !              zeta(i,j,titr) = t
                                              !              zeta(i2,j2,titr) = t

  ! Deallocate seq, profile, angles, ramaseq, iflag
  if(allocated(seq))              deallocate(seq)
  if(allocated(profile))          deallocate(profile) ! issue with this line
  if(allocated(angles))           deallocate(angles)
  if(allocated(aaseq))          deallocate(aaseq)
  if(allocated(ramaseq))          deallocate(ramaseq)
  if(allocated(tmseq))          deallocate(tmseq)
  if ( allocated(seq1) )  deallocate(seq1)
  if ( allocated(ssseq1) )  deallocate(ssseq1)
  if ( allocated(ramaseq1) )  deallocate(ramaseq1)
  if ( allocated(tmseq1) )  deallocate(tmseq1)
  if(allocated(tmseq_unlabel))          deallocate(tmseq_unlabel)
  if(allocated(ssseq))          deallocate(ssseq)
  if(allocated(paradigm))   deallocate(paradigm)
  if(allocated(ground_truth_profiles))          deallocate(ground_truth_profiles)
  ! Deallocate ct, alpha, beta
  if(allocated(ct))               deallocate(ct)
  if(allocated(alpha))            deallocate(alpha)
  if(allocated(beta))             deallocate(beta)

end subroutine run_hmm
        

! TLB 5/25/20
subroutine train_hmm(modelfile, drctfile, epochs) ! epochs = # iterations through training set.
  ! ------------------------------------------------------------------------
  implicit none
  ! ------------------------------------------------------------------------
  character(len=*),intent(inout) :: modelfile
  character(len=*) :: drctfile
  integer :: epochs
  ! ------------------------------------------------------------------------
  integer :: titr, nitr, ei, sitr, dropped_sequences
  integer :: nres
  !---------------------------------------------
  !---------------------------------------------
  ! TLB 5/23/20
  real(8),dimension(:,:),allocatable::newA
  !real(8),dimension(:, :, :),allocatable :: zeta ! (nres, HMMSTR%n, HMMSTR%n) ! outdated
  real(8),dimension(:,:),allocatable::newB_residues
  real(8),dimension(:,:),allocatable::newB_rama
  real(8),dimension(:,:),allocatable::newB_ss
  real(8),dimension(:,:),allocatable::newB_tm
  real(8),dimension(:),allocatable::new_priors
  real,dimension(4775, 2) :: model_probabilities ! P(O|Model), weight
  character(len=150)::ofile,pfile,gfile
  integer :: ios ! check if HMMSTR was written out successfully
  !---------------------------------------------
  ! model file organizer ..! 1 10.12.21
  character(len=100) :: model_path, model_name
  character(len=3) :: mitr, mitrm1
  integer :: model_iteration, dind, slind, gdl
  logical :: dir_e
  integer :: statf
  !!!!!!!!!!!!!!!!!!HMMSTR_RUN_INPUT/OUTPUT!!!!!!!!!!!!!!!
  integer,dimension(4) :: training_flags
  real,dimension(6) :: eval_flags
  real(8),dimension(:, :),allocatable :: tgamma ! (nres, HMMSTR%n)
  real(8),dimension(:,:),allocatable :: zeta
  integer,dimension(:,:),allocatable :: zeta_dict
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! profile sequences
  integer,dimension(:),allocatable :: seq
  real,dimension(:,:),allocatable::profile
  real,dimension(:,:),allocatable::angles
  character(len=1),dimension(:),allocatable::aaseq
  character(len=1),dimension(:),allocatable::ramaseq
  character(len=1),dimension(:),allocatable::ssseq
  character(len=1),dimension(:),allocatable::tmseq
  real,dimension(:,:),allocatable :: mda_weight_accs
  real,dimension(:),allocatable :: mda
  logical :: mda_acc

  slind = index(modelfile, "/", .true.)
  dind = index(modelfile, ".", .true.)
  call model_str(modelfile, model_path, model_name, model_iteration, slind, dind)
  allocate(mda_weight_accs(5000,2))
  allocate(mda(2))
  mda_weight_accs = 0
  training_flags = [1, 1, 0, 1]
                   !AA RAMA SS TM
  ! [ MDA, Viterbi, pt, gamma, .seq]
  eval_flags = [ 0., 0., 0., 0., 0., 0.]  ! DON'T SAVE GAMMA MATRICES
  !eval_flags = [ 0., 0., 0., 1., 0., 0.]  ! SAVE GAMMA MATRICES
  model_probabilities = 0
  HMMSTR = read_hmm(modelfile, .false.) ! READ THE MODEL

  ! allocate newA, newB_residues, newB_rama, new_priors
  ! NOTE: arrays here allocated with HMMSTR%N, during write_new_model,
  ! first k indeces are ignored.... potential speedup possible by allocating to n - k
  allocate(newA(HMMSTR%n, HMMSTR%n))   ; newA = 0
  allocate(newB_residues(HMMSTR%n,m)) ; newB_residues = 0
  allocate(newB_rama(HMMSTR%n,mrama))     ; newB_rama = 0
  allocate(newB_ss(HMMSTR%n,mdssp))     ; newB_ss = 0
  allocate(newB_tm(HMMSTR%n,mtm))     ; newB_tm = 0
  allocate(new_priors(HMMSTR%n)) ; new_priors = 0
  allocate(memobt(HMMSTR%n)); memobt = 0
  do ei=1, epochs
    call get_intrans(intrans,HMMSTR)   ! get indeces of nodes transferring into each node
    call get_outtrans(outtrans,HMMSTR) ! get indeces of nodes transferred to from each node
    mda_acc = .false.
    if ( mod(ei, 10) == 0 ) then
            mda_acc = .true.
    endif
    ! DEBUGGING
    !mda_acc = .true.
    write(mitr,'(I3)') model_iteration + ei
    write(mitrm1,'(I3)') model_iteration + ei - 1

    this_profile = get_profile(drctfile, 0)
    sitr = 1
    dropped_sequences = 0
    model_probabilities = 0 

    do while(this_profile%nres /= -1)
          !write(*,"(a10,i4,a1)") "iteration ", sitr, "."
          if (mod(sitr,100) == 0) then
             write(*,"(a10, i4, a11)") "completed ", sitr, " sequences."
          endif
          
          !write(*,*) this_profile%code
          nres = this_profile%nres-1
          if ( nres < 20  ) then  ! debug loop cycler
          !if ( sitr > 1 ) then
                  sitr = sitr + 1
                  this_profile = get_profile(drctfile, this_profile%last_record)
                  cycle
          endif
          allocate(seq(nres)) ; seq = 0
          allocate(aaseq(nres))
          seq = this_profile%seq_code(:nres)
          ! allocate seq, profile, angles, ramaseq, iflag
          allocate(profile(m,nres)) ; profile = 0
          allocate(ramaseq(nres))
          allocate(angles(3,nres))
          allocate(tmseq(nres))
          allocate(ssseq(nres))

          angles = this_profile%dihs(:,:nres)

          call set_flags(ramaseq,nres,angles)
          aaseq(:nres) = this_profile%seq(:nres)
          ssseq(:nres) = this_profile%ss_seq(:nres)
          tmseq(:nres) = this_profile%tm_seq(:nres)

          do titr=1,nres
              if (all(this_profile%profiles(:,titr) == 0)) then
                  if (seq(titr) == 0 ) then  ! if there is an empty profile, empty seq_code... set to b_ground
                      profile(:,titr) = HMMSTR%AA_ground
                  else
                      profile(seq(titr),titr) = 1.00
                  endif
              else
                  profile(:,titr) = this_profile%profiles(:,titr)
              endif
          enddo
          
          ! ---------------------------------------------------------------------------
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call run_hmm(HMMSTR, tgamma, zeta, zeta_dict, training_flags, eval_flags)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ---------------------------------------------------------------------------

          ! update probabilities scalar
          model_probabilities(sitr, 1) = eval_flags(4)
          model_probabilities(sitr, 2) = nres

          ! Deallocate everything that was allocated ----------------------------------
          !    write(*,*) " Get A"
          call get_A( nres, HMMSTR, newA, outtrans, zeta, zeta_dict) ! 6/22 Added flags ! 10/4 Removed tgamma param
          !    write(*,*) " Get B"
          call get_B(tgamma, nres, HMMSTR, newB_residues, profile, newB_rama, ramaseq, newB_ss, ssseq, newB_tm, tmseq)
          !    write(*,*) " Get priors"
          call get_priors(tgamma, new_priors, HMMSTR%n)

          if(allocated(zeta))             deallocate(zeta)
          if(allocated(tgamma))             deallocate(tgamma)
          if (allocated(seq))   deallocate(seq)
          if (allocated(ramaseq))   deallocate(ramaseq)
          if (allocated(aaseq))   deallocate(aaseq)
          if (allocated(tmseq))   deallocate(tmseq)
          if (allocated(ssseq))   deallocate(ssseq)
          if (allocated(angles))   deallocate(angles)
          if (allocated(profile))   deallocate(profile)
           
          sitr = sitr + 1
          this_profile = get_profile(drctfile, this_profile%last_record)
          ! ---------------------------------------------------------------------------
      enddo
      if ( mda_acc ) then
           call get_mda_acc(mda_weight_accs)
           mda_weight_accs = 0 
           !if (allocated(mda_weight_accs)) deallocate(mda_weight_accs)
      endif
      if ( model_iteration == 1 ) then
          gfile = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitrm1//"GAMMA.drct"
          ofile = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitrm1//"GAMMA/"
          gdl = dind + 9
      else
          gfile = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitrm1//"GAMMA.drct"
          ofile = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitrm1//"GAMMA/"
          gdl = dind + 5
      endif
      write(*,*) gfile
      write(*,*) ofile
      write(*,*) "cat '"//trim(ofile(:gdl))//"'* > '"//trim(gfile(:gdl+4))//"'"
      write(*,*) "rm -r '"//trim(ofile(:gdl))//"'"
      call system("cat '"//trim(ofile(:gdl))//"'* > '"//trim(gfile(:gdl+4))//"'")
      call system("rm -r '"//trim(ofile(:gdl))//"'")

      write(*,*) "TOTAL # SEQUENCES =  ", sitr, "."
      write(*,*) "EPOCH ", ei, " COMPLETED."
      write(*,*) "Dropped ", dropped_sequences, " sequences."
      ! update paramaters in model
      ! SOME OF THESE ERRORS ARE GOING TO BE WRITTEN OUT BECAUSE OF UN-TRAVERSED STATES.
      do nitr=1,HMMSTR%n
          
          if (sum(newA(nitr,:)) == 0.0) then
         !   write(*,*) "OVERWRITING A/B", nitr 
              newA(nitr,:) = HMMSTR%a(nitr,:)
          endif
          if (sum(newB_residues(nitr,:)) == 0.0) then
         !     write(*,*) "OVERWRITING B : AA"
              newB_residues(nitr,:) = HMMSTR%AA(nitr,:)
          endif
          if (sum(newB_rama(nitr,:)) == 0.0) then
         !     write(*,*) "OVERWRITING B : RAMA"
              newB_rama(nitr,:) = HMMSTR%RAMA(nitr,:)
          endif
          if (sum(newB_ss(nitr,:)) == 0.0) then
         !     write(*,*) "OVERWRITING B : SS"
              newB_ss(nitr,:) = HMMSTR%SS(nitr,:)
          endif
          if (sum(newB_tm(nitr,:)) == 0.0) then
         !     write(*,*) "OVERWRITING B : TM"
              newB_tm(nitr,:) = HMMSTR%TM(nitr,:)
          endif
         ! if newPi == 0 then leave it as i
          if (new_priors(nitr) == 0 .or. new_priors(nitr) /= new_priors(nitr)) then
              new_priors(nitr) = HMMSTR%priors(nitr)
          endif
      enddo
      call normalize(newA, HMMSTR%n)
      call normalize(newB_residues, HMMSTR%n)
      call normalize(newB_rama, HMMSTR%n)
      call normalize(newB_ss, HMMSTR%n)
      call normalize(newB_tm, HMMSTR%n)
      call normalize1d(new_priors)
      HMMSTR%a(:HMMSTR%n,:HMMSTR%n)= newA(:HMMSTR%n, :HMMSTR%n)
      HMMSTR%AA(:HMMSTR%n,:)= newB_residues(:HMMSTR%n, :)
      HMMSTR%RAMA(:HMMSTR%n,:)  = newB_rama(:HMMSTR%n, :)
      HMMSTR%SS(:HMMSTR%n,:) = newB_ss(:HMMSTR%n, :)
      HMMSTR%TM(:HMMSTR%n,:) = newB_tm(:HMMSTR%n, :)

      write(*,*) "-----------------------------------------------------------------"
      HMMSTR%priors = new_priors(:)
      
      ! set model ofile and pfile
      write(mitr,'(I3)') model_iteration + ei
      write(mitrm1,'(I3)') model_iteration + ei - 1
      if ( model_iteration == 1 ) then
              ofile = model_path(1:slind)//model_name(:dind-slind-1)//"_"//mitr//".hmm"
              pfile = model_path(1:slind)//model_name(:dind-slind-1)//"_P_"//mitrm1//".txt"
      else
              ofile = model_path(1:slind)//model_name(:dind-slind-5)//"_"//mitr//".hmm"
              pfile = model_path(1:slind)//model_name(:dind-slind-5)//"_P_"//mitrm1//".txt"
      endif
      ofile = trim(ofile)
      pfile = trim(pfile)
      ! write out new data
      write(*,*) "writing model probabilities to :", pfile
      !call write_p(model_probabilities, sitr, pfile)
      write(*,*) "writing new model to : ", ofile
      ios = write_hmm(modelfile,ofile,HMMSTR)
      
      ! Deallocate intrans and outrans
      if(allocated(intrans))    deallocate(intrans)
      if(allocated(outtrans))   deallocate(outtrans)
      ! loop iter
      modelfile = ofile
      HMMSTR = read_hmm(modelfile, .false.) ! READ THE NEW MODEL
      newA = 0
      newB_residues = 0
      newB_rama = 0
      newB_tm = 0
      newB_ss = 0
      new_priors = 0 
  enddo

  !call release_model(HMMSTR)
  ! Deallocate memoization array
  if(allocated(memobt)) deallocate(memobt)
  ! Deallocate newA, newB_residues, newB_rama, new_priors
  if(allocated(newA))             deallocate(newA)
  if(allocated(newB_residues))    deallocate(newB_residues)
  if(allocated(newB_rama))        deallocate(newB_rama)
  if(allocated(newB_tm))        deallocate(newB_tm)
  if(allocated(newB_ss))        deallocate(newB_ss)
  if(allocated(new_priors))       deallocate(new_priors)

  modelfile = ofile

end subroutine train_hmm

  ! TMHMM TM unlabeling scheme, set 
  subroutine unlabel(tmseq, tmseq_unlabel, t)
       implicit none
       character(len=1),dimension(:),intent(in) :: tmseq
       character(len=1),dimension(:),intent(inout) :: tmseq_unlabel
       integer,intent(in) :: t
       integer :: titr, titr2
       integer,dimension(t) :: unlabel_mask
      write(*,*) "STARTED UNLABELING"
       unlabel_mask = 0 
       do titr=1, t
           if (tmseq(titr) /= "U" .and. tmseq(titr) /= "1" .and. tmseq(titr) /= "2" ) then
               if ( tmseq(titr-1) == "U" .or. tmseq(titr-1) == "1" .or. tmseq(titr-1) == "2" ) then
                   do titr2 = max(1,titr-3), min(t,titr+3)
                      unlabel_mask(titr2) = 1
                   enddo
               else if ( tmseq(titr+1) == "U" .or. tmseq(titr+1) == "1" .or. tmseq(titr+1) == "2" ) then
                   do titr2 = max(1, titr-3), min(t, titr+3)
                      unlabel_mask(titr2) = 1
                   enddo
               endif
           endif
       enddo

       do titr=1, t
          if (unlabel_mask(titr) == 1) then
              tmseq_unlabel(titr) = "U"
          else
              tmseq_unlabel(titr) = tmseq(titr)
          endif
       enddo
     write(*,*) "FINISHED UNLABELING"
  end subroutine unlabel

! Zeta functions added 11/18/21 to avoid large, unallocatable 3d array
  subroutine zeta_lookup(i,j,t,zeta,zeta_dict,zv)
     implicit none
     integer,intent(in) :: i,j,t
     real(8),intent(inout) :: zv
     real(8),dimension(:,:),intent(in) :: zeta
     integer,dimension(:,:),intent(in) :: zeta_dict
     integer :: zeta_ind 
     zeta_ind = zeta_dict(i,j)
     zv = zeta(zeta_ind,t)
  end subroutine zeta_lookup

  subroutine zeta_init(zeta,zeta_dict,HMMSTR,nres)
     implicit none
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: nres
     real(8),dimension(:,:),allocatable,intent(inout) :: zeta
     integer,dimension(:,:),allocatable,intent(inout) :: zeta_dict
     integer :: iitr, jitr, zitr, tcount
     
     if (allocated(zeta_dict)) deallocate(zeta_dict)
     allocate(zeta_dict(HMMSTR%n,HMMSTR%n))
     zeta_dict = 0

     zitr = 1
     tcount = 0 
     
     do iitr = 1, HMMSTR%n
        do jitr =1, HMMSTR%n
           if (HMMSTR%a(iitr,jitr) > 0 ) then
               tcount = tcount + 1
               zeta_dict(iitr,jitr) = tcount
           endif
        enddo
     enddo

     if (allocated(zeta)) deallocate(zeta)
     allocate(zeta(tcount,nres))
     zeta = 0
  end subroutine zeta_init

  subroutine zeta_set(i,j,t,zeta,zeta_dict,zv)
      implicit none
      integer,intent(in) :: i,j,t
      real(8),dimension(:,:),intent(inout) :: zeta
      integer,dimension(:,:),intent(in) :: zeta_dict
      real(8),intent(in) :: zv
      integer :: zeta_ind
      zeta_ind = zeta_dict(i,j)
      zeta(zeta_ind,t) = zv
  end subroutine zeta_set


!------------------------------------------------------------------------------------------------------------
  !    write(*,*) " Resolving Ties." ! TLB 12/7
  subroutine untie(HMMSTR, tgamma, zeta, zeta_dict, outtrans, nres) 
          ! if tie states (i,j) exist in B or AB ties: 
      ! for (tiitr, ties i --> j) , update j gamma += sum (i gamma for i --> j )
      ! for i in ( i --> j)
      !      i gamma = j gamma
          ! if tie states(i,j) exist in A or AB ties:
      ! for i,i2 in Atie(tiitr):
      !    for j in i --> j
      !       if j in Atie(tiitr[1,:])
      !         j2 = Atie(tiitr[2,Atie(tiitr[1,:]).index(j)])
      !         for titr in 1,T 
      !              t = sum(zeta(i,j,titr),zeta(i2,j2,titr))
      !              zeta(i,j,titr) = t
      !              zeta(i2,j2,titr) = t
      implicit none
      real(8),dimension(:,:),intent(inout) :: zeta
      integer,dimension(:,:),intent(in) :: zeta_dict
      integer,intent(in) :: nres
      real(8) :: sum_zeta, sum_gamma
      integer :: tiitr, titr,itr, j_ind, nitr_i, nitr_j, nitr_i2, nitr_j2, zind
      type(HMM_), intent(in)::HMMSTR
      type(trans),dimension(:),intent(in)::outtrans
      real(8),dimension(:, :),intent(inout) :: tgamma ! (nres, HMMSTR%n)

      ! state tracker B ties is HMMSTR.u
      !write(*,*) " BTIE"
      do tiitr = 1, HMMSTR%u
        do titr = 1, nres
           do itr = 1, HMMSTR%btie_lookup(tiitr,1)
              nitr_j = HMMSTR%btie_lookup(tiitr,itr)
      !        write(*,*) nitr_j
              tgamma(titr,tiitr) = tgamma(titr, tiitr) + tgamma(titr, nitr_j)
           enddo 
           do itr = 1, HMMSTR%btie_lookup(tiitr,1)
              nitr_j = HMMSTR%btie_lookup(tiitr,itr)
              tgamma(titr,nitr_j) = tgamma(titr, tiitr)
           enddo 
        enddo
      enddo
      !write(*,*) " ATIE"
      ! loop over all Aties
      do tiitr = 1, HMMSTR%nta
         nitr_i = HMMSTR%aties(tiitr,1)
         nitr_i2 = HMMSTR%aties(tiitr,2)
         do itr=1, outtrans(nitr_i)%c(0)
              nitr_j = outtrans(nitr_i)%c(itr)
              j_ind = HMMSTR%atie_lookup(nitr_j)   ! look in atie
              if (j_ind /= 0) then
                  do titr = 1, nres - 1
                      nitr_j2 = HMMSTR%aties(j_ind,2)
                      zind = zeta_dict(nitr_i,nitr_j)
                      sum_zeta = zeta(zind,titr)
                      zind = zeta_dict(nitr_i2,nitr_j2)
                      sum_zeta = sum_zeta + zeta(zind,titr)
                      call zeta_set(nitr_i,nitr_j,titr,zeta,zeta_dict,sum_zeta)
                      call zeta_set(nitr_i2,nitr_j2,titr,zeta,zeta_dict,sum_zeta)
                  enddo 
              endif
         enddo
      enddo 
      ! loop over all AB ties
      do tiitr = 1, HMMSTR%ntab
         nitr_i = HMMSTR%abties(tiitr,1)
         nitr_i2 = HMMSTR%abties(tiitr,2)
         do titr = 1, nres
             sum_gamma = tgamma(titr, nitr_i) + tgamma(titr, nitr_i2)
             tgamma(titr, nitr_i) = sum_gamma
             tgamma(titr, nitr_i2) = sum_gamma
         enddo
         do itr=1, outtrans(nitr_i)%c(0)
              nitr_j = outtrans(nitr_i)%c(itr)
              j_ind = HMMSTR%abtie_lookup(nitr_j)   ! look in atie
              if (j_ind /= 0) then
                  do titr = 1, nres - 1
                      nitr_j2 = HMMSTR%abties(j_ind,2)
                      zind = zeta_dict(nitr_i,nitr_j)
                      sum_zeta = zeta(zind,titr)
                      zind = zeta_dict(nitr_i2,nitr_j2)
                      sum_zeta = sum_zeta + zeta(zind,titr)
                      call zeta_set(nitr_i,nitr_j,titr,zeta,zeta_dict,sum_zeta)
                      call zeta_set(nitr_i2,nitr_j2,titr,zeta,zeta_dict,sum_zeta)
                  enddo
              endif
         enddo
      enddo 
  end subroutine untie

!------------------------------------------------------------------------------------------------------------
  subroutine model_str(model_file, model_path, model_name, model_iteration, slind, dind)
          implicit none
          character(len=*),intent(in) :: model_file
          integer,intent(in) :: slind, dind
          character(len=slind),intent(inout) :: model_path
          character(len=dind-slind-1),intent(inout):: model_name
          integer, intent(inout) :: model_iteration
          integer :: str_len, ind
          integer :: ios

          str_len = len_trim(model_file)
          model_path(:slind) = model_file(:slind)
          model_name(:dind-slind-1) = model_file(slind+1:dind)
          read(model_file(dind-3:dind-1),*,iostat=ios) model_iteration
          if ( ios/= 0) then
            model_iteration = 1
          endif
          !write(*,*) "MODEL_FILE = " , model_file(:len_trim(model_file))
          !write(*,*) "MODEL_PATH = " , model_path(:len_trim(model_path))
          !write(*,*) "MODEL_NAME = " , model_name
          !write(*,*) "MODEL_ITERATION = " , model_iteration
  end subroutine model_str

!------------------------------------------------------------------------------------------------------------
subroutine normalize1d (var)
   real(8),dimension(:),intent(inout) :: var
   real(8) :: sumpos
   sumpos = sum(var(:))
   !if (sumpos == 0 .or. sumpos /= sumpos) then
   !     write(*,*) "ERROR with new_priors"
   !endif
   var(:) = var(:) / sumpos
end subroutine normalize1d
!-----------------------------------------------------------------------------------------------------------
subroutine normalize (var, n)
   real(8),dimension(:,:),intent(inout) :: var
   integer :: n, nitr
   real :: sumpos
   do nitr=1,n
       sumpos = sum(var(nitr,:))
       if (sumpos == 0 .or. sumpos /= sumpos) then
            cycle
       endif
       var(nitr,:) = var(nitr,:) / sumpos
   enddo
end subroutine normalize
!------------------------------------------------------------------------------------------------------------
! TLB 8/23/20
! calculate new prior values
subroutine get_priors(xgamma, new_priors, n)
   ! ----------
   implicit none
   ! ----------
   real(8),dimension(:,:),intent(in)::xgamma ! (nres, HMMSTR%n)
   real(8),dimension(:),intent(inout)::new_priors   ! nres
   integer, intent(in) :: n
   ! ----------
   integer :: nitr
   ! ----------

   ! new_priors(i) += xgamma(i,1)  1 <= i <= modelR%n
   do nitr = 1, n
        if (xgamma(1,nitr) /=xgamma(1,nitr)) cycle
        if (new_priors(nitr) /= new_priors(nitr)) then
                new_priors(nitr) = xgamma(1,nitr)
        else
                new_priors(nitr) = new_priors(nitr) + xgamma(1,nitr)
        endif
   enddo
end subroutine get_priors

!------------------------------------------------------------------------------------------------------------
! TLB 9/28/21 gamma_zeta_naught_recurse
real recursive function gamma_zeta_naught_recurse(jnd, beta, HMMSTR, outrs, titr, trama, tss, ttm, t_prof, flags) result(a_ik)
     implicit none
     integer :: itr2, nitr_k
     real(8),dimension(:,:),intent(in)::beta   ! dimension N x T
     character(len=1),intent(in) :: trama, tss, ttm
     real,dimension(20),intent(in) :: t_prof
      integer,dimension(4),intent(in) :: flags
     !type(trans),dimension(:),intent(in)::intrs
     type(trans),dimension(:),intent(in)::outrs
     integer,intent(in) :: jnd ! current naught node being used
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: titr ! current timestep being assessed.
     a_ik = 0
     do itr2 =1, outrs(jnd)%c(0) ! all nodes which can be reached from the naught state
         nitr_k = outrs(jnd)%c(itr2)
         if ( nitr_k <= HMMSTR%k + HMMSTR%u) then
             if ( memobt(nitr_k) == 0 ) then
                 memobt(nitr_k) = gamma_zeta_naught_recurse(nitr_k,beta,HMMSTR,outrs,titr,trama,tss,ttm,t_prof,flags)
             endif
             a_ik = a_ik + HMMSTR%a(jnd, nitr_k) * memobt(nitr_k)
         else
             if ( memobt(nitr_k) == 0 ) then
                 memobt(nitr_k) = flagrat_brdc_prod(nitr_k,t_prof,HMMSTR,trama,tss,ttm,flags)
             endif
             a_ik = a_ik + beta(nitr_k, titr+1) * HMMSTR%a(jnd, nitr_k) * memobt(nitr_k)
         endif
     enddo                ! Beta & flagrat in ^
end function gamma_zeta_naught_recurse
! TLB 8/6/21 Removed ct scaling 
subroutine get_gamma_zeta( ct, beta, alpha, nres, HMMSTR, profile, ramaseq, ssseq, tmseq, outtrans, flags, tgamma, zeta, zeta_dict)
      ! ------
      implicit none
      ! -----
      character(len=1),dimension(:),intent(in)::ramaseq
      real(8),dimension(:),intent(in)::ct
      character(len=1),dimension(:),intent(in)::tmseq
      character(len=1),dimension(:),intent(in)::ssseq
      integer,dimension(4),intent(in) :: flags
      type(HMM_), intent(in)::HMMSTR
      real(8),dimension(:,:),intent(in)::alpha   ! NxT
      real(8),dimension(:,:),intent(in)::beta    ! NxT
      integer,intent(in)::nres
      real,dimension(:,:),intent(in)::profile  ! 20XT
      integer :: titr, nitr_i, nitr_j, itr
      type(trans),dimension(:),intent(in)::outtrans
      real(8) :: zeta_sum_i_t, zeta_sum_t
      !real(8),dimension(:, :, :),intent(inout) :: zeta ! (nres, HMMSTR%n, HMMSTR%n)
      real(8),dimension(:, :),intent(inout) :: zeta ! (# transitions, t)
      integer,dimension(:, :),intent(in) :: zeta_dict 
      real(8),dimension(:, :),intent(inout) :: tgamma ! (nres, HMMSTR%n)
      real(8) :: tct, tct1
      character(len=1) :: trama, tss, ttm 
      real,dimension(20) :: t_prof
      real(8) :: sum_tgamma
      integer :: n_wrong
      integer :: n_right
      logical :: deallocate_flag
      real(8) :: oldzv, zv
      deallocate_flag = .false.
      n_wrong = 0
      n_right = 0
      zeta = 0
      tgamma = 0
   if ( allocated(memobt) .eqv. .false.) then
          allocate(memobt(HMMSTR%n))
          deallocate_flag = .true.
   endif

        do titr =1,nres-1
          memobt = 0
          tct = ct(titr)
          tct1 = ct(titr+1)
          trama = ramaseq(titr+1)
          tss = ssseq(titr+1)
          ttm = tmseq(titr+1)
          t_prof = profile(:, titr+1)
          do nitr_i =1, HMMSTR%n
            !if ( nitr_i == 1) then
            !    tgamma(titr, nitr_i) = alpha(nitr_i, titr) * beta(nitr_i, titr) 
            !else
            tgamma(titr, nitr_i) = alpha(nitr_i, titr) * beta(nitr_i, titr) / tct 
            !endif
            do itr=1, outtrans(nitr_i)%c(0)
              nitr_j = outtrans(nitr_i)%c(itr)
              if (nitr_j <= HMMSTR%k + HMMSTR%u) then ! naught state traversal
                  if ( memobt(nitr_j) == 0) then
                       memobt(nitr_j) = gamma_zeta_naught_recurse(nitr_j, beta, HMMSTR, outtrans, titr,&
                               trama, tss, ttm, t_prof, flags)
                  endif
                  zv = memobt(nitr_j) * alpha(nitr_i,titr) * HMMSTR%a(nitr_i,nitr_j)
                  call zeta_set(nitr_i, nitr_j, titr, zeta, zeta_dict, zv)
                  !zeta(titr, nitr_i, nitr_j) = memobt(nitr_j) * alpha(nitr_i,titr) * HMMSTR%a(nitr_i, nitr_j)
              else
                  if ( memobt(nitr_j) == 0) then
                       memobt(nitr_j) = flagrat_brdc_prod(nitr_j,t_prof,HMMSTR,trama,tss,ttm, flags) 
                  endif
                  zv =  alpha(nitr_i,titr) * beta(nitr_j,titr+1) * HMMSTR%a(nitr_i,nitr_j) * memobt(nitr_j)
                  call zeta_set(nitr_i,nitr_j,titr,zeta,zeta_dict,zv)
              endif
            enddo
          enddo
          zeta_sum_t = sum(zeta(:,titr))
          if ( zeta_sum_t == 0 ) then
          !        write(*,*) nitr_i, nitr_j
                   cycle
          endif
          sum_tgamma = sum(tgamma(titr, :))
          do nitr_i=1, HMMSTR%n
            zeta_sum_i_t = 0
            do itr=1, outtrans(nitr_i)%c(0)
              nitr_j = outtrans(nitr_i)%c(itr)
              !zeta(titr,nitr_i,nitr_j) = zeta(titr,nitr_i,nitr_j) / zeta_sum_i_t
              call zeta_lookup(nitr_i,nitr_j,titr,zeta,zeta_dict,oldzv)
              zv = oldzv / zeta_sum_t
              call zeta_set(nitr_i,nitr_j,titr,zeta,zeta_dict,zv)
              zeta_sum_i_t = zeta_sum_i_t + zv
            enddo
            ! tgamma(titr, nitr_i) = sum(zeta(titr, nitr_i, : )
            
            tgamma(titr, nitr_i) = tgamma(titr, nitr_i) / sum_tgamma

            if (abs(tgamma(titr, nitr_i) - zeta_sum_i_t)  > 0.01 )  then
            !if ( .false. ) then
                write(*,*) " T                 = ", titr 
                write(*,*) " Node              = ", nitr_i
                write(*,*) " gamma - sum(zeta) = ", tgamma(titr, nitr_i) - zeta_sum_i_t
                write(*,*) " A(node, 1)        = ", HMMSTR%a(nitr_i, 1) 
                write(*,*) " A(1, node)        = ", HMMSTR%a(1, nitr_i)
                write(*,*) " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
                n_wrong = n_wrong + 1
             else if (tgamma(titr, nitr_i) > 0.00) then
                n_right = n_right + 1
             endif 

          enddo
          if ( sum(tgamma(titr, :)) < 0.001 ) then
                  write(*,*) "GAMMA = 0 @ titr = ", titr
          endif
        !write(*,*) sum(zeta(titr,:,:)) / HMMSTR%n
      enddo
     ! if ( n_wrong > 0 ) then
     if ( .false. ) then
       write(*,*) "INCORRECT = ", n_wrong
        write(*,*) "CORRECT   = ", n_right
        stop
     endif
     if ( deallocate_flag) then
             deallocate(memobt)
     endif

end subroutine

!-------------------------------------------------------------------------------------------------------------
! TLB 6/28/21
! TLB 7/12/21 Removed gamma normalization
! TLB 8/5/21 simplified inner loop to vector-summation, added gamma_normalization
! TLB 10.4.21 Removed tgamma input
! TLB 11/18/21 added new zeta logic
! TLB 11/23/21 added tied-state logic
  subroutine get_A(nres, HMMSTR, newA, outtrans, zeta, zeta_dict)
      ! ------
      implicit none
      ! -----
      type(HMM_), intent(in)::HMMSTR
      integer,intent(in)::nres
      real(8),dimension(:,:),intent(inout)::newA
      integer :: nitr_i, nitr_j, itr, nitr_i2, nitr_j2
      type(trans),dimension(:),intent(in)::outtrans
      real(8),dimension(:, :),intent(in) :: zeta 
      integer,dimension(:, :),intent(in) :: zeta_dict 
      real(8) :: B
      integer :: zind 
      do nitr_i = 1,HMMSTR%n
         do itr=1, outtrans(nitr_i)%c(0)
           nitr_j = outtrans(nitr_i)%c(itr)
           zind = zeta_dict(nitr_i,nitr_j)
           B = sum(zeta(zind,:nres-1)) 
           if ( B /= B) cycle 
           if ( HMMSTR%a(nitr_i, nitr_j) > 0 ) then
              newA(nitr_i, nitr_j) = newA(nitr_i, nitr_j) + B
           endif 
         enddo
         B = 0
      enddo
      return

  end subroutine get_A

!-------------------------------------------------------------------------------------------------------------
! TLB 5/21/20
! This recalculates B for residues and ramatypes according to the datapoints given in the profile and ramaseq.
! Can use gamma here because we are comparing aplha and beta in same timestep t, where get_a uses beta at t = t + 1
  subroutine get_B(xgamma, nres, HMMSTR, newB_residues, profile, newB_rama, ramaseq, newB_ss, ssseq, newB_tm, tmseq)
      ! ------
      implicit none
      ! -----
      type(HMM_),intent(in)::HMMSTR
      real(8),dimension(:,:),intent(in)::xgamma
      integer,intent(in)::nres
      real(8),dimension(:,:),intent(inout)::newB_residues
      real(8),dimension(:,:),intent(inout)::newB_rama
      real(8),dimension(:,:),intent(inout)::newB_ss
      real(8),dimension(:,:),intent(inout)::newB_tm
      real,dimension(:,:),intent(in)::profile
      character(len=1),dimension(:),intent(in)::ramaseq
      character(len=1),dimension(:),intent(in)::ssseq
      character(len=1),dimension(:),intent(in)::tmseq
      real(8) :: sab, bpsab, ab, rpsab, sspsab, tmpsab, bval, rval
      integer :: titr, nitr, bitr, tm_bitr
      ! in here, recalculate the value for A.

      !do nitr=1, modelR%n   ! loop over all nodes(N)
      !do bitr=1, 20         ! loop over all residues(B)
      !do titr=1, nres       ! loop over all times(T)


      do nitr=1, HMMSTR%n ! loop over all nodes(N)

          do bitr=1, 20 ! loop across all residues(B)

              bpsab = 0 ! sum of gamma where correct residue is predicted across t 
              sab = 0   ! sum of gamma across t
              rpsab = 0 ! sum of gamma where correct rama is predicted across t
              sspsab = 0 ! sum of gamma where correct ss is predicted across t
              tmpsab = 0 ! sum of gamma where correct tm is predicted across t
              
              do titr=1, nres
                  
                  ab = xgamma(titr,nitr) ! TLB 7/5 changes xgamma from (nitr, titr) to (titr, nitr) 

                  ! see get_gamma_zeta subroutine for details
                  sab = sab + ab
                  bpsab = bpsab + ab * profile(bitr,titr) 

                  if ( bitr <= mrama ) then ! putting rama in same loop as residues cuts down runtime
                      if (bitr == index(rama_string, ramaseq(titr))) then
                          rpsab = rpsab + ab  ! rama will not have prob. dist like aa profile... just add gamma * 1 
                      end if
                  end if
                  if (bitr <= mtm .and. tmseq(titr) /= 'U') then ! TLB added flag to skip undefined TM ground truth titr
                     if (bitr == index(tm_string, tmseq(titr))) then
                          tmpsab = tmpsab + ab
                     endif
                  endif
                  if (bitr <= mdssp) then
                     if (bitr == index(ss_string, ssseq(titr))) then
                          sspsab = sspsab + ab
                     endif
                  endif
              enddo

              ! TLB 10.19 Trying to weight by inverse of ground
              bval = bpsab !/ sab ! see eqn. 40c ! AA
              if (newB_residues(nitr,bitr) == newB_residues(nitr,bitr) .and. bval == bval &
                      .and. HMMSTR%AA_ground(bitr) /= 0) then
                   !newB_residues(nitr,bitr) = newB_residues(nitr,bitr) + ( bval / HMMSTR%AA_ground(bitr))
                   newB_residues(nitr,bitr) = newB_residues(nitr,bitr) + bval
              endif

              if ( bitr <= mrama ) then ! RAMA
                 rval = rpsab !/ sab ! see eqn. 40c
                 if (newB_rama(nitr,bitr) == newB_rama(nitr,bitr) .and. rval == rval &
                         .and. HMMSTR%RAMA_ground(bitr) /= 0 ) then
                     !newB_rama(nitr,bitr) = newB_rama(nitr,bitr) + ( rval / HMMSTR%RAMA_ground(bitr))
                     newB_rama(nitr,bitr) = newB_rama(nitr,bitr) + rval
                 endif
              end if
              if ( bitr <= mtm ) then ! TM
                 rval = tmpsab !/ sab ! see eqn. 40c
                 tm_bitr = bitr
                 if ( tm_bitr == 2 ) then ! ground truth = 2, sway towards 1. 1 supercedes 2, always 
                         tm_bitr = 1
                 endif
                 if (newB_tm(nitr,tm_bitr) == newB_tm(nitr,tm_bitr) .and. rval == rval .and. HMMSTR%TM_ground(tm_bitr) /= 0 ) then
                     !newB_tm(nitr,bitr) = newB_tm(nitr,bitr) + ( rval / HMMSTR%TM_ground(bitr) )
                     newB_tm(nitr,tm_bitr) = newB_tm(nitr,tm_bitr) + rval
                 endif
              end if
              if ( bitr <= mdssp ) then ! SS
                 rval = sspsab !/ sab ! see eqn. 40c
                 if (newB_ss(nitr,bitr) == newB_ss(nitr,bitr) .and. rval == rval .and. HMMSTR%SS_ground(bitr) /= 0) then
                     !newB_ss(nitr,bitr) = newB_ss(nitr,bitr) + ( rval / HMMSTR%SS_ground(bitr) )
                     newB_ss(nitr,bitr) = newB_ss(nitr,bitr) + rval
                 endif
              end if

          enddo

      enddo
  end subroutine get_B

!-------------------------------------------------------------------------------------------------------------
! Gets complex logsum variable expressing space of AA probabilistic distribution at time (t) and i = iq
  real(8) function bprof_bg(iq,profile,HMMSTR)
   !---------------------------
   implicit none
   !---------------------------
   integer,intent(in)::iq
   real,dimension(:),intent(in)::profile
   type(HMM_),intent(in)::HMMSTR
   !---------------------------
   integer::a
   real(8)::logsum, tmp
   !---------------------------
   logsum=0
   do a=1,m
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !tmp = HMMSTR%AA(iq,a)
       tmp = log(((HMMSTR%AA(iq,a)+(HMMSTR%AA_EPS*HMMSTR%AA_ground(a)))/((1.+HMMSTR%AA_EPS)*HMMSTR%AA_ground(a)))) 
       if ( tmp /= tmp .or. tmp / 2 == tmp) cycle
       logsum = logsum + tmp * profile(a)
     !  write(*,*) "STATE EMISSION = ", tmp, "GROUND TRUTH = ", profile(a)
   enddo
   !write(*,*) "LOGSUM = ", logsum
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   bprof_bg = logsum
   
   if ( .false. ) then 
           write(*,*) "bprof_bg = ", bprof_bg
           write(*,*) "profile = ", profile(:)
           write(*,*) "emission = ", HMMSTR%AA(iq,:)
           stop
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!
 end function bprof_bg
!-------------------------------------------------------------------------------------------------------------
  real(8) function flagrat_brdc_prod(iq,profile,HMMSTR,rama_c,dssp_c,tm_c,flags)
   !-----------------------------------------
   implicit none
   !-----------------------------------------
   character(len=1),intent(in)::rama_c,dssp_c,tm_c
   integer,dimension(4),intent(in) :: flags
   integer,intent(in)::iq
   real,dimension(:),intent(in)::profile
   type(HMM_),intent(in)::HMMSTR
   !-----------------------------------------
   integer::k
   real(8)::AA_p, tmp, RAMA_p, SS_p, TM_p
     !-----------------------------------------
   !AA_p=1; RAMA_p=1; SS_p=1; TM_p=1
   AA_p = 0; RAMA_p = 0; SS_p = 0; TM_p = 0
   if(iq==1) then ! Naught state, return 1 
     flagrat_brdc_prod=AA_p
     return
   endif
   ! get AA_p
   if (flags(1) == 1) then
        AA_p = bprof_bg(iq,profile,HMMSTR)
   endif
   ! get RAMA_p
   if (flags(2) == 1 .and. rama_c /= '?') then ! TLB 9.25 code review.. if rama_c == "?" then rama_p == 1 always. 
           do k=1,mrama
              if(rama_string(k:k)==rama_c ) exit
           enddo
           if( k < mrama+1) then
              if (HMMSTR%RAMA_ground(k) /= 0) then
                  RAMA_p = log((HMMSTR%RAMA(iq,k) + HMMSTR%RAMA_EPS*HMMSTR%RAMA_ground(k)) & 
                          /((1.+HMMSTR%RAMA_EPS)*HMMSTR%RAMA_ground(k)))
              endif
           endif
   endif
   ! get SS_p
   if (flags(3) ==1) then
           do k=1,mdssp
              if(ss_string(k:k)== dssp_c) exit
           enddo
           if (k<mdssp+1) then
              if (HMMSTR%SS_ground(k) /= 0) then
                   SS_p = log((HMMSTR%SS(iq,k) + HMMSTR%RAMA_EPS*HMMSTR%SS_ground(k))/((1.+HMMSTR%RAMA_EPS)*HMMSTR%SS_ground(k)))
              endif
           endif
   endif
   ! get TM_p
   if (flags(4) == 1 .and. tm_c /= 'U') then
           do k=1,mtm
              if(tm_string(k:k)== tm_c .or. (tm_c == '2' .and. tm_string(k:k) == '1')) exit
           enddo
           if (k<mtm+1) then
              if (HMMSTR%TM_ground(k) /= 0) then
                   TM_p = log((HMMSTR%TM(iq,k) + HMMSTR%TM_EPS*HMMSTR%TM_ground(k))/((1 + HMMSTR%TM_EPS)*HMMSTR%TM_ground(k)))
              endif
           endif
   endif
  if ( AA_p - 100 == AA_p .or. AA_p /= AA_p ) then
          AA_p = 0
  endif
  if ( SS_p - 100 == SS_p .or. SS_p /= SS_p) then
          SS_p = 0
  endif
  if ( TM_p - 100 == TM_p .or. TM_p /= TM_p) then
         TM_p = 0
  endif
  if ( RAMA_p - 100 == RAMA_p .or. RAMA_p /= RAMA_p) then
          RAMA_p = 0
  endif
  tmp = exp( AA_p + RAMA_p + SS_p + TM_p )
  !flagrat_brdc_prod = ( tmp / ( 1 + tmp ) ) ! Experimental, ensures returned value > 0, < 1
  flagrat_brdc_prod = tmp
  !write(*,*) flagrat_brdc_prod
  if ( .false. ) then
      c = c + 1
      write(*,*) "AAP * RAMAP * TMP = ret", AA_p, RAMA_p, TM_p,  "=", tmp
      write(*,*) "PROFILE = ", profile
      write(*,*) "RAMA Char = ", rama_c
      write(*,*) "TM Char = ", tm_c
      write(*,*) "AA EMISSION = ", HMMSTR%AA(iq,:)
      write(*,*) "RAMA EMISSION = ", HMMSTR%RAMA(iq,:)
      write(*,*) "TM EMISSION = ", HMMSTR%TM(iq,:)
  endif
end function flagrat_brdc_prod
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
  subroutine get_intrans(intrs,HMMSTR)
   ! intrs(i) = ncount, 1...ncount node indeces transition from
   !---------------------------------------
   implicit none
   !---------------------------------------
   type(HMM_),intent(in)::HMMSTR
   type(trans),dimension(:),allocatable,intent(out)::intrs
   !---------------------------------------
   integer::nq
   integer::ncount
   integer::i
   integer::j
   !---------------------------------------
   nq=HMMSTR%n
   allocate(intrs(nq))
   do j=1,nq
     ncount=0
     do i=1,nq
       if(HMMSTR%a(i,j)/=0) ncount = ncount + 1
     enddo
     ! ncount = # transitions into node j 
     allocate(intrs(j)%c(0:ncount))
     ncount=0
     do i=1,nq
       if(HMMSTR%a(i,j)/=0) then
         ncount = ncount + 1
         intrs(j)%c(ncount) = i
       endif
       intrs(j)%c(0) = ncount
     enddo
   enddo
  end subroutine get_intrans
!-------------------------------------------------------------------------------------------------------------
  subroutine get_outtrans(outtrs,HMMSTR)
   ! outrs(i) = ncount, 1...ncount node indeces transition to
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   type(HMM_),intent(in)::HMMSTR
   type(trans),dimension(:),allocatable,intent(out)::outtrs
   !--------------------------------------------
   integer::nq
   integer::ncount
   integer::i
   integer::j
   !--------------------------------------------
   nq=HMMSTR%n
   allocate(outtrs(nq))
   do i=1,nq
     ncount=0
     do j=1,nq
       if(HMMSTR%a(i,j)/=0) ncount = ncount + 1
     enddo
     allocate(outtrs(i)%c(0:ncount))
     ncount=0
     do j=1,nq
       if(HMMSTR%a(i,j)/=0) then
         ncount = ncount + 1
         outtrs(i)%c(ncount) = j
       endif
       outtrs(i)%c(0) = ncount
     enddo
   enddo
  end subroutine get_outtrans
!-------------------------------------------------------------------------------------------------------------
  subroutine set_flags(ramaseq,nres,angles)    !----------------------------------------
   implicit none
   !----------------------------------------
   character(len=1),dimension(:),intent(out)::ramaseq
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::angles
   !----------------------------------------
   integer::i
   !----------------------------------------
   ramaseq="?"
   do i=1,nres
      if(all(angles(:,i)/=999.) .and. all(angles(:,i)/=0.00000000)) then
         ramaseq(i)=ramatype(angles(1,i),angles(2,i),angles(3,i))
      endif
   enddo
  end subroutine set_flags
!-------------------------------
! Given phi,psi,omega angles... return the discretized character representing the ramachandran space type
  character(len=1) function ramatype(ph,ps,om)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   real,intent(in)::ph
   real,intent(in)::ps
   real,intent(in)::om
   !--------------------------------------------
   integer::i
   integer::icen
   real::d
   real::dmin
   real::ds
   real::df
   real,dimension(11,2)::ramacen
   !--------------------------------------------
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
   if(ph>=998.)then; ramatype="?"; return; endif
   if(ps>=998.)then; ramatype="?"; return; endif
   if(om<90.and.om>-90)then; ramatype="c"; return; endif
   dmin=999
   do i=1,mrama-1
     df=abs(ph-ramacen(i,1)); if(df>180.) df= 360.-df
     ds=abs(ps-ramacen(i,2)); if(ds>180.) ds= 360.-ds
     d=sqrt(ds*ds+df*df)
     if(d<dmin) then
       icen=i; dmin=d
     endif
   enddo
   ! rama_string="HGBEdbeLlxc"
   ramatype=rama_string(icen:icen)
  end function ramatype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  
real recursive function alpha_naught_recurse(ind, titr, alpha, HMMSTR, intrs) result(a_ik)
     implicit none
     !real :: a_ik
     integer :: ind2, kitr
     !integer,intent(in) :: nitr ! node of alpha being calculated
     real(8),dimension(:,:),intent(in)::alpha   ! dimension N x T
     type(trans),dimension(:),intent(in)::intrs
     integer,intent(in) :: ind ! current naught node being used
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: titr ! current timestep being assessed.
     a_ik = 0
     do kitr = 1, intrs(ind)%c(0)
        ind2 = intrs(ind)%c(kitr)
        if (ind2 <= HMMSTR%k + HMMSTR%u) then
            if ( memobt(ind2) == 0 ) then
                    memobt(ind2) = alpha_naught_recurse(ind2, titr, alpha, HMMSTR, intrs)
            endif
            a_ik = a_ik + HMMSTR%a(ind2, ind) * memobt(ind2)
        endif
        a_ik = a_ik + HMMSTR%a(ind2, ind) * alpha(ind2, titr - 1)
     enddo
end function alpha_naught_recurse

  !------------------------------------------------------------------------------------------------------------
  subroutine getalpha(ct,alpha,seqres_len,HMMSTR,profile,intrs,outrs,ramaseq,ssseq,tmseq,flags,start_t,end_t)
   !--------------------------------------------------
   implicit none
   !--------------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   character(len=1),dimension(:),intent(in)::ssseq
   character(len=1),dimension(:),intent(in)::tmseq
   integer,dimension(4),intent(in) :: flags
   integer,intent(in)::seqres_len
   real(8),dimension(:,:),intent(out)::alpha   ! dimension N x T
   real,dimension(:,:),intent(in)::profile ! 20 X T  
   real(8),dimension(:),intent(inout)::ct
   type(HMM_),intent(in):: HMMSTR
   type(trans),dimension(:),intent(in)::intrs
   type(trans),dimension(:),intent(in)::outrs
   character(len=1) :: trama, ttm, tss
   integer,intent(in) :: start_t, end_t
   !--------------------------------------------------
   integer :: nitr   ! node iterator ( K + 1 ... N ) (i)
   integer :: titr   ! T iterator ( 2 ... seqres_len )
   integer :: kitr   ! intrans iterator for naught nodde ( 1 ... intrs(iitr)%c(0) ) where iitr <= k
   integer :: iitr   ! intrans iterator ( 1 ... intrs(nitr)%c(0) )
   integer :: ind    ! node which can transfer into node j  
   integer :: ind2   ! node which can transfer into naught node k 
   integer :: n      ! N
   integer :: u      ! U
   integer :: k      ! K
   real(8)    :: b_t1   ! bj(Ot+1)
   real(8)    :: a_ij   ! sum over i ( aij + sum over k ( aik * akj ) )
   real(8)    :: a_ik   ! sum of over i ( aik ) 
   real,dimension(20) :: t_prof
   logical :: naught_scaling
   logical :: deallocate_flag
   deallocate_flag = .false. ! sometimes we allocate & deallocate memobt within this function
   naught_scaling = .false.
   !--------------------------------------------------
   if ( allocated(memobt) .eqv. .false.) then
          allocate(memobt(HMMSTR%n))
          deallocate_flag = .true.
   endif
   !!alpha(nq,nres)  !!intrans(nr)%c(nonzero)! profile(a20,nres)
   ! initialize alpha, ct !
   alpha=0;ct=0
   ! set starting and ending points
   n = HMMSTR%n
   u = HMMSTR%u
   k = HMMSTR%k
   ! Initialize alpha 1 (i) for all i  k+1 <= i <=  n
   ! alpha(1) can either be reached through priors or through the naught node
   do titr = start_t, end_t
    !write(*,*) titr
     memobt = 0
     trama = ramaseq(titr)
     tss = ssseq(titr)
     ttm = tmseq(titr)
     t_prof = profile(:, titr)
     do nitr = u+k+1, n
    !    write(*,*) nitr
        b_t1 = flagrat_brdc_prod( nitr, t_prof, HMMSTR, trama, tss, ttm, flags)
        if ( titr == start_t ) then ! t = 1 uses priors, a special case to enter model
           ! TLB 9/7/2021 removed this naught-initialization step because it 
           ! doesn't follow precendence of how we use the naught state.
           if ( HMMSTR%a(1, nitr) > 0.000) then
              a_ij = HMMSTR%a(1, nitr)
              alpha(nitr, titr) = a_ij * b_t1 * HMMSTR%priors(1)
           endif
           a_ij = HMMSTR%priors(nitr)
           alpha(nitr, titr) = alpha(nitr, titr) + a_ij * b_t1
        else ! t > 1
          a_ij = 0
          ! loop over all intrans....
          do iitr = 1, intrs(nitr)%c(0)
              ind = intrs(nitr)%c(iitr)
              ! if ind refers to a naught node then loop over intrans into naught node ...
              if (ind <= u+k) then  ! note if i == 1 then intrs(1) will be empty (priors)
                 if ( memobt(ind) == 0 ) then
                         memobt(ind) = alpha_naught_recurse(ind, titr, alpha, HMMSTR, intrs) ! returns a_ik
                 endif
                 a_ij = a_ij + HMMSTR%a(ind, nitr) * memobt(ind)
              else
                 a_ij = a_ij + HMMSTR%a(ind, nitr) * alpha(ind, titr - 1)
              endif
           enddo
         alpha(nitr, titr) = a_ij * b_t1
        endif
        ct(titr) = ct(titr) + alpha(nitr, titr)
      enddo
      
      if (naught_scaling) then
              do iitr =1, intrs(1)%c(0) ! loop over all states incoming to naught state
                ind = intrs(1)%c(iitr)
                alpha(1, titr) = alpha(1, titr) + alpha(ind,titr) * HMMSTR%a(ind, 1)
              enddo
              ct(titr) = ct(titr) + alpha(1, titr) ! naught alpha used to calculate scale factor
              ct(titr) = 1. / ct(titr)
              do nitr = 1, n   ! scale all states
                 alpha(nitr, titr) = alpha(nitr, titr) * ct(titr)
              enddo
      else
              
             ct(titr) = 1. / ct(titr)
             do nitr = u+k+1, n  ! scale all states except the naught state
                alpha(nitr, titr) = alpha(nitr, titr) * ct(titr) 
             enddo

             if ( titr > 1 ) then
                     do nitr = u+1, k+u
                             do iitr =1, intrs(nitr)%c(0) ! loop over all states incoming to naught state
                                ind = intrs(nitr)%c(iitr)                      ! TLB 9/10/21
                                if ( ind <= HMMSTR%k + HMMSTR%u) then ! TLB 9/23/21 alpha recurse logic
                                        alpha(nitr,titr) = alpha(nitr, titr) + HMMSTR%a(ind, nitr) * & 
                                                alpha_naught_recurse(ind, titr, alpha, HMMSTR, intrs) 
                                else
                                        alpha(nitr, titr) = alpha(nitr, titr) + HMMSTR%a(ind, nitr) * alpha(ind,titr-1) 
                                endif
                             enddo
                     enddo
             endif
       endif
   enddo 
   do nitr=1,n
       do titr=start_t,end_t
         if (alpha(nitr,titr) /= alpha(nitr,titr)) then ! check for null value
            kitr = 0
                 write(*,*) "ERROR IN getalpha, NaN value"        
         else if (alpha(nitr,titr)-1 == alpha(nitr,titr)) then ! check for infinity
            kitr = 0 
                 write(*,*) "ERROR IN getalpha, infinity value"        
         endif
      enddo
   enddo
    if ( deallocate_flag ) then
            deallocate(memobt)
    endif

  end subroutine

!-------------------------------------------------------------------------------------------------------------
real recursive function beta_naught_recurse(jnd, titr, beta, HMMSTR, outrs, t_prof, trama, tss, ttm, flags) result(a_ik)
     implicit none
     integer :: jnd2, kitr
     real(8),dimension(:,:),intent(in)::beta  ! Dimension N X T
     type(trans),dimension(:),intent(in)::outrs
     integer,intent(in) :: jnd ! current naught node being used
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: titr ! current timestep being assessed.
     character(len=1),intent(in) :: trama, tss, ttm
     real,dimension(20),intent(in) :: t_prof
     integer,dimension(4),intent(in)::flags
     a_ik = 0
     do kitr = 1, outrs(jnd)%c(0)
         jnd2 = outrs(jnd)%c(kitr)
         if ( jnd2 <= HMMSTR%k + HMMSTR%u) then
                 if ( memobt(jnd2) == 0) then
                     memobt(jnd2) = beta_naught_recurse(jnd2, titr, beta, HMMSTR, outrs, t_prof, trama, tss, ttm, flags)
                 endif
                 a_ik = a_ik + HMMSTR%a(jnd,jnd2) * memobt(jnd2)
         else
                 if ( memobt(jnd2) == 0) then
                     memobt(jnd2) = flagrat_brdc_prod( jnd2, t_prof, HMMSTR, trama, tss, ttm, flags)
                 endif
                 a_ik = a_ik + HMMSTR%a(jnd, jnd2) * beta(jnd2, titr+1) * memobt(jnd2)
         endif
     enddo
end function beta_naught_recurse

!-------------------------------------------------------------------------------------------------------------
subroutine getbeta(ct,beta,seqres_len,HMMSTR,profile,intrs,outrs,ramaseq,ssseq,tmseq,flags,start_t,end_t)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   character(len=1),dimension(:),intent(in)::ssseq
   character(len=1),dimension(:),intent(in)::tmseq
   integer,dimension(4),intent(in)::flags
   integer,intent(in)::seqres_len
   integer,intent(in)::start_t,end_t
   real(8),dimension(:),intent(inout)::ct
   real(8),dimension(:,:),intent(out)::beta  ! Dimension N X T
   real,dimension(:,:),intent(in)::profile
   type(HMM_),intent(in)::HMMSTR
   type(trans),dimension(:),intent(in)::outrs
   type(trans),dimension(:),intent(in)::intrs
   !--------------------------------------------
   integer :: nitr   ! node iterator ( K + 1 ... N ) (j)
   integer :: titr   ! T iterator ( seqres_len ... 1 )
   integer :: kitr   ! outrans iterator for naught nodde ( 1 ... ontrs(oitr)%c(0) ) where oitr <= k
   integer :: oitr   ! outrans iterator ( 1 ... otrs(nitr)%c(0) )
   integer :: jnd    ! node which can be transferred to from nitr  
   integer :: jnd2   ! node which can be transferred to from jnd
   integer :: n      ! N
   integer :: u      ! N
   integer :: k      ! K
   real(8)    :: b_t1   ! bj(Ot+1)
   real(8)    :: a_ij   ! sum over i ( aij + sum over k ( aik * akj ) )
   character(len=1) :: trama, tss, ttm 
   real(8) :: tct
   real,dimension(20) :: t_prof
   !--------------------------------------------
   integer :: iitr, ind, ind2 
   logical :: deallocate_flag
   deallocate_flag = .false. ! sometimes we allocate & deallocate memobt within this function
   n = HMMSTR%n
   u = HMMSTR%u
   k = HMMSTR%k  
   beta = 0

   !ct_prod = ct(seqres_len)
   ! initialize beta T  ... see eqn (24) 
   do nitr = 1, n
      beta(nitr, end_t) = ct(end_t) ! Fixed w/ scaled value ... see eqn (94)
   enddo  
      if ( allocated(memobt) .eqv. .false.) then
          allocate(memobt(HMMSTR%n))
          deallocate_flag = .true.
   endif

   ! Inductive step ... see eqn (25)
   do titr = end_t -1, start_t, -1
       memobt = 0
       trama = ramaseq(titr+1)
       tss = ssseq(titr+1)
       ttm = tmseq(titr+1)
       tct = ct(titr)
       t_prof = profile(:, titr+1)
       !ct_prod = ct_prod * tct
       do nitr = 1, n  ! beta(i, t) = sum over j ( beta(j,t+1) * aij * bj(Ot) )
                           ! and ....
                           ! aij = sum over j ( aij + sum over k ( aik * akj ) ) 
                           ! where k is a naught node
           a_ij=0
           do oitr=1,outrs(nitr)%c(0)
              jnd = outrs(nitr)%c(oitr)
              ! if jnd <= k then find all nodes that can be reached from naught node jnd
              if (jnd <= k + u) then
                 if (memobt(jnd) == 0 ) then
                      memobt(jnd) = beta_naught_recurse(jnd,titr,beta,HMMSTR,outrs,t_prof,trama,tss,ttm,flags) 
                 endif
                 a_ij = a_ij + HMMSTR%a(nitr, jnd) * memobt(jnd)
              else
                   if ( memobt(jnd) == 0) then
                         memobt(jnd) = flagrat_brdc_prod( jnd, t_prof, HMMSTR, trama, tss, ttm, flags)
                   endif
                 a_ij = a_ij + HMMSTR%a(nitr, jnd) * beta(jnd, titr + 1) * memobt(jnd) 
              endif
           enddo

           !if (nitr > 1) then
           beta(nitr, titr) =  a_ij * tct  ! beta(i,t) (hat) = ct(t) * beta(i,t)   (94)
           !else
           !       beta(nitr,titr) = a_ij ! only for naught state
           !endif
       enddo
    enddo
    if ( deallocate_flag ) then
            deallocate(memobt)
    endif

  end subroutine
real function get_P(n,t,alpha,ct) result(p)
   implicit none
   !--------------------------------------------
   real,dimension(:,:),intent(in)::alpha
   real,dimension(:),intent(in)::ct
   integer,intent(in) :: n,t
   integer nitr
   p = 0.0
   do nitr = 1, n
      p = p + alpha(nitr, t) / ct(t)
   enddo
end function get_P
subroutine write_p(p_arr, p, pfile)
        implicit none
        real, dimension(:,:), intent(in) :: p_arr ! probability array
        integer, intent(in) :: p ! number of probabilites
        character(len=*),intent(in) :: pfile
        integer :: ios, pitr
        open(3, file = trim(pfile), iostat = ios)
        if (ios /= 0) then
           write(*,*) "IMPROPERLY formatted P file."
           write(*,*) ios
        endif
        do pitr = 1 , p 
           write(3,*) p_arr(pitr,1), p_arr(pitr,2)
        enddo
        close(3)
end subroutine

subroutine gamma_voting(t, HMMSTR, xgamma1, seq1, ssseq1, ramaseq1, tmseq1)
     implicit none
     integer,intent(in) :: t
     Type(HMM_),intent(in) :: HMMSTR
     real(8),dimension(:,:),intent(in) :: xgamma1
     character(len=1),dimension(:),intent(inout) :: seq1
     character(len=1),dimension(:),intent(inout) :: ssseq1
     character(len=1),dimension(:),intent(inout) :: ramaseq1
     character(len=1),dimension(:),intent(inout) :: tmseq1
     real(8),dimension(20) :: aa_prof
     real(8),dimension(6) :: ss_prof
     real(8),dimension(11) :: rama_prof
     real(8),dimension(10) :: tm_prof
     integer :: titr, nitr, ind

     do titr=1,t
        aa_prof = 0
        ss_prof = 0
        tm_prof = 0
        rama_prof = 0
        do nitr=HMMSTR%k + 1,HMMSTR%n
            aa_prof    = aa_prof  + xgamma1(titr,nitr)*HMMSTR%AA(nitr,:)
            ss_prof    = ss_prof  + xgamma1(titr,nitr)*HMMSTR%SS(nitr,:)
            tm_prof    = tm_prof  + xgamma1(titr,nitr)*HMMSTR%TM(nitr,:)
            rama_prof  = rama_prof  + xgamma1(titr,nitr)*HMMSTR%RAMA(nitr,:)
        enddo
        call voting(HMMSTR, seq1, tmseq1, ssseq1, ramaseq1, aa_prof, tm_prof, ss_prof, rama_prof, titr)
     enddo
end subroutine gamma_voting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COPIED FROM DIAGNOSTICS DIRECTORY 1/7/22!!!!!!!
subroutine voting(HMMSTR, aaseq1, tmseq1, ssseq1, ramaseq1, aa_prof, tm_prof, ss_prof, rama_prof, titr)
   implicit none
    type(HMM_),intent(in) :: HMMSTR
    character(len=1),dimension(:),intent(inout) :: aaseq1, tmseq1, ssseq1, ramaseq1
    real(8),dimension(20),intent(in) :: aa_prof
    real(8),dimension(6),intent(inout) :: ss_prof
    real(8),dimension(10),intent(in) :: tm_prof
    real(8),dimension(11),intent(inout) :: rama_prof
    real(8),dimension(6) :: rama_parties
    real(8),dimension(3) :: ss_parties
    real(8),dimension(8) :: tm_parties
    integer, intent(in) :: titr
    integer :: ind
    character(len=1) :: chr

    rama_parties = 0
    ! AA
    ind = maxloc(aa_prof,1)
    chr = aa(ind:ind)
    aaseq1(titr) = chr
    ! TM
    tm_parties(1) = tm_prof(1) ! 1
    tm_parties(2) = tm_prof(3) ! B
    tm_parties(3) = tm_prof(4) ! H
    tm_parties(4) = tm_prof(5) ! C
    tm_parties(5) = tm_prof(6) ! I 
    tm_parties(6) = tm_prof(7) ! L 
    tm_parties(7) = tm_prof(8) ! F
    tm_parties(8) = tm_prof(9) ! U
    if ( any(HMMSTR%TM_reweight /= 0 ) ) then
            tm_parties(1) = tm_parties(1) * HMMSTR%TM_reweight(1) ! 1
            tm_parties(2) = tm_parties(2) * HMMSTR%TM_reweight(2) ! B
            tm_parties(3) = tm_parties(3) * HMMSTR%TM_reweight(3) ! H
            tm_parties(4) = tm_parties(4) * HMMSTR%TM_reweight(4) ! C
            tm_parties(5) = tm_parties(5) * HMMSTR%TM_reweight(5) ! I 
            tm_parties(6) = tm_parties(6) * HMMSTR%TM_reweight(6) ! L 
            tm_parties(7) = tm_parties(7) * HMMSTR%TM_reweight(7) ! F
            tm_parties(8) = tm_parties(8) * HMMSTR%TM_reweight(8) ! U
    endif
    ind = maxloc(tm_parties,1)
    if (ind > 1) then
        tmseq1(titr) = tm_string(ind+1:ind+1)
    else
        tmseq1(titr) = tm_string(ind:ind)
    endif
    ! RAMA ! PRIMARIES:
    ! character(len=mrama),parameter :: rama_string="HGBEdbeLlxc"
    rama_parties(1) = rama_prof(1) + rama_prof(2) ! ( H, G )
    rama_parties(2) = rama_prof(3) + rama_prof(6) ! ( B, b ) 
    rama_parties(3) = rama_prof(4) + rama_prof(5) + rama_prof(7) ! ( E, d, e ) 
    rama_parties(4) = rama_prof(8) + rama_prof(9) ! ( L, l ) 
    rama_parties(5) = rama_prof(10)  ! ( x ) 
    rama_parties(6) = rama_prof(11) ! ( c )

    ! HMMSTR (untrained) reweighting
    if ( any(HMMSTR%RAMA_reweight /= 0) ) then
            rama_parties(1) = rama_parties(1)  * HMMSTR%RAMA_reweight(1)   ! HG
            rama_parties(2) = rama_parties(2)  * HMMSTR%RAMA_reweight(2)   ! Bb
            rama_parties(3) = rama_parties(3)  * HMMSTR%RAMA_reweight(3)   ! Ede
            rama_parties(4) = rama_parties(4)  * HMMSTR%RAMA_reweight(4)   ! Ll
            rama_parties(5) = rama_parties(5)  * HMMSTR%RAMA_reweight(5)   ! x 
            rama_parties(6) = rama_parties(6)  * HMMSTR%RAMA_reweight(6)   ! c
    endif
    ind = maxloc(rama_parties,1)
    !character(len=mrama),parameter :: rama_string="HGBEdbeLlxc"
    select case ( ind )  ! RUNOFF ELECTIONS:
        case ( 1 ) ! (H,G)  
             rama_prof(3:) = 0
             call normalize1d(rama_prof)
             ind = maxloc(rama_prof,1)
             ramaseq1(titr) = rama_string(ind:ind)
        case ( 2 ) ! ( B, b )
             rama_prof(1:2) = 0
             rama_prof(7:11) = 0
             rama_prof(4:5) = 0
             rama_prof(3) = rama_prof(3)
             rama_prof(6) = rama_prof(6)
             call normalize1d(rama_prof)
             ind = maxloc(rama_prof,1)
             ramaseq1(titr) = rama_string(ind:ind)
        case ( 3 ) ! ( E, d, e ) 
             rama_prof(1:3) = 0
             rama_prof(6) = 0
             rama_prof(8:11) = 0
             rama_prof(4) = rama_prof(4)
             rama_prof(5) = rama_prof(5)
             rama_prof(7) = rama_prof(7)
             call normalize1d(rama_prof)
             ind = maxloc(rama_prof,1)
             ramaseq1(titr) = rama_string(ind:ind)
        case ( 4) ! ( L, l) 
             rama_prof(1:7) = 0
             rama_prof(10:11) = 0
             call normalize1d(rama_prof)
             ind = maxloc(rama_prof,1)
             ramaseq1(titr) = rama_string(ind:ind)
        case ( 5 ) ! (x)  
            rama_prof(1:9) = 0
            rama_prof(11:11) = 0
            call normalize1d(rama_prof)
            ind = maxloc(rama_prof,1)
            ramaseq1(titr) = rama_string(ind:ind)
        case ( 6 ) ! ( c )  
            rama_prof(1:10) = 0
            call normalize1d(rama_prof)
            ind = maxloc(rama_prof,1)
            ramaseq1(titr) = rama_string(ind:ind)
    end select
    ! SS PRIMARIES:
    ss_parties(1) = ss_prof(1) + ss_prof(3) ! ( H, G )
    ss_parties(2) = ss_prof(4) + ss_prof(5) + ss_prof(6) ! ( S, T, _ )
    ss_parties(3) = ss_prof(2) ! ( E ) 
    ! HMMSTR (untrained) reweighting

    if ( any(HMMSTR%SS_reweight /= 0)) then
            ss_parties(1) = ss_parties(1)  * HMMSTR%SS_reweight(1) ! HG
            ss_parties(2) = ss_parties(2)  * HMMSTR%SS_reweight(2) ! ST_
            ss_parties(3) = ss_parties(3)  * HMMSTR%SS_reweight(3) ! E
    endif

    ind = maxloc(ss_parties,1)
    !character(len=mdssp),parameter :: ss_string="HEGST_"
    select case(ind)
        case ( 1 ) ! ( H, G ) 
            ss_prof(2:2) = 0
            ss_prof(4:6) = 0
            call normalize1d(ss_prof)
            ind = maxloc(ss_prof, 1)
            ssseq1(titr) = ss_string(ind:ind)
        case ( 2 )  ! ( S, T, _ ) 
            ss_prof(1:3) = 0
            call normalize1d(ss_prof)
            ind = maxloc(ss_prof, 1)
            ssseq1(titr) = ss_string(ind:ind)
        case ( 3 ) ! ( E ) 
            ss_prof(1:1) = 0
            ss_prof(3:) = 0
            ind = maxloc(ss_prof,1)
            ssseq1(titr) = ss_string(ind:ind)
    end select

end subroutine voting

! VITERBI NAUGHT RECURSE CODE 9/26
real(8) recursive function Viterbi_naught_recurse(ind, titr, delta, HMMSTR, intrs) result(a_ik)
     implicit none
     !real :: a_ik
     integer :: ind2, kitr
     real(8) ::  comp
     !integer,intent(in) :: nitr ! node of alpha being calculated
     real(8),dimension(:,:),intent(in)::delta   ! dimension N x T
     type(trans),dimension(:),intent(in)::intrs
     integer,intent(in) :: ind ! current inc-not-naught node being used
     !integer,intent(in) :: jnd ! current out-not-naught node being used
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: titr ! current timestep being assessed.
     a_ik = 0
     !write(*,*) ind, titr
     do kitr = 1, intrs(ind)%c(0)
        ind2 = intrs(ind)%c(kitr)
        if (ind2 <= HMMSTR%k + HMMSTR%u) then
            comp = HMMSTR%a(ind2, ind) * Viterbi_naught_recurse(ind2, titr, delta, HMMSTR, intrs)
        endif
        comp = HMMSTR%a(ind2, ind) * delta(ind2, titr - 1)
        if ( comp > a_ik ) then
                a_ik = comp
        endif
     enddo
end function Viterbi_naught_recurse
real(8) recursive function Viterbi_naught_traverse(ind, jnd, titr, HMMSTR, intrs) result(a_ik)
     implicit none
     !real :: a_ik
     integer :: ind2, kitr
     !integer,intent(in) :: nitr ! node of alpha being calculated
     !real,dimension(:,:),intent(in)::delta   ! dimension N x T
     type(trans),dimension(:),intent(in)::intrs
     !type(trans),dimension(:),intent(in)::outrs
     integer,intent(in) :: ind ! current inc-not-naught node being used
     integer,intent(in) :: jnd ! current out-not-naught node being used
     type(HMM_),intent(in) :: HMMSTR
     integer,intent(in) :: titr ! current timestep being assessed.
     a_ik = 0
     !write(*,*) ind, titr
     do kitr = 1, intrs(jnd)%c(0) 
        ind2 = intrs(jnd)%c(kitr)
        if (ind2 <= HMMSTR%u + HMMSTR%k) then
            a_ik = a_ik + HMMSTR%a(ind2, jnd) * Viterbi_naught_traverse(ind, ind2, titr, HMMSTR, intrs)
        else if ( ind2 == ind ) then
            a_ik = a_ik + HMMSTR%a(ind2, jnd)
        else if ( ind2 > ind .and. ind2 > HMMSTR%k ) then
             exit
        endif
     enddo
end function Viterbi_naught_traverse

    ! TLB 9/8/2021
  !------------------------------------------------------------------------------------------------------------
  subroutine viterbi(seqres_len, HMMSTR, profile, intrs, outrs, seq, ramaseq, ssseq, tmseq, flags, paradigm, QUIET)
   !--------------------------------------------------
   implicit none
   !--------------------------------------------------
   character(len=1),dimension(:),intent(in)::seq
   character(len=1),dimension(:),intent(in)::ramaseq
   character(len=1),dimension(:),intent(in)::ssseq
   character(len=1),dimension(:),intent(in)::tmseq
   integer,dimension(4),intent(in) :: flags
   integer,intent(in)::seqres_len
   real,dimension(:,:),intent(in)::profile ! 20 X T  
   type(HMM_),intent(in):: HMMSTR
   type(trans),dimension(:),intent(in)::intrs
   type(trans),dimension(:),intent(in)::outrs
   character(len=1) :: trama, ttm, tss
   real(8),dimension(20) :: aa_prof
   real(8),dimension(6) :: ss_prof
   real(8),dimension(10) :: tm_prof
   real(8),dimension(11) :: rama_prof
   integer :: loc
   real(8)    :: max_loc, val
   real(8)    :: b_t1   ! bj(Ot+1)
   real(8)    :: a_ij   ! sum over i ( aij + sum over k ( aik * akj ) )
   integer :: titr, nitr, nitr_2, iitr, n, k, ind, jitr, nnitr, u
   real,dimension(20) :: t_prof
   !---------------------------------------
   real(8),dimension(:,:),allocatable :: delta
   integer,dimension(:,:),allocatable :: upsilon
   integer,dimension(:),allocatable ::  Q
   character(len=1),dimension(:),allocatable :: out_seq, out_ssseq, out_ramaseq, out_tmseq
   character(len=1),dimension(:,:),allocatable,intent(inout) :: paradigm
   logical :: QUIET

   n = HMMSTR%n
   k = HMMSTR%k
   u = HMMSTR%u
   allocate(delta(HMMSTR%n, seqres_len))
   allocate(upsilon(HMMSTR%n, seqres_len))
   allocate(Q(seqres_len))
   if(allocated(paradigm))   deallocate(paradigm)

     allocate(paradigm(4,seqres_len))
   delta = 0
   upsilon = 0
   Q = 0
   do titr = 1, seqres_len
     trama = ramaseq(titr)
     tss = ssseq(titr)
     ttm = tmseq(titr)
     t_prof = profile(:, titr)
     do nitr = k+u+1, n
        b_t1 = flagrat_brdc_prod( nitr, t_prof, HMMSTR, trama, tss, ttm, flags)
        if ( titr == 1 ) then ! t = 1 uses priors, a special case to enter model
           do nnitr = u+1,u+k
                   if ( HMMSTR%a(nnitr,nitr) > 0.000) then
                      a_ij = HMMSTR%a(nnitr, nitr)
                      delta(nitr, titr) = delta(nitr,titr) + a_ij * b_t1 * HMMSTR%priors(nnitr)
                   endif
           enddo
           a_ij = HMMSTR%priors(nitr)
           delta(nitr, titr) = delta(nitr, titr) + a_ij * b_t1
        else ! t > 1
           max_loc = 0
           val = 0
           do nitr_2 = 1, HMMSTR%n
               if ( nitr_2 /= nitr ) then
                   a_ij = Viterbi_naught_traverse(nitr_2, nitr, titr, HMMSTR, intrs)
                   val = delta(nitr_2, titr-1) * a_ij
               endif
               if ( val > max_loc ) then
                   max_loc = val
                   loc = nitr_2
               endif
           enddo
           upsilon(nitr,titr) = loc
           delta(nitr,titr) = max_loc * b_t1
        endif
     enddo
     if ( .false. ) then !titr > 1 ) then
             ! FIX THE NAUGHT STATE:
             max_loc = 0
             do nitr = u,u+k ! TLB 9/26 Viterbi Naught Recurse
                     do iitr =1, intrs(nitr)%c(0) ! loop over all states incoming to naught state
                         ind = intrs(nitr)%c(iitr)
                         if ( ind < u+k ) then
                             val = Viterbi_naught_recurse(ind, titr, delta, HMMSTR, intrs) * HMMSTR%a(ind, nitr)
                         else
                             val = delta(ind, titr-1) * HMMSTR%a(ind,nitr) 
                             ! states incoming are at t-1 relative to naught
                         endif
                         if ( val > max_loc) then
                                 max_loc = val
                                 loc = ind
                         endif
                     enddo
                     delta(nitr, titr) = max_loc
                     upsilon(nitr, titr) = loc
             enddo
      endif
   enddo

   allocate(out_seq(seqres_len))
   allocate(out_ssseq(seqres_len))
   allocate(out_ramaseq(seqres_len))
   allocate(out_tmseq(seqres_len))
   ! get path through state sequence backtrack
   do titr=seqres_len,1,-1
      !write(*,*) titr
      if ( titr == seqres_len ) then
           max_loc = 0
           do nitr = u+k+1, n ! fix so we can never choose to end in the naught-state
              if ( delta(nitr,titr) > max_loc ) then
                   max_loc = delta(nitr,titr)
                   Q(titr) = nitr
              endif
           enddo
      else

          Q(titr) = upsilon(Q(titr+1),titr+1)
      endif
      rama_prof = HMMSTR%RAMA(Q(titr),:)
      ss_prof   = HMMSTR%SS(Q(titr),:)
      aa_prof   = HMMSTR%AA(Q(titr),:)
      tm_prof   = HMMSTR%TM(Q(titr),:)
      call voting(HMMSTR, out_seq, out_tmseq, out_ssseq, out_ramaseq, aa_prof, tm_prof, ss_prof, rama_prof, titr)
   enddo
   
   if ( QUIET .eqv. .false. ) then
           write(*,*) "GROUND TRUTH"
           write(*,*) seq(1:seqres_len)
           write(*,*) ssseq(1:seqres_len)
           write(*,*) ramaseq(1:seqres_len)
           write(*,*) tmseq(1:seqres_len)
           write(*,*) "VITERBI PREDICTIONS"
           write(*,*) out_seq(1:seqres_len)
           write(*,*) out_ssseq(1:seqres_len)
           write(*,*) out_ramaseq(1:seqres_len)
           write(*,*) out_tmseq(1:seqres_len)
           !write(*,*) Q(1:seqres_len)
   endif
    paradigm(1,:) = out_seq(:)
    paradigm(2,:) = out_ssseq(:)
    paradigm(3,:) = out_ramaseq(:)
    paradigm(4,:) = out_tmseq(:)

   if(allocated(out_seq)) deallocate(out_seq)
   if(allocated(out_ssseq)) deallocate(out_ssseq)
   if(allocated(out_ramaseq)) deallocate(out_ramaseq)
   if(allocated(out_tmseq)) deallocate(out_tmseq)

   !stop

  end subroutine viterbi
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TLB 9/8/2021 moved from hmm_gamma.f95 to hmmstrtm.f95 !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_mda_acc(mda_weight_accs)
     implicit none
     real,dimension(:,:),intent(in) :: mda_weight_accs
     real :: sum_weighted_accs
     integer :: oitr
     real :: acc
     sum_weighted_accs = 0
     do oitr=1,4999
         !write(*,*) mda_weight_accs(oitr,:)
         if ( sum(mda_weight_accs(oitr,:)) < 0.0001 ) exit
         sum_weighted_accs = sum_weighted_accs + mda_weight_accs(oitr,1) * mda_weight_accs(oitr,2)
     enddo
     acc = sum_weighted_accs / sum(mda_weight_accs(:,1))
     write(*,*) "Database mda accuracy = ", acc
  end subroutine get_mda_acc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TLB 9/8/2021 moved from hmm_gamma.f95 to hmmstrtm.f95 !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine get_mda(profile, paradigm, mda, t, t_mda)
    implicit none
    character(len=1),dimension(:,:),intent(in)   :: profile
    character(len=1),dimension(:,:),intent(in)  :: paradigm ! ( aa_seq, ss_seq, rama_seq, tm_seq ) 
    integer,dimension(:),allocatable,intent(inout) :: t_mda
    character(len=1),dimension(:),allocatable :: t_mda_str
    character(len=1),dimension(:),allocatable :: align
    integer :: correct ,incorrect
    logical :: weighted
    integer, intent(in) :: t
    integer :: titr, titr2
    logical :: any_incorrect
    character(len=1) :: rama_char, rama_ground
    character(len=11) :: rama_s! ="HGBEdbeLlxc"
    character(len=1),dimension(:),allocatable :: seq1, ssseq1, ramaseq1, tmseq1
    real(8),dimension(2),intent(inout) :: mda

    mda = 0
    !write(*,*) "SHAPE OF PROFILE  = ", shape(profile)
    if(allocated(seq1))   deallocate(seq1)
    if(allocated(ssseq1))   deallocate(ssseq1)
    if(allocated(ramaseq1))   deallocate(ramaseq1)
    if(allocated(tmseq1))   deallocate(tmseq1)
    if(allocated(t_mda))   deallocate(t_mda)
    if(allocated(t_mda_str))   deallocate(t_mda_str)
    if(allocated(align))   deallocate(align)
    allocate(seq1(t))
    allocate(ssseq1(t))
    allocate(ramaseq1(t))
    allocate(tmseq1(t))
    allocate(t_mda(t))
    allocate(t_mda_str(t))
    allocate(align(t))
    seq1     = paradigm(1,:)
    ssseq1   = paradigm(2,:)
    ramaseq1 = paradigm(3,:)
    tmseq1   = paradigm(4,:)
    !write(*,*) "LOADED PARADIGM"
    t_mda = 0
    do titr=1, t-7
        any_incorrect = .false.
        rama_char = ramaseq1(titr)
        rama_ground = profile(3,titr)
        if (rama_comp(rama_char, rama_ground)) then
                align(titr) = "0"
        else
                align(titr) = "1"
        endif
        any_incorrect = ( any_incorrect .or. rama_comp(rama_char, rama_ground) )
        do titr2=titr+1, titr+7
            rama_char = ramaseq1(titr2)
            rama_ground = profile(3,titr2)
            if (titr == t -7) then
                    if ( rama_comp(rama_char, rama_ground) ) then
                            align(titr2) = "0"
                    else
                            align(titr2) = "1"
                    endif
            endif
            any_incorrect = ( any_incorrect .or. rama_comp(rama_char, rama_ground) )
        enddo
        if (any_incorrect) then
            if (t_mda(titr) /= 1) then
                t_mda(titr) = 0
                t_mda_str(titr) = "0"
            endif
        else
             t_mda(titr:titr+7) = 1
             t_mda_str(titr:titr+7) = "1"
        endif
        if ( titr == t-7 ) then
            if ( any_incorrect ) then
                t_mda(titr:titr+7) = 0
                t_mda_str(titr:titr+7) = "0"
            else
                t_mda(titr:titr+7) = 1
                t_mda_str(titr:titr+7) = "1"
            endif
        endif
    enddo
    rama_s ="HGBEdbeLlxc"
    do titr = 1, len(rama_s)
       do titr2 = 1, len(rama_s)
            rama_char = rama_s(titr:titr)
            rama_ground = rama_s(titr2:titr2)
            !write(*,*) rama_s(titr:titr), rama_s(titr2:titr2), rama_comp(rama_char, rama_ground)
       enddo
    enddo
    if ( .false. ) then  ! USE FOR DEBUGGING MDA METRIC
            write(*,*) " AA SEQUENCE = ", profile(1,:t)
            write(*,*) " SS SEQUENCE = ", profile(2,:t)
            write(*,*) "GROUND TRUTH = ", profile(3,:t)
            write(*,*) "PREDICTION   = ", ramaseq1(:t)
            write(*,*) "RAMA_COMP    = ", align(:t)
            write(*,*) "MDA          = ", t_mda_str(:t)
            !write(*,*) "ACCURACY = ", (sum(MDA) / size(MDA))
            write(*,*) t ! weight accuracy by length of this chain
            write(*,*) (sum(t_mda) / size(t_mda))
    endif
    mda(1) =  size(t_mda)
    mda(2) =  real((real(sum(t_mda)) / real(size(t_mda))))
    !write(*,*) mda
    if(allocated(seq1))   deallocate(seq1)
    if(allocated(ssseq1))   deallocate(ssseq1)
    if(allocated(ramaseq1))   deallocate(ramaseq1)
    if(allocated(tmseq1))   deallocate(tmseq1)
    !if(allocated(t_mda))   deallocate(t_mda)
    if(allocated(t_mda_str))   deallocate(t_mda_str)
    if(allocated(align))   deallocate(align)

end subroutine get_mda

  logical function rama_comp(rama_char, rama_ground) result(ret)
        implicit none
        character(len=1), intent(in) :: rama_char, rama_ground
        character(len=11) :: rama_string! ="HGBEdbeLlxc"
        real,dimension(11,2)::ramacen
        integer :: ind1, ind2
        real :: dely, delx
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
        rama_string ="HGBEdbeLlxc"
        ind1 = index(rama_string, rama_char)
        ind2 = index(rama_string, rama_ground)
        if ( ind1 == ind2 .or. ind1 == 0 .or. ind2 == 0) then
                ret = .false.
        else
                dely = ramacen(ind1,2) -ramacen(ind2,2)
                if (dely < -180) dely = dely + 360
                if (dely > 180) dely = dely - 360
                delx = ramacen(ind1,1) -ramacen(ind2,1)
                if (delx < -180) delx = delx + 360
                if (delx > 180) delx = delx - 360
                ret = .true.
                if (((delx ** 2 + dely ** 2) ** 0.5) < 120) then
                    ret = .false.
                endif
        endif
end function rama_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hmm
