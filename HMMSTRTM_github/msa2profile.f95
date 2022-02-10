program msa2profile
  !!
  !! read a MSA in FASTA format and write a 
  !! profile (position-specific amino acid probability distribution)
  !! 
  implicit none
  integer,parameter :: LONGLINE=2000, ounit=12, munit=13
  character(len=LONGLINE) :: aline, msafile, outfile, masterseq
  character(len=LONGLINE),dimension(:),allocatable :: myseq
  character(len=LONGLINE),dimension(:),allocatable :: label
  character(len=20),parameter :: amino="ACDEFGHIKLMNPQRSTVWY"
  character(len=6),parameter :: removed="delete"
  real,parameter :: ranpid=0.25
  real :: x,pid,braa,brbb,jarg,pseudo
  real,dimension(:),allocatable :: age
  real,dimension(:,:),allocatable :: dmatrix
  real,dimension(:),allocatable :: seqweight,profile
  integer :: ntaxa, mtaxa, i,j,aa,bb,n,nid,iosa,nres,ios,nongap,iaa,ires
  real,dimension(20),parameter :: Fi= (/ 0.08279,0.01937,0.05855,0.05992,0.04014, &
      0.08089,0.02275, 0.05552,0.05959,0.08020,0.02070,0.04729,0.04599, &
      0.03728,0.04640, 0.06246,0.05888,0.06866,0.01507,0.03756/)
  !!
  !!--------------------------------------------------------------------------
  pseudo=0.001
  jarg = command_argument_count()
  if (jarg < 2) then
    write(*,*) 'Usage: xmsa2profile msafile output_profile [pseudocounts] > profile'
    write(*,*) 'This program takes in a fasta-fomat MSA and reduces in to a profile'
    write(*,*) 'of the first sequence in the MSA. You must first assure that the'
    write(*,*) 'first sequence is the one you want to make a profile of. Then,'
    write(*,*) 'build a MSA by using BLAST and/or UGENE or any multiple sequence'
    write(*,*) 'tool. Prune the MSA if desired using prunefasta.f90, '
    write(*,*) 'but be sure to eliminate any all-gap columns before running this '
    write(*,*) 'program. If using msa2profile.f90 to prepare a profile for'
    write(*,*) 'hmmstrgca.f95, then extract the sequence from the PDB file'
    write(*,*) 'before aligning, and move it to the top before saving the alignment'
    write(*,*) 'msa2profile.f90 will remove any columns that are gaps in'
    write(*,*) 'the master sequence (sequence 1).'
    stop 'msa2profile v.  Wed May 16 16:52:25 EDT 2018'
  endif
  call get_command_argument(1,msafile)
  call get_command_argument(2,outfile)
  if (jarg >= 3) then
    call get_command_argument(3,aline)
    read(aline,*,iostat=ios) pseudo
    if (ios/=0) stop 'Bad format for pseudocounts. Should be a float 0<x<1.0'
    if (pseudo<0.0) stop 'Bad format for pseudocounts. Should be a float 0<x<1.0'
    write(*,'("Pseudocounts = ",f9.3)') pseudo
    if (pseudo>1.0) write(*,'("Are you sure you want pseudocounts to be that big? ")')
  endif
  open(11,file=msafile,status='old',iostat=ios) 
  if (ios/=0) stop 'ERROR opening msa file'
  open(ounit,file=outfile,status='replace',iostat=ios) 
  if (ios/=0) stop 'ERROR opening output file'
  !! Read MSA in fasta format, count sequences
  mtaxa = 0
  do
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:1)==">") mtaxa = mtaxa + 1
  enddo
  allocate(dmatrix(mtaxa,mtaxa),myseq(mtaxa),age(mtaxa),label(mtaxa),stat=ios)
  if (ios/=0) stop 'msa2profile:: ERROR: could not allocate dmatrix'
  allocate(seqweight(mtaxa),profile(20),stat=ios)
  if (ios/=0) stop 'msa2profile:: ERROR: could not allocate seqweight()'
  profile = 0.0
  seqweight=0.0
  !!--------------------------------------------------------------------------
  !! Read MSA in fasta format
  rewind(11)
  ntaxa  = 0
  nres = 0
  do
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:1)==">") then
      if (ntaxa==1) then
        masterseq = " "
        masterseq = trim(myseq(1))
      endif
      ntaxa = ntaxa + 1
      label(ntaxa) = trim(aline(2:))
      myseq(ntaxa) = " "
      cycle
    endif
    myseq(ntaxa) = trim(myseq(ntaxa)) // trim(aline)
    !! nres is the total number of columns, not the sequence length.
    !! Should be the same for all entries in the MSA. 
    nres = max(nres,len_trim(myseq(ntaxa)))
  enddo
  masterseq = trim(myseq(1))
  close(11)
  !!---------------------------------------------------------------------------
  !! trim the MSA to the master seq
  do i=1,ntaxa
    ires = 0
    do j=1,nres
      iaa = index(amino,masterseq(j:j))
      if (iaa/=0.and.iaa<=20) then
        ires = ires + 1
        myseq(i)(ires:ires) = myseq(i)(j:j)
      endif
    enddo
    myseq(i)(ires+1:nres) = " "
  enddo
  nres = len_trim(myseq(1))
  !!---------------------------------------------------------------------------
  !! write out trimmed MSA
  open(munit,file=trim(outfile)//'.fa',status='replace',iostat=ios) 
  if (ios/=0) stop 'ERROR opening output fasta file'
  do i=1,ntaxa
    write(munit,'(">",a)') adjustl(trim(label(i)))
    write(munit,'(a)') trim(myseq(i))
  enddo
  close(munit)
  !!---------------------------------------------------------------------------
  !do i=1,ntaxa
  !  write(*,*) i,trim(label(i)),trim(myseq(i))
  !enddo
  if ( ntaxa > 1 ) then 
          dmatrix = 0.0
          do i=1,ntaxa
            do j=1,i-1
              pid = getpid(myseq(i),myseq(j))
              if (pid <= ranpid) then
                 pid = ranpid + 0.01
                !! diagnostic
                !write(0,*) "WARNING: %id below random ",i,j
                x = 0!huge(x)
              endif
              !! Convert to Jukes-Cantor
              x = -log((pid-ranpid)/(1-ranpid))
              dmatrix(i,j) = x
              dmatrix(j,i) = x
            enddo
          enddo
          !!--------------------------------------------------------------------------
          !! Get sequence weights using sum-of-distance method
          do i=1,ntaxa
            seqweight(i) = sum(dmatrix(i,:))
          enddo
          seqweight = seqweight/sum(seqweight(:))
  else
        seqweight(1) = 1
  endif

  !!--------------------------------------------------------------------------
  !! Calculate profile with pseudocounts
  write(ounit,'("#PROFILE Col# AA_of_seq_1 Prob_A C D E F G H I K L M N P Q R S T V W Y")')
  ires = 0
  do j=1,nres   !! columns
    nongap=0
    profile = Fi*pseudo
    do i=1,ntaxa
      iaa = index(amino,myseq(i)(j:j))
      if (iaa/=0.and.iaa<=20) then
        nongap=nongap+1
        profile(iaa) = profile(iaa) + seqweight(i)
      endif
    enddo
    if (nongap==0) then
      write(*,'("WARNING: no amino acid characters found for column ",i5)') j
      stop 'View the alignment in UGENE and run edit/remove columns of gaps'
    endif
    ires = ires + 1
    profile = profile/sum(profile)
    !! Write out profile. Col# AA_of_seq_1 Prob_A C D E F G H I K L M N P Q R S T V W Y
    write(ounit,'("PROFILE:",i5,1x,a1,20f8.5)') ires,myseq(1)(j:j),profile
  enddo
  close(ounit)

  !! 
CONTAINS
    !!-----------------------------------------------------
    ! end subroutine findmin
    !!-----------------------------------------------------
    real function getpid(seqaa,seqbb)
      implicit none
      character(len=LONGLINE),intent(in) :: seqaa, seqbb
      character(len=20),parameter :: amino="ACDEFGHIKLMNPQRSTVWY"
      integer :: i,n,nid
      real :: x=0
      n = 0
      nid = 0
      do i=1,len_trim(seqaa)
        if (index(amino,seqaa(i:i))==0) cycle
        if (index(amino,seqbb(i:i))==0) cycle
        if (seqaa(i:i)==seqbb(i:i)) nid = nid + 1
        n = n + 1
     enddo
     if (n>0) x = real(nid)/real(n)
     getpid = x
    end function getpid
end program msa2profile
