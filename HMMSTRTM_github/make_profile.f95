! test_profile_io.f95

program test_profile_io
   use profile_io
    
    character (len=30) :: drctfile
    character (len=4) :: code
    integer :: sitr
    
    iarg = command_argument_count()
    call getarg(1, drctfile)
    call getarg(2, code)

    this_profile = get_profile(drctfile, 0)
    sitr = 1
    do while (this_profile%last_record /= -1 )
        if ( this_profile%code == code ) then 
            write(*,*) "WRITING OUT PROFILE FILE FOR CHAIN = " // this_profile%code // this_profile%chain // " ."
            call write_profile(this_profile, './tmp/'//this_profile%code//this_profile%chain//'.profile') 
        endif
        sitr = sitr + 1
        this_profile = get_profile(drctfile, this_profile%last_record)
    enddo

end program test_profile_io
