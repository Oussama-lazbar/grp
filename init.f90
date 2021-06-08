module INIT
    !!! Module qui contient les parametres physiques et numériques
    !!! Il faut faire appel à init_params avant toute manip dans le main

    implicit NONE
    integer, PARAMETER :: pr = 8
    real(pr), parameter :: pi = acos(1._pr) 
    REAL(8) :: a, dt, dx, Tf
    integer :: ordre, Ne, Nt, numPoints


contains

    !! --- lecture des parametres dans params.txt
    subroutine init_params()
        implicit NONE

        open(unit = 10, file = 'params.txt')
        read(10,*) a, dx, dt, Tf, ordre
        close(10)
        Ne = int(1._pr/dx)
        Nt = int(Tf/dt)
        dx = 1._pr/Ne
        dt = Tf/Nt
        numPoints = 64

    end subroutine init_params

    subroutine print_params()
        implicit none
        character(len = 50) :: temp, command
        print*, "******************* PARAMETRES *************************"
        
        print*, "dx = ", dx
        print*, "dt = ", dt
        print*, "a =  ", a
        print*, "Tf = ", Tf
        print*, "ordre = ", ordre
        print*, "Ne = ", Ne
        print*, "Nt = ", Nt
        
        print*, "*********************************************************"
        write(temp,*) numPoints
        command = 'python3 gauss_leg_coeff.py '//trim(adjustl(temp))
        call system(command)
    end subroutine print_params


    !! --- Condition Initiale
    function u_init(x)
        implicit none

        real(8):: x, u_init
        !u_init = cos(2*pi*x)
        u_init = 1._pr

    end function

    !! --- Condition du bord gauche
    function u_g(t)
        implicit NONE
        real(8) :: u_g, t
        u_g = 1._pr
    end function u_g



    
end module INIT