module INIT
    !!! Module qui contient les parametres physiques et numériques
    !!! Il faut faire appel à init_params avant toute manip dans le main

    implicit NONE
    integer, PARAMETER :: pr = 8
    real(pr), parameter ::  pi = acos(-1._pr)
    real(8) :: a, dt, dx, Tf
    integer :: ordre, Ne, Nt, numPoints, ii, jj
    integer, dimension(:,:), allocatable ::pascal_triangle


contains

    !! --- lecture des parametres dans params.txt
    !! --- + generation des poids et noeuds pour 
    !! --- la formule de quadratures
    subroutine init_params(GAUSS_LEGENDRE)
        implicit NONE
        logical, intent(in) :: GAUSS_LEGENDRE
        character(len = 50) :: temp, command

        open(unit = 10, file = 'params.txt')
        read(10,*) a, dx, dt, Tf, ordre
        close(10)

        !! -- Nombre de mailles d'espace 
        Ne = int(1._pr/dx)
        print*,pi

        !! -- Nombre de mailles de temps
        Nt = int(Tf/dt)

        !!-- les pas
        dx = 1._pr/Ne
        dt = Tf/Nt

        !! -- nombre de points dans la formule de gauss
        !! -- on utilise 64 pour verification des calculs
        !! -- On utlisera le nombre min de points requis ulterieurement
        numPoints = 64

        if (GAUSS_LEGENDRE) then
            write(temp,*) numPoints
            print*, 'Génération du fichier n'//trim(adjustl(temp))//'.txt necessaire pour la formule de quadratures.'
            command = 'python3 gauss_leg_coeff.py '//trim(adjustl(temp))
            call system(command)
        endif

        !! initialisation du triangle de pascal
        allocate(pascal_triangle(0:ordre, 0:ordre))
        pascal_triangle = 0
        pascal_triangle(:,0) = 1
        do jj = 1, ordre
            do ii = 1, ordre
                pascal_triangle(ii,jj) = pascal_triangle(ii-1 ,jj) + pascal_triangle(ii-1,jj-1)
            enddo
        enddo


    end subroutine init_params

    subroutine print_params()
        implicit none
        
        print*, "******************* PARAMETRES *************************"
        
        print*, "dx = ", dx
        print*, "dt = ", dt
        print*, "a =  ", a
        print*, "Tf = ", Tf
        print*, "ordre = ", ordre
        print*, "Ne = ", Ne
        print*, "Nt = ", Nt
        print*, "adt/dx = ", a*dt/dx
        
        print*, "*********************************************************"
        
    end subroutine print_params


    !! --- Condition Initiale
    function u_init(x)
        implicit none

        real(8):: x, u_init
        u_init = cos(2*pi*x)
        !u_init = 1._pr

    end function

    !! --- Condition du bord gauche
    function u_g(t)
        implicit NONE
        real(8) :: u_g, t
        !u_g = t
        u_g = 1._pr
    end function u_g



end module INIT