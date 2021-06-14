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

        !! -- Nombre de mailles de temps
        Nt = int(Tf/dt)

        !!-- les pas
        dx = 1._pr/Ne
        dt = Tf/Nt

        !! -- nombre de points dans la formule de gauss
        !! -- on utilise 64 pour verification des calculs
        !! -- On utlisera le nombre min de points requis ulterieurement
        if (2 * (ordre/2) == ordre) then 
            numPoints = ((ordre**2 + 1)/2) + 2
        else
            if (ordre == 1) then
                numPoints = 4
            else
                numPoints = ((ordre**2 + 1)/2) + 1
            endif
        end if
        !numPoints = 64

        if (GAUSS_LEGENDRE) then
            write(temp,*) numPoints
            print*, 'Génération du fichier n'//trim(adjustl(temp))//'.txt necessaire pour la formule de quadratures.'
            command = 'python3 gauss_leg_coeff.py '//trim(adjustl(temp))
            call system(command)
            print*, "Fichier généré."
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

    function u_exacte(t,x)
        implicit NONE
        real(pr) :: t,x, u_exacte
        if (x  - a*t < 0) then
            u_exacte = 1._pr
        else
            u_exacte = u_init(x  - a*t)
        endif
    endfunction u_exacte

    function u_moy(t,x)
        implicit none
        real(pr) :: t,x,u_moy, temp
        integer :: i
        i = int(x/dx)
        if (x  - a*t < 0) then
            u_moy = 1._pr
        else
            temp = 2._pr*pi*(real(i, pr)* dx - a*t)
            u_moy = (1._pr/(2*pi*dx)) * (sin(temp + 2._pr*pi*dx) - sin(temp) )
            !print* , sin(2*pi*((i+1)*dx - a*t)) - sin(2*pi*(i*dx - a*t))
        endif
    end function u_moy

    function u_moy2(t,i)
        implicit none
        real(pr) :: t,u_moy2, temp
        integer :: i
        if (real(i,pr)*dx  - a*t < 0) then
            u_moy2 = 1._pr
        else
            temp = 2._pr*pi*(real(i, pr)* dx - a*t)
            u_moy2 = (1._pr/(2*pi*dx)) * (sin(temp + 2._pr*pi*dx) - sin(temp) )
        endif
    end function u_moy2




end module INIT