program grp
    use FUNC
    use INIT
    use QUAD
    
    implicit none
    real(pr) :: timing, start_time, end_time, time = 0, error_1=0._pr, error_2 = 0._pr
    real(pr), dimension(:),allocatable :: alpha, alpha1, beta, beta1
    real(pr), dimension(:,:), allocatable :: L1, L2, L3
    integer :: i,n,j

    !!! --- initialisation de params physiques et numériques
    !print*, "Initialisation des parametres."
    call init_params(.TRUE.)
    call print_params()

    !!! --- allocation memoire
    !print*, "Allocation mémoire"
    allocate(alpha(0:Ne * (ordre + 1) - 1), alpha1(0:Ne * (ordre + 1) -1), beta(0 : ordre), beta1(0 : ordre))
    allocate(L1(0 : ordre , 0 : ordre), L2(0 : ordre , 0 : ordre), L3(0 : ordre , 0 : ordre))
    !print*, "Fin d'allocation mémoire."

    ! -- initialisation de alpha à l'instant 0
    call init_alpha(alpha)
    alpha1 = 0._pr

    ! -- remplissage des matrices L_{i}
    call init_L(L1, L2, L3)
    !print*, "Fin de l'initialisation."
    !print*, "CALCUL"

    call cpu_time(start_time)
    
    do n = 0, Nt-1
        time = time + dt
        ! -- initialisation de beta
        call init_beta(n,beta)
        !print*, beta
        
        do i = 0, Ne - 1

            !!-----------mise à jour alpha
            alpha1(i*(ordre+1) : (i+1)*(ordre +1) - 1) =(2._pr/dx)* matmul(L1,beta )
            ! -- On multiplie par l'inverse de M.
            do j =0, ordre
                alpha1(i*(ordre+1) + j) = (j+0.5_pr) * alpha1(i*(ordre+1) + j)
            end do
            
            !!----------mise à jour beta
            beta1 = (2._pr/dt)* (matmul(L3,beta ) + matmul(L2,alpha(i*(ordre+1) : (i+1)*(ordre +1) - 1) ))
            ! -- On multiplie par l'inverse de M.
            do j =0, ordre
                beta1(j) = (j+0.5_pr) * beta1(j)
            end do

            beta = beta1

        end do 
        alpha = alpha1
    end do

    call cpu_time(end_time)
    !print*, "FIN DE CALCUL"
    timing = end_time - start_time 
    print*, "Temps de calcul = ", timing
    !print*, alpha

    !! -- On affiche la valeur moyenne de la solution par maille : elle correspond à (alpha_{i}^{n})_{0}
    open(unit = 11, file = 'sol.txt')
    do i = 0, Ne - 1
        write(11,*) i*dx ,alpha(i*(ordre + 1)), u_moy2(time, i)
    end do
    close(11)

    call Err_1(error_1)
    call Err_2(error_2)
    print*, "dx = ", dx, "   ||   erreur_1 = ",error_1, "   ||   erreur_2 = ",error_2 


    deallocate(alpha, alpha1, beta, beta1, L1, L2, L3, pascal_triangle)



    open(unit=12, file="error_data.dat", position='append')
    write(12,*) dx, error_1, error_2
    close(12)

contains

    function sol_num(i,x)
        implicit none 
        real(pr) :: x, sol_num ,temp , xi
        integer :: i ,j 

        xi = (i + 0.5_pr)*dx
        temp = (x - xi)/(dx/2._pr) 
        sol_num = 0._pr

        do j = 0, ordre
            sol_num = sol_num + alpha(i*(ordre + 1) + j)*legendre(j,temp)
        end do
    end function sol_num

    function e1(p,q,x)
        !! --- p et q sont des params entiers.  
        !! -- Ils ne sont la rien que pour pouvoir 
        !! -- utiliser la fonction integrale
        implicit none 
        real(pr) :: x, e1, temp, xi
        integer :: p,q,i
        e1 = 0._pr
        
        do i = 0, Ne-1
            xi = (i + 0.5_pr)*dx
            temp = (dx/2._pr)*x + xi
            e1 = e1 + abs(sol_num(i,temp) - u_exacte(time, temp))
        enddo
        e1 = (dx/2._pr)*e1
    end function e1 

    function e2(p,q,x)
        !! --- p et q sont des params entiers.  
        !! -- Ils ne sont la rien que pour pouvoir 
        !! -- utiliser la fonction integrale
        implicit none 
        real(pr) :: x, e2, temp, xi
        integer :: p,q , i 
        e2 = 0._pr
        
        do i = 0, Ne-1
            xi = (i + 0.5_pr)*dx
            temp = (dx/2._pr)*x + xi
            e2 = e2 + abs(sol_num(i,temp) - u_exacte(time, temp))**2
        enddo
        e2 = (dx/2._pr)*e2

    end function e2

    subroutine Err_1(erreur_1)
        implicit none 
        real(pr), intent(inout) :: erreur_1
        call integrale(0, 0, e1, numPoints, erreur_1)
    end  subroutine Err_1

    subroutine Err_2(erreur_2)
        implicit none 
        real(pr), intent(inout) :: erreur_2
        call integrale(0, 0, e2, numPoints, erreur_2)
        erreur_2 = sqrt(erreur_2)
    end  subroutine Err_2

    
end program grp