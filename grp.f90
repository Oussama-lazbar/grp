program grp
    use FUNC
    
    implicit none
    real(pr), dimension(:),allocatable :: alpha, alpha1, beta, beta1
    real(pr), dimension(:,:), allocatable :: L1, L2, L3
    integer :: i,n,j

    !!! --- initialisation de params physiques et numériques
    call init_params(.TRUE.)
    call print_params()

    !!! --- allocation memoire
    allocate(alpha(0:Ne * (ordre + 1) - 1), alpha1(0:Ne * (ordre + 1) -1), beta(0 : ordre), beta1(0 : ordre))
    allocate(L1(0 : ordre , 0 : ordre), L2(0 : ordre , 0 : ordre), L3(0 : ordre , 0 : ordre))

    ! -- initialisation de alpha à l'instant 0
    call init_alpha(alpha)
    
    ! -- remplissage des matrices L_{i}
    call init_L(L1, L2, L3)

    !! -- affichage du triangle de pascal
    ! do i = 0, ordre
    !     print*, pascal_triangle(i,:)
    ! end do

    
    do n = 0, Nt-1
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


    
    !! -- On affiche la valeur moyenne de la solution par maille : elle correspond à (alpha_{i}^{n})_{0}
    open(unit = 11, file = 'sol.txt')
    do i = 0, Ne - 1
        write(11,*) i*dx , alpha(i*(ordre + 1))
    end do
    close(11)

    deallocate(alpha, alpha1, beta, beta1, L1, L2, L3, pascal_triangle)
    

end program grp