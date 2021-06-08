program grp
    use FUNC
    
    implicit none
    real(pr), dimension(:),allocatable :: alpha, alpha1, beta, beta1
    real(pr), dimension(:,:), allocatable :: L1, L2, L3
    integer :: i,n,j
    real(pr) :: res

    !!! --- initialisation de params physiques et numériques
    call init_params()
    call print_params()
    !!! --- allocation memoire
    allocate(alpha(0:Ne * (ordre + 1) - 1), alpha1(0:Ne * (ordre + 1) -1), beta(0 : ordre), beta1(0 : ordre))
    allocate(L1(0 : ordre , 0 : ordre), L2(0 : ordre , 0 : ordre), L3(0 : ordre , 0 : ordre))

    !call init_alpha(alpha)
    !print*,alpha(0:ordre+13 )
    !call fill(L1, L2, L3)
    !print*,L1
    ! call init_beta(0,beta)
    ! print*, beta

    ! do n = 0, Nt-1
    !     call init_beta(n,beta)
    !     !print*,beta
    !     do i = 0, Ne - 1
    !         !! mise à jour alpha
    !         alpha1(i*(ordre+1) : (i+1)*(ordre +1) - 1) =(2/dx)* matmul(L1,beta )
    !         do j =0, ordre
    !             alpha1(i*(ordre+1) + j) = ((2*j+1)/2._pr) * alpha1(i*(ordre+1) + j)
    !             !print*, alpha1(j)
    !         end do
    !         !! mise à jour beta
    !         beta1 = (2/dt)* (matmul(L3,beta ) + matmul(L2,alpha(i*(ordre+1) : (i+1)*(ordre +1) - 1) ))
    !         do j =0, ordre
    !             beta1(j) = ((2*j+1)/2._pr) * beta1(j)
    !         end do
    !         beta = beta1
    !     end do 
    !     alpha = alpha1
    ! end do
    
    !print*, alpha(0:ordre+1)
    !! disiplay solution
    ! open(unit = 11, file = 'sol.txt')
    ! do i = 0, Ne - 1
    !     write(11,*) i*dx , alpha(i*(ordre + 1)+1)
    ! end do
    ! close(11)

    deallocate(alpha, alpha1, beta, beta1, L1, L2, L3)
    call integrale(1,14, test1, 64, res)
    print*,res
    ! call integrale(1,2, test1, int(ordre/2) + 2, res)
    ! print*,res
    ! call integrale(8,7, test1, int(ordre/2) + 2, res)
    ! print*,res

    contains
    function test1(p,q,x)
        implicit NONE

        integer :: p,q
        real(pr) :: x, test1
        test1 = legendre(16,x)
    end function test1
end program grp