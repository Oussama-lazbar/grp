program test
    use FUNC
    implicit none
    ! integer(kind = 8) :: test1, n = 22
    ! test1 = factorial2(n)
    ! print*, test1
    call init_params(.FALSE.)
    call test_legendre(14)
    deallocate(pascal_triangle)
contains

    subroutine test_fact(n)
        implicit none 
        integer, intent(in) :: n
        integer :: i
        integer :: table(0:n)
        logical :: test = .TRUE.
        do i = 0, n
            table(i) = factorial(i)
        end do
        
        do i = 1, n
            test = (test .AND. (table(i) == table(i-1) * i))
            if (test) then
                print*, 'test passed succesfully !'
            else
                print*, 'error at n = ', i
            end if
        end do
    end  subroutine test_fact
    
    function factorial2(n)
        implicit none
        integer(kind = 8), intent(in) :: n
        integer(kind = 8) :: factorial2, i, temp
        temp = 1
        if (n == 0) then
            factorial2 = 1
        else
            do i = 1,n
                temp = temp*i
            enddo
            factorial2 = temp
        endif
        
    endfunction factorial2

    function prod_legendre(p,q, x)
        implicit NONE
        integer :: p,q 
        real(pr) :: x, prod_legendre
        prod_legendre = legendre(p,x) * legendre(q,x)
    end function

    subroutine test_legendre(n)
        implicit NONE
        integer, intent(in) :: n
        integer :: p,q 
        real(pr) ::  res 
        do p = 0, n
            do q = 0, n
                call integrale(p,q,prod_legendre, 64, res)
                print*, 'integrale du  produit des pol de legendre ', p,' et ', q, ' = ', res 
                if (p == q) then
                    print*, 'la norme du polynome de legendre de degre ' ,p,' = ', 2._pr/(2._pr * p + 1)
                    print*,''
                endif
            enddo
        enddo

    end subroutine test_legendre
    
        
    
    
end program test