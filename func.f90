module FUNC
    use QUAD
    use INIT
    implicit none

contains

        function factorial(n)
            implicit none

            integer, intent(in) :: n
            integer :: factorial, i

            factorial = 1
            if ( n > 1) then
                do i = 2, n
                    factorial = factorial * i
                enddo
            endif
            
        endfunction factorial


        function binomial(n,k)
            implicit none

            integer, intent(in) :: n,k
            integer :: binomial
            binomial = factorial(n)/(factorial(k)* factorial(n-k))

        endfunction binomial


        !!! ---- nième polynome de legendre
        function Legendre(n,x)
            implicit none

            integer, intent(in) :: n
            real(pr), intent(in) :: x
            real(pr) :: Legendre
            integer :: i

            Legendre = 0
            do i = 0, n
                Legendre = Legendre + (binomial(n,i)**2)*((x-1)**(n-i))*((x+1)**i)
            enddo
            Legendre = Legendre/(2**n)

        end function Legendre
        

        function f1(p,q,x)
            implicit none

            integer :: p,q
            real(pr) :: x
            real(pr) :: f1 ,  temp
            temp = 1._pr - (dx/(a*dt))*(1+x)
            f1 = Legendre(p,x) * Legendre(q,temp)
            !f1 = Legendre(p,x) * Legendre(q,x)

        end function f1


        function f2(p,q,x)
            Implicit none

            integer :: p,q
            real(pr) :: x, f2, temp
            if (x > ((2*dx)/(a*dt)) - 1 ) then
                f2 = 0._pr
            else
                temp = 1 - (a*dt/dx)*(1 + x)
                f2 = Legendre(p,x) * Legendre(q,temp)
            endif

        end function f2


        function f3(p,q,x)
            Implicit none

            integer :: p,q
            real(pr) :: x, f3, temp
            if (x < ((2*dx)/(a*dt)) - 1 ) then
                f3 = 0._pr
            else
                temp = x - ((2*dx)/(a*dt))
                f3 = Legendre(p,x) * Legendre(q,temp)
            endif

        end function f3

        function f_alpha(i,j,x)
            implicit NONE

            integer :: j,i
            real(pr) :: x, temp, f_alpha
            temp = (dx/2) * x + (i + 0.5_pr)*dx
            f_alpha = u_init(temp) * Legendre(j,x) 
        
    
        end function f_alpha

        function f_beta(n,j,x)
            implicit none

            integer :: n,j
            real(pr) :: x, f_beta, temp
            temp = (dt/2._pr)*x + (n + 0.5_pr)*dt
            f_beta = u_g(temp) * Legendre(j,x)

        end function f_beta
    
        !! --- initialisation de alpha
        subroutine init_alpha(alpha)
            implicit none

            real(pr), dimension(0:Ne * (ordre + 1) -1), intent(inout) :: alpha
            integer :: i,j
            real(pr) :: temp
            
            do i = 0 , Ne - 1
                do j = 0 , ordre
                    call integrale(i,j, f_alpha,numPoints, temp)
                    temp = (j + 0.5_pr) * temp
                    alpha(i*(ordre + 1) + j) = (dx/2)*temp
                enddo
            enddo 
    
        end subroutine init_alpha

        subroutine init_beta(n,beta)
            implicit none

            real(pr), dimension(0 : ordre), intent(inout) :: beta
            integer, intent(in) :: n
            integer :: j
            real(pr) :: temp

            do j = 0, ordre
                call integrale(n,j,f_beta, numPoints, temp)
                !print*, temp
                temp = (j + 0.5_pr) * temp
                beta(j) = (dt/2)*temp
            enddo

        end subroutine
    
        subroutine fill_L1(L1)
            implicit none

            real(pr), dimension(0 : ordre , 0 : ordre), intent(inout) :: L1
            integer :: p,q
            real(pr) :: temp

            do q = 0, ordre
                do p = 0, ordre
                    call integrale (p,q,f1,numPoints,temp)
                    L1(p,q) = (dx/2) * temp
                end do
            end do

        end subroutine fill_L1

        subroutine fill_L2(L2)
            implicit none

            real(pr), dimension(0 : ordre , 0 : ordre), intent(inout) :: L2
            integer :: p,q
            real(pr) :: temp

            do q = 0, ordre
                do p = 0, ordre
                    call integrale (p,q,f2,numPoints,temp)
                    L2(p,q) = (dt/2) * temp
                end do
            end do

        end subroutine fill_L2

        subroutine fill_L3(L3)
            implicit none

            real(pr), dimension(0 : ordre , 0 : ordre), intent(inout) :: L3
            integer :: p,q
            real(pr) :: temp

            do q = 0, ordre
                do p = 0, ordre
                    call integrale (p,q,f3,numPoints,temp)
                    L3(p,q) = (dt/2) * temp
                end do
            end do

        end subroutine fill_L3

        subroutine fill(L1, L2, l3)
            implicit NONE
            real(pr), dimension(0 : ordre , 0 : ordre), intent(inout) :: L1, L2, L3

            call fill_L1(L1)
            call fill_L2(L2)
            call fill_L3(L3)

        end subroutine fill


end module