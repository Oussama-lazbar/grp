program analyse
    implicit NONE
    call analyse_erreur("./history/ordre1.dat")
    !call analyse_erreur("ordre1.dat")
    contains

    subroutine analyse_erreur(fileName)
        implicit NONE
        character(len=20), intent(in) :: fileName
        real(8) :: deltax_1, erreur_11, erreur_12, deltax_2, erreur_21, erreur_22, order_1, order_2
        integer :: pos = 0
        open(unit=13, file=fileName)
        do while (pos == 0)
            if (pos < 0) then
                print*, "Oups !something went wrong, check your input!"
            else if (pos > 0) then
                print*, "done"
            else
                read(13,*, iostat=pos) deltax_1, erreur_11, erreur_12, deltax_2, erreur_21, erreur_22
                order_1 = log(erreur_21/erreur_11)/log(deltax_2/deltax_1)
                order_2 = log(erreur_22/erreur_12)/log(deltax_2/deltax_1)
                print*, deltax_1, erreur_11, erreur_12, order_1, order_2
                print*, deltax_2, erreur_21, erreur_22, order_1, order_2
            endif
        enddo
        close(13)
    end subroutine

end program 
            