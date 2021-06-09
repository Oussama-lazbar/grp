MODULE QUAD
    !!! Module pour le calcul d'integrales avec la 
    !!! Formule de quadrature de Gauss-Legendre à n points (n < 65).
    !!! Les noeuds et poids de la formule sont lus dans un ficher .txt
    !!! généré par le script python gauss_leg_coeff.py
    implicit none

CONTAINS

    subroutine integrale(p,q,l,np,result)
        implicit none

        integer, INTENT(IN) :: p,q,np
        real(8), INTENT(INOUT) :: result
        integer :: i
        real(8) :: wi ,xi                    ! --- poids et noeuds de la formule
        CHARACTER(LEN = 50) :: df

        interface

            function l(m,n,x)
                implicit none
                integer:: m,n
                real(8) :: x
                real(8) :: l
            end function l

        end interface

        WRITE(df, *) np

        !! --- Lecture des poids et noeuds à partir d'un fichier
        OPEN(unit = 10, file = 'n'//trim(adjustl(df))//'.txt')
        result = 0
        
        !! -- Formule d'integration
        !! -- NOTE : les variables p et q sont là
        !! -- au cas où la fonction à integrer depend 
        !! -- de paramètres entiers. Ce qui est le cas dans toutes 
        !! -- les integrales dans notre cadre. L'integration se fait par
        !! -- rapport à x.
        DO i = 1, np
            READ(10,*) wi , xi
            !print*, wi , xi
            result = result + wi*l(p,q,xi)
        end DO

        CLOSE(10)

    end subroutine

end MODULE QUAD
