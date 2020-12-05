MODULE APROX_FUNC_INTERPOL
    IMPLICIT NONE
    
CONTAINS
    !Plantea el sistema de ecuaciones, luego llama a un método de SEL para resolverlo.
    !Resolviendo el SEL se obtienen los valores de los coeficientes del polinomio interpolante.
    SUBROUTINE POLINOMIO_APROX(X, Y, A)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: A
        !
        
    END SUBROUTINE
END MODULE
