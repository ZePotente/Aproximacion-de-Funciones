MODULE APROX_FUNC_INTERPOL
    USE SEL_MET
    IMPLICIT NONE
CONTAINS
    !Plantea el sistema de ecuaciones, luego llama a un método de SEL para resolverlo.
    !Resolviendo el SEL se obtienen los valores de los coeficientes del polinomio interpolante.
    !Los coeficientes van a ser el vector C, con C de Coeficientes.
    SUBROUTINE POLINOMIO_APROX(X, Y, C)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: C
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: A !Matriz de coeficientes de las ecuaciones lineales.
        REAL(8), DIMENSION(:), ALLOCATABLE :: XAUX !Vector de términos independientes
        INTEGER :: I, J, N
        
        N = SIZE(X)
        ALLOCATE(C(N), A(N,N), XAUX(N))
        !Armado de A
        A(:,1) = 1.
        A(:,2) = X(:)
        XAUX = X
        DO J = 3, N
            XAUX = XAUX*X
            DO I = 1, N
                A(I,J) = XAUX(I)
            END DO
        END DO
        CALL MET_LU_CROUT(A, Y, C) !Uso LU Crout porque gauss es más complicado de llamar, 
        !y la matriz no es EDD (estrictamente diagonalmente dominante) como para llamar un método indirecto.
    END SUBROUTINE
END MODULE
