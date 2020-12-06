MODULE APROX_FUNC_GRAF
    USE APROX_FUNC_MINCUAD
    IMPLICIT NONE
CONTAINS
    SUBROUTINE GUARDAR_POL(A, X0, XN, PASO)
        REAL(8), INTENT(IN) :: A(:), X0, XN, PASO
        !
        REAL(8) :: XACT
        INTEGER :: I
        
        OPEN(1, FILE = 'Pol Minimos Cuadrados.txt')
        XACT = X0
        DO WHILE(XACT < XN)
            WRITE(1, '(2F25.15)') XACT, POLINOMIO(A, XACT)
            I = I + 1
            XACT = XACT + PASO
        END DO
        WRITE(1, '(2F25.15)') XACT, POLINOMIO(A, XACT)
        CLOSE(1)
    END SUBROUTINE
    
    SUBROUTINE GUARDAR_PUNTOS(X, Y)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        !
        INTEGER :: I, N
        N = SIZE(X)
        
        OPEN(1, FILE = 'PuntosPlotMC.txt')
        DO I = 1, N
            WRITE(1, '(2F25.15)') X(I), Y(I)
        END DO
        CLOSE(1)
    END SUBROUTINE
END MODULE
