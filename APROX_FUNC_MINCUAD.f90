MODULE APROX_FUNC_MINCUAD
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MIN_CUAD(X, Y, GRADO, RES)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT):: RES
        INTEGER, INTENT(IN) :: GRADO
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: AUX, DIAG, B, AP
        INTEGER :: N, M, I, J, L, C, INFO
        M = SIZE(X)
        N = GRADO + 1
        !Calcula la Longitud de AP (Matriz A comprimida)
        L = 0
        DO I = 1, N
            L = L + I
        END DO
        ALLOCATE(AUX(M), DIAG(2*N-1), B(N), RES(N), AP(L))
        
        DIAG(1) = M
        AUX = 1
        
        DO I = 1, N
            DIAG(I+1) = DOT_PRODUCT(AUX,X)
            B(I) = DOT_PRODUCT(AUX, Y)
            AUX = AUX*X
        END DO
        DO I = N+1, 2*N-2
            DIAG(I+1) = DOT_PRODUCT(AUX,X)
            AUX = AUX*X
        END DO
        !Armar matriz A comprimida
        C = 1
        DO I = 1, N
            DO J = 0, I-1
                AP(C) = DIAG(J+I)
                C = C + 1
            END DO
        END DO
        CALL DPPSV('U', N, 1, AP, B, N, INFO)
        RES = B
        DEALLOCATE(AUX, DIAG, B, AP)
    END SUBROUTINE
    
    FUNCTION MIN_CUAD_RMS(X, Y, A)
        REAL(8) :: MIN_CUAD_RMS
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y, A
        !
        REAL(8) :: SUMA
        INTEGER :: M!, GRADO
        
        M = SIZE(X)
        SUMA = SUMATORIA_ERROR(X, Y, A, M)
        MIN_CUAD_RMS = SQRT(SUMA/M)
        !GRADO = SIZE(A)-1
        
    END FUNCTION
    
    FUNCTION SUMATORIA_ERROR(X, Y, A, M)
        REAL(8) :: SUMATORIA_ERROR
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y, A
        INTEGER, INTENT(IN) :: M
        !
        INTEGER :: K
        SUMATORIA_ERROR = 0.
        DO K = 1, M
            SUMATORIA_ERROR = SUMATORIA_ERROR + (Y(K)- POLINOMIO(A, X(K)))**2
        END DO
    END FUNCTION
    FUNCTION POLINOMIO(A, X)
        REAL(8) :: POLINOMIO
        REAL(8), DIMENSION(:), INTENT(IN) :: A
        REAL(8), INTENT(IN) :: X
        !
        REAL(8) :: XPOT
        INTEGER :: I, N
        
        N = SIZE(A)
        POLINOMIO = A(1); XPOT = 1.;
        DO I = 2, N
            XPOT = XPOT*X !Calculo X^i
            POLINOMIO = POLINOMIO + A(I)*XPOT
        END DO
    END FUNCTION
END MODULE
