MODULE APROX_FUNC_MINCUAD
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MIN_CUAD_OPTIMO(X, Y, RES_ACT)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT):: RES_ACT
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: RES_SIG !El resultado del grado siguiente
        REAL(8) :: VC_SIG, VC_ACT
        INTEGER :: I
        
        CALL MIN_CUAD(X, Y, 1, RES_ACT); VC_ACT = MIN_CUAD_VARCUAD(X, Y, RES_ACT);
        CALL MIN_CUAD(X, Y, 2, RES_SIG); VC_SIG = MIN_CUAD_VARCUAD(X, Y, RES_SIG);
        WRITE(*, '(A20, F25.15)') 'Actual: ', VC_ACT
        WRITE(*, '(A20, F25.15)') 'Siguiente: ', VC_SIG
        I = 1
        DO WHILE(VC_ACT > VC_SIG .AND. I < SIZE(X)-2)
            I = I + 1!Ya el óptimo es al menos el actual.
            RES_ACT = RES_SIG; VC_ACT = VC_SIG
            CALL MIN_CUAD(X, Y, I+1, RES_SIG); VC_SIG = MIN_CUAD_VARCUAD(X, Y, RES_SIG);
            WRITE(*, '(A20, F25.15)') 'Actual: ', VC_ACT
            WRITE(*, '(A20, F25.15)') 'Siguiente: ', VC_SIG
        END DO
    END SUBROUTINE

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
        IF (ALLOCATED(RES)) DEALLOCATE(RES)
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
    
    FUNCTION MIN_CUAD_VARCUAD(X, Y, A)
        REAL(8) :: MIN_CUAD_VARCUAD
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y, A
        !
        REAL(8) :: SUMA
        INTEGER :: M, GRADO
        
        M = SIZE(X)
        SUMA = SUMATORIA_ERROR(X, Y, A, M)
        GRADO = SIZE(A)-1
        MIN_CUAD_VARCUAD = SUMA/(M-GRADO-1)
    END FUNCTION
    FUNCTION MIN_CUAD_RMS(X, Y, A)
        REAL(8) :: MIN_CUAD_RMS
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y, A
        !
        REAL(8) :: SUMA
        INTEGER :: M
        
        M = SIZE(X)
        SUMA = SUMATORIA_ERROR(X, Y, A, M)
        MIN_CUAD_RMS = SQRT(SUMA/M)
        
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
