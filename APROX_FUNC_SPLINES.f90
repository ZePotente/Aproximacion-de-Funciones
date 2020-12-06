MODULE APROX_FUNC_SPLINES
    IMPLICIT NONE
CONTAINS
    SUBROUTINE CREAR_H(X, H)
        REAL(8), DIMENSION(:), INTENT(IN) :: X
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: H
        !
        INTEGER :: N, I
        N = SIZE(X)
        IF (ALLOCATED(H)) DEALLOCATE(H)
        ALLOCATE(H(0:N-1))
        DO I = 0, N-1
            H(I) = X(I+1) - X(I)
        END DO
    END SUBROUTINE
    
    SUBROUTINE SPLINES_CUBICOS(A, B, C, D, H)
        REAL(8), DIMENSION(0:), INTENT(IN) :: A, H
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: B, C, D
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: UD, DD, LD
        INTEGER :: N, I, INFO
        N = SIZE(A) !Cantidad de puntos total.
        
        ALLOCATE(B(0:N-1), C(0:N-1), D(0:N-1))
        ALLOCATE(UD(0:N-2), DD(0:N-1), LD(0:N-2))
        
        !Preparo el sistema tridiagonal, inicializando
        UD(0) = 0.; UD(0) = 1.; LD(0) = H(0)
        C(0) = 0.
        !preparando el medio
        DO I = 1, N-2
            UD(I) = H(I)
            DD(I) = 2*(H(I-1) + H(I))
            LD(I) = H(I)
            C(I)  = 3*((A(I+1) - A(I)) / H(I) - (A(I) - A(I-1))/ H(I-1))
        END DO
        !ultimos valores
        DD(N-1) = 1.; LD(N-2) = 0.; C(N-1) = 0.
        
        !Resolver el SEL
        CALL DGTSV(N-1, 1, LD, DD, UD, C, N, INFO)
        !O, si no funciona por problemas con la librería , importar SEL_MET y llamar al método de Thomas vectorial.
        !Coeficientes
        DO I = 0, N-2
            B(I) = (A(I+1) - A(I)) / H(I) - H(I)*(2.*C(I) + C(I+1))/3.
            D(I) = (C(I+1) - C(I)) / (3.*H(I))
        END DO
        
        DEALLOCATE(UD, DD, LD)
    END SUBROUTINE
    
    SUBROUTINE GUARDARSPLINE(X, A, B, C, D, PASO)
        REAL(8), DIMENSION(0:), INTENT(IN) :: X, A, B, C, D
        REAL(8), INTENT(IN) :: PASO
        !
        REAL(8) :: XACT
        INTEGER :: I, J, N
        CHARACTER(*), PARAMETER :: ARCH_SPLINES = 'Splines.txt'
        N = SIZE(X)
        !Guardo puntos de los splines
        OPEN(1, FILE = ARCH_SPLINES)
        XACT = X(0)
        I = 0; J = 0;
        DO WHILE(XACT < X(N-1))
            IF (XACT > X(J+1)) J = J+1 !Supongo que PASO es menor que el intervalo de x más chico.
            
            WRITE(1, '(2F25.15)') XACT, POL(J, A, B, C, D, X, XACT)
            I = I + 1
            XACT = XACT + PASO
        END DO
        !Escribo el punto que sigue para que no se corte la gráfica antes del último punto dato cuando compare ambas
        WRITE(1, '(2F25.15)') XACT, POL(J, A, B, C, D, X, XACT)
        CLOSE(1)
        !No guardo los puntos porque yo ya los leo de un archivo.
    END SUBROUTINE
    
    FUNCTION POL(J, A, B, C , D, X, XACT)
        REAL(8) :: POL
        REAL(8), DIMENSION(0:), INTENT(IN) :: X, A, B, C, D
        REAL(8), INTENT(IN) :: XACT
        INTEGER, INTENT(IN) :: J
        !
        REAL(8) :: H
        H = XACT - X(J)
        POL = A(J) + B(J)*H + C(J)*H**2 + D(J)*H**3
    END FUNCTION
END MODULE
