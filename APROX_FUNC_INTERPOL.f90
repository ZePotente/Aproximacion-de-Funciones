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
    
    !---Lagrange, copypasteado---!
    SUBROUTINE POLINOMIO_LAGRANGE(N, x, f, P)
        !Variables
        REAL(8), DIMENSION(0:N-1) :: x, f, P, L
        REAL(8) denom
        INTEGER N, i, j, k
        !Cuerpo
        P = 0.
        DO k=0, N-1
        L = 0.
        L(0) = 1.
        j = 1 ! Tamaño real del polinomio hasta el momento
        denom = 1.
        DO i=0, N-1
            IF (k /= i) THEN
                CALL MULT_VEC_BIN(N, j-1, L, -x(i)) !Llamo con j-1 que será el indice de la última
                !componente a mover (todas las anteriores se moverán también)
                denom = denom * (x(k) - x(i))
                j = j + 1
            END IF
        END DO
        L = L * f(k) / denom
        P = P + L
        END DO
    END SUBROUTINE
    
    SUBROUTINE MULT_VEC_BIN(N, j, L, a)
    !Variables
        REAL(8), DIMENSION(0:N-1) :: L, VAux
        REAL(8) a
        INTEGER N, i, j

        !Cuerpo
        VAux = L
        VAux = VAux * a
        DO i=j, 0, -1 !Hacemos un corrimiento del vector, yendo de atrás hacia adelante para no pisar ningún dato.
            L(i+1) = L(i)
        END DO
        L(0) = 0.
        L = L + VAux
    END SUBROUTINE
    
    FUNCTION Error_estimado(x, N, Punto_a_est, Derivada)
        !Variables
        REAL(8), DIMENSION(0:N) :: x
        REAL(8) :: Punto_a_est, Derivada, Producto, Error_estimado
        INTEGER :: N, i
        
        Producto = 1.
        DO i=0, N
            Producto = Producto * (Punto_a_est - x(i))
        END DO
        Error_estimado = Producto * Derivada / Factorial(N+1)
    END FUNCTION
    
    FUNCTION FACTORIAL(N)
        REAL(8) :: FACTORIAL
        INTEGER, INTENT(IN) :: N
        !
        INTEGER :: I
        FACTORIAL = 1.
        DO I = 2, N
            FACTORIAL = FACTORIAL * I
        END DO
    END FUNCTION
    !Fin Lagrange copypasteado
    
    FUNCTION PUNTO_LAGRANGE(XLG, YLG, XPUNTO)
        REAL(8) :: PUNTO_LAGRANGE
        REAL(8), DIMENSION(:), INTENT(IN) :: XLG, YLG !Se supone que van de (0:N-1)
        REAL(8), INTENT(IN) :: XPUNTO
        !
        REAL(8) :: LK
        INTEGER :: K, N
        
        N = SIZE(XLG)
        
        PUNTO_LAGRANGE = 0.
        DO K = 0, N-1
            LK = CALC_A(K, XLG, XPUNTO) / CALC_A(K, XLG, XLG(K))
            PUNTO_LAGRANGE = PUNTO_LAGRANGE + LK*YLG(K)
        END DO
    END FUNCTION
    
    FUNCTION CALC_A(K, XLG, XPUNTO)
        REAL(8) :: CALC_A
        REAL(8), INTENT(IN) :: XPUNTO, XLG(:) !(va de 0 a n-1)
        INTEGER, INTENT(IN) :: K
        !
        INTEGER :: I, N
        
        N = SIZE(XLG)
        CALC_A = 1.
        DO I = 0, K-1
            CALC_A = CALC_A * (XLG(I)-XPUNTO)
        END DO
        DO I = K+1, N-1
            CALC_A = CALC_A * (XLG(I)-XPUNTO)
        END DO
    END FUNCTION
    
    !---FIN LAGRANGE---!
    !---Inicio Newton---!
    !Newton equiespaciado (diferencias finitas)
    SUBROUTINE MAT_DIF_FIN(Y, MATDIF) ! genera matriz de diferencias asc y desc, pto equispaciados
        REAL(8), DIMENSION(:), INTENT(IN) :: Y
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: MATDIF(:,:)
        !
        INTEGER :: N, M, I, J

        N = SIZE(Y)
        IF (ALLOCATED(MATDIF)) DEALLOCATE(MATDIF)
        
        ALLOCATE(MATDIF(N,N))
        MATDIF = 0.0 ; MATDIF(:, 1) = Y(:)
        M = N
        DO J = 2, N
            M = M - 1
            DO I = 1, M
                MATDIF(i,j) = MATDIF(i+1, j-1) - MATDIF(i, j-1)
            END DO
        END DO
    END SUBROUTINE
    
    SUBROUTINE VEC_DIF_ASC(MATDIF, V_ASC)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATDIF
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V_ASC
        
        ALLOCATE(V_ASC(SIZE(MATDIF,1)))
        V_ASC = MATDIF(1,:)
    END SUBROUTINE

    SUBROUTINE VEC_DIF_DESC(MATDIF, V_DESC)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATDIF
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V_DESC
        !
        INTEGER :: I, N
        N = SIZE(MATDIF,1)
        ALLOCATE(V_DESC(N))
        
        DO I = 1, N
            V_DESC(I) = MATDIF(N - I+1, I)
        END DO
    END SUBROUTINE
    
    SUBROUTINE NEWTON_ASCENDENTE(N, V_ASC, X0, H, A)
        REAL(8), DIMENSION(0:N-1), INTENT(IN) :: V_ASC
        REAL(8), INTENT(IN) :: X0, H
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: A
        INTEGER, INTENT(IN) :: N
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: VEC_S
        REAL(8) :: H_ACT
        INTEGER :: I
        
        IF (ALLOCATED(A)) DEALLOCATE(A)
        ALLOCATE(A(0:N-1), VEC_S(0:N-1))
    
        !Inicializar A
        A = 0.
        A(0) = V_ASC(0) - X0*V_ASC(1)/H
        A(1) = V_ASC(1)/H
        !Inicializar S
        VEC_S = 0.
        VEC_S(0) = -X0
        VEC_S(1) = 1.
        
        H_ACT = H
        
        DO I = 2, N-1
            H_ACT = H_ACT*H
            CALL MULT_VEC_BIN(N, I-1, VEC_S, -(X0 + (I-1)*H)) !S
            
            A = A + VEC_S*V_ASC(I) / (FACTORIAL(I)*H_ACT) !Coeficiente de A
        END DO
        DEALLOCATE(VEC_S)
    END SUBROUTINE
    
    SUBROUTINE NEWTON_DESCENDENTE(N, V_DESC, X0, H, A)
        REAL(8), DIMENSION(0:N-1), INTENT(IN) :: V_DESC
        REAL(8), INTENT(IN) :: X0, H
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: A
        INTEGER, INTENT(IN) :: N
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: VEC_S
        REAL(8) :: H_ACT
        INTEGER :: I
        
        IF (ALLOCATED(A)) DEALLOCATE(A)
        ALLOCATE(A(0:N-1), VEC_S(0:N-1))
        !Inicializar A
        A = 0.
        A(0) = V_DESC(0) - X0*V_DESC(1)/H
        A(1) = V_DESC(1)/H
        !Inicializar S
        VEC_S = 0.
        VEC_S(0) = -X0
        VEC_S(1) = 1.
        
        H_ACT = H
        DO I = 2, N-1
            H_ACT = H_ACT*H
            CALL MULT_VEC_BIN(N, I-1, VEC_S, -(X0 + (I+1)*H)) !S
            
            A = A + VEC_S*V_DESC(I) / (FACTORIAL(I)*H_ACT) !Coeficiente de A
        END DO
        
        DEALLOCATE(VEC_S)
    END SUBROUTINE

    !Newton con Diferencias Divididas
    SUBROUTINE MAT_DIF_DIV(X, Y, MATDIF)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: MATDIF
        !
        INTEGER :: N, M, I, J
        N = SIZE(X)
        IF (ALLOCATED(MATDIF)) DEALLOCATE(MATDIF)
        ALLOCATE(MATDIF(N,N))
        
        MATDIF = 0.
        MATDIF(:,1) = Y(:)
        
        M = N
        DO J = 2, N
            M = M - 1
            DO I = 1, M
                MATDIF(I,J) = (MATDIF(I+1, J-1) - MATDIF(I, J-1)) / (X(I + J-1) - X(I))
            END DO
        END DO
    END SUBROUTINE
    
    SUBROUTINE VEC_DIF_DIV(MATDIF, VECDIF)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATDIF
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: VECDIF
        
        ALLOCATE(VECDIF(SIZE(MATDIF,1)))
        VECDIF = MATDIF(1,:)
    END SUBROUTINE
    
    SUBROUTINE NEWTON_DIFERENCIAS_DIVIDIDAS(N, V_DIV, X, A)
        REAL(8), DIMENSION(0:N-1), INTENT(IN) :: V_DIV, X
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: A
        INTEGER, INTENT(IN) :: N
        !
        REAL(8), ALLOCATABLE :: VEC_AUX(:)
        INTEGER :: I
        
        IF (ALLOCATED(A)) DEALLOCATE(A)
        ALLOCATE(A(0:N-1), VEC_AUX(0:N-1))
        !A
        A = 0.
        A(0) = V_DIV(0) - X(0)*V_DIV(1)
        A(1) = V_DIV(1)
        !AUX similar a S
        VEC_AUX = 0.
        VEC_AUX(0) = -X(0)
        VEC_AUX(1) = 1.
        
        DO I = 2, N-1
            CALL MULT_VEC_BIN(N, I-1, VEC_AUX, -X(I-1)) !AUX
            A = A + VEC_AUX*V_DIV(I)
        END DO
        
        DEALLOCATE(VEC_AUX)
    END SUBROUTINE
    !---Fin Newton---!
END MODULE
