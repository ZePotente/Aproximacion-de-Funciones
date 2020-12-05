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
    SUBROUTINE DIF_ASC(P, V_Asc, N, x0, h_orig)
        !Variables
        REAL(8), DIMENSION(0:N-1) :: V_Asc ,P, s
        REAL(8) x0, h_orig, h_act
        INTEGER N, i
        
        !Cuerpo
        P = 0.
        P(0) = (-x0 * V_Asc(1) / h_orig) + V_Asc(0) !Preparamos el polinomio cargando los dos
        !primeros términos de la fórmula
        P(1) = V_Asc(1) / h_orig
        h_act = h_orig * h_orig
        s = 0. !Cargamos el s inicial
        s(0) = -x0
        s(1) = 1.
        DO i=2, N-1 !Vamos hasta N-1 que es la cantidad de términos del polinomio
        CALL Mult_Vec_Bin(N, i-1, s, -(x0 + (i-1) * h_orig)) !Vamos multiplicando s 
        !La función llama con el opuesto, es decir, si tenemos X-1 debemos llamar con -1
        P = P + s * V_Asc(i) / (Factorial(i) * h_act)
        h_act = h_act * h_orig
        END DO
    END SUBROUTINE
    
    !Newton con Diferencias Divididas
    !---Fin Newton---!
END MODULE
