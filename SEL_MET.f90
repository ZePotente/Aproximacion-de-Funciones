MODULE SEL_MET
    USE VYM_CALCULOS
    USE VYM_MANIP
    USE VYM_IO
    IMPLICIT NONE
CONTAINS        
    
    !Pivotea toda la matriz agarrando desde la diagonal para abajo.
    !Yendo de izquierda a derecha.
    SUBROUTINE SUBRUT_PIVOTEOMATAMP(MAT, B)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        REAL(8), DIMENSION(:), OPTIONAL :: B
        !
        INTEGER :: J, N
        N = SIZE(MAT,1)
        
        DO J = 1, N-1
            IF (PRESENT(B)) THEN
                CALL PIVOTEOMATAMP(MAT, J, B)
            ELSE
                CALL PIVOTEOMATAMP(MAT, J)
            END IF
        END DO
    END SUBROUTINE
    
    SUBROUTINE PIVOTEOMATAMP(MAT, J, B)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        REAL(8), DIMENSION(:), OPTIONAL :: B
        INTEGER, INTENT(IN) :: J
        !
        REAL(8) :: MAXIMO
        INTEGER :: I, INUEVO, N
        N = SIZE(MAT,1)
        MAXIMO = 0.
        INUEVO = J
        DO I = J+1, N
            IF (ABS(MAT(I,J))>MAXIMO) THEN
                MAXIMO = ABS(MAT(I,J))
                INUEVO = I
            END IF
        END DO
        
        IF (INUEVO /= J) THEN
!            PRINT *, 'INTERCAMBIO'
!            PRINT *, 'ANTES:'
!            CALL MAT_MOSTRAR(MAT, '(F8.4)')
            CALL MAT_INTERFILAS(MAT, INUEVO, J)
            IF (PRESENT(B)) CALL VEC_INTERELEM(B, INUEVO, J)
!            WRITE(*,*) 'Se intercambio la fila ', J, ' por la fila ', INUEVO
!            PRINT *, 'DESPUES:'
!            CALL MAT_MOSTRAR(MAT, '(F8.4)')
        END IF
    END SUBROUTINE
    
    SUBROUTINE SUST_REGRESIVA(MATAMP, X)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATAMP
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: X
        !
        INTEGER :: I, K, N, M
        REAL(8), DIMENSION(:), ALLOCATABLE :: SUMA
        
        N = SIZE(MATAMP,1); M = SIZE(MATAMP,2)
        ALLOCATE(X(N,M-N), SUMA(N))
        X = MATAMP(:,N+1:)
        !Paso 1 dividir la ultima fila para tener el resultado del final.
        X(N,:) = X(N,:)/MATAMP(N,N)
        !Paso 2 sustitucion hacia arriba.
        DO I = N-1, 1, -1
            SUMA(:) = 0.
            DO K = I+1, N
                SUMA = SUMA + MATAMP(I,K)*X(K,:)
            END DO
            X(I,:) = (MATAMP(I,N+1:) - SUMA(:))/MATAMP(I,I)
        END DO
        
        DEALLOCATE(SUMA)
    END SUBROUTINE
    
    SUBROUTINE SUST_PROGRESIVA(MATAMP, X)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MATAMP
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: X
        !
        INTEGER :: I, K, N, M
        REAL(8), DIMENSION(:), ALLOCATABLE :: SUMA
        
        N = SIZE(MATAMP,1); M = SIZE(MATAMP,2)
        ALLOCATE(X(N,M-N), SUMA(N))
        X = MATAMP(:,N+1:)
        !Paso 1 dividir la primera fila para tener el resultado del principio.
        X(1,:) = X(1,:)/MATAMP(1,1)
        !Paso 2 sustitucion hacia abajo
        DO I = 2, N
            SUMA(:) = 0.
            DO K = 1, I-1
                SUMA = MATAMP(I,K)*X(K,:)
            END DO
            !x(i) = b(i) - sumatoria(a(i,k)*x(k)) todo /a(i,i)
            X(I,:) = (MATAMP(I,N+1:) - SUMA(:))/MATAMP(I,I)
        END DO
    END SUBROUTINE

    FUNCTION MATRIZAMPLIADA(A, B)
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATRIZAMPLIADA
        REAL(8), DIMENSION(:,:) :: A, B
        
        INTEGER :: N, M, Q
        
        N = SIZE(A,1); M = SIZE(A,2);
        Q = SIZE(B,2);
        
        ALLOCATE(MATRIZAMPLIADA(N,M+Q))
        
        MATRIZAMPLIADA(:,:M) = A
        MATRIZAMPLIADA(:,M+1:) = B
    END FUNCTION

!----------Métodos Directos----------!
    !Devuelve una matriz ampliada en GAUSS, con los resultados correspondientes a B.
    !PUEDE REALIZAR PIVOTEO.
    SUBROUTINE MET_GAUSS(A, B, GAUSS)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GAUSS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        GAUSS = MATRIZAMPLIADA(A, B)
        DO J = 1, ORDEN
            CALL PIVOTEOMATAMP(GAUSS, J)
            DO FILA = J+1, ORDEN
                GAUSS(FILA,J+1:) = GAUSS(FILA,J+1:) - GAUSS(J,J+1:)*GAUSS(FILA,J)/GAUSS(J,J)
                GAUSS(FILA,J) = 0.
            END DO
        END DO
    END SUBROUTINE
    
    !Devuelve una matriz ampliada en GJ, con los resultados correspondientes a B.
    !PUEDE REALIZAR PIVOTEO.
    SUBROUTINE MET_GAUSSJORDAN(A, B, GJ)
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: GJ
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        !
        INTEGER :: J, FILA, ORDEN
        
        ORDEN = SIZE(A,1)
        
        GJ = MATRIZAMPLIADA(A, B)
        
        DO J = 1, ORDEN
            CALL PIVOTEOMATAMP(GJ, J)
            DO FILA = 1, J-1
                GJ(FILA,J+1:) = GJ(FILA,J+1:) - GJ(J,J+1:)*GJ(FILA,J)/GJ(J,J)
                GJ(FILA,J) = 0.
            END DO
            
            DO FILA = J+1, ORDEN
                GJ(FILA,J+1:) = GJ(FILA,J+1:) - GJ(J,J+1:)*GJ(FILA,J)/GJ(J,J)
                GJ(FILA,J) = 0.
            END DO
            !poner 1 en la diagonal, dividiendo toda la fila por el elemento de la diagonal.
            !equivalente a hacer sustitucion regresiva al final.
            GJ(J,J+1:) = GJ(J,J+1:)/GJ(J,J)
            GJ(J,J) = 1.
        END DO
    END SUBROUTINE
    
    SUBROUTINE MAT_INVERSA(MAT, INV)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: INV
        !
        INTEGER :: N, M
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: MATID, MATAMP
        N = SIZE(MAT,1); M = SIZE(MAT,2)
        
        ALLOCATE(INV(N,M), MATAMP(N,N+M))
        CALL MAT_IDENTIDAD(N, MATID)
        CALL MET_GAUSSJORDAN(MAT, MATID, MATAMP)
        INV = MATAMP(:,N+1:)
        
        DEALLOCATE(MATAMP)
    END SUBROUTINE
    
    !Toma una matriz tridiagonal y la divide en tres vectores para resolver.
    !No verifica que la matriz sea tridiagonal. Si vas a usar Thomas, para qué mandarías una que no lo sea.
    !NO REALIZA PIVOTEO
    SUBROUTINE MET_THOMAS(MAT, B, THOMAS)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT, B
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: THOMAS
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: U, D, L
        INTEGER :: I, N
        
        N = SIZE(MAT,1)
        ALLOCATE(U(N),D(N),L(N))
        U(1) = MAT(1,2); U(N) = 0.;
        L(1) = 0.; L(N) = MAT(N-1,N);
        D(1) = MAT(1,1); D(N) = MAT(N,N);
        DO I = 2, N-1
            U(I) = MAT(I,I+1)
            D(I) = MAT(I,I)
            L(I) = MAT(I+1,I)
        END DO
        
        CALL MET_THOMAS_VEC(U, D, L, B, THOMAS)
        DEALLOCATE(U,D,L)
    END SUBROUTINE
    
    !Toma tres vectores provenientes (espero) de una matriz tridiagonal
    !y un conjunto de vectores de términos independientes y resuelve el sistema.
    SUBROUTINE MET_THOMAS_VEC(U_ORIG, D_ORIG, L_ORIG, B_ORIG, THOMAS)
        REAL(8), DIMENSION(:), INTENT(IN) :: U_ORIG, D_ORIG, L_ORIG
        REAL(8), DIMENSION(:,:), INTENT(IN) :: B_ORIG
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: THOMAS
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: U, D, L
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: B
        INTEGER :: I, N, M
        N = SIZE(B_ORIG,1); M = SIZE(B_ORIG,2)
        ALLOCATE(U(N), D(N),L(N), B(N,M))
        U = U_ORIG
        D = D_ORIG
        L = L_ORIG
        B = B_ORIG
        DO I = 1, N-1
            U(I) = U(I)/D(I)
            B(I,:) = B(I,:)/D(I)
            D(I) = 1.
            !
            D(I+1) = D(I+1) - L(I+1)*U(I)
            B(I+1,:) = B(I+1,:) - L(I+1)*B(I,:)
            L(I+1) = 0.
        END DO
        !Sustitucion inversa hardcodeada.
        !Porque no tengo ganas de pasar a matriz para usar la subrutina
        !Paso 1
        B(N,:) = B(N,:)/D(N)
        !Paso 2
        DO I = N-1, 1, -1
            B(I,:) = B(I,:) - U(I)*B(I+1,:)/D(I)
        END DO
        THOMAS = B
        DEALLOCATE(L, D, U, B)
    END SUBROUTINE
    
    !Devuelve un vector solución CROUT, SOLO ACEPTA VECTOR B, NO MATRIZ.
    !NO REALIZA PIVOTEO
    SUBROUTINE MET_LU_CROUT(A, B, CROUT)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: CROUT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: RED
        REAL(8), DIMENSION(:), ALLOCATABLE :: C
        INTEGER :: I, K, N
        
        CALL LU_REDUCCION(A,RED)
        N = SIZE(A,1)
        ALLOCATE(C(N), CROUT(N))
        !Copiado
        !Sustitucion progresiva para calcular C
        C(1) = B(1)/RED(1,1)
        DO I = 2, N
            C(I) = B(I)
            DO K = 1, I-1
                C(I) = C(I) - RED(I,K)*C(K)
            END DO
            C(I) = C(I)/RED(I,I)
        END DO
        !Sustitucion regresiva para calcular X
        CROUT(N) = C(N)
        DO I = N-1, 1, -1
            CROUT(I) = C(I)
            DO K = I+1, N
                CROUT(I) = CROUT(I) - RED(I,K)*CROUT(K)
            END DO
        END DO
        !Fin copiado
        
        DEALLOCATE(C)
    END SUBROUTINE
    
    !A es la combinacion de L y U1 de MAT.
    SUBROUTINE LU_REDUCCION(MAT, RED)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: RED
        !
        INTEGER :: I, FILA, COL, K, N
        N = SIZE(MAT,1)
        !copio
        RED = MAT
        RED(1,2:) = RED(1,2:)/RED(1,1)
        
        DO I = 2, N
            !"Calcula la columna i"
            DO FILA = I, N
                DO K = 1, I-1
                    RED(FILA,I) = RED(FILA,I) - RED(FILA,K)*RED(K,I)
                END DO
            END DO
            !"Calcula fila i"
            DO COL = I+1, N
                DO K = 1, I-1
                    RED(I,COL) = RED(I,COL) - RED(I,K)*RED(K,COL)
                END DO
                RED(I,COL) = RED(I,COL)/RED(I,I)
            END DO
        END DO
    END SUBROUTINE
    
!    SUBROUTINE REF_ITER(A, B, X, TOL, XREF)
!        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
!        REAL(8), DIMENSION(:), INTENT(IN) :: B, X
!        REAL(8), INTENT(IN) :: TOL
!        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: XREF
!        !
!        REAL(8), DIMENSION(:), ALLOCATABLE :: R, DELX
!        INTEGER :: N, ITER, MAXITER
        
!        N = SIZE(B,1)
!        ALLOCATE(R(N), DELX(N), XREF(N))
        
!        XREF = X
!        R = VEC_RESIDUO(A, B, XREF)
        
!        ITER = 0; MAXITER = 100; !Para que no se trabe
!        PRINT *, 'XREF: '; CALL VEC_MOSTRAR(XREF)
!        DO WHILE(ITER < MAXITER .AND. VEC_NORMAE(R) > TOL)
!            CALL MET_LU_CROUT(A, R, DELX)
!            XREF = XREF - DELX
!            R = VEC_RESIDUO(A, B, XREF)
!            ITER = ITER + 1 !Para que no se trabe
!            WRITE(*,*) 'Paso ', ITER; CALL VEC_MOSTRAR(XREF)
!        END DO
        
!        DEALLOCATE(R, DELX)
!    END SUBROUTINE
    SUBROUTINE REF_ITER(A, B, X, TOL, XREF)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B, X
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: XREF
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: BMAT, XMAT
        REAL(8), DIMENSION(:), ALLOCATABLE :: R, DELX
        INTEGER :: N, ITER, MAXITER
        
        N = SIZE(B,1)
        ALLOCATE(R(N), DELX(N), XREF(N), BMAT(N,1), XMAT(N,N+1))
        
        XREF = X
        BMAT(:,1) = B;
        R = VEC_RESIDUO(A, B, XREF)
        
        ITER = 0; MAXITER = 100; !Para que no se trabe
        PRINT *, 'XREF: '; CALL VEC_MOSTRAR(XREF)
        DO WHILE(ITER < MAXITER .AND. VEC_NORMAE(R) > TOL)
            CALL MET_GAUSSJORDAN(A, BMAT, XMAT)
            XREF = XMAT(:,N+1)
            XREF = XREF - DELX
            R = VEC_RESIDUO(A, B, XREF)
            ITER = ITER + 1 !Para que no se trabe
            WRITE(*,*) 'Paso ', ITER; CALL VEC_MOSTRAR(XREF)
        END DO
        
        DEALLOCATE(R, DELX, BMAT, XMAT)
    END SUBROUTINE
!----------Fin Métodos Directos----------!
    
!----------Sensibilidad----------!
    !Toma una matriz SEL o vector de terminos independientes y sus respectivos matriz o vector perturbados
    !Devuelve la cota para el error relativo de la solución de un SEL
    FUNCTION SENS_COTAPERT_SOL(A, APERT, B, BPERT)
        REAL(8) :: SENS_COTAPERT_SOL
        REAL(8), DIMENSION(:,:) :: A
        REAL(8), DIMENSION(:,:), OPTIONAL :: APERT
        REAL(8), DIMENSION(:), OPTIONAL :: B, BPERT
        !
        REAL(8) :: COTA, PROD1, PROD2
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: AINV
        IF (PRESENT(APERT)) THEN
            IF (PRESENT(B) .AND. PRESENT(BPERT)) THEN
                !Perturbacion doble
                CALL MAT_INVERSA(A, AINV)
                IF (MAT_NORMAM(A-APERT) < 1/MAT_NORMAM(AINV)) THEN
                    PROD1 = SENS_NROCOND(A)/(1.-SENS_NROCOND(A)*SENS_ERROR_REL_MAT(A, APERT))
                    PROD2 = SENS_ERROR_REL_MAT(A, APERT) + SENS_ERROR_REL_VEC(B, BPERT)
                    COTA = PROD1*PROD2
                ELSE
                    PRINT *, 'No se cumple la condicion de sensibilidad de perturbacion doble.'
                END IF
            ELSE
            !Perturbacion en A
            COTA = SENS_NROCOND(A)*SENS_ERROR_REL_MAT(A, APERT)
            END IF
        ELSE
            IF (PRESENT(B) .AND. PRESENT(BPERT)) THEN
            !Perturbacion en B
            COTA = SENS_NROCOND(A)*SENS_ERROR_REL_VEC(B, BPERT)
            END IF
        END IF
        SENS_COTAPERT_SOL = COTA
    END FUNCTION
    
    !En la teoria habla de la norma de la matriz, pero no aclara cuál.
    !Yo puse la norma M por costumbre.
    FUNCTION SENS_NROCOND(MAT)
        REAL(8) :: SENS_NROCOND
        REAL(8), DIMENSION(:,:) :: MAT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: INV
        
        CALL MAT_INVERSA(MAT, INV)
        SENS_NROCOND = MAT_NORMAM(MAT)*MAT_NORMAM(INV)
    END FUNCTION
    
    !Calcula el error relativo de una matriz consigo misma pero perturbada.
    FUNCTION SENS_ERROR_REL_MAT(MAT, MATPERT)
        REAL(8) :: SENS_ERROR_REL_MAT
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT, MATPERT
        !
        REAL(8) :: NORMAMNUM, NORMAMDEN
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: RESTA
        ALLOCATE(RESTA(SIZE(MAT,1),SIZE(MAT,2)))
        RESTA = MAT-MATPERT
        
        NORMAMNUM = MAT_NORMAM(RESTA)
        NORMAMDEN = MAT_NORMAM(MAT)
        SENS_ERROR_REL_MAT = NORMAMNUM/NORMAMDEN
        DEALLOCATE(RESTA)
    END FUNCTION
    
    !Calcula el error relativo de un vector consigo mismo pero perturbado.
    FUNCTION SENS_ERROR_REL_VEC(VEC, VECPERT)
        REAL(8) :: SENS_ERROR_REL_VEC
        REAL(8), DIMENSION(:), INTENT(IN) :: VEC, VECPERT
        !
        REAL(8) :: NORMAMNUM, NORMAMDEN
        REAL(8), DIMENSION(:), ALLOCATABLE :: RESTA
        ALLOCATE(RESTA(SIZE(VEC,1)))
        RESTA = VEC - VECPERT
        
        NORMAMNUM = VEC_NORMAM(RESTA)
        NORMAMDEN = VEC_NORMAM(VEC)
        SENS_ERROR_REL_VEC = NORMAMNUM/NORMAMDEN
        DEALLOCATE(RESTA)
    END FUNCTION

!----------Métodos Indirectos----------!
    !Usa A y B para aproximar X dentro de una TOL.
    !No revisa si el SEL converge o diverge, si se quiere se puede revisar antes de llamar.
    !Se por parametro viene una X inicial que sirve de ""semilla"".
    SUBROUTINE MET_JACOBI(AIN, BIN, TOL, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
        REAL(8), DIMENSION(:), INTENT(IN) :: BIN
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8) :: ERROR
        REAL(8), ALLOCATABLE :: XANT(:), BMAT(:,:), B(:), A(:,:), MATAMP(:,:)
        
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
        ALLOCATE(BMAT(N,1))
        BMAT(:,1) = BIN
        
        MATAMP = MATRIZAMPLIADA(AIN, BMAT)
        
        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        CALL SUBRUT_PIVOTEOMATAMP(MATAMP)
        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        ALLOCATE(A(N,N), B(N))
        A = MATAMP(:N,:N)
        B = MATAMP(:,N+1)
        ALLOCATE(XANT(N))
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
            XANT = X
            DO I = 1, N
                X(I) = B(I)
                DO J = 1, I-1
                    X(I) = X(I) - A(I,J)*XANT(J)
                END DO
                
                DO J = I+1, N
                    X(I) = X(I) - A(I,J)*XANT(J)
                END DO
                X(I) = X(I)/A(I,I)
            END DO
            !En la teoria habla de la variacion proporcional.
            !Pero en el codigo solo hace la diferencia, asi que lo dejo asi.
            ERROR = VEC_NORMAM(X-XANT)
!            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
!            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Jacobi diverge'
        DEALLOCATE(XANT)
    END SUBROUTINE
    
    SUBROUTINE MET_GS(AIN, BIN, TOL, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
        REAL(8), DIMENSION(:), INTENT(IN) :: BIN
        REAL(8), INTENT(IN) :: TOL
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8) :: ERROR
        REAL(8), ALLOCATABLE :: XANT(:), BMAT(:,:), B(:), A(:,:), MATAMP(:,:)
        INTEGER :: I, J, N, ITER
        INTEGER, PARAMETER :: MAXITER = 100 !Como para que no se quede en bucle infinito
        
        
        N = SIZE(AIN,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        
        ALLOCATE(BMAT(N,1))
        BMAT(:,1) = BIN
        
        MATAMP = MATRIZAMPLIADA(AIN, BMAT)
        
        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        CALL SUBRUT_PIVOTEOMATAMP(MATAMP)
        CALL MAT_INTERFILAS(MATAMP, 3, 4)
        PRINT *, 'MATAMP'; CALL MAT_MOSTRAR(MATAMP)
        ALLOCATE(A(N,N), B(N))
        A = MATAMP(:N,:N)
        B = MATAMP(:,N+1)

!        ALLOCATE(XANT(N))
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
!            XANT = X
            DO I = 1, N
                X(I) = B(I)
        
                DO J = 1, I-1
                    X(I) = X(I) - A(I,J)*X(J)
                END DO
                
                DO J = I+1, N
                    X(I) = X(I) - A(I,J)*X(J)
                END DO
        
                X(I) = X(I)/A(I,I)
            END DO
!            !En la teoria habla de la variacion proporcional.
!            !Pero en el codigo solo hace la diferencia, asi que lo dejo asi.
!            ERROR = VEC_NORMAM(X-XANT)
            ERROR = VEC_NORMAM(VEC_RESIDUO(A, B, X))
            WRITE(*,*) 'Error', ERROR
            ITER = ITER + 1
        END DO
        IF (ITER >= MAXITER) PRINT *, 'Gauss-Seidel diverge'
        WRITE(*,*) 'Iteraciones: ', ITER
        WRITE(*,*) 'Error: ', ERROR
!        DEALLOCATE(XANT)
    END SUBROUTINE
    
    FUNCTION VEC_RESIDUO(A, B, X)
        REAL(8), DIMENSION(:), ALLOCATABLE :: VEC_RESIDUO
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B, X
        !
        ALLOCATE(VEC_RESIDUO(SIZE(X)))
        VEC_RESIDUO = MATMUL(A,X) - B
    END FUNCTION
    
    !---no entiendo bien el metodo y no se por que no anda---!
    !Es un método que hace sobre Gauss-Seidel.
    !Si W (omega) > 2 diverge.
    !Si W < 1 es sub-relajacion.
    !Si 1 < w < 2 es sobre-relajacion (el OR de SOR)
    !La idea sería entonces que W sobre-relaje.
    SUBROUTINE MET_SOR(A, B, TOL, W, X, XI)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B
        REAL(8), INTENT(IN) :: TOL, W
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        REAL(8), DIMENSION(:), INTENT(IN), OPTIONAL :: XI
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: R, XANT
        REAL(8) :: ERROR
        INTEGER :: I, N, ITER
        INTEGER, PARAMETER :: MAXITER = 1000 !Como para que no se quede en bucle infinito
        
        N = SIZE(A,1)
        IF (.NOT.PRESENT(XI)) THEN
            ALLOCATE(X(N))
            X = 0.
        ELSE
            X = XI
        END IF
        ALLOCATE(R(N), XANT(N))
        ERROR = 2.*TOL !Valor imposible
        ITER = 0; ERROR = 2.*TOL !Valor imposible
        DO WHILE(ITER < MAXITER .AND. ERROR >= TOL)
            I = 0
            XANT = X
            DO I = 1, N
                !Por alguna razon, en la clase 5 el residuo es ax-b, pero en la clase 6 es b-ax.
                !Voy a seguir la clase 6 asi que lo voy a poner negativo.
                R = VEC_RESIDUO(A, X, B)
                X(I) = XANT(I) + W*R(I)/A(I,I)
                CALL VEC_MOSTRAR(X)
                READ(*,*)
            END DO
            ITER = ITER + 1
        END DO
        DEALLOCATE(R)
    END SUBROUTINE
END MODULE
