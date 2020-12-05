PROGRAM APROX_FUNC
    !Modulo
    USE VYM_IO
    USE APROX_FUNC_INTERPOL
    
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: XY
    REAL(8), DIMENSION(:), ALLOCATABLE :: X, Y, RESN
    REAL(8), DIMENSION(:), ALLOCATABLE :: RESLG
    INTEGER :: BANDERA
    
    PRINT *, 'Leyendo matriz.'
    CALL MAT_LEER(XY, BANDERA, 'Puntos.txt')
    IF (BANDERA == 1) GOTO 10
    PRINT *, 'Cantidad de puntos: ', SIZE(XY,1)
    IF (SIZE(XY,1) < 2) THEN; PRINT *, 'Se necesitan al menos 2 puntos para interpolar.'; GOTO 20; END IF
    
    PRINT *, 'Separo X e Y de la matriz de puntos iniciales.'
    CALL SEPARARXY(XY, X, Y)
    PRINT *, 'Matriz de puntos iniciales:'
    CALL MAT_MOSTRAR(XY)
    PRINT *, 'Vector de valores de X:'
    CALL VEC_MOSTRAR(X)
    PRINT *, 'Vector de valores de Y:'
    CALL VEC_MOSTRAR(Y)
    
    !Polinomio interpolante normal
    PRINT *, 'POLINOMIO INTERPOLANTE POR SISTEMA DE ECUACIONES LINEALES'
    PRINT *, 'Obteniendo coeficientes del polinomio interpolante.'
    CALL POLINOMIO_APROX(X, Y, RESN)
    PRINT *, 'Coeficientes: '
    CALL VEC_MOSTRAR(RESN)
    
    !Polinomio interpolante de lagrange
    PRINT *, 'POLINOMIO INTERPOLANTE POR LAGRANGE'
    CALL MET_LAGRANGE(X, Y, RESLG)
    
    GOTO 20
10  PRINT *, 'Error de lectura de datos.'
20  PRINT *, 'Fin.'
CONTAINS
    SUBROUTINE SEPARARXY(XY, X, Y)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: XY
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X, Y
        !
        INTEGER :: N
        
        N = SIZE(XY,1)
        ALLOCATE(X(N), Y(N))
        
        X(:) = XY(:,1)
        Y(:) = XY(:,2)
    END SUBROUTINE
    
    SUBROUTINE MET_LAGRANGE(X, Y, RESLG)
        REAL(8), DIMENSION(:), INTENT(IN) :: X, Y
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: RESLG
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: XLG, YLG
        INTEGER :: N
        
        N = SIZE(X)
        ALLOCATE(XLG(0:N-1), YLG(0:N-1), RESLG(0:N-1))
        XLG(:) = X(:); YLG(:) = Y(:);
        
        PRINT *, 'Obteniendo coeficientes del polinomio interpolante de Lagrange.'
        CALL POLINOMIO_LAGRANGE(N, XLG, YLG, RESLG)
        PRINT *, 'Coeficientes: '
        CALL VEC_MOSTRAR(RESLG)
        
    END SUBROUTINE
END PROGRAM
