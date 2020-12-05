PROGRAM APROX_FUNC
    !Modulo
    USE VYM_IO
    USE APROX_FUNC_INTERPOL
    
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: XY
    REAL(8), DIMENSION(:), ALLOCATABLE :: X, Y
    INTEGER :: BANDERA
    
    PRINT *, 'Leyendo matriz.'
    CALL MAT_LEER(XY, BANDERA, 'Puntos.txt')
    IF (BANDERA == 1) GOTO 10
    PRINT *, 'Separo X e Y de la matriz de puntos iniciales.'
    CALL SEPARARXY(XY, X, Y)
    PRINT *, 'Matriz de puntos iniciales:'
    CALL MAT_MOSTRAR(XY)
    PRINT *, 'Vector de valores de X:'
    CALL VEC_MOSTRAR(X)
    PRINT *, 'Vector de valores de Y:'
    CALL VEC_MOSTRAR(Y)
    
    
    
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
END PROGRAM
