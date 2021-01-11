MODULE VYM_IO
    IMPLICIT NONE
    CHARACTER(*), PARAMETER :: FORMATO_PREDEF = '(F25.15)'
    
    !En todos los casos, BANDERA devuelve 1 si hubo error al abrir el archivo, la accion no se realiza.
    !De lo contrario, BANDERA devuelve 0.
CONTAINS
    
    !Lee un vector de un archivo cuyos elementos deben estar separados por espacios (o lineas (enter))
    !Si no se ingresa un nombre de ARCHIVO, se lee de 'entrada_vector.txt' por defecto.
    SUBROUTINE VEC_LEER(VEC, BANDERA, ARCHIVO)
        REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: VEC
        INTEGER, INTENT(OUT) :: BANDERA
        CHARACTER(*), INTENT(IN), OPTIONAL :: ARCHIVO
        CHARACTER(:), ALLOCATABLE :: ARCH
        
        INTEGER :: UNIDAD = 1, N, IOS = 0
        BANDERA = 0
        !viendo nombre de archivo
        IF (PRESENT(ARCHIVO)) THEN; ARCH = ARCHIVO; ELSE; ARCH = 'entrada_vector'; END IF
        
        OPEN(UNIDAD, FILE = ARCH, IOSTAT = IOS, ACTION = 'READ')
        IF (IOS /= 0) THEN
            BANDERA = 1
            GOTO 10
        END IF
        READ(UNIDAD,*) N!; IF (IOS /= 0) GOTO 10
        IF (ALLOCATED(VEC)) DEALLOCATE(VEC)
        ALLOCATE(VEC(N))
        READ(UNIDAD,*) VEC(:)!; IF (IOS /= 0) GOTO 10
    10  CLOSE(UNIDAD)
        IF (IOS /= 0) BANDERA = 1
    END SUBROUTINE

    !Imprime un vector en un ARCHIVO de texto, con sus elementos separados por al menos un espacio. 
    !Si no se ingresa un nombre de ARCHIVO, se guarda en 'salida_vector.txt' por defecto.
    !Si no se ingresa un FORMATO se utiliza el formato por defecto del modulo.
    SUBROUTINE VEC_GUARDAR(VEC, BANDERA, ARCHIVO, FORMATO)
        REAL(8), DIMENSION(:), INTENT(IN) :: VEC
        INTEGER, INTENT(OUT) :: BANDERA
        CHARACTER(*), INTENT(IN), OPTIONAL :: ARCHIVO, FORMATO
        
        CHARACTER(:), ALLOCATABLE :: ARCH, F
        INTEGER :: UNIDAD = 1, I, N, IOS = 0
        
        N = SIZE(VEC); BANDERA = 0
        ARCH = 'salida_vector.txt'
        F = FORMATO_PREDEF
        IF (PRESENT(ARCHIVO)) ARCH = ARCHIVO
        IF (PRESENT(FORMATO)) F = FORMATO
        F = F(1:(LEN_TRIM(F)-1))//', A)'
        OPEN(UNIDAD, FILE = ARCH, ACTION = 'WRITE', IOSTAT = IOS)
        IF (IOS /= 0) THEN
            BANDERA = 1
            GOTO 10
        END IF
        WRITE(UNIDAD,*) N
        DO I = 1, N
            WRITE(UNIDAD, F, ADVANCE = 'NO') VEC(I), ' '
        END DO
    10  CLOSE(UNIDAD)
    END SUBROUTINE
    
    !Se muestra un vector por pantalla, con el formato ingresado para sus elementos.
    !Si no se ingresa un FORMATO se utiliza el formato por defecto del modulo.
    SUBROUTINE VEC_MOSTRAR(VEC, FORMATO)
        REAL(8), DIMENSION(:), INTENT(IN) :: VEC
        CHARACTER(*), OPTIONAL, INTENT(IN) :: FORMATO
        
        CHARACTER(:), ALLOCATABLE :: F
        IF (PRESENT(FORMATO)) THEN
            F = FORMATO
        ELSE
            F = FORMATO_PREDEF
        END IF
        
        WRITE(*,F) VEC
    END SUBROUTINE
    
    !Lee una matriz de un archivo cuyos elementos deben estar separados por espacios.
    !Si no se ingresa un nombre de ARCHIVO, se lee de 'entrada_matriz.txt' por defecto.
    SUBROUTINE MAT_LEER(MAT, BANDERA, ARCHIVO)
        REAL(8), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: MAT
        INTEGER, INTENT(OUT) :: BANDERA
        CHARACTER(*), INTENT(IN), OPTIONAL :: ARCHIVO
        CHARACTER(:), ALLOCATABLE :: ARCH
        
        INTEGER :: I, N, M, UNIDAD = 1, IOS = 0
        ARCH = 'salida_matriz.txt'; BANDERA = 0
        IF (PRESENT(ARCHIVO)) ARCH = ARCHIVO
        
        OPEN(UNIDAD, FILE = ARCH, ACTION = 'READ', IOSTAT = IOS)
        IF (IOS /= 0) THEN
            BANDERA = 1
            GOTO 10
        END IF
        READ(UNIDAD,*) N, M
        IF (ALLOCATED(MAT)) DEALLOCATE(MAT)
        ALLOCATE(MAT(N,M))
        DO I = 1, N
            READ(UNIDAD,*) MAT(I,:)
        END DO
        
    10  CLOSE(UNIDAD)
    END SUBROUTINE
    
    !Imprime una matriz en un ARCHIVO de texto, con sus elementos separados por al menos un espacio. 
    !Si no se ingresa un nombre de ARCHIVO, se guarda en 'salida_matriz.txt' por defecto.
    !Si no se ingresa un FORMATO se utiliza el formato por defecto del modulo.
    SUBROUTINE MAT_GUARDAR(MAT, BANDERA, ARCHIVO, FORMATO)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        INTEGER, INTENT(OUT) :: BANDERA
        CHARACTER(*), INTENT(IN), OPTIONAL :: ARCHIVO, FORMATO
        
        CHARACTER(:), ALLOCATABLE :: ARCH, F
        INTEGER :: UNIDAD = 1, I, J, N, M, IOS = 0
        
        N = SIZE(MAT, 1); M = SIZE(MAT, 2); BANDERA = 0
        ARCH = 'salida_matriz.txt'
        F = FORMATO_PREDEF
        IF (PRESENT(ARCHIVO)) ARCH = ARCHIVO
        IF (PRESENT(FORMATO)) F = FORMATO
        F = FORMATEAR(F)
        OPEN(UNIDAD, FILE = ARCH, ACTION = 'WRITE', IOSTAT = IOS)
        IF (IOS /= 0) THEN
            BANDERA = 1
            GOTO 10
        END IF
        WRITE(UNIDAD,*) N, M
        DO I = 1, N
            DO J = 1, M
                WRITE(UNIDAD, F, ADVANCE = 'NO') MAT(I,J), ' '
            END DO
            WRITE(UNIDAD,*)
        END DO
    10  CLOSE(UNIDAD)
    END SUBROUTINE
    
    !Muestra una matriz por pantalla, con el formato que se ingrese para sus elementos.
    !Si no se ingresa un FORMATO se utiliza el formato por defecto del modulo.
    SUBROUTINE MAT_MOSTRAR(MAT, FORMATO)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        CHARACTER(*), OPTIONAL, INTENT(IN) :: FORMATO
        
        CHARACTER(:), ALLOCATABLE :: F
        INTEGER :: I, J, N, M

        IF (PRESENT(FORMATO)) THEN; F = FORMATO; ELSE; F = FORMATO_PREDEF;END IF
        
        N = SIZE(MAT, 1); M = SIZE(MAT, 2)
        F = FORMATEAR(F)
        DO I = 1, N
            DO J = 1, M
                WRITE(*, F, ADVANCE = 'NO') MAT(I,J), ' '
            END DO
            WRITE(*,*)
        END DO
        
    END SUBROUTINE
    
    !Se tiene un formato F para escribir un elemento de un arreglo (ej '(F10.2)') 
    !Se lo transforma en un formato del tipo '(MF10.2, A)' (con M = 1 si no se encuentra COLUMNAS)
    FUNCTION FORMATEAR(F)
        CHARACTER(:), ALLOCATABLE :: FORMATEAR
        CHARACTER(*), INTENT(IN) :: F
        
        INTEGER :: NMENOSUNO
        
        NMENOSUNO = LEN_TRIM(F)-1
        FORMATEAR = F(:NMENOSUNO)//', A)'
        !PRINT *, FORMATEAR
    END FUNCTION
END MODULE

MODULE VYM_MANIP
    IMPLICIT NONE
CONTAINS
    SUBROUTINE VEC_INTERELEM(VEC, I1, I2)
        REAL(8), DIMENSION(:), INTENT(INOUT) :: VEC
        INTEGER, INTENT(IN) :: I1, I2
        
        REAL(8) :: AUX
        
        AUX = VEC(I1)
        VEC(I1) = VEC(I2)
        VEC(I2) = AUX
    END SUBROUTINE
    
    SUBROUTINE MAT_INTERELEM(MAT, I1, J1, I2, J2)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        INTEGER, INTENT(IN) :: I1, I2, J1, J2
        
        REAL(8) :: AUX
        
        AUX = MAT(I1,J1)
        MAT(I1,J1) = MAT(I2,J2)
        MAT(I2,J2) = AUX
    END SUBROUTINE
    
    SUBROUTINE MAT_INTERFILAS(MAT, I1, I2)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        INTEGER, INTENT(IN) :: I1, I2
        
        REAL(8), ALLOCATABLE, DIMENSION(:) :: AUX
        ALLOCATE(AUX(SIZE(MAT,2)))
        
        AUX = MAT(I1,:)
        MAT(I1,:) = MAT(I2,:)
        MAT(I2,:) = AUX
            
        DEALLOCATE(AUX)
    END SUBROUTINE
    
    SUBROUTINE MAT_INTERCOLUMNAS(MAT, J1, J2)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: MAT
        INTEGER, INTENT(IN) :: J1, J2
        
        REAL(8), ALLOCATABLE, DIMENSION(:) :: AUX
        ALLOCATE(AUX(SIZE(MAT,1)))
        
        AUX = MAT(:,J1)
        MAT(:,J1) = MAT(:,J2)
        MAT(:,J2) = AUX
        
        DEALLOCATE(AUX)
    END SUBROUTINE
END MODULE

MODULE VYM_CALCULOS
    IMPLICIT NONE
CONTAINS
    SUBROUTINE MAT_IDENTIDAD(N, MATID)
        INTEGER, INTENT(IN) :: N
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: MATID
        !
        INTEGER :: I
        
        ALLOCATE(MATID(N,N))
        
        MATID = 0.
        DO I = 1, N
            MATID(I,I) = 1.
        END DO
        
    END SUBROUTINE
    
    FUNCTION MAT_ISCUADRADA(MAT)
        LOGICAL :: MAT_ISCUADRADA
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        
        MAT_ISCUADRADA = SIZE(MAT,1) .EQ. SIZE(MAT,2)
    END FUNCTION
        
    !Devuelve un logical con true si la matriz es EDD.
    !EDD: Estrictamente diagonalmente dominante.
    FUNCTION MAT_ISEDD(MAT)
        LOGICAL :: MAT_ISEDD
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        
        MAT_ISEDD = MAT_ISEDD_FILAS(MAT) .OR. MAT_ISEDD_COLS(MAT)
    END FUNCTION
    
    !Es edd por filas si el valor absoluto de c/elemento de la diagonal
    !es mayor que la suma del valor absoluto del resto de su fila.
    FUNCTION MAT_ISEDD_FILAS(MAT)
        LOGICAL MAT_ISEDD_FILAS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), ALLOCATABLE :: MATABS(:,:), SUMAFILAS(:)
        REAL(8) :: SUMASD, DIAG !La SUMA Sin la Diagonal y el elemento de la DIAGonal
        INTEGER :: I, N, M
        N = SIZE(MAT,1); M = SIZE(MAT,2)
        
        ALLOCATE(MATABS(N,M));  MATABS = ABS(MAT)
        ALLOCATE(SUMAFILAS(N)); SUMAFILAS = SUM(MATABS,2)

        I = 1
        !do{algo}while()i
        10  SUMASD = SUMAFILAS(I) - MATABS(I,I)
            DIAG = MATABS(I,I)
            I = I + 1
        IF (I <= N .AND. DIAG > SUMASD) GOTO 10 
        !Si al salir i > N llego hasta el final de la matriz, sino, salio por diag <= sumasd
        !Si i > N recorrio toda la matriz, pero no se si diag > sumasd por eso comparo.
        MAT_ISEDD_FILAS = (I > N .AND. DIAG > SUMASD)
        
        DEALLOCATE(MATABS, SUMAFILAS)
    END FUNCTION
    
    !Es edd por columnas si el valor absoluto de c/elemento de la diagonal
    !es mayor que la suma del valor absoluto del resto de su columna.
    FUNCTION MAT_ISEDD_COLS(MAT)
        LOGICAL MAT_ISEDD_COLS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), ALLOCATABLE :: MATABS(:,:), SUMACOLS(:)
        REAL(8) :: SUMASD, DIAG !La SUMA Sin la Diagonal y el elemento de la DIAGonal
        INTEGER :: I, N, M
        N = SIZE(MAT,1); M = SIZE(MAT,2)
        
        ALLOCATE(MATABS(N,M));  MATABS = ABS(MAT)
        ALLOCATE(SUMACOLS(M)); SUMACOLS = SUM(MATABS,1)

        I = 1
        !do{algo}while()i
        10  SUMASD = SUMACOLS(I) - MATABS(I,I)
            DIAG = MATABS(I,I)
            I = I + 1
        IF (I <= N .AND. DIAG > SUMASD) GOTO 10 
        !Si al salir i > N llego hasta el final de la matriz, sino, salio por diag <= sumasd
        !Si i > N recorrio toda la matriz, pero no se si diag > sumasd por eso comparo.
        MAT_ISEDD_COLS = (I > N .AND. DIAG > SUMASD)
        
        DEALLOCATE(MATABS, SUMACOLS)
    END FUNCTION
    
    FUNCTION VEC_NORMAM(VEC)
        REAL(8) :: VEC_NORMAM
        REAL(8), DIMENSION(:), INTENT(IN) :: VEC
        !
        REAL(8), ALLOCATABLE, DIMENSION(:) :: AUX
        
        ALLOCATE(AUX(SIZE(VEC)))
        AUX = ABS(VEC)
        
        VEC_NORMAM = MAXVAL(AUX)
        DEALLOCATE(AUX)
    END FUNCTION
    
    FUNCTION VEC_NORMAE(VEC)
        REAL(8) :: VEC_NORMAE
        REAL(8), DIMENSION(:), INTENT(IN) :: VEC
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: AUX
        ALLOCATE(AUX(SIZE(VEC)))
        AUX = VEC*VEC
        VEC_NORMAE = SQRT(SUM(AUX))
    END FUNCTION
    
    FUNCTION MAT_NORMAM(MAT)
        REAL(8) :: MAT_NORMAM
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), ALLOCATABLE, DIMENSION(:) :: AUXV
        
        ALLOCATE(AUXV(SIZE(MAT,1)))
        
        AUXV = SUM(ABS(MAT),2)
        MAT_NORMAM = MAXVAL(AUXV)
        
        IF (ALLOCATED(AUXV)) DEALLOCATE(AUXV)
    END FUNCTION
    
    FUNCTION MAT_NORMAL(MAT)
        REAL(8) :: MAT_NORMAL
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), ALLOCATABLE, DIMENSION(:) :: AUXV
        
        ALLOCATE(AUXV(SIZE(MAT,2)))
        
        AUXV = SUM(ABS(MAT),1)
        MAT_NORMAL = MAXVAL(AUXV)
        
        IF (ALLOCATED(AUXV)) DEALLOCATE(AUXV)
    END FUNCTION
    
    FUNCTION MAT_NORMAFROBENIUS(MAT)
        REAL(8) :: MAT_NORMAFROBENIUS
        REAL(8), DIMENSION(:,:), INTENT(IN) :: MAT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: AUX
        
        ALLOCATE(AUX(SIZE(MAT,1),SIZE(MAT,2)))
        AUX = ABS(MAT) !No entiendo por qué el valor absoluto si después se eleva al cuadrado, pero bueno.
        AUX = AUX*AUX
        MAT_NORMAFROBENIUS = SQRT(SUM(SUM(AUX,2)))
        DEALLOCATE(AUX)
    END FUNCTION
END MODULE
