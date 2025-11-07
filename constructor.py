

entrenamiento = 'C:/Users/USUARIO/Documents/Fortran/FIS_automatico/Caso02/Entrenamiento.csv' # ruta a los datos de entrenamiento
verificacion = 'C:/Users/USUARIO/Documents/Fortran/FIS_automatico/Caso02/Verificacion.csv' # ruta a los datos de entrenamiento
validacion = 'C:/Users/USUARIO/Documents/Fortran/FIS_automatico/Caso02/Validacion.csv' # ruta a los datos de entrenamiento
carpeta = ''
universo_discurso_entrada_01 = [0.0, 2.0]
universo_discurso_entrada_02 = [0.0, 0.2]
universo_discurso_salida = [-1.0, 1.0]
referencia_entrada_01 = 1.13
referencia_entrada_02 = 0.08
cantidad_variables = 16
vector_reglas = "[1.0d0, 1.0d0, 3.0d0, 3.0d0, 2.0d0, 1.0d0, 2.0d0, 3.0d0, 1.0d0]"
osa = False
reglas_discretas = False

tamanio_poblacion = 40
iteraciones = 10000
pruebas = 10
minimo_global_deseado = -1.0


read_me = f'''
entrenamiento = {entrenamiento}
verificacion = {verificacion}
validacion = {validacion}
carpeta = {carpeta}
universo_discurso_entrada_01 = {universo_discurso_entrada_01}
universo_discurso_entrada_02 = {universo_discurso_entrada_02}
universo_discurso_salida = {universo_discurso_salida}
referencia_entrada_01 = {referencia_entrada_01}
referencia_entrada_02 = {referencia_entrada_02}
cantidad_variables = {cantidad_variables}
vector_reglas = {vector_reglas}
osa = {osa}
reglas_discretas = {reglas_discretas}

tamanio_poblacion = {tamanio_poblacion}
iteraciones = {iteraciones}
pruebas = {pruebas}
minimo_global_deseado = {minimo_global_deseado}
'''

with open("read_me.txt", "w", encoding="utf-8") as f:
    f.write(read_me)

import pandas as pd

metaheuristicas='''
module metaheuristicas
  implicit none

contains

  function differential_evolution(CostFunction, LimInf, LimSup, NumPop, MaxIter, minimo_global_deseado) result(Solution)
    double precision, dimension(:), allocatable :: Solution

    interface
        double precision function CostFunction(x)
            double precision, dimension(:), intent(in) :: x
        end function CostFunction
    end interface


    integer, intent(in) :: NumPop, MaxIter
    double precision, intent(in) :: minimo_global_deseado
    double precision, dimension(:), intent(in) :: LimInf, LimSup
    double precision :: F, Cr, Fbest, FunZ, RndCr, rand_num, calls
    integer :: i, j, k, n, r1, r2, r3, NumVar, ultima_actualizacion

    double precision, dimension(:, :), allocatable :: X
    double precision, dimension(:), allocatable :: Fit
    double precision, dimension(:), allocatable :: Xbest, y, z

    F=0.5
    Cr=0.2
    Fbest = 1e24

    NumVar = size(LimSup)

    allocate(X(NumPop, NumVar))
    allocate(Solution(NumVar + 2))
    allocate(Fit(NumPop))
    allocate(Xbest(NumVar))
    allocate(y(NumVar))
    allocate(z(NumVar))

    calls = 0
    ultima_actualizacion = NumPop * MaxIter

    do i = 1, NumPop
      do j = 1, NumVar
        call random_number(rand_num)
        X(i, j) = LimInf(j) + (LimSup(j) - LimInf(j)) * rand_num
      end do
      Fit(i) = CostFunction(X(i, :))
      calls = calls + 1
    end do

    do i = 1, MaxIter
        ! print *, "iteracion", i
        do j = 1, NumPop

            call random_number(rand_num)
            r1 = int(rand_num * NumPop) + 1

            do 
                call random_number(rand_num)
                r2 = int(rand_num * NumPop) + 1
                if (r2 /= r1) exit 
            end do 
            do 
                call random_number(rand_num)
                r3 = int(rand_num * NumPop) + 1
                if (r3 /= r1 .and. r3 /= r2) exit 
            end do 

            z = X(j, :)

            do k = 1, NumVar
                call random_number(rand_num)
                if (rand_num < Cr) then 
                    z(k) = X(r1,k) + F * (X(r2,k)-X(r3,k))
                end if
                if (z(k) > LimSup(k) .or. z(k) < LimInf(k)) then 
                    do n = 1, NumVar
                        call random_number(rand_num)
                        z(n) = LimInf(n) + (LimSup(n) - LimInf(n)) * rand_num
                    end do
                    exit 
                end if 
            end do

            FunZ = CostFunction(z)
            calls = calls + 1

            if (FunZ < Fit(j)) then 
                X(j, :) = z 
                Fit(j) = FunZ 
                if (FunZ < Fbest) then 
                    Fbest = FunZ
                    Xbest = z 
                    ! ultima_actualizacion = calls
                end if 
            end if
        end do

        if (Fbest < minimo_global_deseado) then 
            ultima_actualizacion = calls
            exit 
         end if 

    end do

    Solution(1:NumVar) = Xbest
    Solution(NumVar + 1) = Fbest
    Solution(NumVar + 2) = ultima_actualizacion
    deallocate(X, Fit, Xbest, z)

  end function differential_evolution

    function particle_swarm_optimization(funcion_costo, limites_inferiores, limites_superiores, &
                                            tamanio_poblacion, maximo_iteraciones) result(solucion)
        ! === Entradas ===
        interface
            double precision function funcion_costo(x)
                double precision, dimension(:), intent(in) :: x
            end function funcion_costo
        end interface

        double precision, intent(in)  :: limites_inferiores(:), limites_superiores(:)
        integer,  intent(in)          :: tamanio_poblacion, maximo_iteraciones
        ! === Salida ===
        double precision, allocatable :: solucion(:)   ! [mejor_posicion_global, mejor_costo_global]

        ! === Parámetros de control ===
        double precision :: w     = 1.0d0
        double precision :: wdamp = 0.99d0
        double precision :: c1    = 1.5d0
        double precision :: c2    = 2.0d0

        ! === Variables ===
        integer :: numero_variables, iter, pop, j, idx, calls, ultima_actualizacion
        double precision, allocatable :: posicion(:,:), velocidad(:,:)
        double precision, allocatable :: mejor_posicion(:,:), mejor_costo(:), costo(:)
        double precision, allocatable :: mejor_posicion_global(:), velocidad_maxima(:), velocidad_minima(:)
        double precision :: mejor_costo_global
        double precision, allocatable :: r1(:), r2(:)

        ! ---- Inicialización ----
        numero_variables = size(limites_superiores)
        allocate(posicion(tamanio_poblacion, numero_variables))
        allocate(velocidad(tamanio_poblacion, numero_variables))
        allocate(mejor_posicion(tamanio_poblacion, numero_variables))
        allocate(mejor_costo(tamanio_poblacion), costo(tamanio_poblacion))
        allocate(mejor_posicion_global(numero_variables))
        allocate(velocidad_maxima(numero_variables), velocidad_minima(numero_variables))
        allocate(r1(numero_variables), r2(numero_variables))

        call random_seed()

        ! Posiciones iniciales
        calls = 0
        do pop = 1, tamanio_poblacion
            do j = 1, numero_variables
                call random_number(posicion(pop,j))
                posicion(pop,j) = limites_inferiores(j) + posicion(pop,j) * (limites_superiores(j) - limites_inferiores(j))
            end do
        end do
        velocidad = 0.0d0

        ! Evaluación inicial
        do pop = 1, tamanio_poblacion
            costo(pop) = funcion_costo(posicion(pop,:))
            calls = calls + 1
            mejor_posicion(pop,:) = posicion(pop,:)
            mejor_costo(pop) = costo(pop)
        end do

        ! Mejor global inicial
        idx = 1
        mejor_costo_global = costo(1)
        do pop = 2, tamanio_poblacion
            if (costo(pop) < mejor_costo_global) then
                mejor_costo_global = costo(pop)
                idx = pop
            end if
        end do
        mejor_posicion_global = posicion(idx,:)
        ultima_actualizacion = calls

        ! Límites de velocidad
        velocidad_maxima = 0.1d0 * (limites_superiores - limites_inferiores)
        velocidad_minima = -velocidad_maxima

        ! Iteraciones
        do iter = 1, maximo_iteraciones
            do pop = 1, tamanio_poblacion
                call random_number(r1)
                call random_number(r2)

                do j = 1, numero_variables
                    velocidad(pop,j) = w*velocidad(pop,j)                                     &
                        + c1*r1(j)*(mejor_posicion(pop,j) - posicion(pop,j))                  &
                        + c2*r2(j)*(mejor_posicion_global(j) - posicion(pop,j))

                    if (velocidad(pop,j) > velocidad_maxima(j)) velocidad(pop,j) = velocidad_maxima(j)
                    if (velocidad(pop,j) < velocidad_minima(j)) velocidad(pop,j) = velocidad_minima(j)
                end do

                posicion(pop,:) = posicion(pop,:) + velocidad(pop,:)

                do j = 1, numero_variables
                    if (posicion(pop,j) > limites_superiores(j) .or. posicion(pop,j) < limites_inferiores(j)) then
                        velocidad(pop,j) = -velocidad(pop,j)
                    end if
                    if (posicion(pop,j) > limites_superiores(j)) posicion(pop,j) = limites_superiores(j)
                    if (posicion(pop,j) < limites_inferiores(j)) posicion(pop,j) = limites_inferiores(j)
                end do

                costo(pop) = funcion_costo(posicion(pop,:))
                calls = calls + 1

                if (costo(pop) < mejor_costo(pop)) then
                    mejor_posicion(pop,:) = posicion(pop,:)
                    mejor_costo(pop) = costo(pop)

                    if (mejor_costo(pop) < mejor_costo_global) then
                        mejor_costo_global = mejor_costo(pop)
                        mejor_posicion_global = mejor_posicion(pop,:)
                        ultima_actualizacion = calls
                    end if
                end if
            end do
            w = w * wdamp
        end do

        ! ---- Salida ----
        allocate(solucion(numero_variables+2))
        solucion(1:numero_variables) = mejor_posicion_global
        solucion(numero_variables+1) = mejor_costo_global
        solucion(numero_variables+2) = ultima_actualizacion
    end function particle_swarm_optimization

    function harmony_search(funcion_costo,limites_inferiores,limites_superiores, &
                            numero_nuevas_armonias,numero_armonias_en_memoria,   &
                            maximo_iteraciones) result(solucion)
        ! === Entradas (nombres iguales al MATLAB) ===
        interface
            double precision function funcion_costo(x)
                double precision, dimension(:), intent(in) :: x
            end function funcion_costo
        end interface
        double precision, intent(in)      :: limites_inferiores(:), limites_superiores(:)
        integer, intent(in)               :: numero_nuevas_armonias
        integer, intent(in)               :: numero_armonias_en_memoria
        integer, intent(in)               :: maximo_iteraciones

        ! === Salida (igual que MATLAB, pero sólo 'solucion') ===
        double precision, allocatable     :: solucion(:)  ! [mejor_arma, mejor_costo]

        ! ==== Parámetros de control (mismos nombres) ====
        double precision :: hmcr, par, fw_damp, mejor_costo_global
        double precision, allocatable :: fw(:)

        ! ==== Variables ====
        integer :: numero_acordes, i, j, iter, numero_armonia, m, nTot, calls, ultima_actualizacion
        double precision, allocatable :: memoria_armonias(:,:), nuevas_armonias(:,:)
        double precision, allocatable :: costo_memoria(:), costo(:), costos_totales(:)
        double precision, allocatable :: armonias_totales(:,:)
        integer, allocatable          :: idx(:)
        integer :: b, ii

        ! ---- Inicialización de parámetros (mismos valores que tu script) ----
        numero_acordes = size(limites_superiores)
        hmcr   = 0.9d0
        par    = 0.1d0
        allocate(fw(numero_acordes))
        fw     = 0.02d0*(limites_superiores - limites_inferiores)
        fw_damp = 0.995d0
        calls = 0

        ! ---- Matrices iniciales ----
        allocate(costo_memoria(numero_armonias_en_memoria))
        allocate(costo(numero_nuevas_armonias))
        allocate(memoria_armonias(numero_armonias_en_memoria, numero_acordes))

        call random_seed()

        ! memoria_armonias inicial ~ U(lb, ub)
        do i = 1, numero_armonias_en_memoria
            do j = 1, numero_acordes
                call random_number(memoria_armonias(i,j))
                memoria_armonias(i,j) = limites_inferiores(j) + memoria_armonias(i,j) * &
                                         (limites_superiores(j) - limites_inferiores(j))
            end do
        end do

        ! costo_memoria
        do i = 1, numero_armonias_en_memoria
            costo_memoria(i) = funcion_costo(memoria_armonias(i,:))
            calls = calls + 1
        end do

        ! Ordenar memoria por costo ascendente
        allocate(idx(numero_armonias_en_memoria))
        do i = 1, numero_armonias_en_memoria
            idx(i) = i
        end do
        call argsort_asc(costo_memoria, idx)
        mejor_costo_global = costo_memoria(1)
        ultima_actualizacion = calls
        memoria_armonias = memoria_armonias(idx,:)
        costo_memoria    = costo_memoria(idx)
        deallocate(idx)

        ! ---- Bucle principal ----
        do iter = 1, maximo_iteraciones

            ! nuevas_armonias de arranque aleatorio (como en MATLAB)
            allocate(nuevas_armonias(numero_nuevas_armonias, numero_acordes))
            do numero_armonia = 1, numero_nuevas_armonias
                do j = 1, numero_acordes
                    call random_number(nuevas_armonias(numero_armonia,j))
                    nuevas_armonias(numero_armonia,j) = limites_inferiores(j) + nuevas_armonias(numero_armonia,j) * &
                                                        (limites_superiores(j) - limites_inferiores(j))
                end do
            end do

            ! Improvisación con HMCR / PAR / FW
            do numero_armonia = 1, numero_nuevas_armonias
                do j = 1, numero_acordes
                    if (randu() < hmcr) then
                        b = 1 + int( randu() * dble(numero_armonias_en_memoria) )
                        if (b > numero_armonias_en_memoria) b = numero_armonias_en_memoria
                        nuevas_armonias(numero_armonia,j) = memoria_armonias(b,j)
                    end if
                    if (randu() < par) then
                        nuevas_armonias(numero_armonia,j) = nuevas_armonias(numero_armonia,j) + fw(j)*randn0()
                    end if
                end do
                ! Clampeo a [lb, ub]
                do j = 1, numero_acordes
                    if (nuevas_armonias(numero_armonia,j) < limites_inferiores(j)) nuevas_armonias(numero_armonia,j) = limites_inferiores(j)
                    if (nuevas_armonias(numero_armonia,j) > limites_superiores(j)) nuevas_armonias(numero_armonia,j) = limites_superiores(j)
                end do
                costo(numero_armonia) = funcion_costo(nuevas_armonias(numero_armonia,:))
                calls = calls + 1
            end do

            ! Unir nuevas + memoria y quedarnos con las mejores 'numero_armonias_en_memoria'
            m    = numero_armonias_en_memoria
            nTot = numero_nuevas_armonias + m
            allocate(costos_totales(nTot))
            allocate(armonias_totales(nTot, numero_acordes))

            costos_totales(1:numero_nuevas_armonias) = costo
            costos_totales(numero_nuevas_armonias+1:nTot) = costo_memoria

            armonias_totales(1:numero_nuevas_armonias,:) = nuevas_armonias
            armonias_totales(numero_nuevas_armonias+1:nTot,:) = memoria_armonias

            allocate(idx(nTot))
            do ii = 1, nTot
                idx(ii) = ii
            end do
            call argsort_asc(costos_totales, idx)

            memoria_armonias = armonias_totales(idx(1:m), :)
            costo_memoria    = costos_totales(idx(1:m))

            if (costos_totales(1) < mejor_costo_global) then 
                ultima_actualizacion = calls 
                mejor_costo_global = costos_totales(1)
            end if

            deallocate(idx, costos_totales, armonias_totales, nuevas_armonias)

            ! Atenuar el ancho de banda
            fw = fw * fw_damp
        end do

        ! ---- Salida tipo MATLAB ----
        allocate(solucion(numero_acordes+2))
        solucion(1:numero_acordes) = memoria_armonias(1,:)
        solucion(numero_acordes+1) = costo_memoria(1)
        solucion(numero_acordes+2) = ultima_actualizacion

        ! Limpieza
        deallocate(memoria_armonias, costo_memoria, costo, fw)

    end function harmony_search

    !==================== Utilidades ====================

    ! Uniforme (0,1)
    double precision function randu()
        double precision :: r
        call random_number(r)
        randu = r
    end function randu

    ! Normal estándar ~ N(0,1) (Box–Muller, una muestra)
    double precision function randn0()
        double precision :: u1, u2
        call random_number(u1)
        call random_number(u2)
        if (u1 <= 0.d0) u1 = 1.d-16
        randn0 = sqrt(-2.d0*log(u1)) * cos(2.d0*acos(-1.d0)*u2)
    end function randn0

    ! Ordena ascendentemente 'a' y aplica el mismo reorden a 'idx' (inserción)
    subroutine argsort_asc(a, idx)
        double precision, intent(inout) :: a(:)
        integer,          intent(inout) :: idx(:)
        integer :: n, i, j, keyi
        double precision :: key
        n = size(a)
        do i = 2, n
            key  = a(i)
            keyi = idx(i)
            j = i - 1
            do while (j >= 1 .and. a(j) > key)
                a(j+1)   = a(j)
                idx(j+1) = idx(j)
                j = j - 1
            end do
            a(j+1)   = key
            idx(j+1) = keyi
        end do
    end subroutine argsort_asc

end module metaheuristicas

'''


with open("metaheuristicas.f90", "w", encoding="utf-8") as f:
    f.write(metaheuristicas)

df_entrenamiento = pd.read_csv(entrenamiento)
df_entrenamiento.columns = ["x_1", "x_2", "respuesta"]
numero_datos_entrenamiento = df_entrenamiento.shape[0]

df_verificacion = pd.read_csv(verificacion)
df_verificacion.columns = ["x_1", "x_2", "respuesta"]
numero_datos_verificacion = df_verificacion.shape[0]

df_validacion = pd.read_csv(validacion)
df_validacion.columns = ["x_1", "x_2", "respuesta"]
numero_datos_validacion = df_validacion.shape[0]

if cantidad_variables == 16:
    call_reglas = "call base_de_reglas(anint(x(1:9)), base)"
    limites_reglas = """  limite_inferior(1:9) = 1.0d0
  limite_superior(1:9) = 3.0d0"""
    definicion_params = f'''
        params(1) = -1000000.0d0
        params(2) = {universo_discurso_entrada_01[0]}d0
        params(3) = params(2) + x(10) * ({universo_discurso_entrada_01[1]}d0 - params(2))
        params(4) = params(3) + x(11) * ({universo_discurso_entrada_01[1]}d0 - params(3))

        params(5) = params(3)
        params(6) = params(4)
        params(7) = params(6) + x(12) * ({universo_discurso_entrada_01[1]}d0 - params(6))

        params(8) = params(6)
        params(9) = params(7)
        params(10) = {universo_discurso_entrada_01[1]}d0
        params(11) = 1000000.0d0

        params(12) = -1000000.0d0
        params(13) = {universo_discurso_entrada_02[0]}d0
        params(14) = params(13) + x(13) * ({universo_discurso_entrada_02[1]}d0 - params(13))
        params(15) = params(14) + x(14) * ({universo_discurso_entrada_02[1]}d0 - params(14))

        params(16) = params(14)
        params(17) = params(15)
        params(18) = params(17) + x(15) * ({universo_discurso_entrada_02[1]}d0 - params(17))

        params(19) = params(17)
        params(20) = params(18)
        params(21) = {universo_discurso_entrada_02[1]}d0
        params(22) = 1000000.0d0
        
        params(23) = x(16) * ({universo_discurso_salida[0]}d0)
        params(24) = params(23)
        params(25) = 0.0d0
        params(26) = 0.0d0
        params(27) = x(16) * ({universo_discurso_salida[1]}d0)
        params(28) = params(27)
'''
    guardado = f'''
  write(10, '(A)') "x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,x_12,x_13,x_14,x_15,x_16,Costo,Times"

  do iteracion = 1, maxima_iteracion
    !             1        2         3         4         5         6         7         8        9          10        11         12        13       14        15         16       17        18        19        20         21       22        23         24        25        26                              
    write(10, '(F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7)') &
    soluciones(iteracion, 1), soluciones(iteracion, 2), soluciones(iteracion, 3), soluciones(iteracion, 4), &
    soluciones(iteracion, 5), soluciones(iteracion, 6), soluciones(iteracion, 7), soluciones(iteracion, 8), &
    soluciones(iteracion, 9), soluciones(iteracion,10), soluciones(iteracion,11), soluciones(iteracion,12), &
    soluciones(iteracion,13), soluciones(iteracion,14), soluciones(iteracion,15), soluciones(iteracion,16), &
    soluciones(iteracion,17), soluciones(iteracion,18)
  end do
  close(10)
'''

if cantidad_variables == 14:
    call_reglas = "call base_de_reglas(anint(x(1:9)), base)"
    limites_reglas = """  limite_inferior(1:9) = 1.0d0
  limite_superior(1:9) = 3.0d0"""
    definicion_params = f'''
        params(1) = -1000000.0d0
        params(2) = {universo_discurso_entrada_01[0]}d0
        params(3) = params(2) + x(10) * ({referencia_entrada_01}d0 - params(2))
        params(4) = {referencia_entrada_01}d0

        params(5) = params(3)
        params(6) = params(4)
        params(7) = params(6) + x(11) * ({universo_discurso_entrada_01[1]}d0 - params(6))

        params(8) = params(6)
        params(9) = params(7)
        params(10) = {universo_discurso_entrada_01[1]}d0
        params(11) = 1000000.0d0

        params(12) = -1000000.0d0
        params(13) = {universo_discurso_entrada_02[0]}d0
        params(14) = params(13) + x(12) * ({referencia_entrada_02}d0 - params(13))
        params(15) = {referencia_entrada_02}d0

        params(16) = params(14)
        params(17) = params(15)
        params(18) = params(17) + x(13) * ({universo_discurso_entrada_02[1]}d0 - params(17))

        params(19) = params(17)
        params(20) = params(18)
        params(21) = {universo_discurso_entrada_02[1]}d0
        params(22) = 1000000.0d0
        
        params(23) = x(14) * ({universo_discurso_salida[0]}d0)
        params(24) = params(23)
        params(25) = 0.0d0
        params(26) = 0.0d0
        params(27) = x(14) * ({universo_discurso_salida[1]}d0)
        params(28) = params(27)
'''
    guardado = f'''
  write(10, '(A)') "x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,x_12,x_13,x_14,Costo,Times"

  do iteracion = 1, maxima_iteracion
    !             1        2         3         4         5         6         7         8        9          10        11         12        13       14        15         16       17        18        19        20         21       22        23         24        25        26                              
    write(10, '(F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7)') &
    soluciones(iteracion, 1), soluciones(iteracion, 2), soluciones(iteracion, 3), soluciones(iteracion, 4), &
    soluciones(iteracion, 5), soluciones(iteracion, 6), soluciones(iteracion, 7), soluciones(iteracion, 8), &
    soluciones(iteracion, 9), soluciones(iteracion,10), soluciones(iteracion,11), soluciones(iteracion,12), &
    soluciones(iteracion,13), soluciones(iteracion,14), soluciones(iteracion,15), soluciones(iteracion,16)
  end do
  close(10)
'''
if cantidad_variables == 7:
    call_reglas = f"call base_de_reglas(anint({vector_reglas}), base)"
    limites_reglas = ""
    definicion_params = f'''
        params(1) = -1000000.0d0
        params(2) = {universo_discurso_entrada_01[0]}d0
        params(3) = params(2) + x(1) * ({universo_discurso_entrada_01[1]}d0 - params(2))
        params(4) = params(3) + x(2) * ({universo_discurso_entrada_01[1]}d0 - params(3))

        params(5) = params(3)
        params(6) = params(4)
        params(7) = params(6) + x(3) * ({universo_discurso_entrada_01[1]}d0 - params(6))

        params(8) = params(6)
        params(9) = params(7)
        params(10) = {universo_discurso_entrada_01[1]}d0
        params(11) = 1000000.0d0

        params(12) = -1000000.0d0
        params(13) = {universo_discurso_entrada_02[0]}d0
        params(14) = params(13) + x(4) * ({universo_discurso_entrada_02[1]}d0 - params(13))
        params(15) = params(14) + x(5) * ({universo_discurso_entrada_02[1]}d0 - params(14))

        params(16) = params(14)
        params(17) = params(15)
        params(18) = params(17) + x(6) * ({universo_discurso_entrada_02[1]}d0 - params(17))

        params(19) = params(17)
        params(20) = params(18)
        params(21) = {universo_discurso_entrada_02[1]}d0
        params(22) = 1000000.0d0
        
        params(23) = x(7) * ({universo_discurso_salida[0]}d0)
        params(24) = params(23)
        params(25) = 0.0d0
        params(26) = 0.0d0
        params(27) = x(7) * ({universo_discurso_salida[1]}d0)
        params(28) = params(27)
'''
    guardado = f'''
  write(10, '(A)') "x_1,x_2,x_3,x_4,x_5,x_6,x_7,Costo,Times"

  do iteracion = 1, maxima_iteracion
    !             1        2         3         4         5         6         7         8        9                                         
    write(10, '(F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7,",",F15.7)') &
    soluciones(iteracion, 1), soluciones(iteracion, 2), soluciones(iteracion, 3), soluciones(iteracion, 4), &
    soluciones(iteracion, 5), soluciones(iteracion, 6), soluciones(iteracion, 7), soluciones(iteracion, 8), &
    soluciones(iteracion, 9)
  end do
  close(10)
'''
if not reglas_discretas:
    call_reglas = call_reglas.replace('anint(', '').replace("))",")")



tsukamoto = f'''
module tsukamoto
    implicit none !Evita que hayan variables implícitas
contains !Todas las funciones y subrutinas 

    function entrenamiento(x) result(objetivo)
        double precision, dimension(:), intent(in) :: x
        double precision, dimension({numero_datos_entrenamiento}) :: x1, x2, respuesta
        double precision, dimension(3,3,3) :: base
        double precision, dimension(28) :: params
        double precision :: suma_cuadrados, prediccion_0, suma_errores, promedio_respuesta, coeficiente_determinacion, objetivo
        integer :: i

        x1 = {list(df_entrenamiento.x_1)}
        x2 = {list(df_entrenamiento.x_2)}
        respuesta = {list(df_entrenamiento.respuesta)}

        suma_cuadrados = 0.0d0 
        suma_errores = 0.0d0
        promedio_respuesta = {float(df_entrenamiento.respuesta.mean())}d0
        {call_reglas}

{definicion_params}

        do i = 1, {numero_datos_entrenamiento}
            prediccion_0 = x1(i) * (1.0d0 - {float(osa)}d0) + fis_tsukamoto(base, params, x1(i), x2(i))
            suma_cuadrados = suma_cuadrados + (respuesta(i) - prediccion_0)**2
            suma_errores = suma_errores + (respuesta(i) - promedio_respuesta)**2
        end do
        coeficiente_determinacion = 1 - suma_cuadrados / suma_errores
        objetivo = -coeficiente_determinacion

    end function entrenamiento

    function verificador(x) result(objetivo)
        double precision, dimension(:), intent(in) :: x
        double precision :: suma_cuadrados, prediccion_0, prediccion_1, prediccion_2, suma_errores, promedio_respuesta, coeficiente_determinacion, objetivo
        double precision, dimension({numero_datos_verificacion}) :: x1, x2, respuesta, prediccion
        double precision, dimension(3,3,3) :: base
        double precision, dimension(28) :: params
        integer :: i, unit, iteracion
        character(len=100) :: nombre_archivo

        x1 = {list(df_verificacion.x_1)}
        x2 = {list(df_verificacion.x_2)}
        respuesta = {list(df_verificacion.respuesta)}

        suma_cuadrados = 0.0d0 
        suma_errores = 0.0d0 
        promedio_respuesta = {float(df_verificacion.respuesta.mean())}d0

        {call_reglas}

{definicion_params}                                    

        prediccion_1 = x1(1)

        do i = 1, {numero_datos_verificacion}
            prediccion_0 = prediccion_1 * (1.0d0 - {float(osa)}d0) + fis_tsukamoto(base, params, prediccion_1, x2(i))
            prediccion(i) = prediccion_0
            suma_cuadrados = suma_cuadrados + (respuesta(i) - prediccion_0)**2
            suma_errores = suma_errores + (respuesta(i) - promedio_respuesta)**2
            prediccion_1 = prediccion_0
        end do
        coeficiente_determinacion = 1 - suma_cuadrados / suma_errores
        objetivo = -coeficiente_determinacion
        
    end function verificador

    function verificacion(x) result(objetivo)
        double precision, dimension(:), intent(in) :: x
        double precision :: suma_cuadrados, prediccion_0, prediccion_1, prediccion_2, suma_errores, promedio_respuesta, coeficiente_determinacion, objetivo
        double precision, dimension({numero_datos_verificacion}) :: x1, x2, respuesta, prediccion
        double precision, dimension(3,3,3) :: base
        double precision, dimension(28) :: params
        integer :: i, unit, iteracion
        character(len=100) :: nombre_archivo

        x1 = {list(df_verificacion.x_1)}
        x2 = {list(df_verificacion.x_2)}
        respuesta = {list(df_verificacion.respuesta)}

        suma_cuadrados = 0.0d0 
        suma_errores = 0.0d0 
        promedio_respuesta = {float(df_verificacion.respuesta.mean())}d0

        {call_reglas}

{definicion_params}                                    

        prediccion_1 = x1(1)

        do i = 1, {numero_datos_verificacion}
            prediccion_0 = prediccion_1 * (1.0d0 - {float(osa)}d0) + fis_tsukamoto(base, params, prediccion_1, x2(i))
            prediccion(i) = prediccion_0
            suma_cuadrados = suma_cuadrados + (respuesta(i) - prediccion_0)**2
            suma_errores = suma_errores + (respuesta(i) - promedio_respuesta)**2
            prediccion_1 = prediccion_0
        end do
        coeficiente_determinacion = 1 - suma_cuadrados / suma_errores
        objetivo = -coeficiente_determinacion


        nombre_archivo = "verificacion.csv"
        open(unit=10, file=trim(nombre_archivo), status="replace")
        ! Escribe encabezados en el archivo CSV
            write(10, '(A)') "x_1,x_2,respuesta,prediccion"

            do iteracion = 1, {numero_datos_verificacion}
                !             1        2         3         4                                     
                write(10, '(F15.7,",",F15.7,",",F15.7,",",F15.7)') &
                x1(iteracion), x2(iteracion), respuesta(iteracion), prediccion(iteracion)
            end do
        close(10)
        
    end function verificacion

    function validacion(x) result(objetivo)
        double precision, dimension(:), intent(in) :: x
        double precision :: suma_cuadrados, prediccion_0, prediccion_1, prediccion_2, suma_errores, promedio_respuesta, coeficiente_determinacion, objetivo
        double precision, dimension({numero_datos_validacion}) :: x1, x2, respuesta, prediccion
        double precision, dimension(3,3,3) :: base
        double precision, dimension(28) :: params
        integer :: i, unit, iteracion
        character(len=100) :: nombre_archivo

        x1 = {list(df_validacion.x_1)}
        x2 = {list(df_validacion.x_2)}
        respuesta = {list(df_validacion.respuesta)}

        suma_cuadrados = 0.0d0 
        suma_errores = 0.0d0 
        promedio_respuesta = {float(df_validacion.respuesta.mean())}d0

        {call_reglas}

{definicion_params}

        prediccion_1 = x1(1)

        do i = 1, {numero_datos_validacion}
            prediccion_0 = prediccion_1 * (1.0d0 - {float(osa)}d0) + fis_tsukamoto(base, params, prediccion_1, x2(i))
            prediccion(i) = prediccion_0
            suma_cuadrados = suma_cuadrados + (respuesta(i) - prediccion_0)**2
            suma_errores = suma_errores + (respuesta(i) - promedio_respuesta)**2
            prediccion_1 = prediccion_0
        end do
        coeficiente_determinacion = 1 - suma_cuadrados / suma_errores
        objetivo = -coeficiente_determinacion


        nombre_archivo = "validacion.csv"
        open(unit=10, file=trim(nombre_archivo), status="replace")
        ! Escribe encabezados en el archivo CSV
            write(10, '(A)') "x_1,x_2,respuesta,prediccion"

            do iteracion = 1, {numero_datos_validacion}
                !             1        2         3         4                                     
                write(10, '(F15.7,",",F15.7,",",F15.7,",",F15.7)') &
                x1(iteracion), x2(iteracion), respuesta(iteracion), prediccion(iteracion)
            end do
        close(10)
        
    end function validacion
    
    subroutine base_de_reglas(reglas, base)
        double precision, dimension(9), intent(in) :: reglas
        double precision, dimension(3, 3, 3), intent(out) :: base
        integer :: fila, columna, idx

        ! reglas son los params

        ! base(entrada 1, entrada 2, salida)
        ! 1: Negativo
        ! 2: Zero
        ! 3: Positivo

        ! r1 r4 r7
        ! r2 r5 r8
        ! r3 r6 r9

        columna = 1
        fila = 1

        do idx = 1, 9
            base(fila, columna, 1) = trapezoidal([-10.0d0, -9.0d0, 1.0d0, 2.0d0], reglas(idx))
            base(fila, columna, 2) = triangular([1.0d0, 2.0d0, 3.0d0], reglas(idx))
            base(fila, columna, 3) = trapezoidal([2.0d0, 3.0d0, 9.0d0, 10.0d0], reglas(idx))

            fila = fila + 1
            if (fila > 3) then 
                fila = 1
                columna = columna + 1
            end if
        end do
    end subroutine base_de_reglas

    double precision function fis_tsukamoto(base, params, entrada_01, entrada_02)
        double precision, dimension(28), intent(in) :: params
        double precision, dimension(3, 3, 3), intent(in) :: base
        double precision, intent(in) :: entrada_01, entrada_02
        double precision, dimension(3) :: antecedente_01, antecedente_02
        double precision :: implicacion_larsen, numerador, denominador
        integer :: i, j, k

        ! Corregir valores en trapezoidal
        antecedente_01(1) = trapezoidal(params(1:4), entrada_01)
        antecedente_01(2) = triangular(params(5:7), entrada_01)
        antecedente_01(3) = trapezoidal(params(8:11), entrada_01)

        antecedente_02(1) = trapezoidal(params(12:15), entrada_02)
        antecedente_02(2) = triangular(params(16:18), entrada_02)
        antecedente_02(3) = trapezoidal(params(19:22), entrada_02)

        denominador = 0.0d0
        numerador = 0.0d0

        do i = 1, 3
            do j = 1, 3
                implicacion_larsen = antecedente_01(i) * antecedente_02(j)
                denominador = denominador + implicacion_larsen 
                do k = 1, 3
                    select case (k)
                        case (1)
                            numerador = numerador + implicacion_larsen * lineal(params(23), params(24), implicacion_larsen) * base(i, j, k)
                        case (2)
                            numerador = numerador + implicacion_larsen * lineal(params(25), params(26), implicacion_larsen) * base(i, j, k)
                        case (3)
                            numerador = numerador + implicacion_larsen * lineal(params(27), params(28), implicacion_larsen) * base(i, j, k)
                    end select 
                end do 
            end do 
        end do 

        ! Evitar división por un número cercano a cero
        denominador = merge(1.0d-10, denominador, abs(denominador) < 1.0d-10)

        fis_tsukamoto = numerador / denominador
    end function fis_tsukamoto

    double precision function trapezoidal(params, w)
        double precision, dimension(4), intent(in) :: params
        double precision, intent(in) :: w
        double precision :: a, b, c, d
        a = params(1)
        b = params(2)
        c = params(3)
        d = params(4)
        trapezoidal = max(min((w - a) / (b - a), min(1.0, (d - w) / (d - c))), 0.0)
    end function trapezoidal 

    double precision function triangular(params, w)
        double precision, dimension(3), intent(in) :: params
        double precision, intent(in) :: w 
        double precision :: a, b, c
        a = params(1)
        b = params(2)
        c = params(3)
        triangular = max(min((w - a) / (b - a), (c - w) / (c - b)), 0.0)
    end function triangular 


    double precision function lineal(x0, x1, y)
        double precision, intent(in) :: x0, x1, y

        lineal = y * (x1 - x0) + x0

    end function lineal

end module tsukamoto
'''


with open("tsukamoto.f90", "w", encoding="utf-8") as f:
    f.write(tsukamoto)

optimizacion = f'''
program OptimizationDriver
! gfortran -c metaheuristicas.f90
! gfortran -c tsukamoto.f90
! gfortran -c OptimizationDriver.f90
! gfortran metaheuristicas.o tsukamoto.o OptimizationDriver.o -o driver.exe

  use metaheuristicas
  use tsukamoto
  implicit none

  integer :: iteracion, tamanio_poblacion, maxima_iteracion_algoritmo, maxima_iteracion
  double precision, dimension({cantidad_variables}) :: limite_inferior, limite_superior
  double precision, dimension({pruebas},{cantidad_variables+2}) :: soluciones
  double precision, dimension({cantidad_variables + 2}) :: solucion
  character(len=100) :: nombre_archivo

  ! Parámetros
  tamanio_poblacion = {tamanio_poblacion}
  maxima_iteracion_algoritmo = {iteraciones}
  limite_inferior = 0.0d0
  limite_superior = 1.0d0

{limites_reglas}

  maxima_iteracion = {pruebas}

  ! Definir el nombre del archivo como una variable
  nombre_archivo = "resultados.csv"
  
  solucion = 0.0        

  ! Realizar iteraciones y guardar soluciones
  do iteracion = 1, maxima_iteracion

    solucion = differential_evolution(CostFunction=entrenamiento, &
                                      LimInf=limite_inferior, &
                                      LimSup=limite_superior, &
                                      NumPop=tamanio_poblacion, &
                                      MaxIter=maxima_iteracion_algoritmo, &
                                      minimo_global_deseado = {minimo_global_deseado}d0)


    soluciones(iteracion, :) = solucion

    print *, "Resultado de la iteración ", iteracion, ": ", solucion
  end do

  ! Guarda los resultados en un archivo CSV
  open(unit=10, file=trim(nombre_archivo), status="replace")
  ! Escribe encabezados en el archivo CSV

{guardado}

  print *, "Resultados guardados en: ", trim(nombre_archivo)

end program OptimizationDriver
'''

with open("OptimizationDriver.f90", "w", encoding="utf-8") as f:
    f.write(optimizacion)

verificaciones = """
import pandas as pd

df = pd.read_csv("resultados.csv")
numero_resultados = df.shape[0]
matriz = [f'    x({i},:) = {list(df.iloc[i-1,0:"""+str(cantidad_variables)+"""])}' for i in range(1,numero_resultados+1)]

completo = []
completo.append(f'''
program verificaciones

    use tsukamoto
    implicit none

    double precision, dimension({numero_resultados},"""+str(cantidad_variables)+""") :: x 
    double precision, dimension({numero_resultados}) :: rs
    character(len=100) :: nombre_archivo
    integer :: i
                
''')

completo.extend(matriz)

completo.append(
    f'''    
    
    do i = 1, {numero_resultados}
        rs(i) = -verificador(x(i,:))
    end do

    nombre_archivo = "rs.csv"
    open(unit=10, file=trim(nombre_archivo), status="replace")
    ! Escribe encabezados en el archivo CSV
        write(10, '(A)') "r_square"

        do i = 1, {numero_resultados}
            write(10, '(F15.7)') rs(i)
        end do
    close(10)
    
end program verificaciones''')

with open("verificaciones.f90", "w", encoding="utf-8") as f:
    f.write("\\n".join(completo))

"""

with open("verificaciones.py", "w", encoding="utf-8") as f:
    f.write(verificaciones)

unir_resultados = f'''
import pandas as pd
rs = pd.read_csv("rs.csv")
df = pd.read_csv("resultados.csv")
df["rs"] = rs
df.sort_values(by="rs", inplace = True)
df.to_csv('{carpeta + "resultados.csv"}', index=False)
'''

with open("unir_resultados.py", "w", encoding="utf-8") as f:
    f.write(unir_resultados)

datos_vv = '''
import pandas as pd
df = pd.read_csv("''' + carpeta + "resultados.csv" + '''")
x = list(df.iloc[-1,0:'''+str(cantidad_variables)+'''])
verificaciones = f"""
program verificaciones
! gfortran -c tsukamoto.f90
! gfortran -c verificaciones.f90
! gfortran tsukamoto.o verificaciones.o -o verificaciones.exe

    use tsukamoto
    implicit none

    double precision, dimension('''+str(cantidad_variables)+''') :: x 
    double precision :: coeficiente_correlacion

    x = {x}

    coeficiente_correlacion = -verificacion(x)
    print *, coeficiente_correlacion

    coeficiente_correlacion = -validacion(x)
    print *, coeficiente_correlacion
    
end program verificaciones
"""
with open("verificaciones.f90", "w", encoding = "utf-8") as f:
    f.write(verificaciones)
'''

with open("datos_vv.py", "w", encoding="utf-8") as f:
    f.write(datos_vv)

