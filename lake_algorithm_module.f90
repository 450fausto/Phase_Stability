module lake_algorithm_module
    implicit none
    contains

    function lake_algorithm(funcion_costo, limite_inferior, limite_superior, iteracion_maxima) result(solution)
        ! Declaración de variables
        integer, intent(in) :: iteracion_maxima
        double precision, intent(in) :: limite_inferior(:), limite_superior(:)
        interface
            double precision function funcion_costo(x)
                double precision, dimension(:), intent(in) :: x
            end function funcion_costo
        end interface

        integer :: numero_variables, intento_maximo, iter, numero_lagos
        double precision :: paso_inicial, tamanio_paso_maximo, tamanio_paso_minimo, numero_aleatorio
        double precision :: paso_dinamico, radio_inicial, tasa_crecimiento, tasa_reduccion, tasa_aumento
        double precision :: probabilidad_cambio, rango_cambio, costo, minimo_local, minimo_global, Ravoid
        double precision :: probabilidad_desaparicion, probabilidad_valor
        double precision, dimension(:), allocatable :: tamanio_espacio, direccion_referencia, solution, x_minimo_global(:)
        double precision, dimension(:), allocatable :: punto_a, punto_b, direccion, nuevo_lago, vector_aleatorio
        double precision, dimension(:, :), allocatable :: ubicacion_lagos
        integer :: fracaso_mejora, fracaso_reduccion

        ! Inicialización de parámetros
        numero_variables = size(limite_inferior)
        intento_maximo = numero_variables * (ceiling(51.0 * exp(-0.5 * numero_variables)) + 3)
        paso_inicial = abs(sum(limite_superior - limite_inferior) / numero_variables)
        tamanio_paso_maximo = paso_inicial
        tamanio_paso_minimo = 1.0e-5
        radio_inicial = paso_inicial / 3.0
        tasa_crecimiento = radio_inicial / 2.0
        tasa_reduccion = (1.0 / 6.0) * (1.0 / 5.0)**(1.0 / intento_maximo)
        tasa_aumento = 0.7 * tasa_reduccion
        probabilidad_cambio = 2.0 * ((real(numero_variables - 1) / real(numero_variables))**numero_variables)
        probabilidad_desaparicion = 0.2
        probabilidad_valor = 0
        numero_lagos = 0
        rango_cambio = 0.3

        ! Inicialización de valores
        fracaso_mejora = 0
        fracaso_reduccion = 1
        Ravoid = 0.0

        allocate(tamanio_espacio(numero_variables))
        allocate(direccion_referencia(numero_variables))
        allocate(punto_a(numero_variables))
        allocate(punto_b(numero_variables))
        allocate(vector_aleatorio(numero_variables))
        allocate(x_minimo_global(numero_variables))
        allocate(nuevo_lago(numero_variables + 1))
        allocate(solution(numero_variables + 1))
        allocate(ubicacion_lagos(int(iteracion_maxima / intento_maximo), numero_variables + 1))

        tamanio_espacio = (limite_superior - limite_inferior) / maxval(limite_superior - limite_inferior)
        paso_dinamico = paso_inicial
        call random_number(vector_aleatorio)
        punto_a = limite_inferior + (limite_superior - limite_inferior) * vector_aleatorio
        call random_number(vector_aleatorio)
        direccion_referencia = normalizar(vector_aleatorio)
        x_minimo_global = punto_a
        minimo_local = funcion_costo(punto_a)
        minimo_global = minimo_local

        ! Iteraciones
        do iter = 1, iteracion_maxima
            if (fracaso_mejora > 1) then
                paso_dinamico = paso_dinamico * (1.0 - tasa_reduccion)
                paso_dinamico = max(min(paso_dinamico, tamanio_paso_maximo), tamanio_paso_minimo)
                direccion_referencia = cambio_signo_valor(direccion_referencia, probabilidad_cambio, probabilidad_valor)
            else
                paso_dinamico = (1.0 + tasa_aumento) * paso_dinamico
                paso_dinamico = max(min(paso_dinamico, tamanio_paso_maximo), tamanio_paso_minimo)
            end if

            call random_number(vector_aleatorio)
            direccion = direccion_referencia + rango_cambio * (2.0 * vector_aleatorio - 1.0)
            direccion = normalizar(direccion)
            punto_b = punto_a + paso_dinamico * direccion * tamanio_espacio
            punto_b = desplazar_punto(punto_b, limite_superior, limite_inferior, ubicacion_lagos, Ravoid, punto_a)
            punto_b = max(limite_inferior, min(limite_superior, punto_b))

            costo = funcion_costo(punto_b)

            if (costo < minimo_local) then
                punto_a = punto_b
                minimo_local = costo
                direccion_referencia = direccion
                if (costo < minimo_global) then
                    minimo_global = costo
                    x_minimo_global = punto_b
                end if
                fracaso_mejora = 0
                fracaso_reduccion = 0
            else if (costo > minimo_local) then
                fracaso_mejora = fracaso_mejora + 1
                fracaso_reduccion = fracaso_reduccion + 1
            else
                fracaso_mejora = 0
                fracaso_reduccion = fracaso_reduccion + 1
                direccion_referencia = cambio_signo_valor(direccion_referencia, probabilidad_cambio, probabilidad_valor)
            end if

            if (fracaso_reduccion > intento_maximo) then
                if (sqrt(sum((punto_a - x_minimo_global) ** 2)) > radio_inicial) then 
                    nuevo_lago(1: numero_variables) = punto_a 
                    nuevo_lago(numero_variables + 1) = radio_inicial 
                    numero_lagos = numero_lagos + 1
                    ubicacion_lagos(numero_lagos, :) = nuevo_lago
                    ubicacion_lagos(1:numero_lagos, :) = aumentar_radio(ubicacion_lagos(1:numero_lagos, :), tasa_crecimiento)
                    numero_lagos = merge(numero_lagos, numero_lagos - 1, ubicacion_lagos(numero_lagos, numero_variables) > 0)
                    call random_number(numero_aleatorio)
                    if (numero_aleatorio < probabilidad_desaparicion) then 
                        ubicacion_lagos(1:numero_lagos, :) = eliminar_fila(ubicacion_lagos(1:numero_lagos, :))
                        numero_lagos = numero_lagos - 1
                    end if
                end if
                call random_number(vector_aleatorio)
                punto_a = (limite_superior - limite_inferior) * vector_aleatorio + limite_inferior
                paso_dinamico = paso_inicial
                minimo_local = funcion_costo(punto_a)
                fracaso_reduccion = 0
                Ravoid = radio_inicial
                fracaso_mejora = 0
            end if
        end do

        ! Asignar la solución final
        solution(1:numero_variables) = x_minimo_global
        solution(numero_variables + 1) = minimo_global
        
        ! do iter = 1, int(iteracion_maxima / intento_maximo)
        !     print *, ubicacion_lagos(iter, :)
        ! end do 
    end function lake_algorithm

    function cambio_signo_valor(direccion, probabilidad_signo, probabilidad_valor) result(nuevo_vector)
        double precision, intent(in) :: direccion(:)
        double precision, intent(in) :: probabilidad_signo, probabilidad_valor
        double precision :: nuevo_vector(size(direccion))
        integer :: i
        double precision :: aleatorio

        do i = 1, size(direccion)
            call random_number(aleatorio)
            if (aleatorio < probabilidad_signo) then
                nuevo_vector(i) = -direccion(i)
            else
                nuevo_vector(i) = direccion(i)
            end if

            call random_number(aleatorio)
            if (aleatorio < probabilidad_valor) then
                nuevo_vector(i) = SIGN(aleatorio, direccion(i))

            end if
        end do

        nuevo_vector = normalizar(nuevo_vector)
    end function cambio_signo_valor

    function normalizar(vector) result(normalized_vector)
        implicit none
        double precision, intent(in) :: vector(:)
        double precision :: normalized_vector(size(vector))
        double precision :: norm

        norm = sqrt(sum(vector**2))
        if (norm > 0.0) then
            normalized_vector = vector / norm
        else
            normalized_vector = vector
        end if
    end function normalizar

    function aumentar_radio(ubicacion_lagos, tasa_crecimiento) result(nueva_matriz)
        implicit none
        double precision, intent(in) :: ubicacion_lagos(:, :), tasa_crecimiento
        double precision, dimension(:), allocatable :: nueva_lista
        double precision :: radios(size(ubicacion_lagos, 1))
        double precision, allocatable :: ubicaciones(:, :)
        integer :: i, filas, columnas
        double precision :: distancia_interna
        double precision, allocatable :: nueva_matriz(:,:)
        logical :: dentro, conservar_ultima

        filas = size(ubicacion_lagos, 1)
        columnas = size(ubicacion_lagos, 2)

        allocate(ubicaciones(filas, columnas - 1))
        allocate(nueva_lista(columnas))
        allocate(nueva_matriz(filas, columnas))
        nueva_matriz = 0
        ubicaciones = ubicacion_lagos(:, 1:columnas - 1)
        radios = ubicacion_lagos(:, columnas)
        conservar_ultima = .true.

        if (filas > 0) then 
            do i = 1, filas - 1
                distancia_interna = sqrt(sum((ubicaciones(filas, :) - ubicaciones(i, :))**2))
                dentro = distancia_interna <= radios(i)
                radios(i) = radios(i) + tasa_crecimiento * merge(1.0, 0.0, dentro)
                conservar_ultima = merge(.false., .true., conservar_ultima .and. dentro)
                if (radios(i) > 0) then 
                    nueva_lista(1: columnas - 1) = ubicaciones(i, :)
                    nueva_lista(columnas) = radios(i)
                    nueva_matriz(i, :) = nueva_lista
                end if
            end do
            if (conservar_ultima) then 
                nueva_matriz(filas, :) = ubicacion_lagos(filas, :)
            end if
        else 
            nueva_matriz = ubicacion_lagos 
        end if
    end function aumentar_radio

    function desplazar_punto(punto_b, limite_superior, limite_inferior, ubicaciones_radios, radio, punto_a) result(punto_c)
        double precision, intent(in) :: punto_b(:), limite_superior(:), limite_inferior(:)
        double precision, intent(in) :: ubicaciones_radios(:, :), radio
        double precision, intent(in) :: punto_a(:)
        double precision :: punto_c(size(punto_b))
        integer :: i, salir, cantidad_filas, cantidad_columnas, todo_fuera
        double precision :: distancia, vector_ab(size(punto_b)), vector_ab_unitario(size(punto_b)), vector_aleatorio(size(punto_b))
        double precision :: aplazamiento, radio_actual

        cantidad_filas = size(ubicaciones_radios, 1)
        cantidad_columnas = size(ubicaciones_radios, 2)
        punto_c = punto_b
        salir = 0

        if (radio > 0 .and. cantidad_filas > 0) then 
            todo_fuera = 0
            do while (salir < 10 .and. todo_fuera == 0)
                salir = salir + 1
                do i = 1, cantidad_filas
                    radio_actual = ubicaciones_radios(i, cantidad_columnas)
                    distancia = sqrt(sum((ubicaciones_radios(i, 1:size(punto_b)) - punto_c)**2))
                    todo_fuera = 1
                    if (distancia <= radio_actual) then
                        vector_ab = punto_c - punto_a
                        vector_ab_unitario = normalizar(vector_ab)
                        aplazamiento = 2.0 * sqrt(radio_actual**2 - distancia**2)
                        punto_c = punto_c + aplazamiento * vector_ab_unitario
                        todo_fuera = 0
                    end if
                end do
            end do
        end if

        if (any(punto_c < limite_inferior) .or. any(punto_c > limite_superior) .or. salir == 10) then
            call random_number(vector_aleatorio)
            punto_c = limite_inferior + (limite_superior - limite_inferior) * vector_aleatorio
        end if
    end function desplazar_punto

    function eliminar_fila(matriz) result(matriz_actualizada)
        implicit none
        double precision, intent(in) :: matriz(:,:)      
        double precision, allocatable :: matriz_actualizada(:,:)  
        integer :: filas, columnas, fila_eliminar
        double precision :: numero_aleatorio

        filas = size(matriz, 1)
        columnas = size(matriz, 2)

        call random_number(numero_aleatorio)
        fila_eliminar = ceiling(numero_aleatorio * real(filas))

        allocate(matriz_actualizada(filas, columnas))
        matriz_actualizada = 0.0

        if (fila_eliminar == 1) then 
            matriz_actualizada(1: filas-1, :) = matriz(2:filas, :)
        else if (fila_eliminar == filas) then 
            matriz_actualizada(1: filas-1, :) = matriz(1:filas-1, :)
        else
            matriz_actualizada(1:fila_eliminar-1, :) = matriz(1: fila_eliminar-1, :)
            matriz_actualizada(fila_eliminar: filas -1, :) = matriz(fila_eliminar + 1: filas, :)
        end if
    end function eliminar_fila


end module lake_algorithm_module 