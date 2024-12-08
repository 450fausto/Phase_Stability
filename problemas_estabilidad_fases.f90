module problemas_estabilidad_fases
  implicit none

contains
    function cardano(a, b, c, d) result(raices)
        double precision, intent(in) :: a, b, c, d 
        double precision p, q, delta, u, v, tetha
        double precision, dimension(3) :: raices
        double precision, parameter :: pi = 3.141592653589793d0

        raices = 0.0

        p = (3 * a * c - b ** 2) / (3 * a ** 2)
        q = (2 * b ** 3 - 9 * a * b * c + 27 * a ** 2 * d) / (27 * a ** 3)

        delta = q ** 2 + (4 * p ** 3) / 27

        if (delta > 0) then 
            u = sign(abs((-q + sqrt(delta)) / 2)**(1.0d0/3.0d0), (-q + sqrt(delta)) / 2)
            v = sign(abs((-q - sqrt(delta)) / 2)**(1.0d0/3.0d0), (-q - sqrt(delta)) / 2)
            raices = u + v - b / (3 * a)
        else if (abs(delta) < 1.0d-12) then 
            if (p == 0) then 
                raices = 0.0
            else
                raices(1) = 3 * q / p
                raices(2) = -3 * q / (2 * p)
                raices(3) = -3 * q / (2 * p)
            end if
        else 
            tetha = acos((3 * q / (2 * p)) * sqrt(-3 / p))
            raices(1) = 2 * sqrt(-p / 3) * cos(tetha / 3) - b / (3 * a)
            raices(2) = 2 * sqrt(-p / 3) * cos((tetha + 2 * pi) / 3) - b / (3 * a)
            raices(3) = 2 * sqrt(-p / 3) * cos((tetha + 4 * pi) / 3) - b / (3 * a)
        
        end if 
    end function cardano

    function ln_gamma_uniquac(x, r, q, qq, tau) result(ln_gamma)
        double precision, dimension(:), allocatable :: phi, termino_a, termino_b, ln_gamma
        double precision, dimension(:), allocatable :: termino_c, termino_d, termino_e, L, tetha, tetha_p
        double precision, dimension(:), intent(in) :: x, r, q, qq 
        double precision, dimension(:, :), intent(in) :: tau 
        double precision, parameter :: z = 10
        integer :: i, j, cantidad_variables
        double precision :: suma

        cantidad_variables = size(x)

        allocate(ln_gamma(cantidad_variables))
        allocate(phi(cantidad_variables))
        allocate(tetha(cantidad_variables))
        allocate(tetha_p(cantidad_variables))
        allocate(L(cantidad_variables))
        allocate(termino_a(cantidad_variables))
        allocate(termino_b(cantidad_variables))
        allocate(termino_c(cantidad_variables))
        allocate(termino_d(cantidad_variables))
        allocate(termino_e(cantidad_variables))

        phi = r * x / sum(x * r)
        tetha = q * x / sum(q * x)
        tetha_p = qq * x / sum(qq * x)
        L = (z / 2) * (r - q) - (r - 1 )

        termino_a = log(phi / x)
        termino_b = (z / 2) * q * log(tetha / phi) 
        termino_c = - (phi / x) * sum(x * L)

        termino_d = 0
        termino_e = 0
        do i = 1, cantidad_variables
            termino_d(i) = - qq(i) * log(sum(tetha_p * tau(i, :)))
            suma = 0.0
            do j = 1, cantidad_variables 
                suma = suma + tetha_p(j) * tau(j, i) / sum(tetha_p * tau(j, :))
            end do 
            termino_e(i) = - qq(i) * suma
        end do

        !ADVERTENCIA: cuando se usa shape en tau declarada como parameter, se convierte en la matriz
        !transpuesta de tau, por esta razón se usan los índices al revés. ej. en lugar de tau_ij, 
        !se tiene que usar tau_ji

        ln_gamma = termino_a + termino_b + L + termino_c + termino_d + qq + termino_e
    end function ln_gamma_uniquac

    function ln_gamma_nrtl(x, tau, G) result(ln_gamma)
        double precision, dimension(:), allocatable :: ln_gamma
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(:, :), intent(in) :: tau, G
        double precision, parameter :: z = 10
        integer :: i, j, cantidad_variables
        double precision :: suma, sumb, sumc

        cantidad_variables = size(x)
        allocate(ln_gamma(cantidad_variables))

        do i = 1, cantidad_variables
            suma = 0.0
            do j = 1, cantidad_variables
                sumb = sum(G(j, :) * x)
                sumc = sum(tau(j, :) * G(j, :) * x)
                suma = suma + x(j) * (tau(j,i) * G(j, i) * sumb - G(j, i) * sumc) / (sumb ** 2)

            end do
            ln_gamma(i) = sum(tau(i, :) * G(i, :) * x) / sum(G(i, :) * x) + suma
        end do

    end function ln_gamma_nrtl

    function ln_phi_srk(x, Tc, Pc, w, k, P, T) result(ln_phi)
        double precision, dimension(:), allocatable :: ln_phi 
        double precision, intent(in) :: P, T
        double precision, dimension(:), intent(in) :: x, tc, pc, w
        double precision, dimension(:, :), intent(in) :: k
        double precision, dimension(:), allocatable :: alpha, A, B, parte_a, parte_b
        double precision, dimension(3) :: raices
        double precision, dimension(:, :), allocatable :: Am
        integer :: cantidad_variables, i, j
        double precision :: presion, temperatura, Bmix, Amix, z
        double precision :: termino_a, termino_b, termino_c, termino_d

        cantidad_variables = size(x)

        allocate(alpha(cantidad_variables))
        allocate(A(cantidad_variables))
        allocate(Am(cantidad_variables, cantidad_variables))
        allocate(B(cantidad_variables))
        allocate(parte_a(cantidad_variables))
        allocate(parte_b(cantidad_variables))

        alpha = (1.0 + (0.480 + 1.574 * w - 0.176 * w ** 2) * (1 - sqrt(T / Tc))) ** 2
        A = 0.42747 * alpha * (P * Tc ** 2 / (Pc * T ** 2))
        B = 0.08664 * (P * Tc / (Pc * T))
        Bmix = sum(x * B)
        Amix = 0.0

        do i = 1, cantidad_variables
            do j = 1, cantidad_variables
                Am(i, j) = (1 - k(i, j)) * sqrt(A(i) * A(j))
                Amix = Amix + x(i) * x(j) * Am(i, j)
            end do 
        end do 

        termino_a = 1
        termino_b = -1
        termino_c = Amix - Bmix - Bmix ** 2
        termino_d = -Amix * Bmix

        raices = cardano(termino_a, termino_b, termino_c, termino_d)
        z = minval(raices)

        parte_a = (B / Bmix) * (z - 1) - log(z - Bmix) 

        do i = 1, cantidad_variables
            parte_b(i) = (Amix / Bmix) * ((B(i) / Bmix) - 2 * sum(Am(i, :) * x) / Amix) * log(1 + Bmix / Z)
        end do

        ln_phi = parte_a + parte_b

    end function ln_phi_srk

!!!! Sistemas !!!!!!!!!!!!!!!!!!!!

        ! C1 + H2S, SRK
    function sistema_01(y) result(respuesta)
        ! C1 + H2S, SRK
        ! x* = 7.6619630769793529E-002  
        ! f* = -3.9319962731680751E-003

        ! Hua, J. Z., Brennecke, J. F., & Stadtherr, M. A. (1998). 
        ! Reliable computation of phase stability using interval analysis. 
        ! Computers & Chemical Engineering, 22(9), 1207–1214. doi:10.1016/s0098-1354(98)00024-6 

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 
        double precision :: criterio
        integer, parameter :: cantidad_variables = 2
        double precision, parameter :: P = 40.53, T = 190
        double precision, parameter :: z(2) = [0.0187, 0.9813]
        double precision, parameter :: Tc(2) = [373.2, 190.6]
        double precision, parameter :: Pc(2) = [89.4, 46.0]
        double precision, parameter :: w(2) = [0.1, 0.008]
        double precision, parameter :: k(2, 2) = reshape([ &
        0.0, 0.08, &
        0.08, 0.0], shape=[2, 2])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))

        x = [y(1), 1 - y(1)]

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_01

        ! N2 + C1 + C2, SRK
    function sistema_02(y) result(respuesta)
        ! N2 + C1 + C2, SRK
        ! x* = 0.13682228197692597, 6.8392258985721613E-002  
        ! f* = -0.014901225016232623

        ! Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006). 
        ! Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures. 
        ! Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g 

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 3
        double precision :: P = 76, T = 270
        double precision, parameter :: z(3) = [0.3, 0.1, 0.6]
        double precision, parameter :: Tc(3) = [126.2, 190.564, 305.32]
        double precision, parameter :: Pc(3) = [33.9, 45.9, 48.5]
        double precision, parameter :: w(3) = [0.037, 0.011, 0.098]
        double precision, parameter :: k(3, 3) = reshape([ &
        0.0, 0.038, 0.080, &
        0.038, 0.0, 0.021, &
        0.080, 0.021, 0.0], shape=[3, 3])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))

        x = [y(1), y(2), 1 - y(1) - y(2)]

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_02

        ! n-Acetato de butilo + agua NRTL
    function sistema_04(w) result(respuesta)
        ! n-Acetato de butilo + agua NRTL
        ! x* = 0.0042086615362195347  
        ! f* = -0.032638313157465078

        ! McDonald C. M. and Floudas C. A. (1994). 
        ! Global Optimization for the phase stability problem. 
        ! AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 2
        double precision, parameter :: z(2) = [0.5, 0.5]
        double precision, parameter :: tau(2, 2) = reshape([ &
        0.0, 3.00498, &
        4.69071, 0.0], shape=[2, 2])

        double precision, parameter :: G(2, 2) = reshape([ &
        1.0, 0.30800, &
        0.15909, 1.0], shape=[2, 2])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), 1 - w(1)]

        lngx = ln_gamma_nrtl(x, tau, G)
        lngz = ln_gamma_nrtl(z, tau, G)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_04

        ! C1 + C2 + C3 + iC4 + C4 + iC5 + C5 + C6 + iC15, SRK
    function sistema_05(y) result(respuesta)
        ! C1 + C2 + C3 + iC4 + C4 + iC5 + C5 + C6 + iC15, SRK
        ! x* = 0.93301293603872271, 0.030185805930259917, 0.018572154791578286, 0.0054796462078695657, 0.0074918054834172730, 0.00061840983302726218, 0.0038403394484845519, 0.00079713025666419972
        ! x* = 0.946076, 0.043568, 0.007851, 0.000675, 0.001247, 0.000197, 0.000264, 0.000121, 0.000000
        ! f* = -1.4571297453173282

        ! Rangaiah, G. P. (2001). 
        ! Evaluation of genetic algorithms and simulated annealing for phase equilibrium and stability problems. 
        ! Fluid Phase Equilibria, 187-188, 83–109. doi:10.1016/s0378-3812(01)00528-3 

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 9
        double precision :: P = 20.1, T = 314
        double precision, parameter :: z(9) = [0.614, 0.10259, 0.04985, 0.008989, 0.02116, 0.00722, 0.01187, 0.01435, 0.16998]
                                              ! C1   +   C2  +  C3  +  iC4  +  C4  +  iC5 + C5 + C6 + iC15
        double precision, parameter :: Tc(9) = [190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 507.4, 707.0]
        double precision, parameter :: Pc(9) = [46.0, 48.84, 42.46, 36.48, 38.0, 33.84, 33.74, 29.69, 15.2]
        double precision, parameter :: w(9) = [0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.706]
        double precision, parameter :: k(9, 9) = reshape([ &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], shape=[9, 9])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))
        x = 0
        x(1: cantidad_variables -1) = [y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8)]
        x(cantidad_variables) = 1 - sum(x)

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_05

        ! C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, SRK
    function sistema_06(y) result(respuesta)
        ! C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, SRK
        ! x* = 0.61410058753526386, 0.074469074863989237, 0.048200677764320229, 0.043016539266069562, 0.031783299830690830, 0.015123336264596665, 0.033972278680988968, 0.042539383539403126, 0.048644720385782388
        ! f* = -0.0000096609553336807350

        ! Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006). 
        ! Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures. 
        ! Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g 

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 10
        double precision :: P = 191.50, T = 435.35
        double precision, parameter :: z(10) = [0.6436, 0.0752, 0.0474, 0.0412, 0.0297, 0.0138, 0.0303, 0.0371, 0.0415, 0.0402]
                                                 !  C1       C2       C3      C4       C5        C6        C7        C8        C9        C10     
        double precision, parameter :: Tc(10) = [190.564, 305.322, 369.89, 425.125, 469.7000, 507.8200, 540.1300, 568.7400, 594.5500, 617.7000]
        double precision, parameter :: Pc(10) = [45.9920, 48.7220, 42.512, 37.9600, 33.67519, 30.44115, 27.36000, 24.83591, 22.81000, 21.03000]
        double precision, parameter :: w(10) =  [0.01142, 0.09900, 0.1521, 0.20081, 0.251032, 0.300319, 0.349000, 0.397528, 0.443300, 0.488400]
        ! Rosa:
        ! double precision, parameter :: Tc(10) = [190.6, 305.4, 369.8, 425.2, 469.7, 507.5, 540.3, 568.8, 594.6, 617.7]
        ! double precision, parameter :: Pc(10) = [46.0,  48.80, 42.50, 38.00, 33.70, 30.10, 27.40, 24.90, 22.90, 21.20]
        ! double precision, parameter :: w(10) =  [0.008, 0.098, 0.153, 0.199, 0.251, 0.299, 0.349, 0.398, 0.445, 0.489]

        double precision, parameter :: k(10, 10) = reshape([ &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], shape=[10, 10])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))
        x = 0
        x(1: cantidad_variables -1) = [y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), y(9)]
        x(cantidad_variables) = 1 - sum(x)

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_06

        ! Tolueno + agua + anilina, NRTL
    function sistema_07(w) result(respuesta)
        ! Tolueno + agua + anilina, NRTL
        ! x* = 0.000066937752536364210, 0.99686525905361334      
        ! f* = -0.29453555322324276

        ! McDonald C. M. and Floudas C. A. (1994). 
        ! Global Optimization for the phase stability problem. 
        ! AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 3
        double precision, parameter :: z(3) = [0.29989, 0.20006, 0.50005]
        double precision, parameter :: tau(3, 3) = reshape([ &
        0.0, 4.93035, 1.59806, &
        7.77063, 0.0, 4.18462, &
        0.03509, 1.27932, 0.0], shape=[3, 3])
        double precision, parameter :: G(3, 3) = reshape([ &
        1.0, 0.29370, 0.61914, &
        0.14500, 1.0, 0.23984, &
        0.98953, 0.64629, 1.0], shape=[3, 3])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), 1 - w(1) - w(2)]

        lngx = ln_gamma_nrtl(x, tau, G)
        lngz = ln_gamma_nrtl(z, tau, G)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_07

        ! Etilenglicol + dodecanol + nitrometano, UNIQUAC
    function sistema_08(w) result(respuesta)
        ! etilenglicol (1), lauril alcohol (2) y nitrometano (3). 
        ! Etilenglicol + dodecanol + nitrometano, UNIQUAC
        ! x* = 0.75425313547529549, 0.0022192732704503513
        ! f* = -0.11395145263903501

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 3
        double precision, parameter :: z(3) = [0.4, 0.3, 0.3]
        double precision, parameter :: r(3) = [2.4088, 8.8495, 2.0086]
        double precision, parameter :: q(3) = [2.248, 7.372, 1.868]
        double precision, parameter :: qq(3) = [2.248, 7.372, 1.868]
        double precision, parameter :: tau(3, 3) = reshape([ &
        1.0, 0.432589, 0.830749, &
        0.789593, 1.0, 0.354992, &
        0.204736, 0.636678, 1.0], shape=[3, 3])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), 1 - w(1) - w(2)]

        lngx = ln_gamma_uniquac(x, r, q, qq, tau)
        lngz = ln_gamma_uniquac(z, r, q, qq, tau)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_08

        ! n-Propanol + n-butanol + agua, NRTL
    function sistema_10(w) result(respuesta)
        ! n-Propanol + n-butanol + agua, NRTL
        ! x* = 0.009.4074726576420242, 0.019052055523897715  
        ! f* = -0.011609328088991122

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 3
        double precision, parameter :: z(3) = [0.04, 0.16, 0.80]
        double precision, parameter :: tau(3, 3) = reshape([ &
        0.0, -0.61259, -0.07149, &
        0.71640, 0.0, 0.90047, &
        2.7425, 3.51307, 0.0], shape=[3, 3])
        double precision, parameter :: G(3, 3) = reshape([ &
        1.0, 1.2017478, 1.021678, &
        0.8066060, 1.0, 0.6490629, &
        0.4392221, 0.1852084, 1.0], shape=[3, 3])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), 1 - w(1) - w(2)]

        lngx = ln_gamma_nrtl(x, tau, G)
        lngz = ln_gamma_nrtl(z, tau, G)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_10

        ! C1 + C3, SRK
    function sistema_12(y) result(respuesta)
        ! C1 + C3, SRK
        ! x* = 0.19276444431534911      
        ! f* = -0.22262682176770848

        ! Hua, J. Z., Brennecke, J. F., & Stadtherr, M. A. (1998). 
        ! Reliable computation of phase stability using interval analysis. 
        ! Computers & Chemical Engineering, 22(9), 1207–1214. doi:10.1016/s0098-1354(98)00024-6
        
        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 2
        double precision :: P = 50.0, T = 277.6
        double precision, parameter :: z(2) = [0.6, 0.4]
        double precision, parameter :: Tc(2) = [190.6, 369.8]
        double precision, parameter :: Pc(2) = [46.0, 42.5]
        double precision, parameter :: w(2) = [0.008, 0.152]
        double precision, parameter :: k(2, 2) = reshape([ &
        0.0, 0.029, &
        0.029, 0.0], shape=[2, 2])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))

        x = [y(1), 1 - y(1)]

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_12

        ! n-Propanol + n-butanol + benceno + agua, NRTL
    function sistema_13(w) result(respuesta)
        ! n-Propanol + n-butanol + benceno + agua, NRTL
        ! x* = 0.018114780903023944, 0.00061999525113539625, 0.0044844004034368520, 0.42237505492523608      
        ! f* = -0.33982213819701940

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 4
        double precision, parameter :: z(4) = [0.148, 0.052, 0.60, 0.20]
        double precision, parameter :: tau(4, 4) = reshape([ &
        0.0, 2.16486, 0.23689, 0.13060, &
        -1.2007, 0.0, -0.09730, 0.19154, &
        2.01911, 1.73912, 0.0, 4.01932, &
        2.31985, 4.31706, 4.09334, 0.0], shape=[4, 4])
        double precision, parameter :: G(4, 4) = reshape([ &
        1.0, 0.34320, 0.93449, 0.96384, &
        1.80967, 1.0, 1.02932, 0.93623, &
        0.56132, 0.59659, 1.0, 0.32322, &
        0.51986, 0.22649, 0.31656, 1.0], shape=[4, 4])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), w(3), 1 - w(1) - w(2) - w(3)]

        lngx = ln_gamma_nrtl(x, tau, G)
        lngz = ln_gamma_nrtl(z, tau, G)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_13

        ! n-Propanol + n-butanol + benceno + etanol + agua, NRTL
    function sistema_14(w) result(respuesta)
        ! n-Propanol + n-butanol + benceno + etanol + agua, NRTL
        ! x* = 0.024314644551625304, 0.00054503387520804633, 0.0017265770398075587, 0.035514954118328224, 0.57410251482629460      
        ! f* = -0.10430190637039251

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 5
        double precision, parameter :: z(5) = [0.148, 0.052, 0.50, 0.10, 0.20]
        double precision, parameter :: tau(5, 5) = reshape([ &
        0.0, 2.16486, 0.23686, 3.78001, 0.13060, &
        -1.20070, 0.0, -0.09730, -1.15187, -0.20374, &
        2.01911, 1.73912, 0.0, 1.85228, 3.73758, &
        -0.10979, 1.16315, 0.47676, 0.0, -0.14651, &
        2.31985, 5.22337, 6.45226, 2.17820, 0.0], shape=[5, 5])
        double precision, parameter :: G(5, 5) = reshape([ &
        1.0, 0.34320, 0.93450, 0.35902, 0.96384, &
        3.91873, 1.0, 1.02931, 1.41931, 1.06238, &
        0.56132, 0.59670, 1.0, 0.57907, 0.36864, &
        1.03030, 0.70216, 0.86880, 1.0, 1.04035, &
        0.51986, 0.21196, 0.17857, 0.55537, 1.0], shape=[5, 5])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), w(3), w(4), 1 - w(1) - w(2) - w(3) - w(4)]

        lngx = ln_gamma_nrtl(x, tau, G)
        lngz = ln_gamma_nrtl(z, tau, G)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_14

        ! Ácido acético + benceno + furfural + ciclohexano + agua, UNIQUAC
    function sistema_17(w) result(respuesta)
        ! Ácido acético + benceno + furfural + ciclohexano + agua, UNIQUAC
        ! x* =  0.22748543782080269, 0.0033447179085991561, 0.047502524800731014, 0.0015842688961442879 
        ! f* = -0.17765416298801889

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 5
        double precision, parameter :: z(5) = [0.2, 0.2, 0.2, 0.2, 0.2]
        double precision, parameter :: r(5) = [2.2024, 3.1878, 3.1680, 4.0464, 0.9200]
        double precision, parameter :: q(5) = [2.072, 2.400, 2.484, 3.240, 1.400]
        double precision, parameter :: qq(5) = [2.072, 2.400, 2.484, 3.240, 1.400]
        double precision, parameter :: tau(5, 5) = reshape([ &
        1.0, 1.26362, 3.36860, 0.85128, 1.54662, &
        0.99972, 1.0, 1.02041, 0.89333, 0.09441, &
        0.31633, 0.79027, 1.0, 0.96249, 0.60488, &
        0.49739, 1.09619, 0.26222, 1.0, 0.08839, &
        2.44225, 0.13507, 0.69066, 0.19491, 1.0], shape=[5, 5])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), w(3), w(4), 1 - w(1) - w(2) - w(3) - w(4)]

        lngx = ln_gamma_uniquac(x, r, q, qq, tau)
        lngz = ln_gamma_uniquac(z, r, q, qq, tau)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_17

        ! Agua + CO2 + 2-propanol + etanol, SRK
    function sistema_19(y) result(respuesta)
        ! Agua + CO2 + 2-propanol + etanol, SRK
        ! x* = 0.18504982605063439, 0.0023796701744409404, 0.45441226048389355       
        ! f* = -0.012649603564122792

        ! Harding, S. T., & Floudas, C. A. (2000). 
        ! Phase stability with cubic equations of state: Global optimization approach. 
        ! AIChE Journal, 46(7), 1422–1440. doi:10.1002/aic.690460715 

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: ln_phi_x, ln_phi_z 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 4
        double precision :: P = 22.5, T = 350
        double precision, parameter :: z(4) = [0.99758, 0.00003, 0.00013, 0.00226]
        ! https://www.kaylaiacovino.com/Petrology_Tools/Critical_Constants_and_Acentric_Factors.htm

        double precision, parameter :: Tc(4) = [647.3, 304.1, 508.3, 513.9]
        double precision, parameter :: Pc(4) = [221.2, 73.8, 47.6, 61.4]
        double precision, parameter :: w(4) = [0.344, 0.239, 0.665, 0.644]
        double precision, parameter :: k(4, 4) = reshape([ &
        0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0], shape=[4, 4])

        allocate(x(cantidad_variables))
        allocate(ln_phi_x(cantidad_variables))
        allocate(ln_phi_z(cantidad_variables))

        x = [y(1), y(2), y(3), 1 - y(1) - y(2) - y(3)]

        ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T)
        ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + ln_phi_x - log(z) - ln_phi_z))
        end if

    end function sistema_19

        ! Ácido acético + benceno + furfural + ciclohexano, UNIQUAC
    function sistema_20(w) result(respuesta)
        ! Ácido acético + benceno + furfural + ciclohexano, UNIQUAC
        ! x* =  0.017514477936952964, 0.19953729228166880, 0.13392662040944914       
        ! f* = -0.0049317723847621814

        ! Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000). 
        ! Reliable phase stability analysis for excess Gibbs energy models. 
        ! Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x

        double precision :: respuesta
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), allocatable :: x
        double precision, dimension(:), allocatable:: lngx, lngz 

        double precision :: criterio
        integer, parameter :: cantidad_variables = 4
        double precision, parameter :: z(4) = [0.05, 0.20, 0.35, 0.40]
        double precision, parameter :: r(4) = [2.2024, 3.1878, 3.1680, 4.0464]
        double precision, parameter :: q(4) = [2.072, 2.400, 2.484, 3.240]
        double precision, parameter :: qq(4) = [2.072, 2.400, 2.484, 3.240]
        double precision, parameter :: tau(4, 4) = reshape([ &
        1.0, 1.26362, 3.36860, 0.85128, & 
        0.99972, 1.0, 1.02041, 0.89333, &
        0.31633, 0.79027, 1.0, 0.96249, &
        0.49739, 1.09619, 0.26222, 1.0], shape=[4, 4])

        allocate(x(cantidad_variables))
        allocate(lngx(cantidad_variables))
        allocate(lngz(cantidad_variables))

        x = [w(1), w(2), w(3), 1 - w(1) - w(2) - w(3)]

        lngx = ln_gamma_uniquac(x, r, q, qq, tau)
        lngz = ln_gamma_uniquac(z, r, q, qq, tau)

        criterio = sum(pack(x, x < 0))

        if (criterio < 0) then 
            respuesta = - criterio * 10.0E5
        else 
            respuesta = sum(x * (log(x) + lngx - log(z) - lngz))
        end if

    end function sistema_20

end module problemas_estabilidad_fases