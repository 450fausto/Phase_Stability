module deprice_module
  implicit none

contains

  function DEPrice(CostFunction, LimInf, LimSup, NumPop, MaxIter) result(Solution)
    double precision, dimension(:), allocatable :: Solution

    interface
        double precision function CostFunction(x)
            double precision, dimension(:), intent(in) :: x
        end function CostFunction
    end interface


    integer, intent(in) :: NumPop, MaxIter
    double precision, dimension(:), intent(in) :: LimInf, LimSup
    double precision :: F, Cr, Fbest, FunZ, RndCr, rand_num
    integer :: i, j, k, n, iter, r1, r2, r3, NumVar

    double precision, dimension(:, :), allocatable :: X
    double precision, dimension(:), allocatable :: Fit
    double precision, dimension(:), allocatable :: Xbest, y, z

    F=0.5
    Cr=0.2
    Fbest = 1000.0

    NumVar = size(LimSup)

    allocate(X(NumPop, NumVar))
    allocate(Solution(NumVar + 1))
    allocate(Fit(NumPop))
    allocate(Xbest(NumVar))
    allocate(y(NumVar))
    allocate(z(NumVar))


    do i = 1, NumPop
      do j = 1, NumVar
        call random_number(rand_num)
        X(i, j) = LimInf(j) + (LimSup(j) - LimInf(j)) * rand_num
      end do
      Fit(i) = CostFunction(X(i, :))
    end do

    do i = 1, MaxIter
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

            if (FunZ < Fit(j)) then 
                X(j, :) = z 
                Fit(j) = FunZ 
                if (FunZ < Fbest) then 
                    Fbest = FunZ
                    Xbest = z 
                end if 
            end if
        end do

    end do

    Solution(1:NumVar) = Xbest
    Solution(NumVar + 1) = Fbest
    deallocate(X, Fit, Xbest, z)

  end function DEPrice

end module deprice_module
