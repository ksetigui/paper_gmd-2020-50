MODULE radm2_JacobianSP
  PUBLIC
  SAVE
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
       3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, &
       8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 13, 13, &
      13, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, &
      17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, &
      20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, &
      23, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 25, &
      25, 25, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, &
      26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, &
      27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29, 29, &
      30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, &
      31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 32, &
      33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, &
      33, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, &
      35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, &
      36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, &
      38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, &
      39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, &
      39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, &
      41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 44, 44, &
      44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48 /)
  INTEGER, PARAMETER, DIMENSION(299) :: LU_IROW_1 = (/ &
      48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59 /)
  INTEGER, PARAMETER, DIMENSION(659) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 5, 53, 2, 24, 27, 28, 29, 46, 3, 27, 28, &
      29, 37, 38, 40, 41, 44, 45, 46, 49, 50, 52, 55, &
      56, 57, 4, 26, 53, 5, 53, 6, 53, 7, 46, 8, &
      53, 9, 53, 10, 53, 11, 53, 12, 45, 59, 13, 53, &
      58, 14, 51, 53, 15, 54, 59, 16, 53, 17, 28, 29, &
      46, 53, 18, 51, 53, 55, 7, 19, 46, 54, 59, 20, &
      51, 53, 59, 21, 51, 52, 53, 10, 11, 22, 53, 54, &
      23, 53, 55, 59, 24, 46, 53, 54, 15, 22, 25, 30, &
      31, 34, 43, 48, 51, 53, 54, 59, 24, 26, 27, 28, &
      29, 31, 34, 43, 45, 46, 48, 52, 53, 54, 55, 58, &
      27, 46, 53, 54, 28, 46, 53, 54, 29, 46, 53, 54, &
      30, 35, 36, 52, 53, 54, 55, 58, 31, 35, 45, 52, &
      53, 54, 55, 58, 22, 32, 51, 52, 53, 54, 55, 59, &
      16, 29, 33, 41, 44, 46, 47, 49, 52, 53, 54, 55, &
      56, 58, 34, 35, 36, 45, 50, 52, 53, 54, 55, 58, &
      10, 35, 51, 52, 53, 55, 58, 11, 36, 51, 52, 53, &
      55, 58, 27, 28, 37, 46, 51, 52, 53, 54, 55, 58, &
      24, 27, 28, 29, 38, 46, 51, 52, 53, 54, 55, 58, &
       8, 9, 22, 23, 39, 45, 51, 52, 53, 54, 55, 58, &
      59, 24, 40, 46, 51, 52, 53, 54, 55, 58, 8, 41, &
      51, 52, 53, 55, 58, 32, 35, 36, 37, 39, 40, 41, &
      42, 44, 45, 46, 49, 50, 51, 52, 53, 54, 55, 56, &
      57, 58, 59, 16, 21, 23, 24, 27, 28, 29, 31, 32, &
      35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 49, &
      50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 9, 44, &
      51, 52, 53, 55, 58, 12, 22, 30, 35, 36, 45, 51, &
      52, 53, 54, 55, 58, 59, 19, 24, 27, 28, 29, 46, &
      51, 53, 54, 58, 59, 32, 38, 41, 44, 46, 47, 51, &
      52, 53, 54, 55, 56, 58, 59, 16, 27, 28, 29, 37, &
      38, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51 /)
  INTEGER, PARAMETER, DIMENSION(299) :: LU_ICOL_1 = (/ &
      52, 53, 54, 55, 56, 57, 58, 59, 29, 46, 49, 51, &
      52, 53, 54, 55, 58, 59, 33, 41, 44, 46, 47, 49, &
      50, 51, 52, 53, 54, 55, 56, 58, 59, 5, 10, 11, &
      14, 16, 20, 21, 22, 24, 26, 27, 28, 29, 30, 31, &
      32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, &
      45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, &
      57, 58, 59, 17, 18, 21, 27, 28, 29, 32, 35, 36, &
      37, 38, 39, 40, 41, 44, 45, 46, 48, 49, 50, 51, &
      52, 53, 54, 55, 56, 57, 58, 59, 5, 6, 7, 8, &
       9, 10, 11, 13, 14, 16, 17, 18, 20, 21, 22, 23, &
      24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, &
      41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, &
      53, 54, 55, 56, 57, 58, 59, 15, 20, 22, 23, 24, &
      25, 27, 28, 29, 30, 31, 34, 35, 36, 43, 44, 45, &
      46, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, &
      59, 18, 23, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
      41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, &
      55, 56, 57, 58, 59, 16, 42, 44, 45, 46, 47, 49, &
      50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 6, 33, &
      41, 44, 46, 47, 49, 51, 52, 53, 54, 55, 56, 57, &
      58, 59, 13, 19, 35, 36, 37, 38, 39, 40, 41, 44, &
      45, 46, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, &
      59, 12, 15, 19, 20, 23, 25, 30, 31, 32, 34, 35, &
      36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, &
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59 /)
  INTEGER, PARAMETER, DIMENSION(659) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1 /)
  INTEGER, PARAMETER, DIMENSION(60) :: LU_CROW = (/ &
       1, 4, 10, 27, 30, 32, 34, 36, 38, 40, 42, 44, &
      47, 50, 53, 56, 58, 63, 67, 72, 76, 80, 85, 89, &
      93,105,121,125,129,133,141,149,157,171,181,188, &
     195,205,217,230,239,246,268,299,306,319,330,344, &
     369,379,394,436,465,512,542,570,587,603,626,660 /)
  INTEGER, PARAMETER, DIMENSION(60) :: LU_DIAG = (/ &
       1, 4, 10, 27, 30, 32, 34, 36, 38, 40, 42, 44, &
      47, 50, 53, 56, 58, 63, 68, 72, 76, 82, 85, 89, &
      95,106,121,125,129,133,141,150,159,171,182,189, &
     197,209,221,231,240,253,284,300,311,324,335,357, &
     371,385,427,457,505,536,565,583,600,624,659,660 /)
END MODULE radm2_JacobianSP
