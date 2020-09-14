MODULE racm_esrlsorg_JacobianSP
  PUBLIC
  SAVE
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, &
       3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, &
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
       4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, &
       8, 9, 9, 10, 10, 11, 11, 11, 12, 12, 13, 13, &
      13, 13, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, &
      17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 19, 19, &
      20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, &
      22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24, &
      24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, &
      25, 25, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, &
      27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, &
      29, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, &
      32, 32, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, &
      33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 35, &
      35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, &
      37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, &
      37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, &
      39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, &
      40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 44, &
      44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, &
      46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, &
      47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, &
      48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, &
      55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80 /)
  INTEGER, PARAMETER, DIMENSION(6) :: LU_IROW_3 = (/ &
      80, 80, 80, 80, 80, 80 /)
  INTEGER, PARAMETER, DIMENSION(1086) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 5, 73, 2, 37, 73, 3, 14, 30, 38, 39, 41, &
      43, 48, 61, 63, 73, 77, 4, 38, 45, 46, 47, 50, &
      53, 55, 61, 63, 65, 66, 67, 68, 69, 71, 72, 75, &
      76, 77, 78, 5, 73, 6, 77, 7, 50, 73, 75, 8, &
      73, 9, 73, 10, 73, 11, 74, 80, 12, 73, 13, 42, &
      73, 75, 14, 73, 15, 63, 65, 73, 77, 16, 61, 73, &
      17, 63, 65, 73, 77, 18, 34, 73, 19, 73, 75, 80, &
      20, 73, 75, 76, 21, 42, 73, 80, 21, 22, 34, 42, &
      73, 79, 80, 23, 28, 29, 31, 73, 79, 80, 24, 36, &
      38, 39, 41, 63, 65, 73, 75, 77, 12, 25, 43, 69, &
      73, 76, 6, 26, 41, 74, 77, 79, 80, 27, 40, 73, &
      74, 75, 80, 9, 28, 73, 77, 80, 10, 29, 73, 77, &
      80, 30, 73, 74, 77, 31, 40, 73, 77, 80, 32, 40, &
      44, 51, 56, 61, 64, 73, 74, 75, 80, 33, 61, 72, &
      73, 75, 77, 34, 39, 50, 73, 74, 79, 35, 43, 72, &
      73, 74, 77, 80, 36, 73, 74, 77, 14, 18, 30, 34, &
      36, 37, 38, 39, 41, 42, 43, 44, 48, 50, 51, 56, &
      61, 63, 64, 65, 73, 74, 77, 79, 38, 73, 74, 77, &
      39, 73, 74, 77, 27, 28, 29, 31, 40, 73, 74, 75, &
      77, 80, 41, 73, 74, 77, 13, 21, 39, 42, 48, 73, &
      74, 75, 77, 79, 80, 43, 73, 74, 77, 78, 80, 14, &
      44, 52, 57, 58, 61, 71, 72, 73, 74, 76, 77, 79, &
      30, 45, 72, 73, 74, 75, 76, 77, 79, 46, 63, 72, &
      73, 74, 75, 76, 79, 47, 65, 72, 73, 74, 75, 76, &
      79, 7, 38, 39, 41, 48, 50, 54, 72, 73, 74, 75, &
      76, 77, 79, 8, 16, 36, 46, 47, 49, 53, 55, 59, &
      60, 61, 63, 65, 66, 67, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 79, 39, 41, 50, 72, 73, 74, 75, 76, &
      77, 79, 22, 25, 34, 39, 42, 43, 48, 50, 51, 52, &
      54, 57, 58, 61, 69, 71, 72, 73, 74, 75, 76, 77 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      78, 79, 80, 31, 40, 52, 72, 73, 74, 75, 76, 77, &
      79, 80, 12, 53, 72, 73, 74, 75, 76, 79, 38, 54, &
      72, 73, 74, 75, 76, 77, 79, 8, 55, 72, 73, 74, &
      75, 76, 79, 14, 18, 20, 25, 30, 33, 34, 35, 38, &
      39, 41, 42, 43, 44, 45, 46, 47, 48, 50, 52, 53, &
      54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, &
      68, 69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      28, 57, 72, 73, 74, 75, 76, 77, 79, 80, 29, 58, &
      72, 73, 74, 75, 76, 77, 79, 80, 36, 59, 72, 73, &
      74, 75, 76, 77, 79, 27, 40, 43, 53, 54, 55, 57, &
      58, 59, 60, 66, 67, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 57, 58, 61, 72, 73, 74, 75, 76, 77, &
      79, 80, 9, 10, 33, 35, 40, 41, 43, 53, 55, 61, &
      62, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 38, 41, 50, 63, 72, 73, 74, 75, 76, 77, 79, &
      12, 14, 16, 36, 45, 46, 47, 53, 55, 59, 60, 61, &
      63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 50, 54, 65, 72, 73, 74, &
      75, 76, 77, 79, 30, 36, 38, 41, 63, 65, 66, 67, &
      72, 73, 74, 75, 76, 77, 79, 30, 36, 38, 41, 63, &
      65, 66, 67, 72, 73, 74, 75, 76, 77, 79, 17, 36, &
      38, 49, 53, 55, 59, 60, 61, 63, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      36, 38, 41, 49, 53, 55, 59, 60, 61, 63, 65, 66, &
      67, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 45, 46, 47, 52, 53, 54, 55, 57, 58, 59, 62, &
      63, 65, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, &
      78, 79, 80, 14, 60, 66, 67, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 25, 33, 35, 39, 41, 42, &
      43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
       5, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, &
      19, 20, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, &
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, &
      48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, &
      60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, &
      72, 73, 74, 75, 76, 77, 78, 79, 80, 11, 19, 26, &
      30, 32, 35, 36, 38, 39, 40, 41, 43, 44, 45, 46, &
      47, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, &
      63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 5, 8, 9, 10, 12, 14, &
      16, 19, 20, 22, 24, 25, 27, 28, 29, 30, 31, 33, &
      34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, &
      47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, &
      60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, &
      72, 73, 74, 75, 76, 77, 78, 79, 80, 15, 20, 33, &
      39, 41, 45, 46, 47, 50, 52, 53, 54, 55, 57, 58, &
      59, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, &
      72, 73, 74, 75, 76, 77, 78, 79, 80, 26, 28, 29, &
      30, 31, 36, 38, 39, 40, 41, 43, 48, 50, 54, 61, &
      63, 65, 72, 73, 74, 75, 76, 77, 78, 79, 80, 43, &
      61, 72, 73, 74, 75, 76, 77, 78, 79, 80, 23, 26, &
      28, 29, 31, 40, 41, 42, 45, 46, 47, 48, 50, 52, &
      53, 54, 55, 57, 58, 59, 62, 63, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      11, 18, 19, 21, 23, 26, 27, 28, 29, 31, 32, 34, &
      35, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, &
      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, &
      63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74 /)
  INTEGER, PARAMETER, DIMENSION(6) :: LU_ICOL_3 = (/ &
      75, 76, 77, 78, 79, 80 /)
  INTEGER, PARAMETER, DIMENSION(1086) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)
  INTEGER, PARAMETER, DIMENSION(81) :: LU_CROW = (/ &
       1, 4, 7, 19, 40, 42, 44, 48, 50, 52, 54, 57, &
      59, 63, 65, 70, 73, 78, 81, 85, 89, 93,100,107, &
     117,123,130,136,141,146,150,155,166,172,178,185, &
     189,213,217,221,231,235,246,252,265,274,282,290, &
     304,329,339,364,375,383,392,400,445,455,465,474, &
     496,507,530,541,571,581,596,611,637,662,688,703, &
     745,814,859,922,958,984,995,1033,1087 /)
  INTEGER, PARAMETER, DIMENSION(81) :: LU_DIAG = (/ &
       1, 4, 7, 19, 40, 42, 44, 48, 50, 52, 54, 57, &
      59, 63, 65, 70, 73, 78, 81, 85, 89, 94,100,107, &
     118,124,130,137,142,146,150,155,166,172,178,185, &
     194,213,217,225,231,238,246,253,266,274,282,294, &
     309,331,347,366,376,384,393,423,446,456,466,483, &
     498,517,533,554,573,587,603,624,650,677,693,736, &
     806,852,916,953,980,992,1031,1086,1087 /)
END MODULE racm_esrlsorg_JacobianSP
