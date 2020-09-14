MODULE saprc99_JacobianSP
  PUBLIC
  SAVE
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, &
       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, &
       4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, &
       5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, &
       6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, &
       9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 11, 11, &
      12, 12, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, &
      16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, &
      20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, &
      25, 25, 25, 25, 26, 26, 27, 27, 28, 28, 28, 28, &
      28, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, &
      31, 31, 32, 32, 32, 32, 33, 33, 34, 34, 35, 35, &
      35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 37, 38, &
      38, 38, 38, 38, 39, 39, 39, 40, 40, 40, 40, 40, &
      40, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, &
      42, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, &
      44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 50, 50, 50, 50, 50, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, &
      52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, &
      55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, &
      57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      59, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72 /)
  INTEGER, PARAMETER, DIMENSION(248) :: LU_IROW_2 = (/ &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79 /)
  INTEGER, PARAMETER, DIMENSION(968) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 10, 71, 2, 39, 50, 67, 3, 30, 39, 48, 50, &
      52, 54, 55, 56, 57, 58, 59, 67, 71, 73, 4, 50, &
      56, 58, 67, 68, 69, 76, 77, 78, 5, 52, 54, 56, &
      57, 58, 67, 69, 70, 72, 76, 77, 78, 79, 6, 68, &
      78, 7, 70, 72, 78, 79, 8, 31, 50, 74, 75, 9, &
      31, 41, 50, 63, 67, 71, 74, 75, 10, 71, 11, 67, &
      12, 71, 13, 71, 14, 27, 54, 67, 71, 15, 68, 74, &
      16, 72, 74, 17, 74, 79, 18, 70, 74, 19, 71, 78, &
      20, 71, 21, 71, 22, 74, 75, 23, 71, 73, 24, 71, &
      24, 25, 71, 74, 26, 71, 27, 71, 28, 69, 71, 76, &
      77, 29, 71, 76, 78, 30, 61, 73, 78, 31, 40, 74, &
      75, 78, 32, 71, 74, 78, 33, 71, 34, 71, 27, 34, &
      35, 71, 27, 34, 36, 71, 27, 34, 37, 71, 75, 27, &
      34, 38, 67, 71, 39, 67, 71, 31, 40, 51, 74, 75, &
      78, 27, 34, 41, 58, 67, 71, 75, 42, 69, 71, 77, &
      78, 34, 43, 51, 71, 75, 78, 27, 34, 35, 36, 37, &
      44, 55, 57, 59, 67, 71, 75, 33, 35, 36, 38, 39, &
      44, 45, 48, 49, 50, 52, 54, 55, 56, 57, 58, 59, &
      60, 61, 63, 64, 67, 71, 75, 20, 24, 25, 26, 33, &
      46, 54, 56, 58, 62, 67, 71, 74, 75, 22, 37, 40, &
      41, 43, 44, 47, 49, 51, 55, 57, 58, 59, 60, 61, &
      64, 67, 71, 74, 75, 78, 48, 63, 67, 71, 75, 27, &
      34, 35, 36, 38, 39, 43, 48, 49, 51, 54, 57, 63, &
      67, 71, 75, 78, 50, 63, 67, 71, 75, 37, 43, 51, &
      68, 70, 71, 72, 73, 74, 75, 78, 79, 52, 63, 67, &
      71, 75, 24, 26, 33, 35, 36, 46, 52, 53, 54, 56, &
      58, 59, 62, 63, 65, 66, 67, 68, 69, 70, 71, 72, &
      73, 74, 75, 76, 77, 78, 79, 54, 63, 67, 71, 75, &
      52, 55, 58, 63, 67, 71, 75, 56, 63, 67, 71, 75, &
      52, 57, 58, 63, 67, 71, 75, 58, 63, 67, 71, 75 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      52, 58, 59, 63, 67, 71, 75, 13, 21, 24, 26, 33, &
      48, 50, 56, 57, 58, 60, 62, 63, 64, 65, 66, 67, &
      68, 70, 71, 72, 73, 75, 79, 21, 24, 26, 28, 29, &
      30, 33, 39, 46, 48, 49, 50, 51, 52, 54, 55, 56, &
      57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, &
      71, 72, 73, 74, 75, 76, 77, 78, 79, 25, 54, 56, &
      57, 58, 62, 63, 67, 69, 71, 73, 74, 75, 11, 48, &
      50, 52, 54, 55, 56, 58, 59, 63, 67, 71, 73, 74, &
      75, 20, 24, 26, 33, 35, 36, 38, 42, 48, 50, 54, &
      55, 56, 57, 58, 59, 62, 63, 64, 65, 66, 67, 69, &
      71, 73, 74, 75, 77, 78, 24, 26, 33, 50, 55, 56, &
      57, 58, 59, 62, 63, 65, 66, 67, 69, 71, 73, 74, &
      75, 76, 77, 26, 33, 34, 52, 54, 56, 57, 58, 59, &
      62, 63, 66, 67, 68, 69, 71, 72, 73, 74, 75, 76, &
      77, 79, 38, 39, 48, 50, 52, 54, 55, 56, 57, 58, &
      59, 63, 67, 68, 70, 71, 72, 73, 74, 75, 78, 79, &
      14, 15, 27, 33, 35, 36, 44, 46, 54, 55, 56, 57, &
      58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
      71, 72, 73, 74, 75, 76, 77, 78, 79, 20, 24, 26, &
      27, 33, 34, 50, 52, 54, 56, 57, 58, 59, 62, 63, &
      64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, &
      76, 77, 78, 79, 18, 52, 55, 57, 58, 59, 63, 67, &
      68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      10, 11, 12, 13, 19, 20, 21, 23, 24, 26, 27, 28, &
      29, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, &
      44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, &
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 16, &
      49, 51, 54, 55, 57, 58, 59, 63, 64, 65, 66, 67, &
      68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79 /)
  INTEGER, PARAMETER, DIMENSION(248) :: LU_ICOL_2 = (/ &
      23, 30, 53, 54, 56, 58, 59, 61, 62, 63, 65, 66, &
      67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, &
      79, 15, 16, 17, 18, 22, 23, 25, 30, 31, 32, 40, &
      47, 49, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, &
      73, 74, 75, 76, 77, 78, 79, 22, 32, 37, 40, 41, &
      43, 44, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, &
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 12, &
      25, 29, 33, 46, 48, 50, 52, 54, 56, 58, 59, 60, &
      62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, &
      74, 75, 76, 77, 78, 79, 13, 20, 21, 24, 26, 27, &
      33, 34, 35, 36, 37, 38, 39, 42, 43, 48, 50, 51, &
      52, 54, 55, 56, 57, 58, 59, 62, 63, 64, 65, 66, &
      67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, &
      79, 10, 19, 21, 23, 27, 28, 29, 30, 31, 32, 34, &
      35, 36, 38, 39, 40, 42, 44, 45, 48, 49, 50, 51, &
      52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, &
      77, 78, 79, 17, 41, 58, 63, 67, 68, 69, 70, 71, &
      72, 73, 74, 75, 76, 77, 78, 79 /)
  INTEGER, PARAMETER, DIMENSION(968) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)
  INTEGER, PARAMETER, DIMENSION(80) :: LU_CROW = (/ &
       1, 4, 8, 23, 33, 47, 50, 55, 60, 69, 71, 73, &
      75, 77, 82, 85, 88, 91, 94, 97, 99,101,104,107, &
     109,113,115,117,122,126,130,135,139,141,143,147, &
     151,156,161,164,170,177,182,188,200,224,238,259, &
     264,281,286,298,303,332,337,344,349,356,361,368, &
     392,430,443,458,487,508,531,553,586,617,637,696, &
     721,746,788,828,859,902,952,969 /)
  INTEGER, PARAMETER, DIMENSION(80) :: LU_DIAG = (/ &
       1, 4, 8, 23, 33, 47, 50, 55, 60, 69, 71, 73, &
      75, 77, 82, 85, 88, 91, 94, 97, 99,101,104,107, &
     110,113,115,117,122,126,130,135,139,141,145,149, &
     153,158,161,165,172,177,183,193,206,229,244,259, &
     272,281,288,298,310,332,338,344,350,356,363,378, &
     412,435,452,476,498,519,543,574,606,627,687,713, &
     739,782,823,855,899,950,968,969 /)
END MODULE saprc99_JacobianSP
