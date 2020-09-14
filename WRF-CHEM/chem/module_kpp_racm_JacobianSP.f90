MODULE racm_JacobianSP
  PUBLIC
  SAVE
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, &
       3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, &
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
       4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, &
       9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, &
      13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 16, 16, &
      16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, &
      18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, &
      20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, &
      22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24, &
      24, 24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, &
      26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, &
      27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, &
      28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, &
      29, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, &
      32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 33, 33, &
      33, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, &
      34, 34, 34, 35, 35, 35, 35, 35, 35, 36, 36, 36, &
      36, 36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 38, &
      38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, &
      39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, &
      40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, &
      43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, &
      44, 44, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, &
      48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66 /)
  INTEGER, PARAMETER, DIMENSION(332) :: LU_IROW_2 = (/ &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73 /)
  INTEGER, PARAMETER, DIMENSION(1052) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 5, 68, 2, 28, 68, 3, 12, 24, 30, 33, 35, &
      37, 51, 54, 57, 63, 68, 4, 30, 36, 38, 39, 44, &
      45, 50, 51, 52, 54, 55, 57, 58, 59, 60, 61, 63, &
      64, 65, 66, 67, 5, 68, 6, 63, 7, 68, 8, 68, &
       9, 68, 10, 70, 73, 11, 68, 12, 68, 13, 57, 58, &
      63, 68, 14, 57, 58, 63, 68, 15, 54, 68, 16, 67, &
      68, 73, 17, 65, 67, 68, 18, 21, 22, 25, 68, 72, &
      73, 19, 29, 30, 33, 37, 57, 58, 63, 67, 68, 20, &
      32, 67, 68, 70, 73, 8, 21, 63, 68, 73, 9, 22, &
      63, 68, 73, 11, 23, 35, 51, 64, 65, 68, 24, 63, &
      68, 70, 25, 32, 63, 68, 73, 26, 54, 63, 66, 67, &
      68, 27, 32, 34, 40, 48, 51, 53, 54, 67, 68, 70, &
      73, 12, 24, 28, 29, 30, 33, 34, 35, 37, 40, 48, &
      51, 53, 54, 57, 58, 62, 63, 68, 70, 29, 63, 68, &
      70, 30, 63, 68, 70, 31, 35, 63, 66, 68, 70, 73, &
      20, 21, 22, 25, 32, 63, 67, 68, 70, 73, 33, 63, &
      68, 70, 12, 34, 41, 46, 47, 52, 54, 63, 65, 66, &
      68, 70, 72, 35, 55, 63, 68, 70, 73, 24, 36, 63, &
      65, 66, 67, 68, 70, 72, 37, 62, 63, 68, 70, 38, &
      57, 65, 66, 67, 68, 70, 72, 39, 58, 65, 66, 67, &
      68, 70, 72, 23, 35, 40, 41, 46, 47, 51, 52, 54, &
      55, 63, 64, 65, 66, 68, 70, 72, 73, 25, 32, 41, &
      63, 65, 66, 67, 68, 70, 72, 73, 7, 15, 29, 38, &
      39, 42, 44, 45, 49, 52, 54, 57, 58, 59, 60, 63, &
      64, 65, 66, 67, 68, 69, 70, 71, 72, 30, 43, 63, &
      65, 66, 67, 68, 70, 72, 7, 44, 65, 66, 67, 68, &
      70, 72, 11, 45, 65, 66, 67, 68, 70, 72, 21, 46, &
      63, 65, 66, 67, 68, 70, 72, 73, 22, 47, 63, 65, &
      66, 67, 68, 70, 72, 73, 12, 17, 23, 24, 26, 30, &
      31, 33, 34, 35, 36, 37, 38, 39, 41, 43, 44, 45 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, &
      59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 70, 72, &
      73, 29, 49, 63, 65, 66, 67, 68, 70, 72, 33, 37, &
      50, 62, 63, 65, 66, 67, 68, 70, 72, 30, 33, 37, &
      43, 50, 51, 62, 63, 65, 66, 67, 68, 70, 72, 12, &
      52, 65, 66, 67, 68, 69, 70, 71, 72, 11, 12, 15, &
      29, 36, 38, 39, 44, 45, 49, 51, 52, 53, 54, 57, &
      58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, &
      70, 71, 72, 37, 46, 47, 54, 62, 63, 65, 66, 67, &
      68, 70, 72, 73, 35, 51, 54, 55, 62, 63, 65, 66, &
      67, 68, 70, 72, 73, 8, 9, 26, 31, 32, 33, 35, &
      37, 44, 45, 51, 52, 54, 55, 56, 62, 63, 64, 65, &
      66, 67, 68, 69, 70, 71, 72, 73, 30, 33, 37, 50, &
      57, 62, 63, 65, 66, 67, 68, 70, 72, 43, 50, 58, &
      62, 63, 65, 66, 67, 68, 70, 72, 24, 29, 30, 33, &
      37, 51, 57, 58, 59, 60, 62, 63, 65, 66, 67, 68, &
      70, 72, 24, 29, 30, 33, 37, 57, 58, 59, 60, 62, &
      63, 65, 66, 67, 68, 70, 72, 13, 29, 30, 42, 44, &
      45, 49, 52, 54, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 66, 67, 68, 69, 70, 71, 72, 73, 6, 33, 37, &
      51, 62, 63, 65, 66, 67, 68, 70, 72, 73, 6, 21, &
      22, 24, 25, 29, 30, 32, 33, 35, 37, 51, 54, 55, &
      57, 58, 62, 63, 65, 66, 67, 68, 70, 72, 73, 29, &
      30, 33, 37, 42, 44, 45, 49, 52, 54, 57, 58, 59, &
      60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, &
      73, 14, 17, 26, 33, 36, 37, 38, 39, 41, 43, 44, &
      45, 46, 47, 49, 50, 52, 53, 54, 55, 56, 57, 58, &
      59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
      71, 72, 73, 23, 26, 31, 33, 35, 36, 37, 38, 39, &
      40, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52 /)
  INTEGER, PARAMETER, DIMENSION(332) :: LU_ICOL_2 = (/ &
      53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
      65, 66, 67, 68, 69, 70, 71, 72, 73, 5, 7, 8, &
       9, 11, 12, 15, 16, 17, 19, 20, 21, 22, 23, 24, &
      25, 26, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, &
      39, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, &
      52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, &
      64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 5, 6, &
       7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, &
      33, 34, 35, 37, 40, 41, 42, 44, 45, 46, 47, 48, &
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, &
      73, 20, 32, 35, 43, 44, 45, 46, 47, 49, 50, 52, &
      55, 59, 60, 62, 63, 65, 66, 67, 68, 69, 70, 71, &
      72, 73, 10, 16, 24, 27, 29, 30, 31, 32, 33, 34, &
      35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, &
      48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, &
      60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, &
      72, 73, 36, 38, 39, 41, 43, 44, 45, 46, 47, 49, &
      50, 51, 52, 55, 56, 57, 58, 61, 62, 63, 64, 65, &
      66, 67, 68, 69, 70, 71, 72, 73, 18, 21, 22, 25, &
      32, 36, 38, 39, 41, 43, 44, 45, 46, 47, 49, 50, &
      52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, &
      66, 67, 68, 69, 70, 71, 72, 73, 10, 16, 18, 20, &
      21, 22, 25, 27, 31, 32, 34, 35, 36, 38, 39, 40, &
      41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, &
      54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, &
      66, 67, 68, 69, 70, 71, 72, 73 /)
  INTEGER, PARAMETER, DIMENSION(1052) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)
  INTEGER, PARAMETER, DIMENSION(74) :: LU_CROW = (/ &
       1, 4, 7, 19, 41, 43, 45, 47, 49, 51, 54, 56, &
      58, 63, 68, 71, 75, 79, 86, 96,102,107,112,119, &
     123,128,134,146,166,170,174,181,191,195,208,214, &
     223,228,236,244,262,273,298,307,315,323,333,343, &
     386,395,406,420,430,460,473,486,513,526,537,555, &
     572,598,611,636,662,700,742,803,866,891,939,969, &
     1005,1053 /)
  INTEGER, PARAMETER, DIMENSION(74) :: LU_DIAG = (/ &
       1, 4, 7, 19, 41, 43, 45, 47, 49, 51, 54, 56, &
      58, 63, 68, 71, 75, 79, 86, 96,103,108,113,119, &
     123,128,134,148,166,170,174,185,191,196,208,215, &
     223,228,236,246,264,278,299,308,316,324,334,363, &
     387,397,411,421,442,463,476,500,517,528,545,563, &
     585,602,628,652,691,734,796,860,886,935,966,1003, &
     1052,1053 /)
END MODULE racm_JacobianSP
