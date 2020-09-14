! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : saprc99_mosaic_8bin_vbs2_JacobianSP.f90
! Time                 : Mon Feb 27 11:09:05 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/saprc99_mosaic_8bin_vbs2
! Equation file        : saprc99_mosaic_8bin_vbs2.kpp
! Output root filename : saprc99_mosaic_8bin_vbs2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_mosaic_8bin_vbs2_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  2,  3,  4,  5,  6,  7, &
       7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, &
       7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, &
       9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 11, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, &
      13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, &
      15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, &
      18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 21, 21, &
      21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, &
      23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, &
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, &
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, &
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, &
      24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, &
      25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26, &
      26, 26, 27, 27, 28, 28, 28, 28, 29, 29, 30, 30, &
      31, 31, 32, 32, 32, 33, 33, 34, 34, 34, 35, 35, &
      36, 36, 37, 37, 38, 38, 38, 39, 39, 39, 40, 40, &
      40, 41, 41, 41, 42, 42, 42, 43, 43, 44, 44, 45, &
      45, 45, 45, 45, 45, 46, 46, 46, 47, 47, 47, 48, &
      48, 49, 49, 49, 49, 50, 50, 51, 51, 52, 52, 52, &
      52, 53, 53, 53, 53, 54, 54, 54, 54, 54, 55, 55, &
      55, 55, 56, 56, 56, 56, 56, 57, 57, 58, 58, 58, &
      58, 59, 59, 59, 59, 60, 60, 60, 60, 60, 61, 61, &
      62, 62, 62, 62, 62, 63, 63, 63, 64, 64, 64, 64, &
      64, 64, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, &
      66, 66, 67, 67, 67, 67, 67, 67, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 73, 73, 73, 73, 73, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, &
      76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 78, &
      78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, 80, &
      80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 82, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 84, &
      84, 84, 84, 84, 85, 85, 85, 85, 85, 85, 85, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 92, &
      92, 92, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_2 = (/ &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, &
      93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, &
      94, 94, 94, 94, 94, 94, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, &
      95, 95, 95, 95, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, &
      96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, &
      97, 97, 97, 97, 97, 97, 97, 97, 98, 98, 98, 98, &
      98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, &
      98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99, 99, 99, 99, 99, 99, 99,100,100,100, &
     100,100,100,100,100,100,100,100,100,100,100,100, &
     100,100,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,101,101,101,101, &
     101,101,101,101,101,101,101,101,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,102,102,102,102,102,102,102,102, &
     102,102,102,102,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103 /)
  INTEGER, PARAMETER, DIMENSION(66) :: LU_IROW_3 = (/ &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,103,103, &
     103,103,103,103,103,103,103,103,103,103,104,104, &
     104,104,104,104,104,104,104,104,104,104,104,104, &
     104,104,104,104,104,104,104,104,104,104,104,104, &
     104,104,104,104,104,104 /)
  INTEGER, PARAMETER, DIMENSION(1146) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2, LU_IROW_3 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 36,103,  2, 63, 73, 92,  3,  4,  5,  6,  7, &
      53, 63, 71, 73, 75, 76, 78, 79, 80, 82, 84, 85, &
      92, 93,103,  8, 73, 76, 84, 92, 94, 95, 96,102, &
     104,  9, 75, 76, 79, 80, 82, 84, 92, 94, 96, 97, &
      98,100,102,104, 10, 95, 96, 11, 96, 97, 98,100, &
      12, 54, 73, 99,101, 13, 54, 66, 73, 88, 92, 99, &
     101,103, 14, 77, 93, 94, 99,102,104, 15, 77, 94, &
      96,102,104, 16, 51, 57,103, 17, 51, 57,103, 18, &
      51, 57,103, 19, 51, 57,103, 20, 75,103, 21, 75, &
      79, 80, 92, 99,103, 22, 75, 79, 80, 92, 99,103, &
      23, 75, 79, 80, 92, 99,103, 24, 35, 36, 37, 43, &
      44, 47, 48, 50, 51, 52, 56, 57, 58, 59, 60, 61, &
      62, 63, 65, 66, 67, 68, 69, 70, 71, 72, 73, 75, &
      76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 89, &
      90, 91, 92, 93, 96, 99,101,103,  3,  4,  5,  6, &
      25, 26, 27, 28, 29, 31, 32, 33, 34,103, 26, 27, &
      31,103, 27,103, 28, 29, 33,103, 29,103, 30, 92, &
      31,103, 31, 32,103, 33,103, 33, 34,103, 35,103, &
      36,103, 37,103, 38, 95,101, 39, 97,101, 40,100, &
     101, 41, 98,101, 42, 96,103, 43,103, 44,103, 45, &
      51, 79, 80, 92,103, 46, 99,101, 47, 93,103, 48, &
     103, 48, 49,101,103, 50,103, 51,103, 52, 96,102, &
     103, 53, 86, 93, 96, 54, 64, 96, 99,101, 55, 96, &
     101,103, 56, 94,102,103,104, 57,103, 51, 57, 58, &
     103, 51, 57, 59,103, 51, 57, 60, 99,103, 61,103, &
      51, 57, 62, 92,103, 63, 92,103, 54, 64, 74, 96, &
      99,101, 65, 94, 96,103,104, 51, 57, 66, 84, 92, &
      99,103, 57, 67, 74, 96, 99,103, 51, 57, 58, 59, &
      60, 68, 78, 82, 85, 92, 99,103, 58, 59, 61, 62, &
      63, 68, 69, 71, 73, 75, 76, 78, 79, 80, 81, 82 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      83, 84, 85, 86, 88, 89, 92, 99,103, 46, 60, 64, &
      66, 67, 68, 70, 74, 78, 81, 82, 83, 84, 85, 86, &
      89, 92, 96, 99,101,103, 71, 88, 92, 99,103, 43, &
      48, 49, 50, 61, 72, 76, 79, 80, 84, 87, 92, 99, &
     101,103, 73, 88, 92, 99,103, 60, 67, 74, 93, 95, &
      96, 97, 98, 99,100,101,103, 75, 88, 92, 99,103, &
      76, 88, 92, 99,103, 48, 50, 58, 59, 61, 72, 75, &
      76, 77, 79, 80, 84, 85, 87, 88, 90, 91, 92, 93, &
      94, 95, 96, 97, 98, 99,100,101,102,103,104, 75, &
      78, 84, 88, 92, 99,103, 79, 88, 92, 99,103, 80, &
      88, 92, 99,103, 51, 57, 58, 59, 62, 63, 67, 71, &
      74, 79, 80, 81, 82, 88, 92, 93, 95, 96, 97, 98, &
      99,100,101,103, 75, 82, 84, 88, 92, 99,103, 37, &
      44, 48, 50, 61, 71, 73, 76, 82, 83, 84, 87, 88, &
      89, 90, 91, 92, 93, 95, 97, 98, 99,100,103, 84, &
      88, 92, 99,103, 75, 84, 85, 88, 92, 99,103, 44, &
      48, 50, 52, 53, 56, 61, 63, 71, 72, 73, 75, 76, &
      78, 79, 80, 81, 82, 84, 85, 86, 87, 88, 90, 91, &
      92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103, &
     104, 49, 76, 79, 80, 82, 84, 87, 88, 92, 93, 99, &
     101,103,104, 30, 71, 73, 75, 76, 78, 79, 80, 84, &
      85, 88, 92, 93, 99,101,103, 43, 48, 50, 58, 59, &
      61, 62, 65, 71, 73, 76, 78, 79, 80, 82, 84, 85, &
      87, 88, 89, 90, 91, 92, 93, 94, 96, 99,101,103, &
     104, 48, 50, 61, 73, 76, 78, 82, 84, 85, 87, 88, &
      90, 91, 92, 93, 94, 99,101,102,103,104, 50, 57, &
      61, 75, 76, 79, 80, 82, 84, 85, 87, 88, 91, 92, &
      93, 94, 95, 97, 99,100,101,102,103,104, 62, 63, &
      71, 73, 75, 76, 78, 79, 80, 82, 84, 85, 88, 92, &
      93, 95, 96, 97, 98, 99,100,101,103, 47, 53, 77 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_2 = (/ &
      79, 80, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99,100,101,102,103,104, 37, 43, &
      44, 48, 50, 51, 57, 58, 59, 60, 61, 62, 63, 65, &
      67, 71, 73, 74, 75, 76, 78, 79, 80, 82, 84, 85, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100,101,102,103,104, 38, 45, 51, 58, 59, 61, &
      68, 72, 76, 78, 79, 80, 82, 83, 84, 85, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
     101,102,103,104, 36, 42, 44, 47, 51, 52, 53, 54, &
      55, 56, 57, 58, 59, 62, 63, 64, 65, 68, 69, 71, &
      73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      98, 99,100,101,102,103,104, 39, 78, 79, 80, 81, &
      82, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 98, 99,100,101,102,103,104, 41, 75, 78, 82, &
      84, 85, 88, 92, 93, 94, 95, 96, 97, 98, 99,100, &
     101,102,103,104, 46, 55, 60, 64, 66, 67, 68, 70, &
      71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, &
      96, 97, 98, 99,100,101,102,103,104, 40, 66, 84, &
      88, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102, &
     103,104, 38, 39, 40, 41, 46, 47, 49, 53, 54, 55, &
      64, 70, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, &
      97, 98, 99,100,101,102,103,104, 35, 49, 52, 61, &
      71, 72, 73, 75, 76, 79, 80, 83, 84, 85, 87, 88, &
      89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100, &
     101,102,103,104,  3,  5, 30, 31, 32, 33, 34, 35, &
      36, 37, 42, 43, 44, 47, 48, 50, 51, 52, 55, 56, &
      57, 58, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69 /)
  INTEGER, PARAMETER, DIMENSION(66) :: LU_ICOL_3 = (/ &
      70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
      95, 96, 97, 98, 99,100,101,102,103,104, 43, 48, &
      50, 51, 57, 61, 73, 75, 76, 79, 80, 82, 84, 85, &
      87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, &
      99,100,101,102,103,104 /)
  INTEGER, PARAMETER, DIMENSION(1146) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2, LU_ICOL_3 /)

  INTEGER, PARAMETER, DIMENSION(105) :: LU_CROW = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 92, 96,100,104,107,114,121,128, &
     177,191,195,197,201,203,205,207,210,212,215,217, &
     219,221,224,227,230,233,236,238,240,246,249,252, &
     254,258,260,262,266,270,275,279,284,286,290,294, &
     299,301,306,309,315,320,327,333,345,370,391,396, &
     411,416,428,433,438,468,475,480,485,509,516,540, &
     545,552,590,604,620,650,671,695,718,743,787,821, &
     872,897,917,958,975,1017,1049,1115,1147 /)

  INTEGER, PARAMETER, DIMENSION(105) :: LU_DIAG = (/ &
       1,  4,  8,  9, 10, 11, 12, 28, 38, 53, 56, 61, &
      66, 75, 82, 88, 92, 96,100,104,107,114,121,128, &
     181,191,195,197,201,203,205,208,210,213,215,217, &
     219,221,224,227,230,233,236,238,240,246,249,252, &
     255,258,260,262,266,270,275,279,284,288,292,296, &
     299,303,306,310,315,322,328,338,351,376,391,401, &
     411,418,428,433,446,469,475,480,496,510,525,540, &
     547,572,596,614,639,661,683,708,731,776,811,863, &
     889,910,952,970,1013,1046,1113,1146,1147 /)


END MODULE saprc99_mosaic_8bin_vbs2_JacobianSP

