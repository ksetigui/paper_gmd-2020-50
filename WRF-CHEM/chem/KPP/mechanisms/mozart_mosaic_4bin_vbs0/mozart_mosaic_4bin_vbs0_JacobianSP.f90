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
! File                 : mozart_mosaic_4bin_vbs0_JacobianSP.f90
! Time                 : Mon Feb 27 11:08:41 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/mozart_mosaic_4bin_vbs0
! Equation file        : mozart_mosaic_4bin_vbs0.kpp
! Output root filename : mozart_mosaic_4bin_vbs0
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_mosaic_4bin_vbs0_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4, &
       5,  5,  5,  6,  6,  6,  7,  7,  8,  8,  8,  9, &
       9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, &
      15, 16, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, &
      19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 22, 23, &
      23, 23, 23, 24, 24, 24, 25, 25, 25, 25, 26, 26, &
      26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 29, &
      29, 29, 29, 30, 30, 30, 30, 31, 31, 31, 31, 32, &
      32, 32, 32, 33, 33, 33, 33, 34, 34, 34, 34, 35, &
      35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 37, 37, &
      37, 37, 37, 37, 38, 38, 38, 38, 39, 39, 39, 39, &
      39, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 42, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, &
      43, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 47, 47, &
      47, 47, 48, 48, 48, 48, 48, 49, 49, 49, 49, 50, &
      50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 52, &
      52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, &
      53, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 55, &
      55, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, &
      57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, &
      60, 60, 60, 61, 61, 61, 61, 61, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, &
      63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, &
      66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, &
      69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 71, 71, &
      71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86 /)
  INTEGER, PARAMETER, DIMENSION(138) :: LU_IROW_2 = (/ &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89 /)
  INTEGER, PARAMETER, DIMENSION(858) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 16, 84,  2, 72, 87,  3, 72, 86,  4, 67, 84, &
       5, 54, 84,  6,  7, 84,  7, 84,  8,  9, 84,  9, &
      84, 10, 84, 11, 84, 12, 84, 13, 84, 14, 84, 15, &
      55, 16, 24, 81, 84, 17, 84, 17, 18, 84, 19, 84, &
      86, 20, 39, 87, 21, 81, 88, 18, 22, 84, 88, 23, &
      76, 84, 86, 24, 81, 84, 25, 64, 82, 84, 26, 78, &
      84, 86, 27, 48, 55, 79, 84, 28, 64, 84, 86, 29, &
      66, 84, 86, 30, 72, 84, 86, 13, 31, 84, 87, 32, &
      57, 84, 86, 33, 82, 84, 86, 34, 84, 86, 88, 35, &
      68, 80, 81, 82, 84, 87, 89, 36, 83, 84, 37, 63, &
      69, 81, 84, 87, 38, 52, 84, 86, 36, 39, 83, 84, &
      87, 40, 84, 85, 86, 89, 41, 56, 63, 84, 87, 36, &
      42, 67, 69, 82, 83, 84, 85, 86, 89, 38, 43, 52, &
      53, 62, 84, 86, 87, 44, 71, 84, 86, 21, 24, 45, &
      70, 75, 79, 81, 84, 88, 46, 59, 84, 86, 47, 84, &
      88, 89, 48, 55, 69, 83, 84, 49, 80, 84, 86, 50, &
      84, 85, 88, 51, 64, 72, 76, 78, 80, 82, 84, 17, &
      38, 52, 84, 86, 87, 22, 38, 52, 53, 84, 86, 87, &
      88, 54, 81, 83, 84, 15, 27, 48, 55, 69, 79, 83, &
      84, 56, 63, 84, 86, 32, 41, 56, 57, 63, 84, 86, &
      87, 29, 31, 44, 56, 58, 63, 66, 71, 82, 84, 86, &
      87, 46, 59, 69, 84, 86, 87, 55, 60, 69, 79, 83, &
      84, 86, 88, 61, 67, 81, 86, 87, 20, 39, 62, 76, &
      78, 81, 82, 83, 84, 87, 89, 14, 56, 63, 84, 86, &
      87, 11, 28, 41, 56, 63, 64, 82, 84, 86, 87, 36, &
      43, 52, 53, 62, 65, 67, 68, 69, 70, 74, 75, 76, &
      77, 78, 79, 81, 82, 83, 84, 86, 87, 88, 89, 12, &
      29, 66, 82, 84, 86, 87, 67, 81, 83, 84, 61, 67, &
      68, 78, 80, 81, 83, 84, 86, 87, 67, 69, 77, 81, &
      83, 84, 25, 28, 31, 32, 46, 56, 57, 59, 63, 64 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      66, 69, 70, 77, 81, 82, 83, 84, 86, 87, 44, 54, &
      71, 81, 83, 84, 86, 87, 30, 58, 63, 66, 71, 72, &
      81, 82, 83, 84, 86, 87, 46, 50, 59, 69, 72, 73, &
      76, 77, 78, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      44, 49, 54, 61, 67, 71, 74, 80, 81, 82, 83, 84, &
      86, 87, 89, 37, 38, 52, 53, 63, 69, 72, 73, 74, &
      75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, &
      88, 89, 23, 35, 49, 68, 76, 78, 80, 81, 82, 83, &
      84, 86, 87, 89, 44, 49, 54, 61, 67, 71, 77, 80, &
      81, 82, 83, 84, 86, 87, 89, 26, 74, 77, 78, 80, &
      81, 82, 83, 84, 86, 87, 89, 20, 30, 31, 33, 36, &
      39, 40, 46, 47, 48, 49, 50, 51, 55, 56, 59, 61, &
      62, 63, 64, 66, 67, 68, 69, 72, 73, 74, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      49, 67, 80, 81, 82, 83, 84, 86, 87, 89, 21, 24, &
      34, 45, 47, 50, 54, 61, 67, 68, 69, 70, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 33, 40, 42, 47, 48, 55, 58, 63, 64, 66, 67, &
      69, 70, 71, 72, 76, 77, 78, 79, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89, 36, 54, 60, 67, 69, 74, &
      77, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      10, 11, 12, 13, 14, 16, 17, 18, 19, 23, 24, 25, &
      26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, &
      40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, &
      52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, &
      65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 26, 50, 67, 74, 78, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 17, 19, 20, 22, 24, 25, 26, 27, &
      28, 29, 31, 33, 34, 36, 39, 43, 44, 46, 48, 49 /)
  INTEGER, PARAMETER, DIMENSION(138) :: LU_ICOL_2 = (/ &
      50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, &
      63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 15, 31, 39, 52, 55, 57, 59, 60, 61, &
      63, 64, 66, 67, 69, 71, 72, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 87, 88, 89, 21, 22, 31, &
      34, 37, 39, 45, 47, 50, 52, 54, 57, 59, 60, 61, &
      63, 64, 66, 67, 68, 69, 70, 71, 72, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      30, 32, 40, 41, 47, 53, 56, 57, 58, 63, 66, 70, &
      71, 72, 73, 74, 75, 76, 77, 78, 80, 81, 82, 83, &
      84, 85, 86, 87, 88, 89 /)
  INTEGER, PARAMETER, DIMENSION(858) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)

  INTEGER, PARAMETER, DIMENSION(90) :: LU_CROW = (/ &
       1,  4,  7, 10, 13, 16, 19, 21, 24, 26, 28, 30, &
      32, 34, 36, 38, 42, 44, 47, 50, 53, 56, 60, 64, &
      67, 71, 75, 80, 84, 88, 92, 96,100,104,108,116, &
     119,125,129,134,139,144,154,162,166,175,179,183, &
     188,192,196,204,210,218,222,230,234,242,254,260, &
     268,273,284,290,300,324,331,335,345,351,371,379, &
     391,409,424,447,461,476,488,529,539,566,595,613, &
     686,701,760,790,829,859 /)

  INTEGER, PARAMETER, DIMENSION(90) :: LU_DIAG = (/ &
       1,  4,  7, 10, 13, 16, 19, 21, 24, 26, 28, 30, &
      32, 34, 36, 38, 42, 45, 47, 50, 53, 57, 60, 64, &
      67, 71, 75, 80, 84, 88, 93, 96,100,104,108,116, &
     119,125,130,134,139,145,155,162,168,175,179,183, &
     188,192,196,206,213,218,225,230,237,246,255,261, &
     268,275,286,295,305,326,331,337,346,363,373,384, &
     396,415,433,451,467,479,518,531,557,587,606,680, &
     696,756,787,827,858,859 /)


END MODULE mozart_mosaic_4bin_vbs0_JacobianSP

