MODULE cbmz_mosaic_4bin_vbs9_JacobianSP
  PUBLIC
  SAVE
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, &
       8, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, &
      10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, &
      13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, &
      15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
      18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 21, &
      21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, &
      27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, &
      33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, &
      39, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, &
      44, 45, 45, 46, 46, 47, 47, 47, 47, 48, 48, 49, &
      49, 49, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, &
      52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, &
      54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56, 56, &
      56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, &
      59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, &
      61, 61, 62, 62, 62, 62, 62, 62, 63, 63, 64, 64, &
      64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, &
      70, 70, 71, 71, 72, 72, 72, 72, 73, 73, 73, 74, &
      74, 74, 74, 74, 74, 75, 75, 75, 75, 76, 76, 76, &
      76, 76, 76, 77, 77, 77, 77, 78, 78, 78, 78, 78, &
      78, 79, 79, 79, 79, 80, 80, 80, 80, 80, 80, 81 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      81, 81, 81, 82, 82, 82, 82, 82, 82, 83, 83, 83, &
      83, 84, 84, 84, 84, 85, 85, 85, 85, 85, 85, 86, &
      86, 86, 86, 86, 86, 87, 87, 87, 88, 88, 88, 89, &
      89, 90, 90, 90, 91, 91, 92, 92, 92, 92, 92, 93, &
      93, 93, 93, 93, 94, 94, 94, 94, 95, 95, 95, 95, &
      96, 96, 96, 96, 97, 97, 97, 97, 97, 98, 98, 98, &
      98, 98, 98, 99, 99, 99, 99, 99, 99, 99, 99, 99, &
      99, 99, 99,100,100,100,100,101,101,101,101,102, &
     102,102,102,102,103,103,103,104,104,104,104,104, &
     104,105,105,105,105,105,106,106,106,106,106,106, &
     107,107,107,107,107,107,107,107,107,107,107,108, &
     108,108,108,108,108,108,108,108,108,108,108,108, &
     108,108,109,109,109,109,109,109,109,109,109,109, &
     109,109,109,109,109,109,109,110,110,110,110,110, &
     110,110,111,111,111,111,111,112,112,112,112,112, &
     113,113,113,113,113,114,114,114,114,115,115,115, &
     115,116,116,116,116,116,116,116,116,116,116,116, &
     116,116,116,116,116,116,116,116,117,117,117,117, &
     117,117,117,117,117,117,117,117,117,117,117,117, &
     117,117,117,118,118,118,118,118,118,118,118,118, &
     118,118,118,118,119,119,119,119,120,120,120,120, &
     120,120,120,120,120,121,121,121,121,121,121,121, &
     121,121,121,121,121,121,121,122,122,122,122,122, &
     122,122,122,122,122,122,122,122,123,123,123,123, &
     123,123,123,123,123,123,123,123,123,123,123,123, &
     123,123,123,123,123,124,124,124,124,124,124,124, &
     124,124,124,124,124,124,124,124,124,124,124,125, &
     125,125,125,125,125,125,125,125,126,126,126,126, &
     126,126,126,126,126,126,126,126,126,126,126,126, &
     127,127,127,127,127,127,127,127,127,127,127,127 /)
  INTEGER, PARAMETER, DIMENSION(309) :: LU_IROW_2 = (/ &
     127,127,127,127,128,128,128,128,128,128,128,128, &
     128,128,128,128,128,128,128,128,128,128,129,129, &
     129,129,129,129,129,129,129,129,129,129,129,129, &
     129,129,129,129,129,129,129,129,129,129,129,129, &
     129,129,129,129,129,129,129,130,130,130,130,130, &
     130,130,130,130,130,130,130,130,130,130,130,130, &
     130,130,130,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,131,131,131, &
     131,131,131,131,131,131,131,131,131,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     132,132,132,132,132,132,132,132,132,132,132,132, &
     133,133,133,133,133,133,133,133,133,133,133,133, &
     133,133,133,133,133,133,133,133,133,133,133,133, &
     133,134,134,134,134,134,134,134,134,134,134,134, &
     134,134,134,135,135,135,135,135,135,135,135,135, &
     135,135,135,135,135,136,136,136,136,136,136,136, &
     136,136,136,136,136,136,136,136,136,136,136,136, &
     136,136,136,136,136,136,136,136,136 /)
  INTEGER, PARAMETER, DIMENSION(1029) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 37,131, 2, 3, 4, 5, 6, 7, 8,103,114, &
     127, 9,114,115,119,125,127,130,132, 10, 89, 91, &
     105,131,136, 11, 89, 91,105,131,136, 12,109,131, &
      13,114,119,127,131,136, 14, 94,127,131, 15, 94, &
     127,131, 16, 95,127,131, 17, 95,127,131, 18, 37, &
      48, 88, 89, 91, 93, 96, 98,100,101,102,103,104, &
     105,107,108,109,110,114,115,116,118,119,121,123, &
     125,126,127,129,131,132,133,135,136, 4, 5, 6, &
       7, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, &
      30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, &
      44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, &
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
      81, 82, 83, 84, 85, 86, 94, 95,131, 20,131, 21, &
     131, 22,131, 23,131, 24,131, 25,131, 26,131, 27, &
     131, 28,131, 29,131, 30,131, 31,131, 32,131, 33, &
     131, 34,131, 35,131, 36,127, 37,131, 38,131, 38, &
      39, 60,131, 40,131, 41,131, 42,131, 43,131, 44, &
     131, 45,131, 46,131, 27, 46, 47,131, 48,131, 46, &
      49,131, 26, 45, 47, 49, 50,131, 45, 49, 51,131, &
      25, 44, 50, 51, 52,131, 44, 51, 53,131, 24, 43, &
      52, 53, 54,131, 43, 53, 55,131, 23, 42, 54, 55, &
      56,131, 42, 55, 57,131, 22, 41, 56, 57, 58,131, &
      41, 57, 59,131, 40, 59, 60,131, 20, 38, 60, 61, &
      62,131, 21, 40, 58, 59, 62,131, 63,131, 63, 64, &
      84,131, 65,131, 66,131, 67,131, 68,131, 69,131, &
      70,131, 71,131, 35, 71, 72,131, 71, 73,131, 34, &
      70, 72, 73, 74,131, 70, 73, 75,131, 33, 69, 74, &
      75, 76,131, 69, 75, 77,131, 32, 68, 76, 77, 78, &
     131, 68, 77, 79,131, 31, 67, 78, 79, 80,131, 67 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      79, 81,131, 30, 66, 80, 81, 82,131, 66, 81, 83, &
     131, 65, 83, 84,131, 28, 63, 84, 85, 86,131, 29, &
      65, 82, 83, 86,131, 87,129,130, 88,131,132, 89, &
     131, 90,129,136, 91,131, 92,105,129,131,136, 93, &
     114,119,127,131, 94,127,131,136, 95,127,131,136, &
      96,129,131,132, 89, 91, 97,131,133, 98,114,119, &
     122,127,131, 99,109,114,119,120,126,127,128,131, &
     133,135,136,100,124,131,132,101,122,131,132,102, &
     129,131,132,133,103,127,131,104,114,119,124,127, &
     131, 89, 91,105,131,136, 36,106,127,129,133,136, &
      90,105,107,116,121,123,125,129,131,132,136,103, &
     108,110,112,114,115,116,119,121,123,125,127,131, &
     133,136, 91, 99,109,111,113,114,119,120,125,126, &
     127,128,131,132,133,135,136, 97,105,110,127,131, &
     133,136,111,115,132,133,136,112,125,131,132,133, &
     113,115,131,132,133,114,127,131,136,115,127,131, &
     136,100,103,104,110,112,113,114,115,116,119,120, &
     124,125,127,131,132,133,134,136, 89, 91,103,105, &
     110,114,115,117,119,121,125,126,127,128,131,132, &
     133,135,136,112,118,119,125,126,127,128,131,132, &
     133,134,135,136,119,127,131,136,114,119,120,126, &
     127,131,132,133,136, 91,110,112,114,119,121,125, &
     127,131,132,133,134,135,136, 98,101,114,119,122, &
     126,127,128,131,132,133,135,136, 48,101,103,110, &
     111,112,114,115,119,120,122,123,125,126,127,128, &
     131,132,133,135,136, 93,100,114,118,119,123,124, &
     125,126,127,128,130,131,132,133,134,135,136,111, &
     113,115,125,127,131,132,133,136, 92, 97,105,111, &
     113,115,120,125,126,127,128,129,131,132,133,136, &
      94, 95,103,106,110,114,115,119,125,127,129,130 /)
  INTEGER, PARAMETER, DIMENSION(309) :: LU_ICOL_2 = (/ &
     131,132,133,136,109,111,113,114,115,119,120,125, &
     126,127,128,129,130,131,132,133,135,136, 87, 90, &
      92, 96, 97,102,105,106,107,111,112,113,115,116, &
     117,119,120,121,122,123,124,125,126,127,128,129, &
     130,131,132,133,134,135,136, 87,110,114,115,118, &
     119,121,123,125,126,127,128,129,130,131,132,133, &
     134,135,136, 4, 5, 6, 7, 20, 21, 22, 23, 24, &
      25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, &
      37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, &
      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
      61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, &
      73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 88, 89, 91, 93, 94, 95, 96, 98,100,101, &
     102,103,104,105,107,108,109,110,111,112,113,114, &
     115,116,118,119,120,121,122,123,124,125,126,127, &
     128,129,130,131,132,133,134,135,136, 37, 48, 88, &
      89, 91, 96, 97,100,101,103,104,105,108,110,111, &
     112,113,114,115,116,117,119,120,121,122,123,124, &
     125,126,127,128,129,130,131,132,133,134,135,136, &
      97,102,106,111,112,113,115,117,119,120,121,122, &
     124,125,126,127,128,129,130,131,132,133,134,135, &
     136,118,119,125,126,127,128,129,130,131,132,133, &
     134,135,136,112,113,115,125,127,128,129,130,131, &
     132,133,134,135,136, 90, 94, 95,105,106,107,114, &
     115,116,117,119,120,121,122,123,124,125,126,127, &
     128,129,130,131,132,133,134,135,136 /)
  INTEGER, PARAMETER, DIMENSION(1029) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)
  INTEGER, PARAMETER, DIMENSION(137) :: LU_CROW = (/ &
       1, 4, 5, 6, 7, 8, 9, 10, 14, 22, 28, 34, &
      37, 43, 47, 51, 55, 59, 94,166,168,170,172,174, &
     176,178,180,182,184,186,188,190,192,194,196,198, &
     200,202,204,208,210,212,214,216,218,220,222,226, &
     228,231,237,241,247,251,257,261,267,271,277,281, &
     285,291,297,299,303,305,307,309,311,313,315,317, &
     321,324,330,334,340,344,350,354,360,364,370,374, &
     378,384,390,393,396,398,401,403,408,413,417,421, &
     425,430,436,448,452,456,461,464,470,475,481,492, &
     507,524,531,536,541,546,550,554,573,592,605,609, &
     618,632,645,666,684,693,709,725,743,776,796,910, &
     949,974,988,1002,1030 /)
  INTEGER, PARAMETER, DIMENSION(137) :: LU_DIAG = (/ &
       1, 4, 5, 6, 7, 8, 9, 10, 14, 22, 28, 34, &
      37, 43, 47, 51, 55, 59, 98,166,168,170,172,174, &
     176,178,180,182,184,186,188,190,192,194,196,198, &
     200,202,205,208,210,212,214,216,218,220,224,226, &
     229,235,239,245,249,255,259,265,269,275,279,283, &
     288,295,297,300,303,305,307,309,311,313,315,319, &
     322,328,332,338,342,348,352,358,362,368,372,376, &
     381,388,390,393,396,398,401,403,408,413,417,421, &
     427,430,436,448,452,456,461,464,472,476,483,493, &
     509,526,531,536,541,546,550,562,580,593,605,611, &
     623,636,656,672,687,701,718,735,768,789,904,944, &
     970,985,1000,1029,1030 /)
END MODULE cbmz_mosaic_4bin_vbs9_JacobianSP
