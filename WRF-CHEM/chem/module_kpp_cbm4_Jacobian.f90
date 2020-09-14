MODULE cbm4_Jacobian
  USE cbm4_Parameters
  USE cbm4_JacobianSP
  IMPLICIT NONE
CONTAINS
SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: UV(NVAR)
  REAL(kind=dp) :: JUV(NVAR)
  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(25)
  JUV(2) = JVS(3)*UV(2)+JVS(4)*UV(27)+JVS(5)*UV(28)
  JUV(3) = JVS(6)*UV(3)+JVS(7)*UV(26)+JVS(8)*UV(32)
  JUV(4) = JVS(9)*UV(4)+JVS(10)*UV(14)+JVS(11)*UV(26)+JVS(12)*UV(27)+JVS(13)*UV(30)
  JUV(5) = JVS(14)*UV(5)+JVS(15)*UV(27)
  JUV(6) = JVS(16)*UV(6)+JVS(17)*UV(26)+JVS(18)*UV(30)
  JUV(7) = JVS(19)*UV(7)+JVS(20)*UV(27)
  JUV(8) = JVS(21)*UV(8)+JVS(22)*UV(13)+JVS(23)*UV(20)+JVS(24)*UV(22)+JVS(25)*UV(23)+JVS(26)*UV(27)+JVS(27)*UV(29)&
             &+JVS(28)*UV(30)+JVS(29)*UV(31)
  JUV(9) = JVS(30)*UV(9)+JVS(31)*UV(26)+JVS(32)*UV(27)+JVS(33)*UV(31)
  JUV(10) = JVS(34)*UV(10)+JVS(35)*UV(26)+JVS(36)*UV(27)+JVS(37)*UV(28)
  JUV(11) = JVS(38)*UV(5)+JVS(39)*UV(7)+JVS(40)*UV(11)+JVS(41)*UV(27)+JVS(42)*UV(31)
  JUV(12) = JVS(43)*UV(6)+JVS(44)*UV(12)+JVS(45)*UV(14)+JVS(46)*UV(21)+JVS(47)*UV(24)+JVS(48)*UV(26)+JVS(49)*UV(27)&
              &+JVS(50)*UV(30)
  JUV(13) = JVS(51)*UV(13)+JVS(52)*UV(20)+JVS(53)*UV(26)+JVS(54)*UV(27)
  JUV(14) = JVS(55)*UV(5)+JVS(56)*UV(7)+JVS(57)*UV(11)+JVS(58)*UV(14)+JVS(59)*UV(27)+JVS(60)*UV(30)
  JUV(15) = JVS(62)*UV(7)+JVS(63)*UV(15)+JVS(64)*UV(19)+JVS(65)*UV(22)+JVS(66)*UV(25)+JVS(67)*UV(27)
  JUV(16) = JVS(68)*UV(15)+JVS(69)*UV(16)+JVS(70)*UV(17)+JVS(71)*UV(19)+JVS(72)*UV(21)+JVS(73)*UV(22)+JVS(74)*UV(23)&
              &+JVS(75)*UV(24)+JVS(76)*UV(25)+JVS(77)*UV(27)+JVS(78)*UV(29)+JVS(79)*UV(30)
  JUV(17) = JVS(80)*UV(17)+JVS(81)*UV(22)+JVS(82)*UV(25)+JVS(83)*UV(27)+JVS(84)*UV(29)
  JUV(18) = JVS(85)*UV(5)+JVS(86)*UV(7)+JVS(87)*UV(13)+JVS(88)*UV(14)+JVS(89)*UV(15)+JVS(90)*UV(17)+JVS(91)*UV(18)&
              &+JVS(92)*UV(19)+JVS(93)*UV(20)+JVS(94)*UV(22)+JVS(95)*UV(23)+JVS(96)*UV(24)+JVS(97)*UV(25)+JVS(99)*UV(27)&
              &+JVS(100)*UV(28)+JVS(101)*UV(29)+JVS(102)*UV(30)+JVS(103)*UV(31)+JVS(104)*UV(32)
  JUV(19) = JVS(105)*UV(11)+JVS(106)*UV(14)+JVS(107)*UV(19)+JVS(108)*UV(25)+JVS(109)*UV(27)+JVS(111)*UV(31)
  JUV(20) = JVS(112)*UV(7)+JVS(113)*UV(13)+JVS(114)*UV(20)+JVS(115)*UV(22)+JVS(116)*UV(23)+JVS(117)*UV(25)+JVS(119)&
              &*UV(27)+JVS(120)*UV(29)+JVS(121)*UV(30)
  JUV(21) = JVS(122)*UV(17)+JVS(123)*UV(19)+JVS(124)*UV(21)+JVS(125)*UV(22)+JVS(126)*UV(23)+JVS(127)*UV(24)+JVS(128)&
              &*UV(25)+JVS(129)*UV(27)+JVS(130)*UV(28)+JVS(131)*UV(29)+JVS(132)*UV(30)+JVS(133)*UV(31)+JVS(134)*UV(32)
  JUV(22) = JVS(135)*UV(22)+JVS(136)*UV(25)+JVS(137)*UV(27)+JVS(138)*UV(29)+JVS(139)*UV(30)
  JUV(23) = JVS(140)*UV(22)+JVS(141)*UV(23)+JVS(142)*UV(25)+JVS(143)*UV(27)+JVS(144)*UV(29)+JVS(145)*UV(30)
  JUV(24) = JVS(146)*UV(13)+JVS(147)*UV(17)+JVS(148)*UV(19)+JVS(149)*UV(20)+JVS(150)*UV(22)+JVS(151)*UV(23)+JVS(152)&
              &*UV(24)+JVS(153)*UV(25)+JVS(155)*UV(27)+JVS(156)*UV(29)+JVS(157)*UV(30)
  JUV(25) = JVS(159)*UV(17)+JVS(160)*UV(19)+JVS(161)*UV(22)+JVS(162)*UV(23)+JVS(163)*UV(25)+JVS(164)*UV(26)+JVS(165)&
              &*UV(27)+JVS(166)*UV(28)+JVS(167)*UV(29)+JVS(169)*UV(31)
  JUV(26) = JVS(170)*UV(3)+JVS(171)*UV(4)+JVS(172)*UV(6)+JVS(173)*UV(9)+JVS(174)*UV(10)+JVS(175)*UV(11)+JVS(176)*UV(13)&
              &+JVS(178)*UV(18)+JVS(182)*UV(23)+JVS(184)*UV(25)+JVS(185)*UV(26)+JVS(186)*UV(27)+JVS(187)*UV(28)+JVS(188)&
              &*UV(29)+JVS(189)*UV(30)+JVS(190)*UV(31)+JVS(191)*UV(32)
  JUV(27) = JVS(192)*UV(1)+JVS(193)*UV(2)+JVS(194)*UV(5)+JVS(195)*UV(7)+JVS(196)*UV(9)+JVS(197)*UV(10)+JVS(198)*UV(12)&
              &+JVS(199)*UV(14)+JVS(200)*UV(15)+JVS(201)*UV(16)+JVS(202)*UV(17)+JVS(203)*UV(19)+JVS(204)*UV(20)+JVS(205)&
              &*UV(21)+JVS(206)*UV(22)+JVS(207)*UV(23)+JVS(208)*UV(24)+JVS(209)*UV(25)+JVS(210)*UV(26)+JVS(211)*UV(27)&
              &+JVS(212)*UV(28)+JVS(213)*UV(29)+JVS(215)*UV(31)+JVS(216)*UV(32)
  JUV(28) = JVS(217)*UV(2)+JVS(218)*UV(5)+JVS(219)*UV(7)+JVS(220)*UV(10)+JVS(221)*UV(11)+JVS(222)*UV(13)+JVS(223)*UV(14)&
              &+JVS(224)*UV(15)+JVS(225)*UV(16)+JVS(226)*UV(17)+JVS(227)*UV(19)+JVS(228)*UV(20)+JVS(229)*UV(21)+JVS(230)&
              &*UV(22)+JVS(231)*UV(23)+JVS(232)*UV(24)+JVS(233)*UV(25)+JVS(234)*UV(26)+JVS(235)*UV(27)+JVS(236)*UV(28)&
              &+JVS(237)*UV(29)+JVS(238)*UV(30)+JVS(239)*UV(31)+JVS(240)*UV(32)
  JUV(29) = JVS(241)*UV(1)+JVS(242)*UV(17)+JVS(243)*UV(21)+JVS(244)*UV(22)+JVS(245)*UV(23)+JVS(246)*UV(24)+JVS(247)&
              &*UV(25)+JVS(248)*UV(26)+JVS(251)*UV(29)+JVS(252)*UV(30)+JVS(253)*UV(31)
  JUV(30) = JVS(255)*UV(6)+JVS(256)*UV(12)+JVS(257)*UV(14)+JVS(258)*UV(21)+JVS(259)*UV(22)+JVS(260)*UV(23)+JVS(261)&
              &*UV(24)+JVS(262)*UV(25)+JVS(263)*UV(26)+JVS(264)*UV(27)+JVS(266)*UV(29)+JVS(267)*UV(30)+JVS(268)*UV(31)
  JUV(31) = JVS(270)*UV(8)+JVS(271)*UV(9)+JVS(272)*UV(11)+JVS(274)*UV(18)+JVS(280)*UV(25)+JVS(281)*UV(26)+JVS(282)&
              &*UV(27)+JVS(283)*UV(28)+JVS(284)*UV(29)+JVS(285)*UV(30)+JVS(286)*UV(31)+JVS(287)*UV(32)
  JUV(32) = JVS(288)*UV(3)+JVS(289)*UV(15)+JVS(290)*UV(19)+JVS(291)*UV(22)+JVS(292)*UV(24)+JVS(293)*UV(25)+JVS(294)&
              &*UV(26)+JVS(295)*UV(27)+JVS(296)*UV(28)+JVS(297)*UV(29)+JVS(298)*UV(30)+JVS(299)*UV(31)+JVS(300)*UV(32)
END SUBROUTINE Jac_SP_Vec
SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: UV(NVAR)
  REAL(kind=dp) :: JTUV(NVAR)
  JTUV(1) = JVS(1)*UV(1)+JVS(192)*UV(27)+JVS(241)*UV(29)
  JTUV(2) = JVS(3)*UV(2)+JVS(193)*UV(27)+JVS(217)*UV(28)
  JTUV(3) = JVS(6)*UV(3)+JVS(170)*UV(26)+JVS(288)*UV(32)
  JTUV(4) = JVS(9)*UV(4)+JVS(171)*UV(26)
  JTUV(5) = JVS(14)*UV(5)+JVS(38)*UV(11)+JVS(55)*UV(14)+JVS(85)*UV(18)+JVS(194)*UV(27)+JVS(218)*UV(28)
  JTUV(6) = JVS(16)*UV(6)+JVS(43)*UV(12)+JVS(172)*UV(26)+JVS(255)*UV(30)
  JTUV(7) = JVS(19)*UV(7)+JVS(39)*UV(11)+JVS(56)*UV(14)+JVS(62)*UV(15)+JVS(86)*UV(18)+JVS(112)*UV(20)+JVS(195)*UV(27)&
              &+JVS(219)*UV(28)
  JTUV(8) = JVS(21)*UV(8)+JVS(270)*UV(31)
  JTUV(9) = JVS(30)*UV(9)+JVS(173)*UV(26)+JVS(196)*UV(27)+JVS(271)*UV(31)
  JTUV(10) = JVS(34)*UV(10)+JVS(174)*UV(26)+JVS(197)*UV(27)+JVS(220)*UV(28)
  JTUV(11) = JVS(40)*UV(11)+JVS(57)*UV(14)+JVS(105)*UV(19)+JVS(175)*UV(26)+JVS(221)*UV(28)+JVS(272)*UV(31)
  JTUV(12) = JVS(44)*UV(12)+JVS(198)*UV(27)+JVS(256)*UV(30)
  JTUV(13) = JVS(22)*UV(8)+JVS(51)*UV(13)+JVS(87)*UV(18)+JVS(113)*UV(20)+JVS(146)*UV(24)+JVS(176)*UV(26)+JVS(222)*UV(28)
  JTUV(14) = JVS(10)*UV(4)+JVS(45)*UV(12)+JVS(58)*UV(14)+JVS(88)*UV(18)+JVS(106)*UV(19)+JVS(199)*UV(27)+JVS(223)*UV(28)&
               &+JVS(257)*UV(30)
  JTUV(15) = JVS(63)*UV(15)+JVS(68)*UV(16)+JVS(89)*UV(18)+JVS(200)*UV(27)+JVS(224)*UV(28)+JVS(289)*UV(32)
  JTUV(16) = JVS(69)*UV(16)+JVS(201)*UV(27)+JVS(225)*UV(28)
  JTUV(17) = JVS(70)*UV(16)+JVS(80)*UV(17)+JVS(90)*UV(18)+JVS(122)*UV(21)+JVS(147)*UV(24)+JVS(159)*UV(25)+JVS(202)&
               &*UV(27)+JVS(226)*UV(28)+JVS(242)*UV(29)
  JTUV(18) = JVS(91)*UV(18)+JVS(178)*UV(26)+JVS(274)*UV(31)
  JTUV(19) = JVS(64)*UV(15)+JVS(71)*UV(16)+JVS(92)*UV(18)+JVS(107)*UV(19)+JVS(123)*UV(21)+JVS(148)*UV(24)+JVS(160)&
               &*UV(25)+JVS(203)*UV(27)+JVS(227)*UV(28)+JVS(290)*UV(32)
  JTUV(20) = JVS(23)*UV(8)+JVS(52)*UV(13)+JVS(93)*UV(18)+JVS(114)*UV(20)+JVS(149)*UV(24)+JVS(204)*UV(27)+JVS(228)*UV(28)
  JTUV(21) = JVS(46)*UV(12)+JVS(72)*UV(16)+JVS(124)*UV(21)+JVS(205)*UV(27)+JVS(229)*UV(28)+JVS(243)*UV(29)+JVS(258)&
               &*UV(30)
  JTUV(22) = JVS(24)*UV(8)+JVS(65)*UV(15)+JVS(73)*UV(16)+JVS(81)*UV(17)+JVS(94)*UV(18)+JVS(115)*UV(20)+JVS(125)*UV(21)&
               &+JVS(135)*UV(22)+JVS(140)*UV(23)+JVS(150)*UV(24)+JVS(161)*UV(25)+JVS(206)*UV(27)+JVS(230)*UV(28)+JVS(244)&
               &*UV(29)+JVS(259)*UV(30)+JVS(291)*UV(32)
  JTUV(23) = JVS(25)*UV(8)+JVS(74)*UV(16)+JVS(95)*UV(18)+JVS(116)*UV(20)+JVS(126)*UV(21)+JVS(141)*UV(23)+JVS(151)*UV(24)&
               &+JVS(162)*UV(25)+JVS(182)*UV(26)+JVS(207)*UV(27)+JVS(231)*UV(28)+JVS(245)*UV(29)+JVS(260)*UV(30)
  JTUV(24) = JVS(47)*UV(12)+JVS(75)*UV(16)+JVS(96)*UV(18)+JVS(127)*UV(21)+JVS(152)*UV(24)+JVS(208)*UV(27)+JVS(232)&
               &*UV(28)+JVS(246)*UV(29)+JVS(261)*UV(30)+JVS(292)*UV(32)
  JTUV(25) = JVS(2)*UV(1)+JVS(66)*UV(15)+JVS(76)*UV(16)+JVS(82)*UV(17)+JVS(97)*UV(18)+JVS(108)*UV(19)+JVS(117)*UV(20)&
               &+JVS(128)*UV(21)+JVS(136)*UV(22)+JVS(142)*UV(23)+JVS(153)*UV(24)+JVS(163)*UV(25)+JVS(184)*UV(26)+JVS(209)&
               &*UV(27)+JVS(233)*UV(28)+JVS(247)*UV(29)+JVS(262)*UV(30)+JVS(280)*UV(31)+JVS(293)*UV(32)
  JTUV(26) = JVS(7)*UV(3)+JVS(11)*UV(4)+JVS(17)*UV(6)+JVS(31)*UV(9)+JVS(35)*UV(10)+JVS(48)*UV(12)+JVS(53)*UV(13)&
               &+JVS(164)*UV(25)+JVS(185)*UV(26)+JVS(210)*UV(27)+JVS(234)*UV(28)+JVS(248)*UV(29)+JVS(263)*UV(30)+JVS(281)&
               &*UV(31)+JVS(294)*UV(32)
  JTUV(27) = JVS(4)*UV(2)+JVS(12)*UV(4)+JVS(15)*UV(5)+JVS(20)*UV(7)+JVS(26)*UV(8)+JVS(32)*UV(9)+JVS(36)*UV(10)+JVS(41)&
               &*UV(11)+JVS(49)*UV(12)+JVS(54)*UV(13)+JVS(59)*UV(14)+JVS(67)*UV(15)+JVS(77)*UV(16)+JVS(83)*UV(17)+JVS(99)&
               &*UV(18)+JVS(109)*UV(19)+JVS(119)*UV(20)+JVS(129)*UV(21)+JVS(137)*UV(22)+JVS(143)*UV(23)+JVS(155)*UV(24)&
               &+JVS(165)*UV(25)+JVS(186)*UV(26)+JVS(211)*UV(27)+JVS(235)*UV(28)+JVS(264)*UV(30)+JVS(282)*UV(31)+JVS(295)&
               &*UV(32)
  JTUV(28) = JVS(5)*UV(2)+JVS(37)*UV(10)+JVS(100)*UV(18)+JVS(130)*UV(21)+JVS(166)*UV(25)+JVS(187)*UV(26)+JVS(212)*UV(27)&
               &+JVS(236)*UV(28)+JVS(283)*UV(31)+JVS(296)*UV(32)
  JTUV(29) = JVS(27)*UV(8)+JVS(78)*UV(16)+JVS(84)*UV(17)+JVS(101)*UV(18)+JVS(120)*UV(20)+JVS(131)*UV(21)+JVS(138)*UV(22)&
               &+JVS(144)*UV(23)+JVS(156)*UV(24)+JVS(167)*UV(25)+JVS(188)*UV(26)+JVS(213)*UV(27)+JVS(237)*UV(28)+JVS(251)&
               &*UV(29)+JVS(266)*UV(30)+JVS(284)*UV(31)+JVS(297)*UV(32)
  JTUV(30) = JVS(13)*UV(4)+JVS(18)*UV(6)+JVS(28)*UV(8)+JVS(50)*UV(12)+JVS(60)*UV(14)+JVS(79)*UV(16)+JVS(102)*UV(18)&
               &+JVS(121)*UV(20)+JVS(132)*UV(21)+JVS(139)*UV(22)+JVS(145)*UV(23)+JVS(157)*UV(24)+JVS(189)*UV(26)+JVS(238)&
               &*UV(28)+JVS(252)*UV(29)+JVS(267)*UV(30)+JVS(285)*UV(31)+JVS(298)*UV(32)
  JTUV(31) = JVS(29)*UV(8)+JVS(33)*UV(9)+JVS(42)*UV(11)+JVS(103)*UV(18)+JVS(111)*UV(19)+JVS(133)*UV(21)+JVS(169)*UV(25)&
               &+JVS(190)*UV(26)+JVS(215)*UV(27)+JVS(239)*UV(28)+JVS(253)*UV(29)+JVS(268)*UV(30)+JVS(286)*UV(31)+JVS(299)&
               &*UV(32)
  JTUV(32) = JVS(8)*UV(3)+JVS(104)*UV(18)+JVS(134)*UV(21)+JVS(191)*UV(26)+JVS(216)*UV(27)+JVS(240)*UV(28)+JVS(287)&
               &*UV(31)+JVS(300)*UV(32)
END SUBROUTINE JacTR_SP_Vec
END MODULE cbm4_Jacobian
