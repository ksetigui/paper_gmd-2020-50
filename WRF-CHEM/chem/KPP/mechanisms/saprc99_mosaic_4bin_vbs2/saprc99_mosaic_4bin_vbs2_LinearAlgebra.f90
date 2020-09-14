! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Linear Algebra Data and Routines File
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
! File                 : saprc99_mosaic_4bin_vbs2_LinearAlgebra.f90
! Time                 : Wed Jun 14 17:35:30 2017
! Working directory    : /data/ksetigui/setigui/WRFnew/chem/KPP/mechanisms/saprc99_mosaic_4bin_vbs2
! Equation file        : saprc99_mosaic_4bin_vbs2.kpp
! Output root filename : saprc99_mosaic_4bin_vbs2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_mosaic_4bin_vbs2_LinearAlgebra

  USE saprc99_mosaic_4bin_vbs2_Parameters
  USE saprc99_mosaic_4bin_vbs2_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! KppSolveTR - sparse, transposed back substitution
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      X         - Vector for variables
!      XX        - Vector for output variables
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE KppSolveTR ( JVS, X, XX )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! X - Vector for variables
  REAL(kind=dp) :: X(NVAR)
! XX - Vector for output variables
  REAL(kind=dp) :: XX(NVAR)

  XX(1) = X(1)/JVS(1)
  XX(2) = X(2)/JVS(4)
  XX(3) = X(3)/JVS(8)
  XX(4) = X(4)/JVS(9)
  XX(5) = X(5)/JVS(10)
  XX(6) = X(6)/JVS(11)
  XX(7) = X(7)/JVS(12)
  XX(8) = X(8)/JVS(28)
  XX(9) = X(9)/JVS(38)
  XX(10) = X(10)/JVS(53)
  XX(11) = X(11)/JVS(56)
  XX(12) = X(12)/JVS(61)
  XX(13) = X(13)/JVS(66)
  XX(14) = X(14)/JVS(75)
  XX(15) = X(15)/JVS(82)
  XX(16) = X(16)/JVS(88)
  XX(17) = X(17)/JVS(94)
  XX(18) = X(18)/JVS(98)
  XX(19) = X(19)/JVS(101)
  XX(20) = X(20)/JVS(105)
  XX(21) = X(21)/JVS(158)
  XX(22) = (X(22)-JVS(159)*XX(21))/(JVS(168))
  XX(23) = (X(23)-JVS(160)*XX(21)-JVS(169)*XX(22))/(JVS(172))
  XX(24) = (X(24)-JVS(161)*XX(21))/(JVS(174))
  XX(25) = (X(25)-JVS(162)*XX(21)-JVS(175)*XX(24))/(JVS(178))
  XX(26) = X(26)/JVS(180)
  XX(27) = (X(27)-JVS(163)*XX(21)-JVS(170)*XX(22))/(JVS(182))
  XX(28) = (X(28)-JVS(164)*XX(21))/(JVS(185))
  XX(29) = (X(29)-JVS(165)*XX(21)-JVS(176)*XX(24))/(JVS(187))
  XX(30) = (X(30)-JVS(166)*XX(21))/(JVS(190))
  XX(31) = (X(31)-JVS(106)*XX(20))/(JVS(192))
  XX(32) = (X(32)-JVS(2)*XX(1)-JVS(107)*XX(20))/(JVS(194))
  XX(33) = (X(33)-JVS(108)*XX(20))/(JVS(196))
  XX(34) = X(34)/JVS(198)
  XX(35) = X(35)/JVS(201)
  XX(36) = X(36)/JVS(204)
  XX(37) = X(37)/JVS(207)
  XX(38) = X(38)/JVS(210)
  XX(39) = (X(39)-JVS(109)*XX(20))/(JVS(213))
  XX(40) = X(40)/JVS(215)
  XX(41) = (X(41)-JVS(110)*XX(20))/(JVS(221))
  XX(42) = X(42)/JVS(223)
  XX(43) = (X(43)-JVS(111)*XX(20))/(JVS(226))
  XX(44) = (X(44)-JVS(112)*XX(20))/(JVS(229))
  XX(45) = X(45)/JVS(232)
  XX(46) = (X(46)-JVS(89)*XX(16)-JVS(113)*XX(20))/(JVS(235))
  XX(47) = (X(47)-JVS(114)*XX(20))/(JVS(237))
  XX(48) = (X(48)-JVS(13)*XX(7))/(JVS(241))
  XX(49) = (X(49)-JVS(62)*XX(12)-JVS(67)*XX(13))/(JVS(245))
  XX(50) = (X(50)-JVS(95)*XX(17)-JVS(115)*XX(20)-JVS(216)*XX(40))/(JVS(250))
  XX(51) = X(51)/JVS(252)
  XX(52) = (X(52)-JVS(116)*XX(20))/(JVS(256))
  XX(53) = (X(53)-JVS(90)*XX(16)-JVS(117)*XX(20))/(JVS(261))
  XX(54) = (X(54)-JVS(96)*XX(17)-JVS(118)*XX(20))/(JVS(263))
  XX(55) = (X(55)-JVS(119)*XX(20))/(JVS(267))
  XX(56) = (X(56)-JVS(120)*XX(20))/(JVS(271))
  XX(57) = (X(57)-JVS(121)*XX(20))/(JVS(275))
  XX(58) = (X(58)-JVS(122)*XX(20))/(JVS(280))
  XX(59) = (X(59)-JVS(5)*XX(2)-JVS(14)*XX(7)-JVS(123)*XX(20))/(JVS(283))
  XX(60) = (X(60)-JVS(124)*XX(20))/(JVS(286))
  XX(61) = (X(61)-JVS(246)*XX(49))/(JVS(292))
  XX(62) = (X(62)-JVS(68)*XX(13)-JVS(125)*XX(20))/(JVS(299))
  XX(63) = (X(63)-JVS(126)*XX(20))/(JVS(305))
  XX(64) = (X(64)-JVS(127)*XX(20))/(JVS(315))
  XX(65) = (X(65)-JVS(128)*XX(20))/(JVS(328))
  XX(66) = (X(66)-JVS(129)*XX(20))/(JVS(353))
  XX(67) = (X(67)-JVS(15)*XX(7)-JVS(130)*XX(20)-JVS(329)*XX(65))/(JVS(368))
  XX(68) = (X(68)-JVS(131)*XX(20))/(JVS(378))
  XX(69) = (X(69)-JVS(6)*XX(2)-JVS(16)*XX(7)-JVS(29)*XX(8)-JVS(63)*XX(12)-JVS(69)*XX(13)-JVS(132)*XX(20)-JVS(330)&
             &*XX(65))/(JVS(388))
  XX(70) = (X(70)-JVS(293)*XX(61)-JVS(306)*XX(63)-JVS(354)*XX(66))/(JVS(395))
  XX(71) = (X(71)-JVS(17)*XX(7)-JVS(39)*XX(9)-JVS(99)*XX(18)-JVS(133)*XX(20)-JVS(331)*XX(65))/(JVS(405))
  XX(72) = (X(72)-JVS(18)*XX(7)-JVS(30)*XX(8)-JVS(40)*XX(9)-JVS(91)*XX(16)-JVS(134)*XX(20)-JVS(332)*XX(65)-JVS(379)&
             &*XX(68))/(JVS(410))
  XX(73) = (X(73)-JVS(76)*XX(14)-JVS(83)*XX(15))/(JVS(423))
  XX(74) = (X(74)-JVS(19)*XX(7)-JVS(135)*XX(20)-JVS(316)*XX(64)-JVS(333)*XX(65)-JVS(355)*XX(66))/(JVS(446))
  XX(75) = (X(75)-JVS(20)*XX(7)-JVS(41)*XX(9)-JVS(102)*XX(19)-JVS(136)*XX(20)-JVS(217)*XX(40)-JVS(334)*XX(65)-JVS(380)&
             &*XX(68)-JVS(424)*XX(73))/(JVS(452))
  XX(76) = (X(76)-JVS(21)*XX(7)-JVS(42)*XX(9)-JVS(103)*XX(19)-JVS(137)*XX(20)-JVS(218)*XX(40)-JVS(335)*XX(65)-JVS(381)&
             &*XX(68)-JVS(425)*XX(73))/(JVS(457))
  XX(77) = (X(77)-JVS(138)*XX(20)-JVS(336)*XX(65)-JVS(356)*XX(66))/(JVS(473))
  XX(78) = (X(78)-JVS(22)*XX(7)-JVS(43)*XX(9)-JVS(139)*XX(20)-JVS(317)*XX(64)-JVS(337)*XX(65)-JVS(357)*XX(66)-JVS(474)&
             &*XX(77))/(JVS(487))
  XX(79) = (X(79)-JVS(140)*XX(20)-JVS(338)*XX(65)-JVS(358)*XX(66))/(JVS(502))
  XX(80) = (X(80)-JVS(23)*XX(7)-JVS(31)*XX(8)-JVS(44)*XX(9)-JVS(92)*XX(16)-JVS(141)*XX(20)-JVS(300)*XX(62)-JVS(339)&
             &*XX(65)-JVS(359)*XX(66)-JVS(382)*XX(68)-JVS(426)*XX(73)-JVS(447)*XX(74)-JVS(488)*XX(78)-JVS(503)*XX(79))&
             &/(JVS(517))
  XX(81) = (X(81)-JVS(24)*XX(7)-JVS(142)*XX(20)-JVS(318)*XX(64)-JVS(340)*XX(65)-JVS(360)*XX(66)-JVS(427)*XX(73))&
             &/(JVS(524))
  XX(82) = (X(82)-JVS(143)*XX(20)-JVS(242)*XX(48)-JVS(341)*XX(65)-JVS(361)*XX(66))/(JVS(549))
  XX(83) = (X(83)-JVS(144)*XX(20)-JVS(383)*XX(68)-JVS(428)*XX(73)-JVS(504)*XX(79)-JVS(550)*XX(82))/(JVS(573))
  XX(84) = (X(84)-JVS(70)*XX(13)-JVS(342)*XX(65)-JVS(369)*XX(67)-JVS(389)*XX(69)-JVS(406)*XX(71)-JVS(411)*XX(72)&
             &-JVS(429)*XX(73)-JVS(448)*XX(74)-JVS(453)*XX(75)-JVS(458)*XX(76)-JVS(475)*XX(77)-JVS(489)*XX(78)-JVS(505)&
             &*XX(79)-JVS(518)*XX(80)-JVS(525)*XX(81)-JVS(551)*XX(82)-JVS(574)*XX(83))/(JVS(591))
  XX(85) = (X(85)-JVS(145)*XX(20)-JVS(343)*XX(65)-JVS(362)*XX(66)-JVS(506)*XX(79))/(JVS(616))
  XX(86) = (X(86)-JVS(146)*XX(20)-JVS(430)*XX(73)-JVS(507)*XX(79)-JVS(552)*XX(82)-JVS(617)*XX(85))/(JVS(638))
  XX(87) = (X(87)-JVS(147)*XX(20)-JVS(431)*XX(73)-JVS(508)*XX(79)-JVS(553)*XX(82)-JVS(618)*XX(85)-JVS(639)*XX(86))&
             &/(JVS(660))
  XX(88) = (X(88)-JVS(7)*XX(2)-JVS(25)*XX(7)-JVS(32)*XX(8)-JVS(45)*XX(9)-JVS(71)*XX(13)-JVS(148)*XX(20)-JVS(181)*XX(26)&
             &-JVS(219)*XX(40)-JVS(281)*XX(58)-JVS(284)*XX(59)-JVS(301)*XX(62)-JVS(319)*XX(64)-JVS(344)*XX(65)-JVS(363)&
             &*XX(66)-JVS(370)*XX(67)-JVS(384)*XX(68)-JVS(390)*XX(69)-JVS(407)*XX(71)-JVS(412)*XX(72)-JVS(432)*XX(73)&
             &-JVS(449)*XX(74)-JVS(454)*XX(75)-JVS(459)*XX(76)-JVS(476)*XX(77)-JVS(490)*XX(78)-JVS(509)*XX(79)-JVS(519)&
             &*XX(80)-JVS(526)*XX(81)-JVS(554)*XX(82)-JVS(575)*XX(83)-JVS(592)*XX(84)-JVS(619)*XX(85)-JVS(640)*XX(86)&
             &-JVS(661)*XX(87))/(JVS(685))
  XX(89) = (X(89)-JVS(33)*XX(8)-JVS(54)*XX(10)-JVS(199)*XX(34)-JVS(396)*XX(70)-JVS(433)*XX(73)-JVS(477)*XX(77)-JVS(510)&
             &*XX(79)-JVS(555)*XX(82)-JVS(662)*XX(87)-JVS(686)*XX(88))/(JVS(717))
  XX(90) = (X(90)-JVS(34)*XX(8)-JVS(46)*XX(9)-JVS(77)*XX(14)-JVS(84)*XX(15)-JVS(257)*XX(52)-JVS(287)*XX(60)-JVS(434)&
             &*XX(73)-JVS(556)*XX(82)-JVS(620)*XX(85)-JVS(641)*XX(86)-JVS(663)*XX(87)-JVS(718)*XX(89))/(JVS(762))
  XX(91) = (X(91)-JVS(35)*XX(8)-JVS(47)*XX(9)-JVS(78)*XX(14)-JVS(85)*XX(15)-JVS(258)*XX(52)-JVS(288)*XX(60)-JVS(435)&
             &*XX(73)-JVS(557)*XX(82)-JVS(576)*XX(83)-JVS(621)*XX(85)-JVS(642)*XX(86)-JVS(664)*XX(87)-JVS(719)*XX(89)&
             &-JVS(763)*XX(90))/(JVS(795))
  XX(92) = (X(92)-JVS(64)*XX(12)-JVS(72)*XX(13)-JVS(149)*XX(20)-JVS(200)*XX(34)-JVS(202)*XX(35)-JVS(205)*XX(36)-JVS(208)&
             &*XX(37)-JVS(224)*XX(42)-JVS(233)*XX(45)-JVS(247)*XX(49)-JVS(253)*XX(51)-JVS(294)*XX(61)-JVS(364)*XX(66)&
             &-JVS(385)*XX(68)-JVS(397)*XX(70)-JVS(436)*XX(73)-JVS(478)*XX(77)-JVS(558)*XX(82)-JVS(577)*XX(83)-JVS(593)&
             &*XX(84)-JVS(622)*XX(85)-JVS(643)*XX(86)-JVS(665)*XX(87)-JVS(687)*XX(88)-JVS(720)*XX(89)-JVS(764)*XX(90)&
             &-JVS(796)*XX(91))/(JVS(838))
  XX(93) = (X(93)-JVS(65)*XX(12)-JVS(73)*XX(13)-JVS(79)*XX(14)-JVS(150)*XX(20)-JVS(225)*XX(42)-JVS(248)*XX(49)-JVS(276)&
             &*XX(57)-JVS(295)*XX(61)-JVS(302)*XX(62)-JVS(307)*XX(63)-JVS(320)*XX(64)-JVS(345)*XX(65)-JVS(365)*XX(66)&
             &-JVS(371)*XX(67)-JVS(386)*XX(68)-JVS(391)*XX(69)-JVS(398)*XX(70)-JVS(408)*XX(71)-JVS(413)*XX(72)-JVS(437)&
             &*XX(73)-JVS(450)*XX(74)-JVS(455)*XX(75)-JVS(460)*XX(76)-JVS(479)*XX(77)-JVS(491)*XX(78)-JVS(511)*XX(79)&
             &-JVS(520)*XX(80)-JVS(527)*XX(81)-JVS(559)*XX(82)-JVS(578)*XX(83)-JVS(594)*XX(84)-JVS(623)*XX(85)-JVS(644)&
             &*XX(86)-JVS(666)*XX(87)-JVS(688)*XX(88)-JVS(721)*XX(89)-JVS(765)*XX(90)-JVS(797)*XX(91)-JVS(839)*XX(92))&
             &/(JVS(880))
  XX(94) = (X(94)-JVS(36)*XX(8)-JVS(48)*XX(9)-JVS(80)*XX(14)-JVS(86)*XX(15)-JVS(238)*XX(47)-JVS(259)*XX(52)-JVS(438)&
             &*XX(73)-JVS(560)*XX(82)-JVS(645)*XX(86)-JVS(667)*XX(87)-JVS(722)*XX(89)-JVS(766)*XX(90)-JVS(798)*XX(91)&
             &-JVS(840)*XX(92)-JVS(881)*XX(93))/(JVS(913))
  XX(95) = (X(95)-JVS(49)*XX(9)-JVS(57)*XX(11)-JVS(206)*XX(36)-JVS(399)*XX(70)-JVS(439)*XX(73)-JVS(480)*XX(77)-JVS(512)&
             &*XX(79)-JVS(561)*XX(82)-JVS(668)*XX(87)-JVS(689)*XX(88)-JVS(723)*XX(89)-JVS(767)*XX(90)-JVS(799)*XX(91)&
             &-JVS(841)*XX(92)-JVS(882)*XX(93)-JVS(914)*XX(94))/(JVS(931))
  XX(96) = (X(96)-JVS(50)*XX(9)-JVS(58)*XX(11)-JVS(209)*XX(37)-JVS(400)*XX(70)-JVS(440)*XX(73)-JVS(481)*XX(77)-JVS(513)&
             &*XX(79)-JVS(562)*XX(82)-JVS(690)*XX(88)-JVS(724)*XX(89)-JVS(768)*XX(90)-JVS(800)*XX(91)-JVS(842)*XX(92)&
             &-JVS(883)*XX(93)-JVS(915)*XX(94)-JVS(932)*XX(95))/(JVS(952))
  XX(97) = (X(97)-JVS(37)*XX(8)-JVS(51)*XX(9)-JVS(55)*XX(10)-JVS(59)*XX(11)-JVS(87)*XX(15)-JVS(151)*XX(20)-JVS(211)&
             &*XX(38)-JVS(239)*XX(47)-JVS(243)*XX(48)-JVS(249)*XX(49)-JVS(254)*XX(51)-JVS(289)*XX(60)-JVS(296)*XX(61)&
             &-JVS(308)*XX(63)-JVS(366)*XX(66)-JVS(401)*XX(70)-JVS(441)*XX(73)-JVS(482)*XX(77)-JVS(563)*XX(82)-JVS(624)&
             &*XX(85)-JVS(691)*XX(88)-JVS(725)*XX(89)-JVS(769)*XX(90)-JVS(801)*XX(91)-JVS(843)*XX(92)-JVS(884)*XX(93)&
             &-JVS(916)*XX(94)-JVS(933)*XX(95)-JVS(953)*XX(96))/(JVS(1004))
  XX(98) = (X(98)-JVS(3)*XX(1)-JVS(26)*XX(7)-JVS(74)*XX(13)-JVS(93)*XX(16)-JVS(97)*XX(17)-JVS(100)*XX(18)-JVS(104)&
             &*XX(19)-JVS(152)*XX(20)-JVS(167)*XX(21)-JVS(171)*XX(22)-JVS(173)*XX(23)-JVS(177)*XX(24)-JVS(179)*XX(25)&
             &-JVS(183)*XX(27)-JVS(186)*XX(28)-JVS(188)*XX(29)-JVS(191)*XX(30)-JVS(193)*XX(31)-JVS(195)*XX(32)-JVS(197)&
             &*XX(33)-JVS(212)*XX(38)-JVS(214)*XX(39)-JVS(220)*XX(40)-JVS(222)*XX(41)-JVS(227)*XX(43)-JVS(230)*XX(44)&
             &-JVS(234)*XX(45)-JVS(236)*XX(46)-JVS(240)*XX(47)-JVS(251)*XX(50)-JVS(255)*XX(51)-JVS(260)*XX(52)-JVS(262)&
             &*XX(53)-JVS(264)*XX(54)-JVS(268)*XX(55)-JVS(272)*XX(56)-JVS(277)*XX(57)-JVS(282)*XX(58)-JVS(285)*XX(59)&
             &-JVS(290)*XX(60)-JVS(303)*XX(62)-JVS(309)*XX(63)-JVS(321)*XX(64)-JVS(346)*XX(65)-JVS(367)*XX(66)-JVS(372)&
             &*XX(67)-JVS(387)*XX(68)-JVS(392)*XX(69)-JVS(402)*XX(70)-JVS(409)*XX(71)-JVS(414)*XX(72)-JVS(442)*XX(73)&
             &-JVS(451)*XX(74)-JVS(456)*XX(75)-JVS(461)*XX(76)-JVS(483)*XX(77)-JVS(492)*XX(78)-JVS(514)*XX(79)-JVS(521)&
             &*XX(80)-JVS(528)*XX(81)-JVS(564)*XX(82)-JVS(579)*XX(83)-JVS(595)*XX(84)-JVS(625)*XX(85)-JVS(646)*XX(86)&
             &-JVS(669)*XX(87)-JVS(692)*XX(88)-JVS(726)*XX(89)-JVS(770)*XX(90)-JVS(802)*XX(91)-JVS(844)*XX(92)-JVS(885)&
             &*XX(93)-JVS(917)*XX(94)-JVS(934)*XX(95)-JVS(954)*XX(96)-JVS(1005)*XX(97))/(JVS(1071))
  XX(99) = (X(99)-JVS(52)*XX(9)-JVS(60)*XX(11)-JVS(203)*XX(35)-JVS(403)*XX(70)-JVS(443)*XX(73)-JVS(484)*XX(77)-JVS(515)&
             &*XX(79)-JVS(565)*XX(82)-JVS(670)*XX(87)-JVS(693)*XX(88)-JVS(727)*XX(89)-JVS(771)*XX(90)-JVS(803)*XX(91)&
             &-JVS(845)*XX(92)-JVS(886)*XX(93)-JVS(918)*XX(94)-JVS(935)*XX(95)-JVS(955)*XX(96)-JVS(1006)*XX(97)-JVS(1072)&
             &*XX(98))/(JVS(1097))
  XX(100) = (X(100)-JVS(27)*XX(7)-JVS(81)*XX(14)-JVS(153)*XX(20)-JVS(228)*XX(43)-JVS(244)*XX(48)-JVS(404)*XX(70)&
              &-JVS(444)*XX(73)-JVS(485)*XX(77)-JVS(516)*XX(79)-JVS(566)*XX(82)-JVS(580)*XX(83)-JVS(596)*XX(84)-JVS(626)&
              &*XX(85)-JVS(647)*XX(86)-JVS(671)*XX(87)-JVS(694)*XX(88)-JVS(728)*XX(89)-JVS(772)*XX(90)-JVS(804)*XX(91)&
              &-JVS(846)*XX(92)-JVS(887)*XX(93)-JVS(919)*XX(94)-JVS(936)*XX(95)-JVS(956)*XX(96)-JVS(1007)*XX(97)-JVS(1073)&
              &*XX(98)-JVS(1098)*XX(99))/(JVS(1123))
  XX(100) = XX(100)
  XX(99) = XX(99)-JVS(1122)*XX(100)
  XX(98) = XX(98)-JVS(1096)*XX(99)-JVS(1121)*XX(100)
  XX(97) = XX(97)-JVS(1070)*XX(98)-JVS(1095)*XX(99)-JVS(1120)*XX(100)
  XX(96) = XX(96)-JVS(1003)*XX(97)-JVS(1069)*XX(98)-JVS(1094)*XX(99)-JVS(1119)*XX(100)
  XX(95) = XX(95)-JVS(951)*XX(96)-JVS(1002)*XX(97)-JVS(1068)*XX(98)-JVS(1093)*XX(99)-JVS(1118)*XX(100)
  XX(94) = XX(94)-JVS(930)*XX(95)-JVS(950)*XX(96)-JVS(1001)*XX(97)-JVS(1067)*XX(98)-JVS(1092)*XX(99)-JVS(1117)*XX(100)
  XX(93) = XX(93)-JVS(912)*XX(94)-JVS(929)*XX(95)-JVS(949)*XX(96)-JVS(1000)*XX(97)-JVS(1066)*XX(98)-JVS(1091)*XX(99)&
             &-JVS(1116)*XX(100)
  XX(92) = XX(92)-JVS(879)*XX(93)-JVS(911)*XX(94)-JVS(928)*XX(95)-JVS(948)*XX(96)-JVS(999)*XX(97)-JVS(1065)*XX(98)&
             &-JVS(1090)*XX(99)-JVS(1115)*XX(100)
  XX(91) = XX(91)-JVS(837)*XX(92)-JVS(878)*XX(93)-JVS(910)*XX(94)-JVS(927)*XX(95)-JVS(947)*XX(96)-JVS(998)*XX(97)&
             &-JVS(1064)*XX(98)-JVS(1089)*XX(99)-JVS(1114)*XX(100)
  XX(90) = XX(90)-JVS(794)*XX(91)-JVS(836)*XX(92)-JVS(877)*XX(93)-JVS(909)*XX(94)-JVS(926)*XX(95)-JVS(946)*XX(96)&
             &-JVS(997)*XX(97)-JVS(1063)*XX(98)-JVS(1088)*XX(99)-JVS(1113)*XX(100)
  XX(89) = XX(89)-JVS(761)*XX(90)-JVS(793)*XX(91)-JVS(835)*XX(92)-JVS(876)*XX(93)-JVS(908)*XX(94)-JVS(925)*XX(95)&
             &-JVS(945)*XX(96)-JVS(996)*XX(97)-JVS(1062)*XX(98)-JVS(1087)*XX(99)-JVS(1112)*XX(100)
  XX(88) = XX(88)-JVS(716)*XX(89)-JVS(760)*XX(90)-JVS(792)*XX(91)-JVS(834)*XX(92)-JVS(875)*XX(93)-JVS(907)*XX(94)&
             &-JVS(924)*XX(95)-JVS(944)*XX(96)-JVS(995)*XX(97)-JVS(1061)*XX(98)-JVS(1086)*XX(99)-JVS(1111)*XX(100)
  XX(87) = XX(87)-JVS(715)*XX(89)-JVS(759)*XX(90)-JVS(791)*XX(91)-JVS(833)*XX(92)-JVS(874)*XX(93)-JVS(906)*XX(94)&
             &-JVS(994)*XX(97)-JVS(1060)*XX(98)-JVS(1085)*XX(99)-JVS(1110)*XX(100)
  XX(86) = XX(86)-JVS(714)*XX(89)-JVS(758)*XX(90)-JVS(790)*XX(91)-JVS(832)*XX(92)-JVS(873)*XX(93)-JVS(905)*XX(94)&
             &-JVS(993)*XX(97)-JVS(1059)*XX(98)-JVS(1084)*XX(99)-JVS(1109)*XX(100)
  XX(85) = XX(85)-JVS(713)*XX(89)-JVS(757)*XX(90)-JVS(789)*XX(91)-JVS(831)*XX(92)-JVS(872)*XX(93)-JVS(904)*XX(94)&
             &-JVS(992)*XX(97)-JVS(1058)*XX(98)-JVS(1083)*XX(99)
  XX(84) = XX(84)-JVS(615)*XX(85)-JVS(637)*XX(86)-JVS(659)*XX(87)-JVS(684)*XX(88)-JVS(712)*XX(89)-JVS(756)*XX(90)&
             &-JVS(788)*XX(91)-JVS(830)*XX(92)-JVS(871)*XX(93)-JVS(903)*XX(94)-JVS(923)*XX(95)-JVS(943)*XX(96)-JVS(991)&
             &*XX(97)-JVS(1057)*XX(98)-JVS(1082)*XX(99)-JVS(1108)*XX(100)
  XX(83) = XX(83)-JVS(614)*XX(85)-JVS(636)*XX(86)-JVS(658)*XX(87)-JVS(711)*XX(89)-JVS(755)*XX(90)-JVS(787)*XX(91)&
             &-JVS(829)*XX(92)-JVS(870)*XX(93)-JVS(902)*XX(94)-JVS(990)*XX(97)-JVS(1056)*XX(98)-JVS(1107)*XX(100)
  XX(82) = XX(82)-JVS(828)*XX(92)-JVS(869)*XX(93)-JVS(989)*XX(97)-JVS(1055)*XX(98)-JVS(1106)*XX(100)
  XX(81) = XX(81)-JVS(548)*XX(82)-JVS(590)*XX(84)-JVS(613)*XX(85)-JVS(635)*XX(86)-JVS(657)*XX(87)-JVS(683)*XX(88)&
             &-JVS(710)*XX(89)-JVS(754)*XX(90)-JVS(786)*XX(91)-JVS(827)*XX(92)-JVS(868)*XX(93)-JVS(901)*XX(94)-JVS(942)&
             &*XX(96)-JVS(988)*XX(97)-JVS(1054)*XX(98)-JVS(1081)*XX(99)-JVS(1105)*XX(100)
  XX(80) = XX(80)-JVS(523)*XX(81)-JVS(547)*XX(82)-JVS(572)*XX(83)-JVS(589)*XX(84)-JVS(612)*XX(85)-JVS(634)*XX(86)&
             &-JVS(656)*XX(87)-JVS(682)*XX(88)-JVS(709)*XX(89)-JVS(753)*XX(90)-JVS(785)*XX(91)-JVS(826)*XX(92)-JVS(867)&
             &*XX(93)-JVS(900)*XX(94)-JVS(922)*XX(95)-JVS(941)*XX(96)-JVS(987)*XX(97)-JVS(1053)*XX(98)-JVS(1080)*XX(99)&
             &-JVS(1104)*XX(100)
  XX(79) = XX(79)-JVS(708)*XX(89)-JVS(825)*XX(92)-JVS(866)*XX(93)-JVS(899)*XX(94)-JVS(986)*XX(97)-JVS(1052)*XX(98)
  XX(78) = XX(78)-JVS(501)*XX(79)-JVS(546)*XX(82)-JVS(571)*XX(83)-JVS(611)*XX(85)-JVS(633)*XX(86)-JVS(655)*XX(87)&
             &-JVS(681)*XX(88)-JVS(707)*XX(89)-JVS(752)*XX(90)-JVS(784)*XX(91)-JVS(824)*XX(92)-JVS(865)*XX(93)-JVS(940)&
             &*XX(96)-JVS(985)*XX(97)-JVS(1051)*XX(98)-JVS(1079)*XX(99)
  XX(77) = XX(77)-JVS(545)*XX(82)-JVS(823)*XX(92)-JVS(864)*XX(93)-JVS(984)*XX(97)-JVS(1050)*XX(98)-JVS(1078)*XX(99)
  XX(76) = XX(76)-JVS(472)*XX(77)-JVS(544)*XX(82)-JVS(570)*XX(83)-JVS(588)*XX(84)-JVS(610)*XX(85)-JVS(654)*XX(87)&
             &-JVS(680)*XX(88)-JVS(706)*XX(89)-JVS(751)*XX(90)-JVS(783)*XX(91)-JVS(822)*XX(92)-JVS(863)*XX(93)-JVS(898)&
             &*XX(94)-JVS(983)*XX(97)-JVS(1049)*XX(98)-JVS(1077)*XX(99)-JVS(1103)*XX(100)
  XX(75) = XX(75)-JVS(471)*XX(77)-JVS(543)*XX(82)-JVS(569)*XX(83)-JVS(587)*XX(84)-JVS(609)*XX(85)-JVS(653)*XX(87)&
             &-JVS(679)*XX(88)-JVS(705)*XX(89)-JVS(750)*XX(90)-JVS(782)*XX(91)-JVS(821)*XX(92)-JVS(862)*XX(93)-JVS(897)&
             &*XX(94)-JVS(982)*XX(97)-JVS(1048)*XX(98)-JVS(1076)*XX(99)-JVS(1102)*XX(100)
  XX(74) = XX(74)-JVS(542)*XX(82)-JVS(586)*XX(84)-JVS(608)*XX(85)-JVS(632)*XX(86)-JVS(678)*XX(88)-JVS(704)*XX(89)&
             &-JVS(749)*XX(90)-JVS(820)*XX(92)-JVS(861)*XX(93)-JVS(939)*XX(96)-JVS(981)*XX(97)-JVS(1047)*XX(98)-JVS(1075)&
             &*XX(99)
  XX(73) = XX(73)-JVS(819)*XX(92)-JVS(860)*XX(93)-JVS(1101)*XX(100)
  XX(72) = XX(72)-JVS(422)*XX(73)-JVS(500)*XX(79)-JVS(541)*XX(82)-JVS(568)*XX(83)-JVS(585)*XX(84)-JVS(607)*XX(85)&
             &-JVS(631)*XX(86)-JVS(652)*XX(87)-JVS(677)*XX(88)-JVS(703)*XX(89)-JVS(748)*XX(90)-JVS(781)*XX(91)-JVS(859)&
             &*XX(93)-JVS(896)*XX(94)-JVS(980)*XX(97)-JVS(1046)*XX(98)
  XX(71) = XX(71)-JVS(421)*XX(73)-JVS(445)*XX(74)-JVS(486)*XX(78)-JVS(522)*XX(81)-JVS(540)*XX(82)-JVS(584)*XX(84)&
             &-JVS(651)*XX(87)-JVS(676)*XX(88)-JVS(747)*XX(90)-JVS(780)*XX(91)-JVS(818)*XX(92)-JVS(858)*XX(93)-JVS(895)&
             &*XX(94)-JVS(938)*XX(96)-JVS(979)*XX(97)-JVS(1045)*XX(98)
  XX(70) = XX(70)-JVS(470)*XX(77)-JVS(746)*XX(90)-JVS(817)*XX(92)-JVS(857)*XX(93)-JVS(978)*XX(97)-JVS(1044)*XX(98)
  XX(69) = XX(69)-JVS(499)*XX(79)-JVS(539)*XX(82)-JVS(583)*XX(84)-JVS(606)*XX(85)-JVS(630)*XX(86)-JVS(675)*XX(88)&
             &-JVS(745)*XX(90)-JVS(779)*XX(91)-JVS(856)*XX(93)-JVS(894)*XX(94)-JVS(977)*XX(97)-JVS(1043)*XX(98)
  XX(68) = XX(68)-JVS(420)*XX(73)-JVS(538)*XX(82)-JVS(702)*XX(89)-JVS(893)*XX(94)-JVS(1042)*XX(98)
  XX(67) = XX(67)-JVS(469)*XX(77)-JVS(498)*XX(79)-JVS(537)*XX(82)-JVS(582)*XX(84)-JVS(605)*XX(85)-JVS(674)*XX(88)&
             &-JVS(744)*XX(90)-JVS(855)*XX(93)-JVS(892)*XX(94)-JVS(976)*XX(97)-JVS(1041)*XX(98)
  XX(66) = XX(66)-JVS(816)*XX(92)-JVS(854)*XX(93)-JVS(1040)*XX(98)
  XX(65) = XX(65)-JVS(975)*XX(97)-JVS(1039)*XX(98)
  XX(64) = XX(64)-JVS(327)*XX(65)-JVS(352)*XX(66)-JVS(701)*XX(89)-JVS(853)*XX(93)-JVS(974)*XX(97)-JVS(1038)*XX(98)
  XX(63) = XX(63)-JVS(351)*XX(66)-JVS(394)*XX(70)-JVS(468)*XX(77)-JVS(743)*XX(90)-JVS(852)*XX(93)-JVS(1037)*XX(98)
  XX(62) = XX(62)-JVS(350)*XX(66)-JVS(851)*XX(93)-JVS(921)*XX(95)-JVS(1036)*XX(98)
  XX(61) = XX(61)-JVS(349)*XX(66)-JVS(815)*XX(92)-JVS(850)*XX(93)-JVS(973)*XX(97)
  XX(60) = XX(60)-JVS(604)*XX(85)-JVS(742)*XX(90)-JVS(972)*XX(97)-JVS(1035)*XX(98)
  XX(59) = XX(59)-JVS(326)*XX(65)-JVS(467)*XX(77)-JVS(536)*XX(82)-JVS(673)*XX(88)-JVS(741)*XX(90)-JVS(971)*XX(97)&
             &-JVS(1034)*XX(98)
  XX(58) = XX(58)-JVS(325)*XX(65)-JVS(466)*XX(77)-JVS(603)*XX(85)-JVS(672)*XX(88)-JVS(740)*XX(90)-JVS(970)*XX(97)&
             &-JVS(1033)*XX(98)
  XX(57) = XX(57)-JVS(314)*XX(64)-JVS(348)*XX(66)-JVS(393)*XX(70)-JVS(739)*XX(90)-JVS(849)*XX(93)-JVS(1032)*XX(98)
  XX(56) = XX(56)-JVS(313)*XX(64)-JVS(324)*XX(65)-JVS(419)*XX(73)-JVS(465)*XX(77)-JVS(602)*XX(85)-JVS(700)*XX(89)&
             &-JVS(738)*XX(90)-JVS(969)*XX(97)-JVS(1031)*XX(98)
  XX(55) = XX(55)-JVS(312)*XX(64)-JVS(323)*XX(65)-JVS(418)*XX(73)-JVS(464)*XX(77)-JVS(601)*XX(85)-JVS(699)*XX(89)&
             &-JVS(737)*XX(90)-JVS(968)*XX(97)-JVS(1030)*XX(98)
  XX(54) = XX(54)-JVS(266)*XX(55)-JVS(270)*XX(56)-JVS(274)*XX(57)-JVS(279)*XX(58)-JVS(298)*XX(62)-JVS(304)*XX(63)&
             &-JVS(311)*XX(64)-JVS(463)*XX(77)-JVS(650)*XX(87)-JVS(736)*XX(90)-JVS(778)*XX(91)-JVS(967)*XX(97)-JVS(1029)&
             &*XX(98)
  XX(53) = XX(53)-JVS(322)*XX(65)-JVS(377)*XX(68)-JVS(417)*XX(73)-JVS(497)*XX(79)-JVS(535)*XX(82)-JVS(600)*XX(85)&
             &-JVS(629)*XX(86)-JVS(649)*XX(87)-JVS(698)*XX(89)-JVS(735)*XX(90)-JVS(777)*XX(91)-JVS(891)*XX(94)-JVS(1028)&
             &*XX(98)
  XX(52) = XX(52)-JVS(534)*XX(82)-JVS(966)*XX(97)-JVS(1027)*XX(98)
  XX(51) = XX(51)-JVS(814)*XX(92)-JVS(848)*XX(93)-JVS(965)*XX(97)-JVS(1026)*XX(98)
  XX(50) = XX(50)-JVS(265)*XX(55)-JVS(269)*XX(56)-JVS(273)*XX(57)-JVS(278)*XX(58)-JVS(297)*XX(62)-JVS(310)*XX(64)&
             &-JVS(462)*XX(77)-JVS(697)*XX(89)-JVS(734)*XX(90)-JVS(776)*XX(91)-JVS(964)*XX(97)-JVS(1025)*XX(98)
  XX(49) = XX(49)-JVS(291)*XX(61)-JVS(813)*XX(92)-JVS(963)*XX(97)
  XX(48) = XX(48)-JVS(533)*XX(82)-JVS(812)*XX(92)-JVS(962)*XX(97)-JVS(1100)*XX(100)
  XX(47) = XX(47)-JVS(532)*XX(82)-JVS(890)*XX(94)-JVS(961)*XX(97)-JVS(1024)*XX(98)
  XX(46) = XX(46)-JVS(376)*XX(68)-JVS(416)*XX(73)-JVS(496)*XX(79)-JVS(531)*XX(82)-JVS(599)*XX(85)-JVS(628)*XX(86)&
             &-JVS(648)*XX(87)-JVS(733)*XX(90)-JVS(775)*XX(91)-JVS(1023)*XX(98)
  XX(45) = XX(45)-JVS(375)*XX(68)-JVS(567)*XX(83)-JVS(811)*XX(92)-JVS(889)*XX(94)
  XX(44) = XX(44)-JVS(231)*XX(45)-JVS(374)*XX(68)-JVS(415)*XX(73)-JVS(495)*XX(79)-JVS(530)*XX(82)-JVS(598)*XX(85)&
             &-JVS(627)*XX(86)-JVS(732)*XX(90)-JVS(774)*XX(91)-JVS(1022)*XX(98)
  XX(43) = XX(43)-JVS(810)*XX(92)-JVS(960)*XX(97)-JVS(1021)*XX(98)-JVS(1099)*XX(100)
  XX(42) = XX(42)-JVS(347)*XX(66)-JVS(809)*XX(92)-JVS(847)*XX(93)
  XX(41) = XX(41)-JVS(494)*XX(79)-JVS(529)*XX(82)-JVS(731)*XX(90)-JVS(959)*XX(97)-JVS(1020)*XX(98)
  XX(40) = XX(40)-JVS(696)*XX(89)
  XX(39) = XX(39)-JVS(373)*XX(68)-JVS(597)*XX(85)-JVS(730)*XX(90)-JVS(773)*XX(91)-JVS(1019)*XX(98)
  XX(38) = XX(38)-JVS(958)*XX(97)-JVS(1018)*XX(98)
  XX(37) = XX(37)-JVS(808)*XX(92)-JVS(937)*XX(96)
  XX(36) = XX(36)-JVS(807)*XX(92)-JVS(920)*XX(95)
  XX(35) = XX(35)-JVS(806)*XX(92)-JVS(1074)*XX(99)
  XX(34) = XX(34)-JVS(695)*XX(89)-JVS(805)*XX(92)
  XX(33) = XX(33)-JVS(493)*XX(79)-JVS(729)*XX(90)-JVS(1017)*XX(98)
  XX(32) = XX(32)-JVS(957)*XX(97)-JVS(1016)*XX(98)
  XX(31) = XX(31)-JVS(888)*XX(94)-JVS(1015)*XX(98)
  XX(30) = XX(30)-JVS(1014)*XX(98)
  XX(29) = XX(29)-JVS(189)*XX(30)-JVS(1013)*XX(98)
  XX(28) = XX(28)-JVS(1012)*XX(98)
  XX(27) = XX(27)-JVS(184)*XX(28)-JVS(1011)*XX(98)
  XX(26) = XX(26)-JVS(581)*XX(84)-JVS(1010)*XX(98)
  XX(25) = XX(25)
  XX(24) = XX(24)
  XX(23) = XX(23)
  XX(22) = XX(22)
  XX(21) = XX(21)
  XX(20) = XX(20)
  XX(19) = XX(19)
  XX(18) = XX(18)
  XX(17) = XX(17)
  XX(16) = XX(16)
  XX(15) = XX(15)
  XX(14) = XX(14)
  XX(13) = XX(13)
  XX(12) = XX(12)
  XX(11) = XX(11)
  XX(10) = XX(10)
  XX(9) = XX(9)
  XX(8) = XX(8)
  XX(7) = XX(7)
  XX(6) = XX(6)-JVS(157)*XX(21)
  XX(5) = XX(5)-JVS(156)*XX(21)-JVS(1009)*XX(98)
  XX(4) = XX(4)-JVS(155)*XX(21)
  XX(3) = XX(3)-JVS(154)*XX(21)-JVS(1008)*XX(98)
  XX(2) = XX(2)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE saprc99_mosaic_4bin_vbs2_LinearAlgebra

