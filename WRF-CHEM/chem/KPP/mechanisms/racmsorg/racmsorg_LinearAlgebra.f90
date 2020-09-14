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
! File                 : racmsorg_LinearAlgebra.f90
! Time                 : Mon Feb 27 11:08:56 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/racmsorg
! Equation file        : racmsorg.kpp
! Output root filename : racmsorg
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racmsorg_LinearAlgebra

  USE racmsorg_Parameters
  USE racmsorg_JacobianSP

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
  XX(3) = X(3)/JVS(7)
  XX(4) = X(4)/JVS(19)
  XX(5) = (X(5)-JVS(2)*XX(1))/(JVS(41))
  XX(6) = X(6)/JVS(43)
  XX(7) = X(7)/JVS(45)
  XX(8) = X(8)/JVS(47)
  XX(9) = X(9)/JVS(49)
  XX(10) = X(10)/JVS(51)
  XX(11) = X(11)/JVS(54)
  XX(12) = (X(12)-JVS(8)*XX(3))/(JVS(56))
  XX(13) = X(13)/JVS(58)
  XX(14) = X(14)/JVS(63)
  XX(15) = X(15)/JVS(68)
  XX(16) = X(16)/JVS(71)
  XX(17) = X(17)/JVS(75)
  XX(18) = X(18)/JVS(79)
  XX(19) = X(19)/JVS(86)
  XX(20) = X(20)/JVS(96)
  XX(21) = (X(21)-JVS(80)*XX(18))/(JVS(103))
  XX(22) = (X(22)-JVS(81)*XX(18))/(JVS(108))
  XX(23) = X(23)/JVS(113)
  XX(24) = (X(24)-JVS(9)*XX(3))/(JVS(119))
  XX(25) = (X(25)-JVS(82)*XX(18))/(JVS(123))
  XX(26) = X(26)/JVS(128)
  XX(27) = X(27)/JVS(134)
  XX(28) = (X(28)-JVS(5)*XX(2))/(JVS(148))
  XX(29) = (X(29)-JVS(87)*XX(19)-JVS(149)*XX(28))/(JVS(166))
  XX(30) = (X(30)-JVS(10)*XX(3)-JVS(20)*XX(4)-JVS(88)*XX(19)-JVS(150)*XX(28))/(JVS(170))
  XX(31) = X(31)/JVS(174)
  XX(32) = (X(32)-JVS(97)*XX(20)-JVS(124)*XX(25)-JVS(135)*XX(27))/(JVS(185))
  XX(33) = (X(33)-JVS(11)*XX(3)-JVS(89)*XX(19)-JVS(151)*XX(28))/(JVS(191))
  XX(34) = (X(34)-JVS(136)*XX(27)-JVS(152)*XX(28))/(JVS(196))
  XX(35) = (X(35)-JVS(12)*XX(3)-JVS(114)*XX(23)-JVS(153)*XX(28)-JVS(175)*XX(31))/(JVS(208))
  XX(36) = (X(36)-JVS(21)*XX(4))/(JVS(215))
  XX(37) = (X(37)-JVS(13)*XX(3)-JVS(90)*XX(19)-JVS(154)*XX(28))/(JVS(223))
  XX(38) = (X(38)-JVS(22)*XX(4))/(JVS(228))
  XX(39) = (X(39)-JVS(23)*XX(4))/(JVS(236))
  XX(40) = (X(40)-JVS(137)*XX(27)-JVS(155)*XX(28))/(JVS(246))
  XX(41) = (X(41)-JVS(197)*XX(34)-JVS(247)*XX(40))/(JVS(264))
  XX(42) = X(42)/JVS(278)
  XX(43) = X(43)/JVS(299)
  XX(44) = (X(44)-JVS(24)*XX(4)-JVS(279)*XX(42))/(JVS(308))
  XX(45) = (X(45)-JVS(25)*XX(4)-JVS(280)*XX(42))/(JVS(316))
  XX(46) = (X(46)-JVS(198)*XX(34)-JVS(248)*XX(40))/(JVS(324))
  XX(47) = (X(47)-JVS(199)*XX(34)-JVS(249)*XX(40))/(JVS(334))
  XX(48) = (X(48)-JVS(138)*XX(27)-JVS(156)*XX(28))/(JVS(363))
  XX(49) = (X(49)-JVS(281)*XX(42)-JVS(364)*XX(48))/(JVS(387))
  XX(50) = (X(50)-JVS(26)*XX(4)-JVS(365)*XX(48))/(JVS(397))
  XX(51) = (X(51)-JVS(14)*XX(3)-JVS(27)*XX(4)-JVS(115)*XX(23)-JVS(139)*XX(27)-JVS(157)*XX(28)-JVS(250)*XX(40)-JVS(366)&
             &*XX(48))/(JVS(411))
  XX(52) = (X(52)-JVS(28)*XX(4)-JVS(200)*XX(34)-JVS(251)*XX(40)-JVS(282)*XX(42)-JVS(367)*XX(48))/(JVS(421))
  XX(53) = (X(53)-JVS(140)*XX(27)-JVS(158)*XX(28))/(JVS(442))
  XX(54) = (X(54)-JVS(15)*XX(3)-JVS(29)*XX(4)-JVS(69)*XX(15)-JVS(129)*XX(26)-JVS(141)*XX(27)-JVS(159)*XX(28)-JVS(201)&
             &*XX(34)-JVS(252)*XX(40)-JVS(283)*XX(42)-JVS(368)*XX(48)-JVS(443)*XX(53))/(JVS(463))
  XX(55) = (X(55)-JVS(30)*XX(4)-JVS(209)*XX(35)-JVS(253)*XX(40)-JVS(369)*XX(48))/(JVS(476))
  XX(56) = (X(56)-JVS(370)*XX(48))/(JVS(500))
  XX(57) = (X(57)-JVS(16)*XX(3)-JVS(31)*XX(4)-JVS(59)*XX(13)-JVS(64)*XX(14)-JVS(91)*XX(19)-JVS(160)*XX(28)-JVS(229)&
             &*XX(38)-JVS(284)*XX(42)-JVS(371)*XX(48)-JVS(444)*XX(53))/(JVS(517))
  XX(58) = (X(58)-JVS(32)*XX(4)-JVS(60)*XX(13)-JVS(65)*XX(14)-JVS(92)*XX(19)-JVS(161)*XX(28)-JVS(237)*XX(39)-JVS(285)&
             &*XX(42)-JVS(372)*XX(48)-JVS(445)*XX(53))/(JVS(528))
  XX(59) = (X(59)-JVS(33)*XX(4)-JVS(286)*XX(42)-JVS(373)*XX(48)-JVS(446)*XX(53))/(JVS(545))
  XX(60) = (X(60)-JVS(34)*XX(4)-JVS(287)*XX(42)-JVS(374)*XX(48)-JVS(447)*XX(53)-JVS(546)*XX(59))/(JVS(563))
  XX(61) = (X(61)-JVS(35)*XX(4)-JVS(375)*XX(48)-JVS(448)*XX(53))/(JVS(585))
  XX(62) = (X(62)-JVS(162)*XX(28)-JVS(224)*XX(37)-JVS(376)*XX(48)-JVS(398)*XX(50)-JVS(412)*XX(51)-JVS(449)*XX(53)&
             &-JVS(464)*XX(54)-JVS(477)*XX(55)-JVS(501)*XX(56)-JVS(518)*XX(57)-JVS(529)*XX(58)-JVS(547)*XX(59)-JVS(564)&
             &*XX(60)-JVS(586)*XX(61))/(JVS(602))
  XX(63) = (X(63)-JVS(17)*XX(3)-JVS(36)*XX(4)-JVS(44)*XX(6)-JVS(61)*XX(13)-JVS(66)*XX(14)-JVS(93)*XX(19)-JVS(104)*XX(21)&
             &-JVS(109)*XX(22)-JVS(120)*XX(24)-JVS(125)*XX(25)-JVS(130)*XX(26)-JVS(163)*XX(28)-JVS(167)*XX(29)-JVS(171)&
             &*XX(30)-JVS(176)*XX(31)-JVS(186)*XX(32)-JVS(192)*XX(33)-JVS(202)*XX(34)-JVS(210)*XX(35)-JVS(216)*XX(36)&
             &-JVS(225)*XX(37)-JVS(254)*XX(40)-JVS(265)*XX(41)-JVS(288)*XX(42)-JVS(300)*XX(43)-JVS(325)*XX(46)-JVS(335)&
             &*XX(47)-JVS(377)*XX(48)-JVS(388)*XX(49)-JVS(399)*XX(50)-JVS(413)*XX(51)-JVS(450)*XX(53)-JVS(465)*XX(54)&
             &-JVS(478)*XX(55)-JVS(502)*XX(56)-JVS(519)*XX(57)-JVS(530)*XX(58)-JVS(548)*XX(59)-JVS(565)*XX(60)-JVS(587)&
             &*XX(61)-JVS(603)*XX(62))/(JVS(628))
  XX(64) = (X(64)-JVS(37)*XX(4)-JVS(116)*XX(23)-JVS(255)*XX(40)-JVS(289)*XX(42)-JVS(378)*XX(48)-JVS(451)*XX(53)-JVS(503)&
             &*XX(56)-JVS(588)*XX(61))/(JVS(652))
  XX(65) = (X(65)-JVS(38)*XX(4)-JVS(76)*XX(17)-JVS(117)*XX(23)-JVS(203)*XX(34)-JVS(217)*XX(36)-JVS(230)*XX(38)-JVS(238)&
             &*XX(39)-JVS(256)*XX(40)-JVS(266)*XX(41)-JVS(290)*XX(42)-JVS(301)*XX(43)-JVS(309)*XX(44)-JVS(317)*XX(45)&
             &-JVS(326)*XX(46)-JVS(336)*XX(47)-JVS(379)*XX(48)-JVS(389)*XX(49)-JVS(400)*XX(50)-JVS(414)*XX(51)-JVS(422)&
             &*XX(52)-JVS(452)*XX(53)-JVS(466)*XX(54)-JVS(479)*XX(55)-JVS(504)*XX(56)-JVS(520)*XX(57)-JVS(531)*XX(58)&
             &-JVS(549)*XX(59)-JVS(566)*XX(60)-JVS(589)*XX(61)-JVS(604)*XX(62)-JVS(629)*XX(63)-JVS(653)*XX(64))/(JVS(691))
  XX(66) = (X(66)-JVS(39)*XX(4)-JVS(131)*XX(26)-JVS(177)*XX(31)-JVS(204)*XX(34)-JVS(218)*XX(36)-JVS(231)*XX(38)-JVS(239)&
             &*XX(39)-JVS(257)*XX(40)-JVS(267)*XX(41)-JVS(291)*XX(42)-JVS(302)*XX(43)-JVS(310)*XX(44)-JVS(318)*XX(45)&
             &-JVS(327)*XX(46)-JVS(337)*XX(47)-JVS(380)*XX(48)-JVS(390)*XX(49)-JVS(401)*XX(50)-JVS(415)*XX(51)-JVS(423)&
             &*XX(52)-JVS(453)*XX(53)-JVS(467)*XX(54)-JVS(480)*XX(55)-JVS(505)*XX(56)-JVS(521)*XX(57)-JVS(532)*XX(58)&
             &-JVS(550)*XX(59)-JVS(567)*XX(60)-JVS(590)*XX(61)-JVS(605)*XX(62)-JVS(630)*XX(63)-JVS(654)*XX(64)-JVS(692)&
             &*XX(65))/(JVS(734))
  XX(67) = (X(67)-JVS(40)*XX(4)-JVS(72)*XX(16)-JVS(77)*XX(17)-JVS(94)*XX(19)-JVS(98)*XX(20)-JVS(132)*XX(26)-JVS(142)&
             &*XX(27)-JVS(187)*XX(32)-JVS(219)*XX(36)-JVS(232)*XX(38)-JVS(240)*XX(39)-JVS(268)*XX(41)-JVS(292)*XX(42)&
             &-JVS(303)*XX(43)-JVS(311)*XX(44)-JVS(319)*XX(45)-JVS(328)*XX(46)-JVS(338)*XX(47)-JVS(381)*XX(48)-JVS(391)&
             &*XX(49)-JVS(402)*XX(50)-JVS(416)*XX(51)-JVS(424)*XX(52)-JVS(454)*XX(53)-JVS(468)*XX(54)-JVS(481)*XX(55)&
             &-JVS(506)*XX(56)-JVS(522)*XX(57)-JVS(533)*XX(58)-JVS(551)*XX(59)-JVS(568)*XX(60)-JVS(591)*XX(61)-JVS(606)&
             &*XX(62)-JVS(631)*XX(63)-JVS(655)*XX(64)-JVS(693)*XX(65)-JVS(735)*XX(66))/(JVS(796))
  XX(68) = (X(68)-JVS(3)*XX(1)-JVS(6)*XX(2)-JVS(18)*XX(3)-JVS(42)*XX(5)-JVS(46)*XX(7)-JVS(48)*XX(8)-JVS(50)*XX(9)&
             &-JVS(55)*XX(11)-JVS(57)*XX(12)-JVS(62)*XX(13)-JVS(67)*XX(14)-JVS(70)*XX(15)-JVS(73)*XX(16)-JVS(78)*XX(17)&
             &-JVS(83)*XX(18)-JVS(95)*XX(19)-JVS(99)*XX(20)-JVS(105)*XX(21)-JVS(110)*XX(22)-JVS(118)*XX(23)-JVS(121)*XX(24)&
             &-JVS(126)*XX(25)-JVS(133)*XX(26)-JVS(143)*XX(27)-JVS(164)*XX(28)-JVS(168)*XX(29)-JVS(172)*XX(30)-JVS(178)&
             &*XX(31)-JVS(188)*XX(32)-JVS(193)*XX(33)-JVS(205)*XX(34)-JVS(211)*XX(35)-JVS(220)*XX(36)-JVS(226)*XX(37)&
             &-JVS(233)*XX(38)-JVS(241)*XX(39)-JVS(258)*XX(40)-JVS(269)*XX(41)-JVS(293)*XX(42)-JVS(304)*XX(43)-JVS(312)&
             &*XX(44)-JVS(320)*XX(45)-JVS(329)*XX(46)-JVS(339)*XX(47)-JVS(382)*XX(48)-JVS(392)*XX(49)-JVS(403)*XX(50)&
             &-JVS(417)*XX(51)-JVS(425)*XX(52)-JVS(455)*XX(53)-JVS(469)*XX(54)-JVS(482)*XX(55)-JVS(507)*XX(56)-JVS(523)&
             &*XX(57)-JVS(534)*XX(58)-JVS(552)*XX(59)-JVS(569)*XX(60)-JVS(592)*XX(61)-JVS(607)*XX(62)-JVS(632)*XX(63)&
             &-JVS(656)*XX(64)-JVS(694)*XX(65)-JVS(736)*XX(66)-JVS(797)*XX(67))/(JVS(860))
  XX(69) = (X(69)-JVS(294)*XX(42)-JVS(426)*XX(52)-JVS(456)*XX(53)-JVS(508)*XX(56)-JVS(593)*XX(61)-JVS(657)*XX(64)&
             &-JVS(695)*XX(65)-JVS(737)*XX(66)-JVS(798)*XX(67)-JVS(861)*XX(68))/(JVS(886))
  XX(70) = (X(70)-JVS(52)*XX(10)-JVS(100)*XX(20)-JVS(122)*XX(24)-JVS(144)*XX(27)-JVS(165)*XX(28)-JVS(169)*XX(29)&
             &-JVS(173)*XX(30)-JVS(179)*XX(31)-JVS(189)*XX(32)-JVS(194)*XX(33)-JVS(206)*XX(34)-JVS(212)*XX(35)-JVS(221)&
             &*XX(36)-JVS(227)*XX(37)-JVS(234)*XX(38)-JVS(242)*XX(39)-JVS(259)*XX(40)-JVS(270)*XX(41)-JVS(295)*XX(42)&
             &-JVS(305)*XX(43)-JVS(313)*XX(44)-JVS(321)*XX(45)-JVS(330)*XX(46)-JVS(340)*XX(47)-JVS(383)*XX(48)-JVS(393)&
             &*XX(49)-JVS(404)*XX(50)-JVS(418)*XX(51)-JVS(427)*XX(52)-JVS(457)*XX(53)-JVS(470)*XX(54)-JVS(483)*XX(55)&
             &-JVS(509)*XX(56)-JVS(524)*XX(57)-JVS(535)*XX(58)-JVS(553)*XX(59)-JVS(570)*XX(60)-JVS(594)*XX(61)-JVS(608)&
             &*XX(62)-JVS(633)*XX(63)-JVS(658)*XX(64)-JVS(696)*XX(65)-JVS(738)*XX(66)-JVS(799)*XX(67)-JVS(862)*XX(68)&
             &-JVS(887)*XX(69))/(JVS(935))
  XX(71) = (X(71)-JVS(296)*XX(42)-JVS(428)*XX(52)-JVS(458)*XX(53)-JVS(510)*XX(56)-JVS(595)*XX(61)-JVS(659)*XX(64)&
             &-JVS(697)*XX(65)-JVS(739)*XX(66)-JVS(800)*XX(67)-JVS(863)*XX(68)-JVS(888)*XX(69)-JVS(936)*XX(70))/(JVS(966))
  XX(72) = (X(72)-JVS(84)*XX(18)-JVS(207)*XX(34)-JVS(222)*XX(36)-JVS(235)*XX(38)-JVS(243)*XX(39)-JVS(260)*XX(40)&
             &-JVS(271)*XX(41)-JVS(297)*XX(42)-JVS(306)*XX(43)-JVS(314)*XX(44)-JVS(322)*XX(45)-JVS(331)*XX(46)-JVS(341)&
             &*XX(47)-JVS(384)*XX(48)-JVS(394)*XX(49)-JVS(405)*XX(50)-JVS(419)*XX(51)-JVS(429)*XX(52)-JVS(459)*XX(53)&
             &-JVS(471)*XX(54)-JVS(484)*XX(55)-JVS(511)*XX(56)-JVS(525)*XX(57)-JVS(536)*XX(58)-JVS(554)*XX(59)-JVS(571)&
             &*XX(60)-JVS(596)*XX(61)-JVS(609)*XX(62)-JVS(634)*XX(63)-JVS(660)*XX(64)-JVS(698)*XX(65)-JVS(740)*XX(66)&
             &-JVS(801)*XX(67)-JVS(864)*XX(68)-JVS(889)*XX(69)-JVS(937)*XX(70)-JVS(967)*XX(71))/(JVS(1003))
  XX(73) = (X(73)-JVS(53)*XX(10)-JVS(74)*XX(16)-JVS(85)*XX(18)-JVS(101)*XX(20)-JVS(106)*XX(21)-JVS(111)*XX(22)-JVS(127)&
             &*XX(25)-JVS(145)*XX(27)-JVS(180)*XX(31)-JVS(190)*XX(32)-JVS(213)*XX(35)-JVS(261)*XX(40)-JVS(272)*XX(41)&
             &-JVS(332)*XX(46)-JVS(342)*XX(47)-JVS(385)*XX(48)-JVS(472)*XX(54)-JVS(485)*XX(55)-JVS(512)*XX(56)-JVS(597)&
             &*XX(61)-JVS(610)*XX(62)-JVS(635)*XX(63)-JVS(661)*XX(64)-JVS(699)*XX(65)-JVS(741)*XX(66)-JVS(802)*XX(67)&
             &-JVS(865)*XX(68)-JVS(890)*XX(69)-JVS(938)*XX(70)-JVS(968)*XX(71)-JVS(1004)*XX(72))/(JVS(1052))
  XX(73) = XX(73)
  XX(72) = XX(72)-JVS(1051)*XX(73)
  XX(71) = XX(71)-JVS(1002)*XX(72)-JVS(1050)*XX(73)
  XX(70) = XX(70)-JVS(965)*XX(71)-JVS(1001)*XX(72)-JVS(1049)*XX(73)
  XX(69) = XX(69)-JVS(934)*XX(70)-JVS(964)*XX(71)-JVS(1000)*XX(72)-JVS(1048)*XX(73)
  XX(68) = XX(68)-JVS(885)*XX(69)-JVS(933)*XX(70)-JVS(963)*XX(71)-JVS(999)*XX(72)-JVS(1047)*XX(73)
  XX(67) = XX(67)-JVS(859)*XX(68)-JVS(884)*XX(69)-JVS(932)*XX(70)-JVS(962)*XX(71)-JVS(998)*XX(72)-JVS(1046)*XX(73)
  XX(66) = XX(66)-JVS(795)*XX(67)-JVS(858)*XX(68)-JVS(883)*XX(69)-JVS(931)*XX(70)-JVS(961)*XX(71)-JVS(997)*XX(72)&
             &-JVS(1045)*XX(73)
  XX(65) = XX(65)-JVS(733)*XX(66)-JVS(794)*XX(67)-JVS(857)*XX(68)-JVS(882)*XX(69)-JVS(930)*XX(70)-JVS(960)*XX(71)&
             &-JVS(996)*XX(72)-JVS(1044)*XX(73)
  XX(64) = XX(64)-JVS(690)*XX(65)-JVS(732)*XX(66)-JVS(793)*XX(67)-JVS(856)*XX(68)-JVS(929)*XX(70)-JVS(959)*XX(71)&
             &-JVS(995)*XX(72)-JVS(1043)*XX(73)
  XX(63) = XX(63)-JVS(651)*XX(64)-JVS(689)*XX(65)-JVS(731)*XX(66)-JVS(792)*XX(67)-JVS(855)*XX(68)-JVS(881)*XX(69)&
             &-JVS(928)*XX(70)-JVS(958)*XX(71)-JVS(994)*XX(72)-JVS(1042)*XX(73)
  XX(62) = XX(62)-JVS(627)*XX(63)-JVS(650)*XX(64)-JVS(688)*XX(65)-JVS(730)*XX(66)-JVS(791)*XX(67)-JVS(854)*XX(68)&
             &-JVS(880)*XX(69)-JVS(927)*XX(70)-JVS(957)*XX(71)-JVS(993)*XX(72)-JVS(1041)*XX(73)
  XX(61) = XX(61)-JVS(687)*XX(65)-JVS(729)*XX(66)-JVS(790)*XX(67)-JVS(853)*XX(68)-JVS(926)*XX(70)-JVS(956)*XX(71)&
             &-JVS(992)*XX(72)-JVS(1040)*XX(73)
  XX(60) = XX(60)-JVS(584)*XX(61)-JVS(649)*XX(64)-JVS(686)*XX(65)-JVS(728)*XX(66)-JVS(789)*XX(67)-JVS(852)*XX(68)&
             &-JVS(879)*XX(69)-JVS(925)*XX(70)-JVS(991)*XX(72)-JVS(1039)*XX(73)
  XX(59) = XX(59)-JVS(562)*XX(60)-JVS(583)*XX(61)-JVS(648)*XX(64)-JVS(685)*XX(65)-JVS(727)*XX(66)-JVS(788)*XX(67)&
             &-JVS(851)*XX(68)-JVS(878)*XX(69)-JVS(924)*XX(70)-JVS(990)*XX(72)-JVS(1038)*XX(73)
  XX(58) = XX(58)-JVS(544)*XX(59)-JVS(561)*XX(60)-JVS(582)*XX(61)-JVS(626)*XX(63)-JVS(647)*XX(64)-JVS(684)*XX(65)&
             &-JVS(726)*XX(66)-JVS(787)*XX(67)-JVS(850)*XX(68)-JVS(923)*XX(70)-JVS(955)*XX(71)-JVS(989)*XX(72)-JVS(1037)&
             &*XX(73)
  XX(57) = XX(57)-JVS(543)*XX(59)-JVS(560)*XX(60)-JVS(581)*XX(61)-JVS(625)*XX(63)-JVS(646)*XX(64)-JVS(683)*XX(65)&
             &-JVS(725)*XX(66)-JVS(786)*XX(67)-JVS(849)*XX(68)-JVS(922)*XX(70)-JVS(954)*XX(71)-JVS(988)*XX(72)-JVS(1036)&
             &*XX(73)
  XX(56) = XX(56)-JVS(682)*XX(65)-JVS(724)*XX(66)-JVS(785)*XX(67)-JVS(848)*XX(68)-JVS(921)*XX(70)-JVS(953)*XX(71)&
             &-JVS(987)*XX(72)-JVS(1035)*XX(73)
  XX(55) = XX(55)-JVS(499)*XX(56)-JVS(624)*XX(63)-JVS(681)*XX(65)-JVS(723)*XX(66)-JVS(784)*XX(67)-JVS(847)*XX(68)&
             &-JVS(877)*XX(69)-JVS(920)*XX(70)-JVS(952)*XX(71)-JVS(986)*XX(72)-JVS(1034)*XX(73)
  XX(54) = XX(54)-JVS(475)*XX(55)-JVS(498)*XX(56)-JVS(580)*XX(61)-JVS(623)*XX(63)-JVS(645)*XX(64)-JVS(680)*XX(65)&
             &-JVS(722)*XX(66)-JVS(783)*XX(67)-JVS(846)*XX(68)-JVS(919)*XX(70)-JVS(1033)*XX(73)
  XX(53) = XX(53)-JVS(679)*XX(65)-JVS(721)*XX(66)-JVS(782)*XX(67)-JVS(845)*XX(68)-JVS(918)*XX(70)-JVS(1032)*XX(73)
  XX(52) = XX(52)-JVS(441)*XX(53)-JVS(497)*XX(56)-JVS(579)*XX(61)-JVS(644)*XX(64)-JVS(678)*XX(65)-JVS(720)*XX(66)&
             &-JVS(781)*XX(67)-JVS(844)*XX(68)-JVS(876)*XX(69)-JVS(917)*XX(70)-JVS(951)*XX(71)-JVS(985)*XX(72)-JVS(1031)&
             &*XX(73)
  XX(51) = XX(51)-JVS(440)*XX(53)-JVS(474)*XX(55)-JVS(496)*XX(56)-JVS(542)*XX(59)-JVS(601)*XX(62)-JVS(622)*XX(63)&
             &-JVS(719)*XX(66)-JVS(780)*XX(67)-JVS(843)*XX(68)-JVS(916)*XX(70)-JVS(950)*XX(71)-JVS(1030)*XX(73)
  XX(50) = XX(50)-JVS(410)*XX(51)-JVS(516)*XX(57)-JVS(527)*XX(58)-JVS(677)*XX(65)-JVS(718)*XX(66)-JVS(779)*XX(67)&
             &-JVS(842)*XX(68)-JVS(875)*XX(69)-JVS(915)*XX(70)-JVS(949)*XX(71)-JVS(984)*XX(72)-JVS(1029)*XX(73)
  XX(49) = XX(49)-JVS(439)*XX(53)-JVS(578)*XX(61)-JVS(643)*XX(64)-JVS(676)*XX(65)-JVS(717)*XX(66)-JVS(778)*XX(67)&
             &-JVS(841)*XX(68)-JVS(874)*XX(69)-JVS(914)*XX(70)-JVS(948)*XX(71)-JVS(983)*XX(72)-JVS(1028)*XX(73)
  XX(48) = XX(48)-JVS(777)*XX(67)-JVS(840)*XX(68)-JVS(913)*XX(70)-JVS(1027)*XX(73)
  XX(47) = XX(47)-JVS(362)*XX(48)-JVS(462)*XX(54)-JVS(675)*XX(65)-JVS(716)*XX(66)-JVS(776)*XX(67)-JVS(839)*XX(68)&
             &-JVS(873)*XX(69)-JVS(912)*XX(70)-JVS(947)*XX(71)-JVS(982)*XX(72)-JVS(1026)*XX(73)
  XX(46) = XX(46)-JVS(361)*XX(48)-JVS(461)*XX(54)-JVS(674)*XX(65)-JVS(715)*XX(66)-JVS(775)*XX(67)-JVS(838)*XX(68)&
             &-JVS(872)*XX(69)-JVS(911)*XX(70)-JVS(946)*XX(71)-JVS(981)*XX(72)-JVS(1025)*XX(73)
  XX(45) = XX(45)-JVS(360)*XX(48)-JVS(438)*XX(53)-JVS(495)*XX(56)-JVS(577)*XX(61)-JVS(642)*XX(64)-JVS(673)*XX(65)&
             &-JVS(714)*XX(66)-JVS(774)*XX(67)-JVS(837)*XX(68)-JVS(871)*XX(69)-JVS(910)*XX(70)-JVS(945)*XX(71)-JVS(980)&
             &*XX(72)-JVS(1024)*XX(73)
  XX(44) = XX(44)-JVS(359)*XX(48)-JVS(437)*XX(53)-JVS(494)*XX(56)-JVS(576)*XX(61)-JVS(641)*XX(64)-JVS(672)*XX(65)&
             &-JVS(713)*XX(66)-JVS(773)*XX(67)-JVS(836)*XX(68)-JVS(870)*XX(69)-JVS(909)*XX(70)-JVS(944)*XX(71)-JVS(979)&
             &*XX(72)-JVS(1023)*XX(73)
  XX(43) = XX(43)-JVS(358)*XX(48)-JVS(409)*XX(51)-JVS(526)*XX(58)-JVS(671)*XX(65)-JVS(712)*XX(66)-JVS(772)*XX(67)&
             &-JVS(869)*XX(69)-JVS(908)*XX(70)-JVS(943)*XX(71)-JVS(978)*XX(72)-JVS(1022)*XX(73)
  XX(42) = XX(42)-JVS(575)*XX(61)-JVS(640)*XX(64)-JVS(711)*XX(66)-JVS(835)*XX(68)
  XX(41) = XX(41)-JVS(357)*XX(48)-JVS(670)*XX(65)-JVS(710)*XX(66)-JVS(771)*XX(67)-JVS(834)*XX(68)-JVS(907)*XX(70)&
             &-JVS(942)*XX(71)-JVS(977)*XX(72)-JVS(1021)*XX(73)
  XX(40) = XX(40)-JVS(709)*XX(66)-JVS(770)*XX(67)-JVS(833)*XX(68)-JVS(906)*XX(70)-JVS(1020)*XX(73)
  XX(39) = XX(39)-JVS(277)*XX(42)-JVS(356)*XX(48)-JVS(436)*XX(53)-JVS(669)*XX(65)-JVS(708)*XX(66)-JVS(769)*XX(67)&
             &-JVS(905)*XX(70)-JVS(941)*XX(71)-JVS(976)*XX(72)-JVS(1019)*XX(73)
  XX(38) = XX(38)-JVS(276)*XX(42)-JVS(355)*XX(48)-JVS(435)*XX(53)-JVS(668)*XX(65)-JVS(707)*XX(66)-JVS(768)*XX(67)&
             &-JVS(904)*XX(70)-JVS(940)*XX(71)-JVS(975)*XX(72)-JVS(1018)*XX(73)
  XX(37) = XX(37)-JVS(354)*XX(48)-JVS(396)*XX(50)-JVS(408)*XX(51)-JVS(460)*XX(54)-JVS(493)*XX(56)-JVS(515)*XX(57)&
             &-JVS(541)*XX(59)-JVS(559)*XX(60)-JVS(600)*XX(62)-JVS(621)*XX(63)-JVS(639)*XX(64)-JVS(667)*XX(65)-JVS(706)&
             &*XX(66)-JVS(767)*XX(67)-JVS(832)*XX(68)-JVS(903)*XX(70)
  XX(36) = XX(36)-JVS(353)*XX(48)-JVS(434)*XX(53)-JVS(666)*XX(65)-JVS(705)*XX(66)-JVS(766)*XX(67)-JVS(902)*XX(70)&
             &-JVS(939)*XX(71)-JVS(974)*XX(72)-JVS(1017)*XX(73)
  XX(35) = XX(35)-JVS(245)*XX(40)-JVS(352)*XX(48)-JVS(473)*XX(55)-JVS(492)*XX(56)-JVS(620)*XX(63)-JVS(704)*XX(66)&
             &-JVS(765)*XX(67)-JVS(831)*XX(68)-JVS(868)*XX(69)-JVS(901)*XX(70)-JVS(1016)*XX(73)
  XX(34) = XX(34)-JVS(351)*XX(48)-JVS(764)*XX(67)-JVS(830)*XX(68)-JVS(900)*XX(70)-JVS(1015)*XX(73)
  XX(33) = XX(33)-JVS(350)*XX(48)-JVS(395)*XX(50)-JVS(407)*XX(51)-JVS(491)*XX(56)-JVS(514)*XX(57)-JVS(540)*XX(59)&
             &-JVS(558)*XX(60)-JVS(599)*XX(62)-JVS(619)*XX(63)-JVS(638)*XX(64)-JVS(665)*XX(65)-JVS(703)*XX(66)-JVS(763)&
             &*XX(67)-JVS(829)*XX(68)-JVS(899)*XX(70)
  XX(32) = XX(32)-JVS(263)*XX(41)-JVS(490)*XX(56)-JVS(618)*XX(63)-JVS(762)*XX(67)-JVS(828)*XX(68)-JVS(867)*XX(69)&
             &-JVS(898)*XX(70)-JVS(973)*XX(72)-JVS(1014)*XX(73)
  XX(31) = XX(31)-JVS(349)*XX(48)-JVS(489)*XX(56)-JVS(702)*XX(66)-JVS(827)*XX(68)-JVS(897)*XX(70)-JVS(1013)*XX(73)
  XX(30) = XX(30)-JVS(298)*XX(43)-JVS(348)*XX(48)-JVS(406)*XX(51)-JVS(513)*XX(57)-JVS(539)*XX(59)-JVS(557)*XX(60)&
             &-JVS(574)*XX(61)-JVS(617)*XX(63)-JVS(637)*XX(64)-JVS(761)*XX(67)-JVS(826)*XX(68)-JVS(896)*XX(70)
  XX(29) = XX(29)-JVS(275)*XX(42)-JVS(386)*XX(49)-JVS(433)*XX(53)-JVS(538)*XX(59)-JVS(556)*XX(60)-JVS(573)*XX(61)&
             &-JVS(616)*XX(63)-JVS(636)*XX(64)-JVS(760)*XX(67)-JVS(825)*XX(68)-JVS(895)*XX(70)
  XX(28) = XX(28)-JVS(759)*XX(67)-JVS(824)*XX(68)
  XX(27) = XX(27)-JVS(823)*XX(68)-JVS(894)*XX(70)-JVS(1012)*XX(73)
  XX(26) = XX(26)-JVS(347)*XX(48)-JVS(488)*XX(56)-JVS(664)*XX(65)-JVS(701)*XX(66)-JVS(758)*XX(67)-JVS(822)*XX(68)
  XX(25) = XX(25)-JVS(184)*XX(32)-JVS(262)*XX(41)-JVS(615)*XX(63)-JVS(757)*XX(67)-JVS(821)*XX(68)-JVS(972)*XX(72)&
             &-JVS(1011)*XX(73)
  XX(24) = XX(24)-JVS(147)*XX(28)-JVS(214)*XX(36)-JVS(346)*XX(48)-JVS(537)*XX(59)-JVS(555)*XX(60)-JVS(614)*XX(63)&
             &-JVS(756)*XX(67)-JVS(820)*XX(68)-JVS(893)*XX(70)
  XX(23) = XX(23)-JVS(244)*XX(40)-JVS(345)*XX(48)-JVS(700)*XX(66)-JVS(755)*XX(67)-JVS(819)*XX(68)
  XX(22) = XX(22)-JVS(183)*XX(32)-JVS(333)*XX(47)-JVS(613)*XX(63)-JVS(754)*XX(67)-JVS(818)*XX(68)-JVS(971)*XX(72)&
             &-JVS(1010)*XX(73)
  XX(21) = XX(21)-JVS(182)*XX(32)-JVS(323)*XX(46)-JVS(612)*XX(63)-JVS(753)*XX(67)-JVS(817)*XX(68)-JVS(970)*XX(72)&
             &-JVS(1009)*XX(73)
  XX(20) = XX(20)-JVS(181)*XX(32)-JVS(752)*XX(67)-JVS(866)*XX(69)-JVS(1008)*XX(73)
  XX(19) = XX(19)-JVS(751)*XX(67)-JVS(816)*XX(68)
  XX(18) = XX(18)-JVS(815)*XX(68)-JVS(969)*XX(72)-JVS(1007)*XX(73)
  XX(17) = XX(17)-JVS(344)*XX(48)-JVS(663)*XX(65)-JVS(750)*XX(67)-JVS(814)*XX(68)
  XX(16) = XX(16)-JVS(749)*XX(67)-JVS(813)*XX(68)-JVS(892)*XX(70)-JVS(1006)*XX(73)
  XX(15) = XX(15)-JVS(274)*XX(42)-JVS(432)*XX(53)-JVS(748)*XX(67)-JVS(812)*XX(68)
  XX(14) = XX(14)-JVS(662)*XX(65)-JVS(811)*XX(68)
  XX(13) = XX(13)-JVS(572)*XX(61)-JVS(810)*XX(68)
  XX(12) = XX(12)-JVS(146)*XX(28)-JVS(195)*XX(34)-JVS(343)*XX(48)-JVS(420)*XX(52)-JVS(431)*XX(53)-JVS(747)*XX(67)&
             &-JVS(809)*XX(68)
  XX(11) = XX(11)-JVS(112)*XX(23)-JVS(315)*XX(45)-JVS(430)*XX(53)-JVS(746)*XX(67)-JVS(808)*XX(68)
  XX(10) = XX(10)-JVS(891)*XX(70)-JVS(1005)*XX(73)
  XX(9) = XX(9)-JVS(107)*XX(22)-JVS(487)*XX(56)-JVS(745)*XX(67)-JVS(807)*XX(68)
  XX(8) = XX(8)-JVS(102)*XX(21)-JVS(486)*XX(56)-JVS(744)*XX(67)-JVS(806)*XX(68)
  XX(7) = XX(7)-JVS(273)*XX(42)-JVS(307)*XX(44)-JVS(743)*XX(67)-JVS(805)*XX(68)
  XX(6) = XX(6)-JVS(598)*XX(62)-JVS(611)*XX(63)-JVS(804)*XX(68)
  XX(5) = XX(5)-JVS(742)*XX(67)-JVS(803)*XX(68)
  XX(4) = XX(4)
  XX(3) = XX(3)
  XX(2) = XX(2)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE racmsorg_LinearAlgebra

