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
! File                 : racm_mim_LinearAlgebra.f90
! Time                 : Wed Jun 14 17:35:26 2017
! Working directory    : /data/ksetigui/setigui/WRFnew/chem/KPP/mechanisms/racm_mim
! Equation file        : racm_mim.kpp
! Output root filename : racm_mim
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racm_mim_LinearAlgebra

  USE racm_mim_Parameters
  USE racm_mim_JacobianSP

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
  XX(5) = (X(5)-JVS(2)*XX(1))/(JVS(40))
  XX(6) = X(6)/JVS(42)
  XX(7) = X(7)/JVS(44)
  XX(8) = X(8)/JVS(48)
  XX(9) = X(9)/JVS(50)
  XX(10) = X(10)/JVS(52)
  XX(11) = X(11)/JVS(54)
  XX(12) = X(12)/JVS(57)
  XX(13) = X(13)/JVS(59)
  XX(14) = (X(14)-JVS(8)*XX(3))/(JVS(63))
  XX(15) = X(15)/JVS(65)
  XX(16) = X(16)/JVS(70)
  XX(17) = X(17)/JVS(73)
  XX(18) = X(18)/JVS(78)
  XX(19) = X(19)/JVS(81)
  XX(20) = X(20)/JVS(85)
  XX(21) = X(21)/JVS(89)
  XX(22) = X(22)/JVS(94)
  XX(23) = X(23)/JVS(100)
  XX(24) = X(24)/JVS(107)
  XX(25) = X(25)/JVS(118)
  XX(26) = X(26)/JVS(124)
  XX(27) = X(27)/JVS(130)
  XX(28) = (X(28)-JVS(101)*XX(23))/(JVS(137))
  XX(29) = (X(29)-JVS(102)*XX(23))/(JVS(142))
  XX(30) = (X(30)-JVS(9)*XX(3))/(JVS(146))
  XX(31) = (X(31)-JVS(103)*XX(23))/(JVS(150))
  XX(32) = X(32)/JVS(155)
  XX(33) = X(33)/JVS(166)
  XX(34) = (X(34)-JVS(79)*XX(18)-JVS(95)*XX(22))/(JVS(172))
  XX(35) = X(35)/JVS(178)
  XX(36) = (X(36)-JVS(108)*XX(24))/(JVS(185))
  XX(37) = (X(37)-JVS(5)*XX(2))/(JVS(194))
  XX(38) = (X(38)-JVS(10)*XX(3)-JVS(20)*XX(4)-JVS(109)*XX(24)-JVS(195)*XX(37))/(JVS(213))
  XX(39) = (X(39)-JVS(11)*XX(3)-JVS(110)*XX(24)-JVS(173)*XX(34)-JVS(196)*XX(37))/(JVS(217))
  XX(40) = (X(40)-JVS(131)*XX(27)-JVS(151)*XX(31)-JVS(156)*XX(32))/(JVS(225))
  XX(41) = (X(41)-JVS(12)*XX(3)-JVS(111)*XX(24)-JVS(125)*XX(26)-JVS(197)*XX(37))/(JVS(231))
  XX(42) = (X(42)-JVS(60)*XX(13)-JVS(90)*XX(21)-JVS(96)*XX(22)-JVS(198)*XX(37))/(JVS(238))
  XX(43) = (X(43)-JVS(13)*XX(3)-JVS(119)*XX(25)-JVS(179)*XX(35)-JVS(199)*XX(37))/(JVS(246))
  XX(44) = (X(44)-JVS(157)*XX(32)-JVS(200)*XX(37))/(JVS(253))
  XX(45) = (X(45)-JVS(21)*XX(4))/(JVS(266))
  XX(46) = (X(46)-JVS(22)*XX(4))/(JVS(274))
  XX(47) = (X(47)-JVS(23)*XX(4))/(JVS(282))
  XX(48) = (X(48)-JVS(14)*XX(3)-JVS(201)*XX(37)-JVS(239)*XX(42))/(JVS(294))
  XX(49) = X(49)/JVS(309)
  XX(50) = (X(50)-JVS(24)*XX(4)-JVS(45)*XX(7)-JVS(174)*XX(34)-JVS(202)*XX(37)-JVS(295)*XX(48))/(JVS(331))
  XX(51) = (X(51)-JVS(158)*XX(32)-JVS(203)*XX(37))/(JVS(347))
  XX(52) = (X(52)-JVS(254)*XX(44)-JVS(348)*XX(51))/(JVS(366))
  XX(53) = (X(53)-JVS(25)*XX(4)-JVS(310)*XX(49))/(JVS(376))
  XX(54) = (X(54)-JVS(296)*XX(48)-JVS(349)*XX(51))/(JVS(384))
  XX(55) = (X(55)-JVS(26)*XX(4)-JVS(311)*XX(49))/(JVS(393))
  XX(56) = (X(56)-JVS(159)*XX(32)-JVS(204)*XX(37))/(JVS(423))
  XX(57) = (X(57)-JVS(255)*XX(44)-JVS(350)*XX(51)-JVS(424)*XX(56))/(JVS(446))
  XX(58) = (X(58)-JVS(256)*XX(44)-JVS(351)*XX(51)-JVS(425)*XX(56))/(JVS(456))
  XX(59) = (X(59)-JVS(312)*XX(49)-JVS(426)*XX(56))/(JVS(466))
  XX(60) = (X(60)-JVS(313)*XX(49))/(JVS(483))
  XX(61) = (X(61)-JVS(15)*XX(3)-JVS(27)*XX(4)-JVS(71)*XX(16)-JVS(160)*XX(32)-JVS(167)*XX(33)-JVS(205)*XX(37)-JVS(257)&
             &*XX(44)-JVS(314)*XX(49)-JVS(352)*XX(51)-JVS(427)*XX(56))/(JVS(498))
  XX(62) = (X(62)-JVS(428)*XX(56))/(JVS(517))
  XX(63) = (X(63)-JVS(16)*XX(3)-JVS(28)*XX(4)-JVS(66)*XX(15)-JVS(74)*XX(17)-JVS(112)*XX(24)-JVS(206)*XX(37)-JVS(275)&
             &*XX(46)-JVS(315)*XX(49)-JVS(429)*XX(56))/(JVS(533))
  XX(64) = (X(64)-JVS(161)*XX(32)-JVS(207)*XX(37))/(JVS(554))
  XX(65) = (X(65)-JVS(29)*XX(4)-JVS(67)*XX(15)-JVS(75)*XX(17)-JVS(113)*XX(24)-JVS(208)*XX(37)-JVS(283)*XX(47)-JVS(316)&
             &*XX(49)-JVS(430)*XX(56)-JVS(555)*XX(64))/(JVS(573))
  XX(66) = (X(66)-JVS(30)*XX(4)-JVS(317)*XX(49)-JVS(431)*XX(56)-JVS(484)*XX(60)-JVS(556)*XX(64))/(JVS(587))
  XX(67) = (X(67)-JVS(31)*XX(4)-JVS(318)*XX(49)-JVS(432)*XX(56)-JVS(485)*XX(60)-JVS(557)*XX(64)-JVS(588)*XX(66))&
             &/(JVS(603))
  XX(68) = (X(68)-JVS(32)*XX(4)-JVS(433)*XX(56)-JVS(558)*XX(64))/(JVS(624))
  XX(69) = (X(69)-JVS(33)*XX(4)-JVS(120)*XX(25)-JVS(319)*XX(49)-JVS(353)*XX(51)-JVS(434)*XX(56)-JVS(518)*XX(62)-JVS(559)&
             &*XX(64)-JVS(625)*XX(68))/(JVS(650))
  XX(70) = (X(70)-JVS(320)*XX(49)-JVS(519)*XX(62)-JVS(560)*XX(64)-JVS(626)*XX(68)-JVS(651)*XX(69))/(JVS(677))
  XX(71) = (X(71)-JVS(34)*XX(4)-JVS(258)*XX(44)-JVS(321)*XX(49)-JVS(354)*XX(51)-JVS(435)*XX(56)-JVS(486)*XX(60)-JVS(520)&
             &*XX(62)-JVS(561)*XX(64)-JVS(627)*XX(68)-JVS(652)*XX(69)-JVS(678)*XX(70))/(JVS(693))
  XX(72) = (X(72)-JVS(35)*XX(4)-JVS(168)*XX(33)-JVS(180)*XX(35)-JVS(259)*XX(44)-JVS(267)*XX(45)-JVS(276)*XX(46)-JVS(284)&
             &*XX(47)-JVS(297)*XX(48)-JVS(322)*XX(49)-JVS(332)*XX(50)-JVS(355)*XX(51)-JVS(367)*XX(52)-JVS(377)*XX(53)&
             &-JVS(385)*XX(54)-JVS(394)*XX(55)-JVS(436)*XX(56)-JVS(447)*XX(57)-JVS(457)*XX(58)-JVS(467)*XX(59)-JVS(487)&
             &*XX(60)-JVS(499)*XX(61)-JVS(521)*XX(62)-JVS(534)*XX(63)-JVS(562)*XX(64)-JVS(574)*XX(65)-JVS(589)*XX(66)&
             &-JVS(604)*XX(67)-JVS(628)*XX(68)-JVS(653)*XX(69)-JVS(679)*XX(70)-JVS(694)*XX(71))/(JVS(736))
  XX(73) = (X(73)-JVS(3)*XX(1)-JVS(6)*XX(2)-JVS(17)*XX(3)-JVS(41)*XX(5)-JVS(46)*XX(7)-JVS(49)*XX(8)-JVS(51)*XX(9)&
             &-JVS(53)*XX(10)-JVS(58)*XX(12)-JVS(61)*XX(13)-JVS(64)*XX(14)-JVS(68)*XX(15)-JVS(72)*XX(16)-JVS(76)*XX(17)&
             &-JVS(80)*XX(18)-JVS(82)*XX(19)-JVS(86)*XX(20)-JVS(91)*XX(21)-JVS(97)*XX(22)-JVS(104)*XX(23)-JVS(114)*XX(24)&
             &-JVS(121)*XX(25)-JVS(132)*XX(27)-JVS(138)*XX(28)-JVS(143)*XX(29)-JVS(147)*XX(30)-JVS(152)*XX(31)-JVS(162)&
             &*XX(32)-JVS(169)*XX(33)-JVS(175)*XX(34)-JVS(181)*XX(35)-JVS(186)*XX(36)-JVS(209)*XX(37)-JVS(214)*XX(38)&
             &-JVS(218)*XX(39)-JVS(226)*XX(40)-JVS(232)*XX(41)-JVS(240)*XX(42)-JVS(247)*XX(43)-JVS(260)*XX(44)-JVS(268)&
             &*XX(45)-JVS(277)*XX(46)-JVS(285)*XX(47)-JVS(298)*XX(48)-JVS(323)*XX(49)-JVS(333)*XX(50)-JVS(356)*XX(51)&
             &-JVS(368)*XX(52)-JVS(378)*XX(53)-JVS(386)*XX(54)-JVS(395)*XX(55)-JVS(437)*XX(56)-JVS(448)*XX(57)-JVS(458)&
             &*XX(58)-JVS(468)*XX(59)-JVS(488)*XX(60)-JVS(500)*XX(61)-JVS(522)*XX(62)-JVS(535)*XX(63)-JVS(563)*XX(64)&
             &-JVS(575)*XX(65)-JVS(590)*XX(66)-JVS(605)*XX(67)-JVS(629)*XX(68)-JVS(654)*XX(69)-JVS(680)*XX(70)-JVS(695)&
             &*XX(71)-JVS(737)*XX(72))/(JVS(806))
  XX(74) = (X(74)-JVS(55)*XX(11)-JVS(126)*XX(26)-JVS(133)*XX(27)-JVS(148)*XX(30)-JVS(163)*XX(32)-JVS(176)*XX(34)&
             &-JVS(182)*XX(35)-JVS(187)*XX(36)-JVS(210)*XX(37)-JVS(215)*XX(38)-JVS(219)*XX(39)-JVS(227)*XX(40)-JVS(233)&
             &*XX(41)-JVS(241)*XX(42)-JVS(248)*XX(43)-JVS(261)*XX(44)-JVS(269)*XX(45)-JVS(278)*XX(46)-JVS(286)*XX(47)&
             &-JVS(299)*XX(48)-JVS(324)*XX(49)-JVS(334)*XX(50)-JVS(357)*XX(51)-JVS(369)*XX(52)-JVS(379)*XX(53)-JVS(387)&
             &*XX(54)-JVS(396)*XX(55)-JVS(438)*XX(56)-JVS(449)*XX(57)-JVS(459)*XX(58)-JVS(469)*XX(59)-JVS(489)*XX(60)&
             &-JVS(501)*XX(61)-JVS(523)*XX(62)-JVS(536)*XX(63)-JVS(564)*XX(64)-JVS(576)*XX(65)-JVS(591)*XX(66)-JVS(606)&
             &*XX(67)-JVS(630)*XX(68)-JVS(655)*XX(69)-JVS(681)*XX(70)-JVS(696)*XX(71)-JVS(738)*XX(72)-JVS(807)*XX(73))&
             &/(JVS(852))
  XX(75) = (X(75)-JVS(36)*XX(4)-JVS(47)*XX(7)-JVS(62)*XX(13)-JVS(83)*XX(19)-JVS(87)*XX(20)-JVS(115)*XX(24)-JVS(134)&
             &*XX(27)-JVS(164)*XX(32)-JVS(170)*XX(33)-JVS(228)*XX(40)-JVS(242)*XX(42)-JVS(270)*XX(45)-JVS(279)*XX(46)&
             &-JVS(287)*XX(47)-JVS(300)*XX(48)-JVS(325)*XX(49)-JVS(335)*XX(50)-JVS(358)*XX(51)-JVS(370)*XX(52)-JVS(380)&
             &*XX(53)-JVS(388)*XX(54)-JVS(397)*XX(55)-JVS(439)*XX(56)-JVS(450)*XX(57)-JVS(460)*XX(58)-JVS(470)*XX(59)&
             &-JVS(490)*XX(60)-JVS(502)*XX(61)-JVS(524)*XX(62)-JVS(537)*XX(63)-JVS(565)*XX(64)-JVS(577)*XX(65)-JVS(592)&
             &*XX(66)-JVS(607)*XX(67)-JVS(631)*XX(68)-JVS(656)*XX(69)-JVS(682)*XX(70)-JVS(697)*XX(71)-JVS(739)*XX(72)&
             &-JVS(808)*XX(73)-JVS(853)*XX(74))/(JVS(916))
  XX(76) = (X(76)-JVS(37)*XX(4)-JVS(88)*XX(20)-JVS(122)*XX(25)-JVS(262)*XX(44)-JVS(271)*XX(45)-JVS(280)*XX(46)-JVS(288)&
             &*XX(47)-JVS(301)*XX(48)-JVS(326)*XX(49)-JVS(336)*XX(50)-JVS(359)*XX(51)-JVS(371)*XX(52)-JVS(381)*XX(53)&
             &-JVS(389)*XX(54)-JVS(398)*XX(55)-JVS(440)*XX(56)-JVS(451)*XX(57)-JVS(461)*XX(58)-JVS(471)*XX(59)-JVS(491)&
             &*XX(60)-JVS(503)*XX(61)-JVS(525)*XX(62)-JVS(538)*XX(63)-JVS(566)*XX(64)-JVS(578)*XX(65)-JVS(593)*XX(66)&
             &-JVS(608)*XX(67)-JVS(632)*XX(68)-JVS(657)*XX(69)-JVS(683)*XX(70)-JVS(698)*XX(71)-JVS(740)*XX(72)-JVS(809)&
             &*XX(73)-JVS(854)*XX(74)-JVS(917)*XX(75))/(JVS(953))
  XX(77) = (X(77)-JVS(18)*XX(3)-JVS(38)*XX(4)-JVS(43)*XX(6)-JVS(69)*XX(15)-JVS(77)*XX(17)-JVS(116)*XX(24)-JVS(127)&
             &*XX(26)-JVS(139)*XX(28)-JVS(144)*XX(29)-JVS(149)*XX(30)-JVS(153)*XX(31)-JVS(171)*XX(33)-JVS(183)*XX(35)&
             &-JVS(188)*XX(36)-JVS(211)*XX(37)-JVS(216)*XX(38)-JVS(220)*XX(39)-JVS(229)*XX(40)-JVS(234)*XX(41)-JVS(243)&
             &*XX(42)-JVS(249)*XX(43)-JVS(263)*XX(44)-JVS(272)*XX(45)-JVS(302)*XX(48)-JVS(327)*XX(49)-JVS(337)*XX(50)&
             &-JVS(360)*XX(51)-JVS(372)*XX(52)-JVS(390)*XX(54)-JVS(441)*XX(56)-JVS(452)*XX(57)-JVS(462)*XX(58)-JVS(472)&
             &*XX(59)-JVS(492)*XX(60)-JVS(504)*XX(61)-JVS(526)*XX(62)-JVS(539)*XX(63)-JVS(567)*XX(64)-JVS(579)*XX(65)&
             &-JVS(594)*XX(66)-JVS(609)*XX(67)-JVS(633)*XX(68)-JVS(658)*XX(69)-JVS(684)*XX(70)-JVS(699)*XX(71)-JVS(741)&
             &*XX(72)-JVS(810)*XX(73)-JVS(855)*XX(74)-JVS(918)*XX(75)-JVS(954)*XX(76))/(JVS(981))
  XX(78) = (X(78)-JVS(39)*XX(4)-JVS(250)*XX(43)-JVS(361)*XX(51)-JVS(442)*XX(56)-JVS(493)*XX(60)-JVS(527)*XX(62)-JVS(568)&
             &*XX(64)-JVS(634)*XX(68)-JVS(659)*XX(69)-JVS(685)*XX(70)-JVS(700)*XX(71)-JVS(742)*XX(72)-JVS(811)*XX(73)&
             &-JVS(856)*XX(74)-JVS(919)*XX(75)-JVS(955)*XX(76)-JVS(982)*XX(77))/(JVS(993))
  XX(79) = (X(79)-JVS(98)*XX(22)-JVS(105)*XX(23)-JVS(128)*XX(26)-JVS(177)*XX(34)-JVS(212)*XX(37)-JVS(244)*XX(42)&
             &-JVS(264)*XX(44)-JVS(273)*XX(45)-JVS(281)*XX(46)-JVS(289)*XX(47)-JVS(303)*XX(48)-JVS(328)*XX(49)-JVS(338)&
             &*XX(50)-JVS(362)*XX(51)-JVS(373)*XX(52)-JVS(382)*XX(53)-JVS(391)*XX(54)-JVS(399)*XX(55)-JVS(443)*XX(56)&
             &-JVS(453)*XX(57)-JVS(463)*XX(58)-JVS(473)*XX(59)-JVS(494)*XX(60)-JVS(505)*XX(61)-JVS(528)*XX(62)-JVS(540)&
             &*XX(63)-JVS(569)*XX(64)-JVS(580)*XX(65)-JVS(595)*XX(66)-JVS(610)*XX(67)-JVS(635)*XX(68)-JVS(660)*XX(69)&
             &-JVS(686)*XX(70)-JVS(701)*XX(71)-JVS(743)*XX(72)-JVS(812)*XX(73)-JVS(857)*XX(74)-JVS(920)*XX(75)-JVS(956)&
             &*XX(76)-JVS(983)*XX(77)-JVS(994)*XX(78))/(JVS(1032))
  XX(80) = (X(80)-JVS(56)*XX(11)-JVS(84)*XX(19)-JVS(92)*XX(21)-JVS(99)*XX(22)-JVS(106)*XX(23)-JVS(129)*XX(26)-JVS(135)&
             &*XX(27)-JVS(140)*XX(28)-JVS(145)*XX(29)-JVS(154)*XX(31)-JVS(165)*XX(32)-JVS(184)*XX(35)-JVS(230)*XX(40)&
             &-JVS(245)*XX(42)-JVS(251)*XX(43)-JVS(363)*XX(51)-JVS(374)*XX(52)-JVS(444)*XX(56)-JVS(454)*XX(57)-JVS(464)&
             &*XX(58)-JVS(495)*XX(60)-JVS(506)*XX(61)-JVS(529)*XX(62)-JVS(570)*XX(64)-JVS(636)*XX(68)-JVS(661)*XX(69)&
             &-JVS(687)*XX(70)-JVS(702)*XX(71)-JVS(744)*XX(72)-JVS(813)*XX(73)-JVS(858)*XX(74)-JVS(921)*XX(75)-JVS(957)&
             &*XX(76)-JVS(984)*XX(77)-JVS(995)*XX(78)-JVS(1033)*XX(79))/(JVS(1087))
  XX(80) = XX(80)
  XX(79) = XX(79)-JVS(1086)*XX(80)
  XX(78) = XX(78)-JVS(1031)*XX(79)-JVS(1085)*XX(80)
  XX(77) = XX(77)-JVS(992)*XX(78)-JVS(1030)*XX(79)-JVS(1084)*XX(80)
  XX(76) = XX(76)-JVS(980)*XX(77)-JVS(991)*XX(78)-JVS(1029)*XX(79)-JVS(1083)*XX(80)
  XX(75) = XX(75)-JVS(952)*XX(76)-JVS(979)*XX(77)-JVS(990)*XX(78)-JVS(1028)*XX(79)-JVS(1082)*XX(80)
  XX(74) = XX(74)-JVS(915)*XX(75)-JVS(951)*XX(76)-JVS(978)*XX(77)-JVS(989)*XX(78)-JVS(1027)*XX(79)-JVS(1081)*XX(80)
  XX(73) = XX(73)-JVS(851)*XX(74)-JVS(914)*XX(75)-JVS(950)*XX(76)-JVS(977)*XX(77)-JVS(988)*XX(78)-JVS(1026)*XX(79)&
             &-JVS(1080)*XX(80)
  XX(72) = XX(72)-JVS(805)*XX(73)-JVS(850)*XX(74)-JVS(913)*XX(75)-JVS(949)*XX(76)-JVS(976)*XX(77)-JVS(987)*XX(78)&
             &-JVS(1025)*XX(79)-JVS(1079)*XX(80)
  XX(71) = XX(71)-JVS(735)*XX(72)-JVS(804)*XX(73)-JVS(849)*XX(74)-JVS(912)*XX(75)-JVS(948)*XX(76)-JVS(1024)*XX(79)&
             &-JVS(1078)*XX(80)
  XX(70) = XX(70)-JVS(692)*XX(71)-JVS(734)*XX(72)-JVS(803)*XX(73)-JVS(848)*XX(74)-JVS(911)*XX(75)-JVS(947)*XX(76)&
             &-JVS(1023)*XX(79)-JVS(1077)*XX(80)
  XX(69) = XX(69)-JVS(676)*XX(70)-JVS(733)*XX(72)-JVS(802)*XX(73)-JVS(847)*XX(74)-JVS(910)*XX(75)-JVS(946)*XX(76)&
             &-JVS(1022)*XX(79)-JVS(1076)*XX(80)
  XX(68) = XX(68)-JVS(675)*XX(70)-JVS(732)*XX(72)-JVS(801)*XX(73)-JVS(846)*XX(74)-JVS(909)*XX(75)-JVS(945)*XX(76)&
             &-JVS(1021)*XX(79)-JVS(1075)*XX(80)
  XX(67) = XX(67)-JVS(623)*XX(68)-JVS(649)*XX(69)-JVS(691)*XX(71)-JVS(731)*XX(72)-JVS(800)*XX(73)-JVS(845)*XX(74)&
             &-JVS(908)*XX(75)-JVS(944)*XX(76)-JVS(1020)*XX(79)-JVS(1074)*XX(80)
  XX(66) = XX(66)-JVS(602)*XX(67)-JVS(622)*XX(68)-JVS(648)*XX(69)-JVS(690)*XX(71)-JVS(730)*XX(72)-JVS(799)*XX(73)&
             &-JVS(844)*XX(74)-JVS(907)*XX(75)-JVS(943)*XX(76)-JVS(1019)*XX(79)-JVS(1073)*XX(80)
  XX(65) = XX(65)-JVS(586)*XX(66)-JVS(601)*XX(67)-JVS(621)*XX(68)-JVS(647)*XX(69)-JVS(674)*XX(70)-JVS(729)*XX(72)&
             &-JVS(798)*XX(73)-JVS(843)*XX(74)-JVS(906)*XX(75)-JVS(942)*XX(76)-JVS(975)*XX(77)-JVS(1018)*XX(79)-JVS(1072)&
             &*XX(80)
  XX(64) = XX(64)-JVS(728)*XX(72)-JVS(797)*XX(73)-JVS(842)*XX(74)-JVS(905)*XX(75)-JVS(941)*XX(76)-JVS(1071)*XX(80)
  XX(63) = XX(63)-JVS(553)*XX(64)-JVS(585)*XX(66)-JVS(600)*XX(67)-JVS(620)*XX(68)-JVS(646)*XX(69)-JVS(673)*XX(70)&
             &-JVS(727)*XX(72)-JVS(796)*XX(73)-JVS(841)*XX(74)-JVS(904)*XX(75)-JVS(940)*XX(76)-JVS(974)*XX(77)-JVS(1017)&
             &*XX(79)-JVS(1070)*XX(80)
  XX(62) = XX(62)-JVS(672)*XX(70)-JVS(726)*XX(72)-JVS(795)*XX(73)-JVS(840)*XX(74)-JVS(903)*XX(75)-JVS(939)*XX(76)&
             &-JVS(1016)*XX(79)-JVS(1069)*XX(80)
  XX(61) = XX(61)-JVS(516)*XX(62)-JVS(552)*XX(64)-JVS(619)*XX(68)-JVS(645)*XX(69)-JVS(725)*XX(72)-JVS(794)*XX(73)&
             &-JVS(839)*XX(74)-JVS(902)*XX(75)-JVS(938)*XX(76)-JVS(973)*XX(77)-JVS(986)*XX(78)-JVS(1068)*XX(80)
  XX(60) = XX(60)-JVS(551)*XX(64)-JVS(618)*XX(68)-JVS(644)*XX(69)-JVS(689)*XX(71)-JVS(724)*XX(72)-JVS(793)*XX(73)&
             &-JVS(901)*XX(75)-JVS(1067)*XX(80)
  XX(59) = XX(59)-JVS(482)*XX(60)-JVS(550)*XX(64)-JVS(617)*XX(68)-JVS(643)*XX(69)-JVS(671)*XX(70)-JVS(723)*XX(72)&
             &-JVS(792)*XX(73)-JVS(838)*XX(74)-JVS(900)*XX(75)-JVS(937)*XX(76)-JVS(1015)*XX(79)-JVS(1066)*XX(80)
  XX(58) = XX(58)-JVS(481)*XX(60)-JVS(497)*XX(61)-JVS(670)*XX(70)-JVS(722)*XX(72)-JVS(791)*XX(73)-JVS(837)*XX(74)&
             &-JVS(899)*XX(75)-JVS(936)*XX(76)-JVS(1014)*XX(79)-JVS(1065)*XX(80)
  XX(57) = XX(57)-JVS(480)*XX(60)-JVS(496)*XX(61)-JVS(669)*XX(70)-JVS(721)*XX(72)-JVS(790)*XX(73)-JVS(836)*XX(74)&
             &-JVS(898)*XX(75)-JVS(935)*XX(76)-JVS(1013)*XX(79)-JVS(1064)*XX(80)
  XX(56) = XX(56)-JVS(789)*XX(73)-JVS(835)*XX(74)-JVS(897)*XX(75)-JVS(1063)*XX(80)
  XX(55) = XX(55)-JVS(422)*XX(56)-JVS(479)*XX(60)-JVS(515)*XX(62)-JVS(549)*XX(64)-JVS(616)*XX(68)-JVS(642)*XX(69)&
             &-JVS(668)*XX(70)-JVS(720)*XX(72)-JVS(788)*XX(73)-JVS(834)*XX(74)-JVS(896)*XX(75)-JVS(934)*XX(76)-JVS(1012)&
             &*XX(79)-JVS(1062)*XX(80)
  XX(54) = XX(54)-JVS(421)*XX(56)-JVS(478)*XX(60)-JVS(572)*XX(65)-JVS(667)*XX(70)-JVS(719)*XX(72)-JVS(787)*XX(73)&
             &-JVS(833)*XX(74)-JVS(895)*XX(75)-JVS(933)*XX(76)-JVS(972)*XX(77)-JVS(1011)*XX(79)-JVS(1061)*XX(80)
  XX(53) = XX(53)-JVS(420)*XX(56)-JVS(477)*XX(60)-JVS(514)*XX(62)-JVS(548)*XX(64)-JVS(615)*XX(68)-JVS(641)*XX(69)&
             &-JVS(666)*XX(70)-JVS(718)*XX(72)-JVS(786)*XX(73)-JVS(832)*XX(74)-JVS(894)*XX(75)-JVS(932)*XX(76)-JVS(1010)&
             &*XX(79)-JVS(1060)*XX(80)
  XX(52) = XX(52)-JVS(419)*XX(56)-JVS(665)*XX(70)-JVS(717)*XX(72)-JVS(785)*XX(73)-JVS(831)*XX(74)-JVS(893)*XX(75)&
             &-JVS(931)*XX(76)-JVS(1009)*XX(79)-JVS(1059)*XX(80)
  XX(51) = XX(51)-JVS(716)*XX(72)-JVS(784)*XX(73)-JVS(830)*XX(74)-JVS(892)*XX(75)-JVS(1058)*XX(80)
  XX(50) = XX(50)-JVS(346)*XX(51)-JVS(418)*XX(56)-JVS(532)*XX(63)-JVS(571)*XX(65)-JVS(715)*XX(72)-JVS(783)*XX(73)&
             &-JVS(891)*XX(75)-JVS(930)*XX(76)-JVS(971)*XX(77)-JVS(1008)*XX(79)-JVS(1057)*XX(80)
  XX(49) = XX(49)-JVS(614)*XX(68)-JVS(640)*XX(69)-JVS(714)*XX(72)-JVS(782)*XX(73)
  XX(48) = XX(48)-JVS(345)*XX(51)-JVS(417)*XX(56)-JVS(713)*XX(72)-JVS(781)*XX(73)-JVS(890)*XX(75)-JVS(970)*XX(77)&
             &-JVS(1007)*XX(79)-JVS(1056)*XX(80)
  XX(47) = XX(47)-JVS(308)*XX(49)-JVS(416)*XX(56)-JVS(547)*XX(64)-JVS(664)*XX(70)-JVS(712)*XX(72)-JVS(829)*XX(74)&
             &-JVS(889)*XX(75)-JVS(929)*XX(76)-JVS(1006)*XX(79)-JVS(1055)*XX(80)
  XX(46) = XX(46)-JVS(307)*XX(49)-JVS(415)*XX(56)-JVS(546)*XX(64)-JVS(663)*XX(70)-JVS(711)*XX(72)-JVS(828)*XX(74)&
             &-JVS(888)*XX(75)-JVS(928)*XX(76)-JVS(1005)*XX(79)-JVS(1054)*XX(80)
  XX(45) = XX(45)-JVS(414)*XX(56)-JVS(545)*XX(64)-JVS(662)*XX(70)-JVS(710)*XX(72)-JVS(827)*XX(74)-JVS(887)*XX(75)&
             &-JVS(927)*XX(76)-JVS(1004)*XX(79)-JVS(1053)*XX(80)
  XX(44) = XX(44)-JVS(413)*XX(56)-JVS(780)*XX(73)-JVS(826)*XX(74)-JVS(886)*XX(75)-JVS(1052)*XX(80)
  XX(43) = XX(43)-JVS(344)*XX(51)-JVS(412)*XX(56)-JVS(476)*XX(60)-JVS(513)*XX(62)-JVS(709)*XX(72)-JVS(779)*XX(73)&
             &-JVS(825)*XX(74)-JVS(885)*XX(75)-JVS(969)*XX(77)-JVS(985)*XX(78)-JVS(1051)*XX(80)
  XX(42) = XX(42)-JVS(343)*XX(51)-JVS(411)*XX(56)-JVS(708)*XX(72)-JVS(778)*XX(73)-JVS(884)*XX(75)-JVS(1003)*XX(79)&
             &-JVS(1050)*XX(80)
  XX(41) = XX(41)-JVS(293)*XX(48)-JVS(330)*XX(50)-JVS(410)*XX(56)-JVS(512)*XX(62)-JVS(531)*XX(63)-JVS(584)*XX(66)&
             &-JVS(599)*XX(67)-JVS(639)*XX(69)-JVS(707)*XX(72)-JVS(777)*XX(73)-JVS(824)*XX(74)-JVS(883)*XX(75)-JVS(926)&
             &*XX(76)-JVS(968)*XX(77)-JVS(1002)*XX(79)-JVS(1049)*XX(80)
  XX(40) = XX(40)-JVS(365)*XX(52)-JVS(475)*XX(60)-JVS(511)*XX(62)-JVS(776)*XX(73)-JVS(823)*XX(74)-JVS(882)*XX(75)&
             &-JVS(967)*XX(77)-JVS(1001)*XX(79)-JVS(1048)*XX(80)
  XX(39) = XX(39)-JVS(237)*XX(42)-JVS(292)*XX(48)-JVS(329)*XX(50)-JVS(342)*XX(51)-JVS(409)*XX(56)-JVS(706)*XX(72)&
             &-JVS(775)*XX(73)-JVS(822)*XX(74)-JVS(881)*XX(75)-JVS(925)*XX(76)-JVS(966)*XX(77)-JVS(1047)*XX(80)
  XX(38) = XX(38)-JVS(291)*XX(48)-JVS(383)*XX(54)-JVS(408)*XX(56)-JVS(530)*XX(63)-JVS(583)*XX(66)-JVS(598)*XX(67)&
             &-JVS(613)*XX(68)-JVS(638)*XX(69)-JVS(774)*XX(73)-JVS(821)*XX(74)-JVS(880)*XX(75)-JVS(965)*XX(77)
  XX(37) = XX(37)-JVS(773)*XX(73)-JVS(879)*XX(75)
  XX(36) = XX(36)-JVS(193)*XX(37)-JVS(306)*XX(49)-JVS(465)*XX(59)-JVS(544)*XX(64)-JVS(582)*XX(66)-JVS(597)*XX(67)&
             &-JVS(612)*XX(68)-JVS(637)*XX(69)-JVS(772)*XX(73)-JVS(820)*XX(74)-JVS(878)*XX(75)-JVS(964)*XX(77)
  XX(35) = XX(35)-JVS(407)*XX(56)-JVS(510)*XX(62)-JVS(705)*XX(72)-JVS(771)*XX(73)-JVS(819)*XX(74)-JVS(1046)*XX(80)
  XX(34) = XX(34)-JVS(192)*XX(37)-JVS(341)*XX(51)-JVS(406)*XX(56)-JVS(770)*XX(73)-JVS(877)*XX(75)-JVS(1045)*XX(80)
  XX(33) = XX(33)-JVS(405)*XX(56)-JVS(509)*XX(62)-JVS(704)*XX(72)-JVS(769)*XX(73)-JVS(876)*XX(75)-JVS(924)*XX(76)
  XX(32) = XX(32)-JVS(768)*XX(73)-JVS(818)*XX(74)-JVS(1044)*XX(80)
  XX(31) = XX(31)-JVS(224)*XX(40)-JVS(364)*XX(52)-JVS(767)*XX(73)-JVS(875)*XX(75)-JVS(963)*XX(77)-JVS(1000)*XX(79)&
             &-JVS(1043)*XX(80)
  XX(30) = XX(30)-JVS(191)*XX(37)-JVS(265)*XX(45)-JVS(404)*XX(56)-JVS(581)*XX(66)-JVS(596)*XX(67)-JVS(766)*XX(73)&
             &-JVS(817)*XX(74)-JVS(874)*XX(75)-JVS(962)*XX(77)
  XX(29) = XX(29)-JVS(223)*XX(40)-JVS(455)*XX(58)-JVS(765)*XX(73)-JVS(873)*XX(75)-JVS(961)*XX(77)-JVS(999)*XX(79)&
             &-JVS(1042)*XX(80)
  XX(28) = XX(28)-JVS(222)*XX(40)-JVS(445)*XX(57)-JVS(764)*XX(73)-JVS(872)*XX(75)-JVS(960)*XX(77)-JVS(998)*XX(79)&
             &-JVS(1041)*XX(80)
  XX(27) = XX(27)-JVS(221)*XX(40)-JVS(474)*XX(60)-JVS(871)*XX(75)-JVS(1040)*XX(80)
  XX(26) = XX(26)-JVS(816)*XX(74)-JVS(959)*XX(77)-JVS(997)*XX(79)-JVS(1039)*XX(80)
  XX(25) = XX(25)-JVS(340)*XX(51)-JVS(403)*XX(56)-JVS(703)*XX(72)-JVS(763)*XX(73)-JVS(870)*XX(75)
  XX(24) = XX(24)-JVS(762)*XX(73)-JVS(869)*XX(75)
  XX(23) = XX(23)-JVS(761)*XX(73)-JVS(996)*XX(79)-JVS(1038)*XX(80)
  XX(22) = XX(22)-JVS(339)*XX(51)-JVS(760)*XX(73)-JVS(868)*XX(75)
  XX(21) = XX(21)-JVS(93)*XX(22)-JVS(236)*XX(42)-JVS(759)*XX(73)-JVS(1037)*XX(80)
  XX(20) = XX(20)-JVS(402)*XX(56)-JVS(758)*XX(73)-JVS(867)*XX(75)-JVS(923)*XX(76)
  XX(19) = XX(19)-JVS(757)*XX(73)-JVS(815)*XX(74)-JVS(866)*XX(75)-JVS(1036)*XX(80)
  XX(18) = XX(18)-JVS(190)*XX(37)-JVS(401)*XX(56)-JVS(756)*XX(73)-JVS(1035)*XX(80)
  XX(17) = XX(17)-JVS(611)*XX(68)-JVS(755)*XX(73)
  XX(16) = XX(16)-JVS(305)*XX(49)-JVS(543)*XX(64)-JVS(754)*XX(73)-JVS(865)*XX(75)
  XX(15) = XX(15)-JVS(753)*XX(73)-JVS(922)*XX(76)
  XX(14) = XX(14)-JVS(189)*XX(37)-JVS(252)*XX(44)-JVS(400)*XX(56)-JVS(542)*XX(64)-JVS(688)*XX(71)-JVS(752)*XX(73)&
             &-JVS(864)*XX(75)
  XX(13) = XX(13)-JVS(235)*XX(42)-JVS(751)*XX(73)
  XX(12) = XX(12)-JVS(117)*XX(25)-JVS(375)*XX(53)-JVS(541)*XX(64)-JVS(750)*XX(73)-JVS(863)*XX(75)
  XX(11) = XX(11)-JVS(814)*XX(74)-JVS(1034)*XX(80)
  XX(10) = XX(10)-JVS(141)*XX(29)-JVS(508)*XX(62)-JVS(749)*XX(73)-JVS(862)*XX(75)
  XX(9) = XX(9)-JVS(136)*XX(28)-JVS(507)*XX(62)-JVS(748)*XX(73)-JVS(861)*XX(75)
  XX(8) = XX(8)-JVS(304)*XX(49)-JVS(392)*XX(55)-JVS(747)*XX(73)-JVS(860)*XX(75)
  XX(7) = XX(7)-JVS(290)*XX(48)
  XX(6) = XX(6)-JVS(123)*XX(26)-JVS(746)*XX(73)-JVS(958)*XX(77)
  XX(5) = XX(5)-JVS(745)*XX(73)-JVS(859)*XX(75)
  XX(4) = XX(4)
  XX(3) = XX(3)
  XX(2) = XX(2)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE racm_mim_LinearAlgebra

