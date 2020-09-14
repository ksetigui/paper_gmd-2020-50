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
! File                 : cbmz_mosaic_4bin_vbs9_LinearAlgebra.f90
! Time                 : Mon Feb 27 11:08:25 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/cbmz_mosaic_4bin_vbs9
! Equation file        : cbmz_mosaic_4bin_vbs9.kpp
! Output root filename : cbmz_mosaic_4bin_vbs9
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbmz_mosaic_4bin_vbs9_LinearAlgebra

  USE cbmz_mosaic_4bin_vbs9_Parameters
  USE cbmz_mosaic_4bin_vbs9_JacobianSP

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
  XX(3) = X(3)/JVS(5)
  XX(4) = X(4)/JVS(6)
  XX(5) = X(5)/JVS(7)
  XX(6) = X(6)/JVS(8)
  XX(7) = X(7)/JVS(9)
  XX(8) = X(8)/JVS(10)
  XX(9) = X(9)/JVS(14)
  XX(10) = X(10)/JVS(22)
  XX(11) = X(11)/JVS(28)
  XX(12) = X(12)/JVS(34)
  XX(13) = X(13)/JVS(37)
  XX(14) = X(14)/JVS(43)
  XX(15) = X(15)/JVS(47)
  XX(16) = X(16)/JVS(51)
  XX(17) = X(17)/JVS(55)
  XX(18) = X(18)/JVS(59)
  XX(19) = X(19)/JVS(98)
  XX(20) = (X(20)-JVS(99)*XX(19))/(JVS(166))
  XX(21) = (X(21)-JVS(100)*XX(19))/(JVS(168))
  XX(22) = (X(22)-JVS(101)*XX(19))/(JVS(170))
  XX(23) = (X(23)-JVS(102)*XX(19))/(JVS(172))
  XX(24) = (X(24)-JVS(103)*XX(19))/(JVS(174))
  XX(25) = (X(25)-JVS(104)*XX(19))/(JVS(176))
  XX(26) = (X(26)-JVS(105)*XX(19))/(JVS(178))
  XX(27) = (X(27)-JVS(106)*XX(19))/(JVS(180))
  XX(28) = (X(28)-JVS(107)*XX(19))/(JVS(182))
  XX(29) = (X(29)-JVS(108)*XX(19))/(JVS(184))
  XX(30) = (X(30)-JVS(109)*XX(19))/(JVS(186))
  XX(31) = (X(31)-JVS(110)*XX(19))/(JVS(188))
  XX(32) = (X(32)-JVS(111)*XX(19))/(JVS(190))
  XX(33) = (X(33)-JVS(112)*XX(19))/(JVS(192))
  XX(34) = (X(34)-JVS(113)*XX(19))/(JVS(194))
  XX(35) = (X(35)-JVS(114)*XX(19))/(JVS(196))
  XX(36) = X(36)/JVS(198)
  XX(37) = (X(37)-JVS(2)*XX(1)-JVS(60)*XX(18))/(JVS(200))
  XX(38) = (X(38)-JVS(115)*XX(19))/(JVS(202))
  XX(39) = (X(39)-JVS(116)*XX(19))/(JVS(205))
  XX(40) = (X(40)-JVS(117)*XX(19))/(JVS(208))
  XX(41) = (X(41)-JVS(118)*XX(19))/(JVS(210))
  XX(42) = (X(42)-JVS(119)*XX(19))/(JVS(212))
  XX(43) = (X(43)-JVS(120)*XX(19))/(JVS(214))
  XX(44) = (X(44)-JVS(121)*XX(19))/(JVS(216))
  XX(45) = (X(45)-JVS(122)*XX(19))/(JVS(218))
  XX(46) = (X(46)-JVS(123)*XX(19))/(JVS(220))
  XX(47) = (X(47)-JVS(124)*XX(19))/(JVS(224))
  XX(48) = (X(48)-JVS(61)*XX(18))/(JVS(226))
  XX(49) = (X(49)-JVS(125)*XX(19))/(JVS(229))
  XX(50) = (X(50)-JVS(126)*XX(19))/(JVS(235))
  XX(51) = (X(51)-JVS(127)*XX(19))/(JVS(239))
  XX(52) = (X(52)-JVS(128)*XX(19))/(JVS(245))
  XX(53) = (X(53)-JVS(129)*XX(19))/(JVS(249))
  XX(54) = (X(54)-JVS(130)*XX(19))/(JVS(255))
  XX(55) = (X(55)-JVS(131)*XX(19))/(JVS(259))
  XX(56) = (X(56)-JVS(132)*XX(19))/(JVS(265))
  XX(57) = (X(57)-JVS(133)*XX(19))/(JVS(269))
  XX(58) = (X(58)-JVS(134)*XX(19))/(JVS(275))
  XX(59) = (X(59)-JVS(135)*XX(19))/(JVS(279))
  XX(60) = (X(60)-JVS(136)*XX(19)-JVS(206)*XX(39))/(JVS(283))
  XX(61) = (X(61)-JVS(137)*XX(19))/(JVS(288))
  XX(62) = (X(62)-JVS(138)*XX(19)-JVS(289)*XX(61))/(JVS(295))
  XX(63) = (X(63)-JVS(139)*XX(19))/(JVS(297))
  XX(64) = (X(64)-JVS(140)*XX(19))/(JVS(300))
  XX(65) = (X(65)-JVS(141)*XX(19))/(JVS(303))
  XX(66) = (X(66)-JVS(142)*XX(19))/(JVS(305))
  XX(67) = (X(67)-JVS(143)*XX(19))/(JVS(307))
  XX(68) = (X(68)-JVS(144)*XX(19))/(JVS(309))
  XX(69) = (X(69)-JVS(145)*XX(19))/(JVS(311))
  XX(70) = (X(70)-JVS(146)*XX(19))/(JVS(313))
  XX(71) = (X(71)-JVS(147)*XX(19))/(JVS(315))
  XX(72) = (X(72)-JVS(148)*XX(19))/(JVS(319))
  XX(73) = (X(73)-JVS(149)*XX(19))/(JVS(322))
  XX(74) = (X(74)-JVS(150)*XX(19))/(JVS(328))
  XX(75) = (X(75)-JVS(151)*XX(19))/(JVS(332))
  XX(76) = (X(76)-JVS(152)*XX(19))/(JVS(338))
  XX(77) = (X(77)-JVS(153)*XX(19))/(JVS(342))
  XX(78) = (X(78)-JVS(154)*XX(19))/(JVS(348))
  XX(79) = (X(79)-JVS(155)*XX(19))/(JVS(352))
  XX(80) = (X(80)-JVS(156)*XX(19))/(JVS(358))
  XX(81) = (X(81)-JVS(157)*XX(19))/(JVS(362))
  XX(82) = (X(82)-JVS(158)*XX(19))/(JVS(368))
  XX(83) = (X(83)-JVS(159)*XX(19))/(JVS(372))
  XX(84) = (X(84)-JVS(160)*XX(19)-JVS(301)*XX(64))/(JVS(376))
  XX(85) = (X(85)-JVS(161)*XX(19))/(JVS(381))
  XX(86) = (X(86)-JVS(162)*XX(19)-JVS(382)*XX(85))/(JVS(388))
  XX(87) = X(87)/JVS(390)
  XX(88) = (X(88)-JVS(62)*XX(18))/(JVS(393))
  XX(89) = (X(89)-JVS(23)*XX(10)-JVS(29)*XX(11)-JVS(63)*XX(18))/(JVS(396))
  XX(90) = X(90)/JVS(398)
  XX(91) = (X(91)-JVS(24)*XX(10)-JVS(30)*XX(11)-JVS(64)*XX(18))/(JVS(401))
  XX(92) = X(92)/JVS(403)
  XX(93) = (X(93)-JVS(65)*XX(18))/(JVS(408))
  XX(94) = (X(94)-JVS(44)*XX(14)-JVS(48)*XX(15)-JVS(163)*XX(19))/(JVS(413))
  XX(95) = (X(95)-JVS(52)*XX(16)-JVS(56)*XX(17)-JVS(164)*XX(19))/(JVS(417))
  XX(96) = (X(96)-JVS(66)*XX(18))/(JVS(421))
  XX(97) = X(97)/JVS(427)
  XX(98) = (X(98)-JVS(67)*XX(18))/(JVS(430))
  XX(99) = X(99)/JVS(436)
  XX(100) = (X(100)-JVS(68)*XX(18))/(JVS(448))
  XX(101) = (X(101)-JVS(69)*XX(18))/(JVS(452))
  XX(102) = (X(102)-JVS(70)*XX(18))/(JVS(456))
  XX(103) = (X(103)-JVS(11)*XX(8)-JVS(71)*XX(18))/(JVS(461))
  XX(104) = (X(104)-JVS(72)*XX(18))/(JVS(464))
  XX(105) = (X(105)-JVS(25)*XX(10)-JVS(31)*XX(11)-JVS(73)*XX(18)-JVS(404)*XX(92))/(JVS(472))
  XX(106) = X(106)/JVS(476)
  XX(107) = (X(107)-JVS(74)*XX(18))/(JVS(483))
  XX(108) = (X(108)-JVS(75)*XX(18))/(JVS(493))
  XX(109) = (X(109)-JVS(35)*XX(12)-JVS(76)*XX(18)-JVS(437)*XX(99))/(JVS(509))
  XX(110) = (X(110)-JVS(77)*XX(18)-JVS(494)*XX(108))/(JVS(526))
  XX(111) = (X(111)-JVS(510)*XX(109))/(JVS(531))
  XX(112) = (X(112)-JVS(495)*XX(108))/(JVS(536))
  XX(113) = (X(113)-JVS(511)*XX(109))/(JVS(541))
  XX(114) = (X(114)-JVS(12)*XX(8)-JVS(15)*XX(9)-JVS(38)*XX(13)-JVS(78)*XX(18)-JVS(409)*XX(93)-JVS(431)*XX(98)-JVS(438)&
              &*XX(99)-JVS(465)*XX(104)-JVS(496)*XX(108)-JVS(512)*XX(109))/(JVS(546))
  XX(115) = (X(115)-JVS(16)*XX(9)-JVS(79)*XX(18)-JVS(497)*XX(108)-JVS(532)*XX(111)-JVS(542)*XX(113))/(JVS(550))
  XX(116) = (X(116)-JVS(80)*XX(18)-JVS(484)*XX(107)-JVS(498)*XX(108))/(JVS(562))
  XX(117) = X(117)/JVS(580)
  XX(118) = (X(118)-JVS(81)*XX(18))/(JVS(593))
  XX(119) = (X(119)-JVS(17)*XX(9)-JVS(39)*XX(13)-JVS(82)*XX(18)-JVS(410)*XX(93)-JVS(432)*XX(98)-JVS(439)*XX(99)-JVS(466)&
              &*XX(104)-JVS(499)*XX(108)-JVS(513)*XX(109)-JVS(563)*XX(116)-JVS(581)*XX(117)-JVS(594)*XX(118))/(JVS(605))
  XX(120) = (X(120)-JVS(440)*XX(99)-JVS(514)*XX(109)-JVS(564)*XX(116))/(JVS(611))
  XX(121) = (X(121)-JVS(83)*XX(18)-JVS(485)*XX(107)-JVS(500)*XX(108)-JVS(582)*XX(117))/(JVS(623))
  XX(122) = (X(122)-JVS(433)*XX(98)-JVS(453)*XX(101))/(JVS(636))
  XX(123) = (X(123)-JVS(84)*XX(18)-JVS(486)*XX(107)-JVS(501)*XX(108))/(JVS(656))
  XX(124) = (X(124)-JVS(449)*XX(100)-JVS(467)*XX(104)-JVS(565)*XX(116))/(JVS(672))
  XX(125) = (X(125)-JVS(18)*XX(9)-JVS(85)*XX(18)-JVS(487)*XX(107)-JVS(502)*XX(108)-JVS(515)*XX(109)-JVS(537)*XX(112)&
              &-JVS(566)*XX(116)-JVS(583)*XX(117)-JVS(595)*XX(118)-JVS(624)*XX(121)-JVS(657)*XX(123)-JVS(673)*XX(124))&
              &/(JVS(687))
  XX(126) = (X(126)-JVS(86)*XX(18)-JVS(441)*XX(99)-JVS(516)*XX(109)-JVS(584)*XX(117)-JVS(596)*XX(118)-JVS(612)*XX(120)&
              &-JVS(637)*XX(122)-JVS(658)*XX(123)-JVS(674)*XX(124))/(JVS(701))
  XX(127) = (X(127)-JVS(13)*XX(8)-JVS(19)*XX(9)-JVS(40)*XX(13)-JVS(45)*XX(14)-JVS(49)*XX(15)-JVS(53)*XX(16)-JVS(57)&
              &*XX(17)-JVS(87)*XX(18)-JVS(199)*XX(36)-JVS(411)*XX(93)-JVS(414)*XX(94)-JVS(418)*XX(95)-JVS(434)*XX(98)&
              &-JVS(442)*XX(99)-JVS(462)*XX(103)-JVS(468)*XX(104)-JVS(477)*XX(106)-JVS(503)*XX(108)-JVS(517)*XX(109)&
              &-JVS(527)*XX(110)-JVS(547)*XX(114)-JVS(551)*XX(115)-JVS(567)*XX(116)-JVS(585)*XX(117)-JVS(597)*XX(118)&
              &-JVS(606)*XX(119)-JVS(613)*XX(120)-JVS(625)*XX(121)-JVS(638)*XX(122)-JVS(659)*XX(123)-JVS(675)*XX(124)&
              &-JVS(688)*XX(125)-JVS(702)*XX(126))/(JVS(718))
  XX(128) = (X(128)-JVS(443)*XX(99)-JVS(518)*XX(109)-JVS(586)*XX(117)-JVS(598)*XX(118)-JVS(639)*XX(122)-JVS(660)*XX(123)&
              &-JVS(676)*XX(124)-JVS(703)*XX(126))/(JVS(735))
  XX(129) = (X(129)-JVS(88)*XX(18)-JVS(391)*XX(87)-JVS(399)*XX(90)-JVS(405)*XX(92)-JVS(422)*XX(96)-JVS(457)*XX(102)&
              &-JVS(478)*XX(106)-JVS(488)*XX(107)-JVS(704)*XX(126)-JVS(719)*XX(127)-JVS(736)*XX(128))/(JVS(768))
  XX(130) = (X(130)-JVS(20)*XX(9)-JVS(392)*XX(87)-JVS(677)*XX(124)-JVS(720)*XX(127)-JVS(737)*XX(128)-JVS(769)*XX(129))&
              &/(JVS(789))
  XX(131) = (X(131)-JVS(3)*XX(1)-JVS(26)*XX(10)-JVS(32)*XX(11)-JVS(36)*XX(12)-JVS(41)*XX(13)-JVS(46)*XX(14)-JVS(50)&
              &*XX(15)-JVS(54)*XX(16)-JVS(58)*XX(17)-JVS(89)*XX(18)-JVS(165)*XX(19)-JVS(167)*XX(20)-JVS(169)*XX(21)-JVS(171)&
              &*XX(22)-JVS(173)*XX(23)-JVS(175)*XX(24)-JVS(177)*XX(25)-JVS(179)*XX(26)-JVS(181)*XX(27)-JVS(183)*XX(28)&
              &-JVS(185)*XX(29)-JVS(187)*XX(30)-JVS(189)*XX(31)-JVS(191)*XX(32)-JVS(193)*XX(33)-JVS(195)*XX(34)-JVS(197)&
              &*XX(35)-JVS(201)*XX(37)-JVS(203)*XX(38)-JVS(207)*XX(39)-JVS(209)*XX(40)-JVS(211)*XX(41)-JVS(213)*XX(42)&
              &-JVS(215)*XX(43)-JVS(217)*XX(44)-JVS(219)*XX(45)-JVS(221)*XX(46)-JVS(225)*XX(47)-JVS(227)*XX(48)-JVS(230)&
              &*XX(49)-JVS(236)*XX(50)-JVS(240)*XX(51)-JVS(246)*XX(52)-JVS(250)*XX(53)-JVS(256)*XX(54)-JVS(260)*XX(55)&
              &-JVS(266)*XX(56)-JVS(270)*XX(57)-JVS(276)*XX(58)-JVS(280)*XX(59)-JVS(284)*XX(60)-JVS(290)*XX(61)-JVS(296)&
              &*XX(62)-JVS(298)*XX(63)-JVS(302)*XX(64)-JVS(304)*XX(65)-JVS(306)*XX(66)-JVS(308)*XX(67)-JVS(310)*XX(68)&
              &-JVS(312)*XX(69)-JVS(314)*XX(70)-JVS(316)*XX(71)-JVS(320)*XX(72)-JVS(323)*XX(73)-JVS(329)*XX(74)-JVS(333)&
              &*XX(75)-JVS(339)*XX(76)-JVS(343)*XX(77)-JVS(349)*XX(78)-JVS(353)*XX(79)-JVS(359)*XX(80)-JVS(363)*XX(81)&
              &-JVS(369)*XX(82)-JVS(373)*XX(83)-JVS(377)*XX(84)-JVS(383)*XX(85)-JVS(389)*XX(86)-JVS(394)*XX(88)-JVS(397)&
              &*XX(89)-JVS(402)*XX(91)-JVS(406)*XX(92)-JVS(412)*XX(93)-JVS(415)*XX(94)-JVS(419)*XX(95)-JVS(423)*XX(96)&
              &-JVS(428)*XX(97)-JVS(435)*XX(98)-JVS(444)*XX(99)-JVS(450)*XX(100)-JVS(454)*XX(101)-JVS(458)*XX(102)-JVS(463)&
              &*XX(103)-JVS(469)*XX(104)-JVS(473)*XX(105)-JVS(489)*XX(107)-JVS(504)*XX(108)-JVS(519)*XX(109)-JVS(528)&
              &*XX(110)-JVS(538)*XX(112)-JVS(543)*XX(113)-JVS(548)*XX(114)-JVS(552)*XX(115)-JVS(568)*XX(116)-JVS(587)&
              &*XX(117)-JVS(599)*XX(118)-JVS(607)*XX(119)-JVS(614)*XX(120)-JVS(626)*XX(121)-JVS(640)*XX(122)-JVS(661)&
              &*XX(123)-JVS(678)*XX(124)-JVS(689)*XX(125)-JVS(705)*XX(126)-JVS(721)*XX(127)-JVS(738)*XX(128)-JVS(770)&
              &*XX(129)-JVS(790)*XX(130))/(JVS(904))
  XX(132) = (X(132)-JVS(21)*XX(9)-JVS(90)*XX(18)-JVS(395)*XX(88)-JVS(424)*XX(96)-JVS(451)*XX(100)-JVS(455)*XX(101)&
              &-JVS(459)*XX(102)-JVS(490)*XX(107)-JVS(520)*XX(109)-JVS(533)*XX(111)-JVS(539)*XX(112)-JVS(544)*XX(113)&
              &-JVS(569)*XX(116)-JVS(588)*XX(117)-JVS(600)*XX(118)-JVS(615)*XX(120)-JVS(627)*XX(121)-JVS(641)*XX(122)&
              &-JVS(662)*XX(123)-JVS(679)*XX(124)-JVS(690)*XX(125)-JVS(706)*XX(126)-JVS(722)*XX(127)-JVS(739)*XX(128)&
              &-JVS(771)*XX(129)-JVS(791)*XX(130)-JVS(905)*XX(131))/(JVS(944))
  XX(133) = (X(133)-JVS(91)*XX(18)-JVS(429)*XX(97)-JVS(445)*XX(99)-JVS(460)*XX(102)-JVS(479)*XX(106)-JVS(505)*XX(108)&
              &-JVS(521)*XX(109)-JVS(529)*XX(110)-JVS(534)*XX(111)-JVS(540)*XX(112)-JVS(545)*XX(113)-JVS(570)*XX(116)&
              &-JVS(589)*XX(117)-JVS(601)*XX(118)-JVS(616)*XX(120)-JVS(628)*XX(121)-JVS(642)*XX(122)-JVS(663)*XX(123)&
              &-JVS(680)*XX(124)-JVS(691)*XX(125)-JVS(707)*XX(126)-JVS(723)*XX(127)-JVS(740)*XX(128)-JVS(772)*XX(129)&
              &-JVS(792)*XX(130)-JVS(906)*XX(131)-JVS(945)*XX(132))/(JVS(970))
  XX(134) = (X(134)-JVS(571)*XX(116)-JVS(602)*XX(118)-JVS(629)*XX(121)-JVS(681)*XX(124)-JVS(773)*XX(129)-JVS(793)&
              &*XX(130)-JVS(907)*XX(131)-JVS(946)*XX(132)-JVS(971)*XX(133))/(JVS(985))
  XX(135) = (X(135)-JVS(92)*XX(18)-JVS(446)*XX(99)-JVS(522)*XX(109)-JVS(590)*XX(117)-JVS(603)*XX(118)-JVS(630)*XX(121)&
              &-JVS(643)*XX(122)-JVS(664)*XX(123)-JVS(682)*XX(124)-JVS(741)*XX(128)-JVS(774)*XX(129)-JVS(794)*XX(130)&
              &-JVS(908)*XX(131)-JVS(947)*XX(132)-JVS(972)*XX(133)-JVS(986)*XX(134))/(JVS(1000))
  XX(136) = (X(136)-JVS(27)*XX(10)-JVS(33)*XX(11)-JVS(42)*XX(13)-JVS(93)*XX(18)-JVS(400)*XX(90)-JVS(407)*XX(92)-JVS(416)&
              &*XX(94)-JVS(420)*XX(95)-JVS(447)*XX(99)-JVS(474)*XX(105)-JVS(480)*XX(106)-JVS(491)*XX(107)-JVS(506)*XX(108)&
              &-JVS(523)*XX(109)-JVS(530)*XX(110)-JVS(535)*XX(111)-JVS(549)*XX(114)-JVS(553)*XX(115)-JVS(572)*XX(116)&
              &-JVS(591)*XX(117)-JVS(604)*XX(118)-JVS(608)*XX(119)-JVS(617)*XX(120)-JVS(631)*XX(121)-JVS(644)*XX(122)&
              &-JVS(665)*XX(123)-JVS(683)*XX(124)-JVS(692)*XX(125)-JVS(708)*XX(126)-JVS(724)*XX(127)-JVS(742)*XX(128)&
              &-JVS(775)*XX(129)-JVS(795)*XX(130)-JVS(909)*XX(131)-JVS(948)*XX(132)-JVS(973)*XX(133)-JVS(987)*XX(134)&
              &-JVS(1001)*XX(135))/(JVS(1029))
  XX(136) = XX(136)
  XX(135) = XX(135)-JVS(1028)*XX(136)
  XX(134) = XX(134)-JVS(999)*XX(135)-JVS(1027)*XX(136)
  XX(133) = XX(133)-JVS(984)*XX(134)-JVS(998)*XX(135)-JVS(1026)*XX(136)
  XX(132) = XX(132)-JVS(969)*XX(133)-JVS(983)*XX(134)-JVS(997)*XX(135)-JVS(1025)*XX(136)
  XX(131) = XX(131)-JVS(943)*XX(132)-JVS(968)*XX(133)-JVS(982)*XX(134)-JVS(996)*XX(135)-JVS(1024)*XX(136)
  XX(130) = XX(130)-JVS(903)*XX(131)-JVS(942)*XX(132)-JVS(967)*XX(133)-JVS(981)*XX(134)-JVS(995)*XX(135)-JVS(1023)&
              &*XX(136)
  XX(129) = XX(129)-JVS(788)*XX(130)-JVS(902)*XX(131)-JVS(941)*XX(132)-JVS(966)*XX(133)-JVS(980)*XX(134)-JVS(994)&
              &*XX(135)-JVS(1022)*XX(136)
  XX(128) = XX(128)-JVS(767)*XX(129)-JVS(787)*XX(130)-JVS(901)*XX(131)-JVS(940)*XX(132)-JVS(965)*XX(133)-JVS(979)&
              &*XX(134)-JVS(993)*XX(135)-JVS(1021)*XX(136)
  XX(127) = XX(127)-JVS(734)*XX(128)-JVS(766)*XX(129)-JVS(786)*XX(130)-JVS(900)*XX(131)-JVS(939)*XX(132)-JVS(964)&
              &*XX(133)-JVS(978)*XX(134)-JVS(992)*XX(135)-JVS(1020)*XX(136)
  XX(126) = XX(126)-JVS(733)*XX(128)-JVS(765)*XX(129)-JVS(785)*XX(130)-JVS(899)*XX(131)-JVS(938)*XX(132)-JVS(963)&
              &*XX(133)-JVS(977)*XX(134)-JVS(1019)*XX(136)
  XX(125) = XX(125)-JVS(700)*XX(126)-JVS(717)*XX(127)-JVS(732)*XX(128)-JVS(764)*XX(129)-JVS(784)*XX(130)-JVS(898)&
              &*XX(131)-JVS(937)*XX(132)-JVS(962)*XX(133)-JVS(976)*XX(134)-JVS(991)*XX(135)-JVS(1018)*XX(136)
  XX(124) = XX(124)-JVS(763)*XX(129)-JVS(897)*XX(131)-JVS(936)*XX(132)-JVS(961)*XX(133)-JVS(1017)*XX(136)
  XX(123) = XX(123)-JVS(671)*XX(124)-JVS(762)*XX(129)-JVS(783)*XX(130)-JVS(896)*XX(131)-JVS(935)*XX(132)-JVS(1016)&
              &*XX(136)
  XX(122) = XX(122)-JVS(655)*XX(123)-JVS(761)*XX(129)-JVS(895)*XX(131)-JVS(934)*XX(132)-JVS(960)*XX(133)-JVS(1015)&
              &*XX(136)
  XX(121) = XX(121)-JVS(760)*XX(129)-JVS(782)*XX(130)-JVS(894)*XX(131)-JVS(933)*XX(132)-JVS(959)*XX(133)-JVS(1014)&
              &*XX(136)
  XX(120) = XX(120)-JVS(654)*XX(123)-JVS(699)*XX(126)-JVS(731)*XX(128)-JVS(759)*XX(129)-JVS(893)*XX(131)-JVS(932)&
              &*XX(132)-JVS(958)*XX(133)-JVS(1013)*XX(136)
  XX(119) = XX(119)-JVS(610)*XX(120)-JVS(622)*XX(121)-JVS(635)*XX(122)-JVS(653)*XX(123)-JVS(670)*XX(124)-JVS(716)&
              &*XX(127)-JVS(730)*XX(128)-JVS(758)*XX(129)-JVS(781)*XX(130)-JVS(892)*XX(131)-JVS(931)*XX(132)-JVS(957)&
              &*XX(133)-JVS(975)*XX(134)-JVS(1012)*XX(136)
  XX(118) = XX(118)-JVS(669)*XX(124)-JVS(780)*XX(130)-JVS(891)*XX(131)-JVS(974)*XX(134)
  XX(117) = XX(117)-JVS(757)*XX(129)-JVS(930)*XX(132)-JVS(956)*XX(133)-JVS(1011)*XX(136)
  XX(116) = XX(116)-JVS(756)*XX(129)-JVS(890)*XX(131)-JVS(929)*XX(132)-JVS(1010)*XX(136)
  XX(115) = XX(115)-JVS(561)*XX(116)-JVS(579)*XX(117)-JVS(652)*XX(123)-JVS(686)*XX(125)-JVS(698)*XX(126)-JVS(715)&
              &*XX(127)-JVS(729)*XX(128)-JVS(755)*XX(129)-JVS(779)*XX(130)-JVS(889)*XX(131)-JVS(928)*XX(132)-JVS(955)&
              &*XX(133)-JVS(990)*XX(135)-JVS(1009)*XX(136)
  XX(114) = XX(114)-JVS(560)*XX(116)-JVS(578)*XX(117)-JVS(609)*XX(120)-JVS(621)*XX(121)-JVS(634)*XX(122)-JVS(651)&
              &*XX(123)-JVS(668)*XX(124)-JVS(714)*XX(127)-JVS(728)*XX(128)-JVS(778)*XX(130)-JVS(888)*XX(131)-JVS(927)&
              &*XX(132)-JVS(1008)*XX(136)
  XX(113) = XX(113)-JVS(559)*XX(116)-JVS(685)*XX(125)-JVS(697)*XX(126)-JVS(727)*XX(128)-JVS(754)*XX(129)-JVS(887)&
              &*XX(131)-JVS(926)*XX(132)-JVS(954)*XX(133)-JVS(989)*XX(135)
  XX(112) = XX(112)-JVS(558)*XX(116)-JVS(592)*XX(118)-JVS(620)*XX(121)-JVS(650)*XX(123)-JVS(753)*XX(129)-JVS(886)&
              &*XX(131)-JVS(925)*XX(132)-JVS(953)*XX(133)-JVS(988)*XX(135)
  XX(111) = XX(111)-JVS(649)*XX(123)-JVS(684)*XX(125)-JVS(696)*XX(126)-JVS(726)*XX(128)-JVS(752)*XX(129)-JVS(885)&
              &*XX(131)-JVS(924)*XX(132)-JVS(952)*XX(133)
  XX(110) = XX(110)-JVS(557)*XX(116)-JVS(577)*XX(117)-JVS(619)*XX(121)-JVS(648)*XX(123)-JVS(713)*XX(127)-JVS(777)&
              &*XX(130)-JVS(884)*XX(131)-JVS(923)*XX(132)
  XX(109) = XX(109)-JVS(725)*XX(128)-JVS(883)*XX(131)
  XX(108) = XX(108)-JVS(882)*XX(131)-JVS(922)*XX(132)
  XX(107) = XX(107)-JVS(751)*XX(129)-JVS(881)*XX(131)-JVS(1007)*XX(136)
  XX(106) = XX(106)-JVS(712)*XX(127)-JVS(750)*XX(129)-JVS(951)*XX(133)-JVS(1006)*XX(136)
  XX(105) = XX(105)-JVS(482)*XX(107)-JVS(525)*XX(110)-JVS(576)*XX(117)-JVS(695)*XX(126)-JVS(749)*XX(129)-JVS(880)&
              &*XX(131)-JVS(921)*XX(132)-JVS(1005)*XX(136)
  XX(104) = XX(104)-JVS(556)*XX(116)-JVS(879)*XX(131)-JVS(920)*XX(132)
  XX(103) = XX(103)-JVS(492)*XX(108)-JVS(555)*XX(116)-JVS(575)*XX(117)-JVS(647)*XX(123)-JVS(711)*XX(127)-JVS(878)&
              &*XX(131)-JVS(919)*XX(132)
  XX(102) = XX(102)-JVS(748)*XX(129)-JVS(877)*XX(131)-JVS(950)*XX(133)
  XX(101) = XX(101)-JVS(633)*XX(122)-JVS(646)*XX(123)-JVS(876)*XX(131)-JVS(918)*XX(132)
  XX(100) = XX(100)-JVS(554)*XX(116)-JVS(667)*XX(124)-JVS(875)*XX(131)-JVS(917)*XX(132)
  XX(99) = XX(99)-JVS(508)*XX(109)
  XX(98) = XX(98)-JVS(632)*XX(122)-JVS(874)*XX(131)
  XX(97) = XX(97)-JVS(524)*XX(110)-JVS(694)*XX(126)-JVS(747)*XX(129)-JVS(916)*XX(132)-JVS(949)*XX(133)
  XX(96) = XX(96)-JVS(746)*XX(129)-JVS(873)*XX(131)-JVS(915)*XX(132)
  XX(95) = XX(95)-JVS(710)*XX(127)-JVS(872)*XX(131)-JVS(1004)*XX(136)
  XX(94) = XX(94)-JVS(709)*XX(127)-JVS(871)*XX(131)-JVS(1003)*XX(136)
  XX(93) = XX(93)-JVS(666)*XX(124)-JVS(870)*XX(131)
  XX(92) = XX(92)-JVS(693)*XX(126)-JVS(745)*XX(129)
  XX(91) = XX(91)-JVS(426)*XX(97)-JVS(471)*XX(105)-JVS(507)*XX(109)-JVS(574)*XX(117)-JVS(618)*XX(121)-JVS(869)*XX(131)&
             &-JVS(914)*XX(132)
  XX(90) = XX(90)-JVS(481)*XX(107)-JVS(744)*XX(129)-JVS(1002)*XX(136)
  XX(89) = XX(89)-JVS(425)*XX(97)-JVS(470)*XX(105)-JVS(573)*XX(117)-JVS(868)*XX(131)-JVS(913)*XX(132)
  XX(88) = XX(88)-JVS(867)*XX(131)-JVS(912)*XX(132)
  XX(87) = XX(87)-JVS(743)*XX(129)-JVS(776)*XX(130)
  XX(86) = XX(86)-JVS(866)*XX(131)
  XX(85) = XX(85)-JVS(865)*XX(131)
  XX(84) = XX(84)-JVS(380)*XX(85)-JVS(864)*XX(131)
  XX(83) = XX(83)-JVS(375)*XX(84)-JVS(387)*XX(86)-JVS(863)*XX(131)
  XX(82) = XX(82)-JVS(386)*XX(86)-JVS(862)*XX(131)
  XX(81) = XX(81)-JVS(367)*XX(82)-JVS(371)*XX(83)-JVS(861)*XX(131)
  XX(80) = XX(80)-JVS(366)*XX(82)-JVS(860)*XX(131)
  XX(79) = XX(79)-JVS(357)*XX(80)-JVS(361)*XX(81)-JVS(859)*XX(131)
  XX(78) = XX(78)-JVS(356)*XX(80)-JVS(858)*XX(131)
  XX(77) = XX(77)-JVS(347)*XX(78)-JVS(351)*XX(79)-JVS(857)*XX(131)
  XX(76) = XX(76)-JVS(346)*XX(78)-JVS(856)*XX(131)
  XX(75) = XX(75)-JVS(337)*XX(76)-JVS(341)*XX(77)-JVS(855)*XX(131)
  XX(74) = XX(74)-JVS(336)*XX(76)-JVS(854)*XX(131)
  XX(73) = XX(73)-JVS(327)*XX(74)-JVS(331)*XX(75)-JVS(853)*XX(131)
  XX(72) = XX(72)-JVS(326)*XX(74)-JVS(852)*XX(131)
  XX(71) = XX(71)-JVS(318)*XX(72)-JVS(321)*XX(73)-JVS(851)*XX(131)
  XX(70) = XX(70)-JVS(325)*XX(74)-JVS(330)*XX(75)-JVS(850)*XX(131)
  XX(69) = XX(69)-JVS(335)*XX(76)-JVS(340)*XX(77)-JVS(849)*XX(131)
  XX(68) = XX(68)-JVS(345)*XX(78)-JVS(350)*XX(79)-JVS(848)*XX(131)
  XX(67) = XX(67)-JVS(355)*XX(80)-JVS(360)*XX(81)-JVS(847)*XX(131)
  XX(66) = XX(66)-JVS(365)*XX(82)-JVS(370)*XX(83)-JVS(846)*XX(131)
  XX(65) = XX(65)-JVS(374)*XX(84)-JVS(385)*XX(86)-JVS(845)*XX(131)
  XX(64) = XX(64)-JVS(844)*XX(131)
  XX(63) = XX(63)-JVS(299)*XX(64)-JVS(379)*XX(85)-JVS(843)*XX(131)
  XX(62) = XX(62)-JVS(842)*XX(131)
  XX(61) = XX(61)-JVS(841)*XX(131)
  XX(60) = XX(60)-JVS(287)*XX(61)-JVS(840)*XX(131)
  XX(59) = XX(59)-JVS(282)*XX(60)-JVS(294)*XX(62)-JVS(839)*XX(131)
  XX(58) = XX(58)-JVS(293)*XX(62)-JVS(838)*XX(131)
  XX(57) = XX(57)-JVS(274)*XX(58)-JVS(278)*XX(59)-JVS(837)*XX(131)
  XX(56) = XX(56)-JVS(273)*XX(58)-JVS(836)*XX(131)
  XX(55) = XX(55)-JVS(264)*XX(56)-JVS(268)*XX(57)-JVS(835)*XX(131)
  XX(54) = XX(54)-JVS(263)*XX(56)-JVS(834)*XX(131)
  XX(53) = XX(53)-JVS(254)*XX(54)-JVS(258)*XX(55)-JVS(833)*XX(131)
  XX(52) = XX(52)-JVS(253)*XX(54)-JVS(832)*XX(131)
  XX(51) = XX(51)-JVS(244)*XX(52)-JVS(248)*XX(53)-JVS(831)*XX(131)
  XX(50) = XX(50)-JVS(243)*XX(52)-JVS(830)*XX(131)
  XX(49) = XX(49)-JVS(234)*XX(50)-JVS(238)*XX(51)-JVS(829)*XX(131)
  XX(48) = XX(48)-JVS(645)*XX(123)-JVS(828)*XX(131)-JVS(911)*XX(132)
  XX(47) = XX(47)-JVS(233)*XX(50)-JVS(827)*XX(131)
  XX(46) = XX(46)-JVS(223)*XX(47)-JVS(228)*XX(49)-JVS(826)*XX(131)
  XX(45) = XX(45)-JVS(232)*XX(50)-JVS(237)*XX(51)-JVS(825)*XX(131)
  XX(44) = XX(44)-JVS(242)*XX(52)-JVS(247)*XX(53)-JVS(824)*XX(131)
  XX(43) = XX(43)-JVS(252)*XX(54)-JVS(257)*XX(55)-JVS(823)*XX(131)
  XX(42) = XX(42)-JVS(262)*XX(56)-JVS(267)*XX(57)-JVS(822)*XX(131)
  XX(41) = XX(41)-JVS(272)*XX(58)-JVS(277)*XX(59)-JVS(821)*XX(131)
  XX(40) = XX(40)-JVS(281)*XX(60)-JVS(292)*XX(62)-JVS(820)*XX(131)
  XX(39) = XX(39)-JVS(819)*XX(131)
  XX(38) = XX(38)-JVS(204)*XX(39)-JVS(286)*XX(61)-JVS(818)*XX(131)
  XX(37) = XX(37)-JVS(817)*XX(131)-JVS(910)*XX(132)
  XX(36) = XX(36)-JVS(475)*XX(106)-JVS(816)*XX(131)
  XX(35) = XX(35)-JVS(317)*XX(72)-JVS(815)*XX(131)
  XX(34) = XX(34)-JVS(324)*XX(74)-JVS(814)*XX(131)
  XX(33) = XX(33)-JVS(334)*XX(76)-JVS(813)*XX(131)
  XX(32) = XX(32)-JVS(344)*XX(78)-JVS(812)*XX(131)
  XX(31) = XX(31)-JVS(354)*XX(80)-JVS(811)*XX(131)
  XX(30) = XX(30)-JVS(364)*XX(82)-JVS(810)*XX(131)
  XX(29) = XX(29)-JVS(384)*XX(86)-JVS(809)*XX(131)
  XX(28) = XX(28)-JVS(378)*XX(85)-JVS(808)*XX(131)
  XX(27) = XX(27)-JVS(222)*XX(47)-JVS(807)*XX(131)
  XX(26) = XX(26)-JVS(231)*XX(50)-JVS(806)*XX(131)
  XX(25) = XX(25)-JVS(241)*XX(52)-JVS(805)*XX(131)
  XX(24) = XX(24)-JVS(251)*XX(54)-JVS(804)*XX(131)
  XX(23) = XX(23)-JVS(261)*XX(56)-JVS(803)*XX(131)
  XX(22) = XX(22)-JVS(271)*XX(58)-JVS(802)*XX(131)
  XX(21) = XX(21)-JVS(291)*XX(62)-JVS(801)*XX(131)
  XX(20) = XX(20)-JVS(285)*XX(61)-JVS(800)*XX(131)
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
  XX(7) = XX(7)-JVS(97)*XX(19)-JVS(799)*XX(131)
  XX(6) = XX(6)-JVS(96)*XX(19)-JVS(798)*XX(131)
  XX(5) = XX(5)-JVS(95)*XX(19)-JVS(797)*XX(131)
  XX(4) = XX(4)-JVS(94)*XX(19)-JVS(796)*XX(131)
  XX(3) = XX(3)
  XX(2) = XX(2)
  XX(1) = XX(1)
      
END SUBROUTINE KppSolveTR

! End of KppSolveTR function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbmz_mosaic_4bin_vbs9_LinearAlgebra

