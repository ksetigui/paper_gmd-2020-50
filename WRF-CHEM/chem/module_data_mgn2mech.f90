MODULE module_data_mgn2mech
  USE module_state_description
  USE module_data_megan2
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: non_react = 9999
  INTEGER, PARAMETER :: n_megan2cbmz = 173
  INTEGER, DIMENSION (n_megan2cbmz) :: p_of_megan2cbmz, p_of_cbmz
  REAL, DIMENSION (n_megan2cbmz) :: cbmz_per_megan
  DATA p_of_cbmz / n_megan2cbmz*non_react /
  INTEGER, PARAMETER :: n_megan2radm2 = 165
  INTEGER, DIMENSION (n_megan2radm2) :: p_of_megan2radm2, p_of_radm2
  REAL, DIMENSION (n_megan2radm2) :: radm2_per_megan
  DATA p_of_radm2 / n_megan2radm2*non_react /
  INTEGER, PARAMETER :: n_megan2racm = 147
  INTEGER, DIMENSION (n_megan2racm) :: p_of_megan2racm, p_of_racm
  REAL, DIMENSION (n_megan2racm) :: racm_per_megan
  DATA p_of_racm / n_megan2racm*non_react /
  INTEGER, PARAMETER :: n_megan2racmSOA = 178
  INTEGER, DIMENSION (n_megan2racmSOA) :: p_of_megan2racmSOA, p_of_racmSOA
  REAL, DIMENSION (n_megan2racmSOA) :: racmSOA_per_megan
  DATA p_of_racmSOA / n_megan2racmSOA*non_react /
  INTEGER, PARAMETER :: n_megan2mozcart = 142
  INTEGER, DIMENSION (n_megan2mozcart) :: p_of_megan2mozcart, p_of_mozcart
  REAL, DIMENSION (n_megan2mozcart) :: mozcart_per_megan
  DATA p_of_mozcart / n_megan2mozcart*non_react /
  INTEGER, PARAMETER :: n_megan2saprcnov = 138
  INTEGER, DIMENSION (n_megan2saprcnov) :: p_of_megan2saprcnov, p_of_saprcnov
  REAL, DIMENSION (n_megan2saprcnov) :: saprcnov_per_megan
  DATA p_of_saprcnov / n_megan2saprcnov*non_react /
  INTEGER, PARAMETER :: n_megan2crimech = 188
  INTEGER, DIMENSION (n_megan2crimech) :: p_of_megan2crimech, p_of_crimech
  REAL, DIMENSION (n_megan2crimech) :: crimech_per_megan
  DATA p_of_crimech / n_megan2crimech*non_react /
CONTAINS
  SUBROUTINE get_megan2mozcart_table
    p_of_megan2mozcart( 1) = is_isoprene ; p_of_mozcart( 1) = p_isopr ; mozcart_per_megan( 1) = 1.
    p_of_megan2mozcart( 2) = is_myrcene ; p_of_mozcart( 2) = p_c10h16 ; mozcart_per_megan( 2) = 1.
    p_of_megan2mozcart( 3) = is_sabinene ; p_of_mozcart( 3) = p_c10h16 ; mozcart_per_megan( 3) = 1.
    p_of_megan2mozcart( 4) = is_limonene ; p_of_mozcart( 4) = p_c10h16 ; mozcart_per_megan( 4) = 1.
    p_of_megan2mozcart( 5) = is_carene_3 ; p_of_mozcart( 5) = p_c10h16 ; mozcart_per_megan( 5) = 1.
    p_of_megan2mozcart( 6) = is_ocimene_t_b ; p_of_mozcart( 6) = p_c10h16 ; mozcart_per_megan( 6) = 1.
    p_of_megan2mozcart( 7) = is_pinene_b ; p_of_mozcart( 7) = p_c10h16 ; mozcart_per_megan( 7) = 1.
    p_of_megan2mozcart( 8) = is_pinene_a ; p_of_mozcart( 8) = p_c10h16 ; mozcart_per_megan( 8) = 1.
    p_of_megan2mozcart( 9) = is_2met_styrene ; p_of_mozcart( 9) = p_c10h16 ; mozcart_per_megan( 9) = 1.
    p_of_megan2mozcart( 10) = is_cymene_p ; p_of_mozcart( 10) = p_tol ; mozcart_per_megan( 10) = 1.5
    p_of_megan2mozcart( 11) = is_cymene_o ; p_of_mozcart( 11) = p_tol ; mozcart_per_megan( 11) = 1.5
    p_of_megan2mozcart( 12) = is_phellandrene_a ; p_of_mozcart( 12) = p_c10h16 ; mozcart_per_megan( 12) = 1.
    p_of_megan2mozcart( 13) = is_thujene_a ; p_of_mozcart( 13) = p_c10h16 ; mozcart_per_megan( 13) = 1.
    p_of_megan2mozcart( 14) = is_terpinene_a ; p_of_mozcart( 14) = p_c10h16 ; mozcart_per_megan( 14) = 1.
    p_of_megan2mozcart( 15) = is_terpinene_g ; p_of_mozcart( 15) = p_c10h16 ; mozcart_per_megan( 15) = 1.
    p_of_megan2mozcart( 16) = is_terpinolene ; p_of_mozcart( 16) = p_c10h16 ; mozcart_per_megan( 16) = 1.
    p_of_megan2mozcart( 17) = is_phellandrene_b ; p_of_mozcart( 17) = p_c10h16 ; mozcart_per_megan( 17) = 1.
    p_of_megan2mozcart( 18) = is_camphene ; p_of_mozcart( 18) = p_c10h16 ; mozcart_per_megan( 18) = 1.
    p_of_megan2mozcart( 19) = is_bornene ; p_of_mozcart( 19) = p_c10h16 ; mozcart_per_megan( 19) = 1.
    p_of_megan2mozcart( 20) = is_fenchene_a ; p_of_mozcart( 20) = p_c10h16 ; mozcart_per_megan( 20) = 1.
    p_of_megan2mozcart( 21) = is_ocimene_al ; p_of_mozcart( 21) = p_c10h16 ; mozcart_per_megan( 21) = 1.
    p_of_megan2mozcart( 22) = is_ocimene_c_b ; p_of_mozcart( 22) = p_c10h16 ; mozcart_per_megan( 22) = 1.
    p_of_megan2mozcart( 23) = is_tricyclene ; p_of_mozcart( 23) = p_c10h16 ; mozcart_per_megan( 23) = 1.
    p_of_megan2mozcart( 24) = is_estragole ; p_of_mozcart( 24) = p_c10h16 ; mozcart_per_megan( 24) = 1.
    p_of_megan2mozcart( 25) = is_camphor ; p_of_mozcart( 25) = p_bigalk ; mozcart_per_megan( 25) = 2.
    p_of_megan2mozcart( 26) = is_fenchone ; p_of_mozcart( 26) = p_bigalk ; mozcart_per_megan( 26) = 2.
    p_of_megan2mozcart( 27) = is_piperitone ; p_of_mozcart( 27) = p_c10h16 ; mozcart_per_megan( 27) = 1.
    p_of_megan2mozcart( 28) = is_thujone_a ; p_of_mozcart( 28) = p_bigalk ; mozcart_per_megan( 28) = 2.
    p_of_megan2mozcart( 29) = is_thujone_b ; p_of_mozcart( 29) = p_bigalk ; mozcart_per_megan( 29) = 2.
    p_of_megan2mozcart( 30) = is_cineole_1_8 ; p_of_mozcart( 30) = p_bigalk ; mozcart_per_megan( 30) = 2.
    p_of_megan2mozcart( 31) = is_borneol ; p_of_mozcart( 31) = p_bigalk ; mozcart_per_megan( 31) = 2.
    p_of_megan2mozcart( 32) = is_linalool ; p_of_mozcart( 32) = p_c10h16 ; mozcart_per_megan( 32) = 1.
    p_of_megan2mozcart( 33) = is_terpineol_4 ; p_of_mozcart( 33) = p_c10h16 ; mozcart_per_megan( 33) = 1.
    p_of_megan2mozcart( 34) = is_terpineol_a ; p_of_mozcart( 34) = p_c10h16 ; mozcart_per_megan( 34) = 1.
    p_of_megan2mozcart( 35) = is_linalool_oxd_c ; p_of_mozcart( 35) = p_c10h16 ; mozcart_per_megan( 35) = 1.25
    p_of_megan2mozcart( 36) = is_linalool_oxd_t ; p_of_mozcart( 36) = p_c10h16 ; mozcart_per_megan( 36) = 1.25
    p_of_megan2mozcart( 37) = is_ionone_b ; p_of_mozcart( 37) = p_c10h16 ; mozcart_per_megan( 37) = 1.4
    p_of_megan2mozcart( 38) = is_bornyl_act ; p_of_mozcart( 38) = p_bigalk ; mozcart_per_megan( 38) = 2.7
    p_of_megan2mozcart( 39) = is_farnescene_a ; p_of_mozcart( 39) = p_c10h16 ; mozcart_per_megan( 39) = 1.5
    p_of_megan2mozcart( 40) = is_caryophyllene_b ; p_of_mozcart( 40) = p_c10h16 ; mozcart_per_megan( 40) = 1.5
    p_of_megan2mozcart( 41) = is_acoradiene ; p_of_mozcart( 41) = p_c10h16 ; mozcart_per_megan( 41) = 1.5
    p_of_megan2mozcart( 42) = is_aromadendrene ; p_of_mozcart( 42) = p_c10h16 ; mozcart_per_megan( 42) = 1.5
    p_of_megan2mozcart( 43) = is_bergamotene_a ; p_of_mozcart( 43) = p_c10h16 ; mozcart_per_megan( 43) = 1.5
    p_of_megan2mozcart( 44) = is_bergamotene_b ; p_of_mozcart( 44) = p_c10h16 ; mozcart_per_megan( 44) = 1.5
    p_of_megan2mozcart( 45) = is_bisabolene_a ; p_of_mozcart( 45) = p_c10h16 ; mozcart_per_megan( 45) = 1.5
    p_of_megan2mozcart( 46) = is_bisabolene_b ; p_of_mozcart( 46) = p_c10h16 ; mozcart_per_megan( 46) = 1.5
    p_of_megan2mozcart( 47) = is_bourbonene_b ; p_of_mozcart( 47) = p_c10h16 ; mozcart_per_megan( 47) = 1.5
    p_of_megan2mozcart( 48) = is_cadinene_d ; p_of_mozcart( 48) = p_c10h16 ; mozcart_per_megan( 48) = 1.5
    p_of_megan2mozcart( 49) = is_cadinene_g ; p_of_mozcart( 49) = p_c10h16 ; mozcart_per_megan( 49) = 1.5
    p_of_megan2mozcart( 50) = is_cedrene_a ; p_of_mozcart( 50) = p_c10h16 ; mozcart_per_megan( 50) = 1.5
    p_of_megan2mozcart( 51) = is_copaene_a ; p_of_mozcart( 51) = p_c10h16 ; mozcart_per_megan( 51) = 1.5
    p_of_megan2mozcart( 52) = is_cubebene_a ; p_of_mozcart( 52) = p_c10h16 ; mozcart_per_megan( 52) = 1.5
    p_of_megan2mozcart( 53) = is_cubebene_b ; p_of_mozcart( 53) = p_c10h16 ; mozcart_per_megan( 53) = 1.5
    p_of_megan2mozcart( 54) = is_elemene_b ; p_of_mozcart( 54) = p_c10h16 ; mozcart_per_megan( 54) = 1.5
    p_of_megan2mozcart( 55) = is_farnescene_b ; p_of_mozcart( 55) = p_c10h16 ; mozcart_per_megan( 55) = 1.5
    p_of_megan2mozcart( 56) = is_germacrene_B ; p_of_mozcart( 56) = p_c10h16 ; mozcart_per_megan( 56) = 1.5
    p_of_megan2mozcart( 57) = is_germacrene_D ; p_of_mozcart( 57) = p_c10h16 ; mozcart_per_megan( 57) = 1.5
    p_of_megan2mozcart( 58) = is_gurjunene_b ; p_of_mozcart( 58) = p_c10h16 ; mozcart_per_megan( 58) = 1.5
    p_of_megan2mozcart( 59) = is_humulene_a ; p_of_mozcart( 59) = p_c10h16 ; mozcart_per_megan( 59) = 1.5
    p_of_megan2mozcart( 60) = is_humulene_g ; p_of_mozcart( 60) = p_c10h16 ; mozcart_per_megan( 60) = 1.5
    p_of_megan2mozcart( 61) = is_isolongifolene ; p_of_mozcart( 61) = p_c10h16 ; mozcart_per_megan( 61) = 1.5
    p_of_megan2mozcart( 62) = is_longifolene ; p_of_mozcart( 62) = p_c10h16 ; mozcart_per_megan( 62) = 1.5
    p_of_megan2mozcart( 63) = is_longipinene ; p_of_mozcart( 63) = p_c10h16 ; mozcart_per_megan( 63) = 1.5
    p_of_megan2mozcart( 64) = is_muurolene_a ; p_of_mozcart( 64) = p_c10h16 ; mozcart_per_megan( 64) = 1.5
    p_of_megan2mozcart( 65) = is_muurolene_g ; p_of_mozcart( 65) = p_c10h16 ; mozcart_per_megan( 65) = 1.5
    p_of_megan2mozcart( 66) = is_selinene_b ; p_of_mozcart( 66) = p_c10h16 ; mozcart_per_megan( 66) = 1.5
    p_of_megan2mozcart( 67) = is_selinene_d ; p_of_mozcart( 67) = p_c10h16 ; mozcart_per_megan( 67) = 1.5
    p_of_megan2mozcart( 68) = is_nerolidol_c ; p_of_mozcart( 68) = p_c10h16 ; mozcart_per_megan( 68) = 1.5
    p_of_megan2mozcart( 69) = is_nerolidol_t ; p_of_mozcart( 69) = p_c10h16 ; mozcart_per_megan( 69) = 1.5
    p_of_megan2mozcart( 70) = is_cedrol ; p_of_mozcart( 70) = p_bigalk ; mozcart_per_megan( 70) = 3.
    p_of_megan2mozcart( 71) = is_mbo_2m3e2ol ; p_of_mozcart( 71) = p_isopr ; mozcart_per_megan( 71) = 2.4
    p_of_megan2mozcart( 72) = is_methanol ; p_of_mozcart( 72) = p_ch3oh ; mozcart_per_megan( 72) = 1.
    p_of_megan2mozcart( 73) = is_acetone ; p_of_mozcart( 73) = p_acet ; mozcart_per_megan( 73) = 1.
    p_of_megan2mozcart( 74) = is_methane ; p_of_mozcart( 74) = non_react ; mozcart_per_megan( 74) = 1.
    p_of_megan2mozcart( 75) = is_ammonia ; p_of_mozcart( 75) = p_nh3 ; mozcart_per_megan( 75) = 1.
    p_of_megan2mozcart( 76) = is_nitrous_oxd ; p_of_mozcart( 76) = non_react ; mozcart_per_megan( 76) = 1.
    p_of_megan2mozcart( 77) = is_nitric_oxd ; p_of_mozcart( 77) = p_no ; mozcart_per_megan( 77) = 1.
    p_of_megan2mozcart( 78) = is_acetaldehyde ; p_of_mozcart( 78) = p_ald ; mozcart_per_megan( 78) = 1.
    p_of_megan2mozcart( 79) = is_ethanol ; p_of_mozcart( 79) = p_c2h5oh ; mozcart_per_megan( 79) = 1.
    p_of_megan2mozcart( 80) = is_formic_acid ; p_of_mozcart( 80) = non_react ; mozcart_per_megan( 80) = 1.
    p_of_megan2mozcart( 81) = is_formaldehyde ; p_of_mozcart( 81) = p_hcho ; mozcart_per_megan( 81) = 1.
    p_of_megan2mozcart( 82) = is_acetic_acid ; p_of_mozcart( 82) = p_ch3cooh ; mozcart_per_megan( 82) = 1.
    p_of_megan2mozcart( 83) = is_mbo_3m2e1ol ; p_of_mozcart( 83) = p_isopr ; mozcart_per_megan( 83) = 1.25
    p_of_megan2mozcart( 84) = is_benzaldehyde ; p_of_mozcart( 84) = p_tol ; mozcart_per_megan( 84) = 1.1
    p_of_megan2mozcart( 85) = is_butanone_2 ; p_of_mozcart( 85) = p_mek ; mozcart_per_megan( 85) = 1.
    p_of_megan2mozcart( 86) = is_decanal ; p_of_mozcart( 86) = p_bigalk ; mozcart_per_megan( 86) = 3.
    p_of_megan2mozcart( 87) = is_dodecene_1 ; p_of_mozcart( 87) = p_bigene ; mozcart_per_megan( 87) = 2.25
    p_of_megan2mozcart( 88) = is_geranyl_acetone ; p_of_mozcart( 88) = p_c10h16 ; mozcart_per_megan( 88) = 1.4
    p_of_megan2mozcart( 89) = is_heptanal ; p_of_mozcart( 89) = p_bigalk ; mozcart_per_megan( 89) = 2.
    p_of_megan2mozcart( 90) = is_heptane ; p_of_mozcart( 90) = p_bigalk ; mozcart_per_megan( 90) = 2.
    p_of_megan2mozcart( 91) = is_hexane ; p_of_mozcart( 91) = p_bigalk ; mozcart_per_megan( 91) = 1.5
    p_of_megan2mozcart( 92) = is_met_benzoate ; p_of_mozcart( 92) = p_tol ; mozcart_per_megan( 92) = 1.5
    p_of_megan2mozcart( 93) = is_met_heptenone ; p_of_mozcart( 93) = p_bigene ; mozcart_per_megan( 93) = 1.75
    p_of_megan2mozcart( 94) = is_neryl_acetone ; p_of_mozcart( 94) = p_bigene ; mozcart_per_megan( 94) = 2.7
    p_of_megan2mozcart( 95) = is_nonanal ; p_of_mozcart( 95) = p_bigalk ; mozcart_per_megan( 95) = 2.5
    p_of_megan2mozcart( 96) = is_nonenal ; p_of_mozcart( 96) = p_bigene ; mozcart_per_megan( 96) = 2.
    p_of_megan2mozcart( 97) = is_octanal ; p_of_mozcart( 97) = p_bigalk ; mozcart_per_megan( 97) = 2.3
    p_of_megan2mozcart( 98) = is_octanol ; p_of_mozcart( 98) = p_bigalk ; mozcart_per_megan( 98) = 2.3
    p_of_megan2mozcart( 99) = is_octenol_1e3ol ; p_of_mozcart( 99) = p_bigene ; mozcart_per_megan( 99) = 1.7
    p_of_megan2mozcart(100) = is_oxopentanal ; p_of_mozcart(100) = p_mek ; mozcart_per_megan(100) = 1.4
    p_of_megan2mozcart(101) = is_pentane ; p_of_mozcart(101) = p_bigalk ; mozcart_per_megan(101) = 1.25
    p_of_megan2mozcart(102) = is_phenyl_cco ; p_of_mozcart(102) = p_tol ; mozcart_per_megan(102) = 1.3
    p_of_megan2mozcart(103) = is_pyruvic_acid ; p_of_mozcart(103) = non_react ; mozcart_per_megan(103) = 1.
    p_of_megan2mozcart(104) = is_terpinyl_act_a ; p_of_mozcart(104) = p_c10h16 ; mozcart_per_megan(104) = 1.4
    p_of_megan2mozcart(105) = is_tetradecene_1 ; p_of_mozcart(105) = p_bigene ; mozcart_per_megan(105) = 3.5
    p_of_megan2mozcart(106) = is_toluene ; p_of_mozcart(106) = p_tol ; mozcart_per_megan(106) = 1.
    p_of_megan2mozcart(107) = is_carbon_monoxide ; p_of_mozcart(107) = p_co ; mozcart_per_megan(107) = 1.
    p_of_megan2mozcart(108) = is_butene ; p_of_mozcart(108) = p_bigene ; mozcart_per_megan(108) = .8
    p_of_megan2mozcart(109) = is_ethane ; p_of_mozcart(109) = p_c2h6 ; mozcart_per_megan(109) = 1.
    p_of_megan2mozcart(110) = is_ethene ; p_of_mozcart(110) = p_c2h4 ; mozcart_per_megan(110) = 1.
    p_of_megan2mozcart(111) = is_hydrogen_cyanide ; p_of_mozcart(111) = non_react ; mozcart_per_megan(111) = 1.
    p_of_megan2mozcart(112) = is_propane ; p_of_mozcart(112) = p_c3h8 ; mozcart_per_megan(112) = 1.
    p_of_megan2mozcart(113) = is_propene ; p_of_mozcart(113) = p_c3h6 ; mozcart_per_megan(113) = 1.
    p_of_megan2mozcart(114) = is_carbon_2s ; p_of_mozcart(114) = non_react ; mozcart_per_megan(114) = 1.
    p_of_megan2mozcart(115) = is_carbonyl_s ; p_of_mozcart(115) = non_react ; mozcart_per_megan(115) = 1.
    p_of_megan2mozcart(116) = is_diallyl_2s ; p_of_mozcart(116) = p_bigene ; mozcart_per_megan(116) = .66
    p_of_megan2mozcart(117) = is_diallyl_2s ; p_of_mozcart(117) = p_so2 ; mozcart_per_megan(117) = 1.53
    p_of_megan2mozcart(118) = is_2met_2s ; p_of_mozcart(118) = p_c2h6 ; mozcart_per_megan(118) = 1.
    p_of_megan2mozcart(119) = is_2met_2s ; p_of_mozcart(119) = p_so2 ; mozcart_per_megan(119) = 1.
    p_of_megan2mozcart(120) = is_met_chloride ; p_of_mozcart(120) = non_react ; mozcart_per_megan(120) = 1.
    p_of_megan2mozcart(121) = is_met_bromide ; p_of_mozcart(121) = non_react ; mozcart_per_megan(121) = 1.
    p_of_megan2mozcart(122) = is_met_iodide ; p_of_mozcart(122) = non_react ; mozcart_per_megan(122) = 1.
    p_of_megan2mozcart(123) = is_hydrogen_s ; p_of_mozcart(123) = p_so2 ; mozcart_per_megan(123) = .5
    p_of_megan2mozcart(124) = is_met_mercaptan ; p_of_mozcart(124) = p_so2 ; mozcart_per_megan(124) = .75
    p_of_megan2mozcart(125) = is_met_propenyl_2s ; p_of_mozcart(125) = p_c3h6 ; mozcart_per_megan(125) = 2.8
    p_of_megan2mozcart(126) = is_met_propenyl_2s ; p_of_mozcart(126) = p_so2 ; mozcart_per_megan(126) = 1.8
    p_of_megan2mozcart(127) = is_pppp_2s ; p_of_mozcart(127) = p_c3h6 ; mozcart_per_megan(127) = 3.5
    p_of_megan2mozcart(128) = is_pppp_2s ; p_of_mozcart(128) = p_so2 ; mozcart_per_megan(128) = 2.3
    p_of_megan2mozcart(129) = is_2met_nonatriene ; p_of_mozcart(129) = p_c10h16 ; mozcart_per_megan(129) = 1.1
    p_of_megan2mozcart(130) = is_met_salicylate ; p_of_mozcart(130) = p_tol ; mozcart_per_megan(130) = 1.6
    p_of_megan2mozcart(131) = is_indole ; p_of_mozcart(131) = non_react ; mozcart_per_megan(131) = 1.
    p_of_megan2mozcart(132) = is_jasmone ; p_of_mozcart(132) = p_c10h16 ; mozcart_per_megan(132) = 1.2
    p_of_megan2mozcart(133) = is_met_jasmonate ; p_of_mozcart(133) = p_c10h16 ; mozcart_per_megan(133) = 1.6
    p_of_megan2mozcart(134) = is_3met_3dctt ; p_of_mozcart(134) = p_bigene ; mozcart_per_megan(134) = 3.
    p_of_megan2mozcart(135) = is_hexanal ; p_of_mozcart(135) = p_bigalk ; mozcart_per_megan(135) = 1.8
    p_of_megan2mozcart(136) = is_hexanol_1 ; p_of_mozcart(136) = p_bigalk ; mozcart_per_megan(136) = 1.8
    p_of_megan2mozcart(137) = is_hexenal_c3 ; p_of_mozcart(137) = p_bigene ; mozcart_per_megan(137) = 1.4
    p_of_megan2mozcart(138) = is_hexenal_t2 ; p_of_mozcart(138) = p_bigene ; mozcart_per_megan(138) = 1.4
    p_of_megan2mozcart(139) = is_hexenol_c3 ; p_of_mozcart(139) = p_bigene ; mozcart_per_megan(139) = 1.4
    p_of_megan2mozcart(140) = is_hexenyl_act_c3 ; p_of_mozcart(140) = p_bigene ; mozcart_per_megan(140) = 2.
    p_of_megan2mozcart(141) = is_mbo_3m3e1ol ; p_of_mozcart(141) = p_isopr ; mozcart_per_megan(141) = 1.25
    p_of_megan2mozcart(142) = is_2met_s ; p_of_mozcart(142) = p_dms ; mozcart_per_megan(142) = 1.
  END SUBROUTINE get_megan2mozcart_table
  SUBROUTINE get_megan2cbmz_table
    p_of_megan2cbmz( 1) = is_isoprene ; p_of_cbmz( 1) = p_iso ; cbmz_per_megan( 1) = 1.
    p_of_megan2cbmz( 2) = is_myrcene ; p_of_cbmz( 2) = p_iso ; cbmz_per_megan( 2) = 2.
    p_of_megan2cbmz( 3) = is_sabinene ; p_of_cbmz( 3) = p_iso ; cbmz_per_megan( 3) = 2.
    p_of_megan2cbmz( 4) = is_limonene ; p_of_cbmz( 4) = p_iso ; cbmz_per_megan( 4) = 2.
    p_of_megan2cbmz( 5) = is_carene_3 ; p_of_cbmz( 5) = p_iso ; cbmz_per_megan( 5) = 2.
    p_of_megan2cbmz( 6) = is_ocimene_t_b ; p_of_cbmz( 6) = p_iso ; cbmz_per_megan( 6) = 2.
    p_of_megan2cbmz( 7) = is_pinene_b ; p_of_cbmz( 7) = p_iso ; cbmz_per_megan( 7) = 2.
    p_of_megan2cbmz( 8) = is_pinene_a ; p_of_cbmz( 8) = p_iso ; cbmz_per_megan( 8) = 2.
    p_of_megan2cbmz( 9) = is_2met_styrene ; p_of_cbmz( 9) = p_iso ; cbmz_per_megan( 9) = 2.
    p_of_megan2cbmz( 10) = is_cymene_p ; p_of_cbmz( 10) = p_iso ; cbmz_per_megan( 10) = 2.
    p_of_megan2cbmz( 11) = is_cymene_o ; p_of_cbmz( 11) = p_iso ; cbmz_per_megan( 11) = 2.
    p_of_megan2cbmz( 12) = is_phellandrene_a ; p_of_cbmz( 12) = p_iso ; cbmz_per_megan( 12) = 2.
    p_of_megan2cbmz( 13) = is_thujene_a ; p_of_cbmz( 13) = p_iso ; cbmz_per_megan( 13) = 2.
    p_of_megan2cbmz( 14) = is_terpinene_a ; p_of_cbmz( 14) = p_iso ; cbmz_per_megan( 14) = 2.
    p_of_megan2cbmz( 15) = is_terpinene_g ; p_of_cbmz( 15) = p_iso ; cbmz_per_megan( 15) = 2.
    p_of_megan2cbmz( 16) = is_terpinolene ; p_of_cbmz( 16) = p_iso ; cbmz_per_megan( 16) = 2.
    p_of_megan2cbmz( 17) = is_phellandrene_b ; p_of_cbmz( 17) = p_iso ; cbmz_per_megan( 17) = 2.
    p_of_megan2cbmz( 18) = is_camphene ; p_of_cbmz( 18) = p_iso ; cbmz_per_megan( 18) = 2.
    p_of_megan2cbmz( 19) = is_bornene ; p_of_cbmz( 19) = p_iso ; cbmz_per_megan( 19) = 2.
    p_of_megan2cbmz( 20) = is_fenchene_a ; p_of_cbmz( 20) = p_iso ; cbmz_per_megan( 20) = 2.
    p_of_megan2cbmz( 21) = is_ocimene_al ; p_of_cbmz( 21) = p_iso ; cbmz_per_megan( 21) = 2.
    p_of_megan2cbmz( 22) = is_ocimene_c_b ; p_of_cbmz( 22) = p_iso ; cbmz_per_megan( 22) = 2.
    p_of_megan2cbmz( 23) = is_tricyclene ; p_of_cbmz( 23) = p_iso ; cbmz_per_megan( 23) = 2.
    p_of_megan2cbmz( 24) = is_estragole ; p_of_cbmz( 24) = p_iso ; cbmz_per_megan( 24) = 2.
    p_of_megan2cbmz( 25) = is_camphor ; p_of_cbmz( 25) = p_iso ; cbmz_per_megan( 25) = 2.
    p_of_megan2cbmz( 26) = is_fenchone ; p_of_cbmz( 26) = p_iso ; cbmz_per_megan( 26) = 2.
    p_of_megan2cbmz( 27) = is_piperitone ; p_of_cbmz( 27) = p_iso ; cbmz_per_megan( 27) = 2.
    p_of_megan2cbmz( 28) = is_thujone_a ; p_of_cbmz( 28) = p_iso ; cbmz_per_megan( 28) = 2.
    p_of_megan2cbmz( 29) = is_thujone_b ; p_of_cbmz( 29) = p_iso ; cbmz_per_megan( 29) = 2.
    p_of_megan2cbmz( 30) = is_cineole_1_8 ; p_of_cbmz( 30) = p_iso ; cbmz_per_megan( 30) = 2.
    p_of_megan2cbmz( 31) = is_borneol ; p_of_cbmz( 31) = p_iso ; cbmz_per_megan( 31) = 2.
    p_of_megan2cbmz( 32) = is_linalool ; p_of_cbmz( 32) = p_iso ; cbmz_per_megan( 32) = 2.
    p_of_megan2cbmz( 33) = is_terpineol_4 ; p_of_cbmz( 33) = p_iso ; cbmz_per_megan( 33) = 2.
    p_of_megan2cbmz( 34) = is_terpineol_a ; p_of_cbmz( 34) = p_iso ; cbmz_per_megan( 34) = 2.
    p_of_megan2cbmz( 35) = is_linalool_oxd_c ; p_of_cbmz( 35) = p_iso ; cbmz_per_megan( 35) = 2.
    p_of_megan2cbmz( 36) = is_linalool_oxd_t ; p_of_cbmz( 36) = p_iso ; cbmz_per_megan( 36) = 2.
    p_of_megan2cbmz( 37) = is_ionone_b ; p_of_cbmz( 37) = p_iso ; cbmz_per_megan( 37) = 3.
    p_of_megan2cbmz( 38) = is_bornyl_act ; p_of_cbmz( 38) = p_iso ; cbmz_per_megan( 38) = 2.
    p_of_megan2cbmz( 39) = is_farnescene_a ; p_of_cbmz( 39) = p_iso ; cbmz_per_megan( 39) = 3.
    p_of_megan2cbmz( 40) = is_caryophyllene_b ; p_of_cbmz( 40) = p_iso ; cbmz_per_megan( 40) = 3.
    p_of_megan2cbmz( 41) = is_acoradiene ; p_of_cbmz( 41) = p_iso ; cbmz_per_megan( 41) = 3.
    p_of_megan2cbmz( 42) = is_aromadendrene ; p_of_cbmz( 42) = p_iso ; cbmz_per_megan( 42) = 3.
    p_of_megan2cbmz( 43) = is_bergamotene_a ; p_of_cbmz( 43) = p_iso ; cbmz_per_megan( 43) = 3.
    p_of_megan2cbmz( 44) = is_bergamotene_b ; p_of_cbmz( 44) = p_iso ; cbmz_per_megan( 44) = 3.
    p_of_megan2cbmz( 45) = is_bisabolene_a ; p_of_cbmz( 45) = p_iso ; cbmz_per_megan( 45) = 3.
    p_of_megan2cbmz( 46) = is_bisabolene_b ; p_of_cbmz( 46) = p_iso ; cbmz_per_megan( 46) = 3.
    p_of_megan2cbmz( 47) = is_bourbonene_b ; p_of_cbmz( 47) = p_iso ; cbmz_per_megan( 47) = 3.
    p_of_megan2cbmz( 48) = is_cadinene_d ; p_of_cbmz( 48) = p_iso ; cbmz_per_megan( 48) = 3.
    p_of_megan2cbmz( 49) = is_cadinene_g ; p_of_cbmz( 49) = p_iso ; cbmz_per_megan( 49) = 3.
    p_of_megan2cbmz( 50) = is_cedrene_a ; p_of_cbmz( 50) = p_iso ; cbmz_per_megan( 50) = 3.
    p_of_megan2cbmz( 51) = is_copaene_a ; p_of_cbmz( 51) = p_iso ; cbmz_per_megan( 51) = 3.
    p_of_megan2cbmz( 52) = is_cubebene_a ; p_of_cbmz( 52) = p_iso ; cbmz_per_megan( 52) = 3.
    p_of_megan2cbmz( 53) = is_cubebene_b ; p_of_cbmz( 53) = p_iso ; cbmz_per_megan( 53) = 3.
    p_of_megan2cbmz( 54) = is_elemene_b ; p_of_cbmz( 54) = p_iso ; cbmz_per_megan( 54) = 3.
    p_of_megan2cbmz( 55) = is_farnescene_b ; p_of_cbmz( 55) = p_iso ; cbmz_per_megan( 55) = 3.
    p_of_megan2cbmz( 56) = is_germacrene_B ; p_of_cbmz( 56) = p_iso ; cbmz_per_megan( 56) = 3.
    p_of_megan2cbmz( 57) = is_germacrene_D ; p_of_cbmz( 57) = p_iso ; cbmz_per_megan( 57) = 3.
    p_of_megan2cbmz( 58) = is_gurjunene_b ; p_of_cbmz( 58) = p_iso ; cbmz_per_megan( 58) = 3.
    p_of_megan2cbmz( 59) = is_humulene_a ; p_of_cbmz( 59) = p_iso ; cbmz_per_megan( 59) = 3.
    p_of_megan2cbmz( 60) = is_humulene_g ; p_of_cbmz( 60) = p_iso ; cbmz_per_megan( 60) = 3.
    p_of_megan2cbmz( 61) = is_isolongifolene ; p_of_cbmz( 61) = p_iso ; cbmz_per_megan( 61) = 3.
    p_of_megan2cbmz( 62) = is_longifolene ; p_of_cbmz( 62) = p_iso ; cbmz_per_megan( 62) = 3.
    p_of_megan2cbmz( 63) = is_longipinene ; p_of_cbmz( 63) = p_iso ; cbmz_per_megan( 63) = 3.
    p_of_megan2cbmz( 64) = is_muurolene_a ; p_of_cbmz( 64) = p_iso ; cbmz_per_megan( 64) = 3.
    p_of_megan2cbmz( 65) = is_muurolene_g ; p_of_cbmz( 65) = p_iso ; cbmz_per_megan( 65) = 3.
    p_of_megan2cbmz( 66) = is_selinene_b ; p_of_cbmz( 66) = p_iso ; cbmz_per_megan( 66) = 3.
    p_of_megan2cbmz( 67) = is_selinene_d ; p_of_cbmz( 67) = p_iso ; cbmz_per_megan( 67) = 3.
    p_of_megan2cbmz( 68) = is_nerolidol_c ; p_of_cbmz( 68) = p_iso ; cbmz_per_megan( 68) = 3.
    p_of_megan2cbmz( 69) = is_nerolidol_t ; p_of_cbmz( 69) = p_iso ; cbmz_per_megan( 69) = 3.
    p_of_megan2cbmz( 70) = is_cedrol ; p_of_cbmz( 70) = p_iso ; cbmz_per_megan( 70) = 3.
    p_of_megan2cbmz( 71) = is_mbo_2m3e2ol ; p_of_cbmz( 71) = p_olt ; cbmz_per_megan( 71) = 1.
    p_of_megan2cbmz( 72) = is_mbo_2m3e2ol ; p_of_cbmz( 72) = p_par ; cbmz_per_megan( 72) = 3.
    p_of_megan2cbmz( 73) = is_methanol ; p_of_cbmz( 73) = p_ch3oh ; cbmz_per_megan( 73) = 1.
    p_of_megan2cbmz( 74) = is_acetone ; p_of_cbmz( 74) = p_ket ; cbmz_per_megan( 74) = 1.
    p_of_megan2cbmz( 75) = is_methane ; p_of_cbmz( 75) = p_ch4 ; cbmz_per_megan( 75) = 1
    p_of_megan2cbmz( 76) = is_ammonia ; p_of_cbmz( 76) = p_nh3 ; cbmz_per_megan( 76) = 1.
    p_of_megan2cbmz( 77) = is_nitrous_oxd ; p_of_cbmz( 77) = non_react ; cbmz_per_megan( 77) = 1.
    p_of_megan2cbmz( 78) = is_nitric_oxd ; p_of_cbmz( 78) = p_no ; cbmz_per_megan( 78) = 1.
    p_of_megan2cbmz( 79) = is_acetaldehyde ; p_of_cbmz( 79) = p_ald ; cbmz_per_megan( 79) = 1.
    p_of_megan2cbmz( 80) = is_ethanol ; p_of_cbmz( 80) = p_c2h5oh ; cbmz_per_megan( 80) = 1.
    p_of_megan2cbmz( 81) = is_formic_acid ; p_of_cbmz( 81) = p_ora1 ; cbmz_per_megan( 81) = 1.
    p_of_megan2cbmz( 82) = is_formaldehyde ; p_of_cbmz( 82) = p_hcho ; cbmz_per_megan( 82) = 1.
    p_of_megan2cbmz( 83) = is_acetic_acid ; p_of_cbmz( 83) = p_ora2 ; cbmz_per_megan( 83) = 1.
    p_of_megan2cbmz( 84) = is_mbo_3m2e1ol ; p_of_cbmz( 84) = p_ald ; cbmz_per_megan( 84) = 1.
    p_of_megan2cbmz( 85) = is_mbo_3m2e1ol ; p_of_cbmz( 85) = p_par ; cbmz_per_megan( 85) = 3.
    p_of_megan2cbmz( 86) = is_mbo_3m3e1ol ; p_of_cbmz( 86) = p_hcho ; cbmz_per_megan( 86) = 1.
    p_of_megan2cbmz( 87) = is_mbo_3m3e1ol ; p_of_cbmz( 87) = p_par ; cbmz_per_megan( 87) = 4.
    p_of_megan2cbmz( 88) = is_benzaldehyde ; p_of_cbmz( 88) = p_tol ; cbmz_per_megan( 88) = 1.
    p_of_megan2cbmz( 89) = is_butanone_2 ; p_of_cbmz( 89) = p_ket ; cbmz_per_megan( 89) = 1.
    p_of_megan2cbmz( 90) = is_butanone_2 ; p_of_cbmz( 90) = p_par ; cbmz_per_megan( 90) = 1.
    p_of_megan2cbmz( 91) = is_decanal ; p_of_cbmz( 91) = p_ald ; cbmz_per_megan( 91) = 1.
    p_of_megan2cbmz( 92) = is_decanal ; p_of_cbmz( 92) = p_par ; cbmz_per_megan( 92) = 8.
    p_of_megan2cbmz( 93) = is_dodecene_1 ; p_of_cbmz( 93) = p_olt ; cbmz_per_megan( 93) = 1.
    p_of_megan2cbmz( 94) = is_dodecene_1 ; p_of_cbmz( 94) = p_par ; cbmz_per_megan( 94) = 10.
    p_of_megan2cbmz( 95) = is_geranyl_acetone ; p_of_cbmz( 95) = p_iso ; cbmz_per_megan( 95) = 3.
    p_of_megan2cbmz( 96) = is_heptanal ; p_of_cbmz( 96) = p_ald ; cbmz_per_megan( 96) = 1.
    p_of_megan2cbmz( 97) = is_heptanal ; p_of_cbmz( 97) = p_par ; cbmz_per_megan( 97) = 5.
    p_of_megan2cbmz( 98) = is_heptane ; p_of_cbmz( 98) = p_par ; cbmz_per_megan( 98) = 7.
    p_of_megan2cbmz( 99) = is_hexane ; p_of_cbmz( 99) = p_par ; cbmz_per_megan( 99) = 6.
    p_of_megan2cbmz(100) = is_met_benzoate ; p_of_cbmz(100) = p_tol ; cbmz_per_megan(100) = 1.
    p_of_megan2cbmz(101) = is_met_heptenone ; p_of_cbmz(101) = p_ket ; cbmz_per_megan(101) = 1.
    p_of_megan2cbmz(102) = is_met_heptenone ; p_of_cbmz(102) = p_par ; cbmz_per_megan(102) = 3.
    p_of_megan2cbmz(103) = is_met_heptenone ; p_of_cbmz(103) = p_olt ; cbmz_per_megan(103) = 1.
    p_of_megan2cbmz(104) = is_neryl_acetone ; p_of_cbmz(104) = p_ket ; cbmz_per_megan(104) = 1.
    p_of_megan2cbmz(105) = is_neryl_acetone ; p_of_cbmz(105) = p_par ; cbmz_per_megan(105) = 8.
    p_of_megan2cbmz(106) = is_neryl_acetone ; p_of_cbmz(106) = p_oli ; cbmz_per_megan(106) = 2.
    p_of_megan2cbmz(107) = is_nonanal ; p_of_cbmz(107) = p_ald ; cbmz_per_megan(107) = 1.
    p_of_megan2cbmz(108) = is_nonanal ; p_of_cbmz(108) = p_par ; cbmz_per_megan(108) = 7.
    p_of_megan2cbmz(109) = is_nonenal ; p_of_cbmz(109) = p_ald ; cbmz_per_megan(109) = 1.
    p_of_megan2cbmz(110) = is_nonenal ; p_of_cbmz(110) = p_par ; cbmz_per_megan(110) = 6.
    p_of_megan2cbmz(111) = is_nonenal ; p_of_cbmz(111) = p_oli ; cbmz_per_megan(111) = 1.
    p_of_megan2cbmz(112) = is_octanal ; p_of_cbmz(112) = p_ald ; cbmz_per_megan(112) = 1.
    p_of_megan2cbmz(113) = is_octanal ; p_of_cbmz(113) = p_par ; cbmz_per_megan(113) = 6.
    p_of_megan2cbmz(114) = is_octanol ; p_of_cbmz(114) = p_par ; cbmz_per_megan(114) = 8.
    p_of_megan2cbmz(115) = is_octenol_1e3ol ; p_of_cbmz(115) = p_par ; cbmz_per_megan(115) = 6.
    p_of_megan2cbmz(116) = is_octenol_1e3ol ; p_of_cbmz(116) = p_olt ; cbmz_per_megan(116) = 1.
    p_of_megan2cbmz(117) = is_oxopentanal ; p_of_cbmz(117) = p_ald ; cbmz_per_megan(117) = 1.
    p_of_megan2cbmz(118) = is_oxopentanal ; p_of_cbmz(118) = p_par ; cbmz_per_megan(118) = 3.
    p_of_megan2cbmz(119) = is_pentane ; p_of_cbmz(119) = p_par ; cbmz_per_megan(119) = 5.
    p_of_megan2cbmz(120) = is_phenyl_cco ; p_of_cbmz(120) = p_ald ; cbmz_per_megan(120) = 1
    p_of_megan2cbmz(121) = is_phenyl_cco ; p_of_cbmz(121) = p_tol ; cbmz_per_megan(121) = 1.
    p_of_megan2cbmz(122) = is_pyruvic_acid ; p_of_cbmz(122) = p_ora1 ; cbmz_per_megan(122) = 1.
    p_of_megan2cbmz(123) = is_pyruvic_acid ; p_of_cbmz(123) = p_ket ; cbmz_per_megan(123) = 1.
    p_of_megan2cbmz(124) = is_terpinyl_act_a ; p_of_cbmz(124) = p_iso ; cbmz_per_megan(124) = 2.
    p_of_megan2cbmz(125) = is_tetradecene_1 ; p_of_cbmz(125) = p_par ; cbmz_per_megan(125) = 12.
    p_of_megan2cbmz(126) = is_tetradecene_1 ; p_of_cbmz(126) = p_olt ; cbmz_per_megan(126) = 1.
    p_of_megan2cbmz(127) = is_toluene ; p_of_cbmz(127) = p_tol ; cbmz_per_megan(127) = 1.
    p_of_megan2cbmz(128) = is_carbon_monoxide ; p_of_cbmz(128) = p_co ; cbmz_per_megan(128) = 1.
    p_of_megan2cbmz(129) = is_butene ; p_of_cbmz(129) = p_olt ; cbmz_per_megan(129) = 1.
    p_of_megan2cbmz(130) = is_butene ; p_of_cbmz(130) = p_par ; cbmz_per_megan(130) = 2.
    p_of_megan2cbmz(131) = is_ethane ; p_of_cbmz(131) = p_eth ; cbmz_per_megan(131) = 1.
    p_of_megan2cbmz(132) = is_ethene ; p_of_cbmz(132) = p_ol2 ; cbmz_per_megan(132) = 1.
    p_of_megan2cbmz(133) = is_hydrogen_cyanide ; p_of_cbmz(133) = non_react ; cbmz_per_megan(133) = 1.
    p_of_megan2cbmz(134) = is_propane ; p_of_cbmz(134) = p_par ; cbmz_per_megan(134) = 3.
    p_of_megan2cbmz(135) = is_propene ; p_of_cbmz(135) = p_olt ; cbmz_per_megan(135) = 1.
    p_of_megan2cbmz(136) = is_propene ; p_of_cbmz(136) = p_par ; cbmz_per_megan(136) = 1.
    p_of_megan2cbmz(137) = is_carbon_2s ; p_of_cbmz(137) = non_react ; cbmz_per_megan(137) = 1.
    p_of_megan2cbmz(138) = is_carbonyl_s ; p_of_cbmz(138) = non_react ; cbmz_per_megan(138) = 1.
    p_of_megan2cbmz(139) = is_diallyl_2s ; p_of_cbmz(139) = p_dms ; cbmz_per_megan(139) = 1.
    p_of_megan2cbmz(140) = is_diallyl_2s ; p_of_cbmz(140) = p_par ; cbmz_per_megan(140) = 2.
    p_of_megan2cbmz(141) = is_diallyl_2s ; p_of_cbmz(141) = p_olt ; cbmz_per_megan(141) = 2.
    p_of_megan2cbmz(142) = is_2met_2s ; p_of_cbmz(142) = p_dms ; cbmz_per_megan(142) = 1.
    p_of_megan2cbmz(143) = is_2met_s ; p_of_cbmz(143) = p_dms ; cbmz_per_megan(143) = 1.
    p_of_megan2cbmz(144) = is_met_chloride ; p_of_cbmz(144) = non_react ; cbmz_per_megan(144) = 1.
    p_of_megan2cbmz(145) = is_met_bromide ; p_of_cbmz(145) = non_react ; cbmz_per_megan(145) = 1.
    p_of_megan2cbmz(146) = is_met_iodide ; p_of_cbmz(146) = non_react ; cbmz_per_megan(146) = 1.
    p_of_megan2cbmz(147) = is_hydrogen_s ; p_of_cbmz(147) = non_react ; cbmz_per_megan(147) = 1.
    p_of_megan2cbmz(148) = is_met_mercaptan ; p_of_cbmz(148) = p_par ; cbmz_per_megan(148) = 1.
    p_of_megan2cbmz(149) = is_met_propenyl_2s ; p_of_cbmz(149) = p_dms ; cbmz_per_megan(149) = 1.
    p_of_megan2cbmz(150) = is_met_propenyl_2s ; p_of_cbmz(150) = p_oli ; cbmz_per_megan(150) = 1.
    p_of_megan2cbmz(151) = is_pppp_2s ; p_of_cbmz(151) = p_dms ; cbmz_per_megan(151) = 1.
    p_of_megan2cbmz(152) = is_pppp_2s ; p_of_cbmz(152) = p_par ; cbmz_per_megan(152) = 2.
    p_of_megan2cbmz(153) = is_pppp_2s ; p_of_cbmz(153) = p_oli ; cbmz_per_megan(153) = 1.
    p_of_megan2cbmz(154) = is_2met_nonatriene ; p_of_cbmz(154) = p_iso ; cbmz_per_megan(154) = 2.
    p_of_megan2cbmz(155) = is_met_salicylate ; p_of_cbmz(155) = p_tol ; cbmz_per_megan(155) = 1.
    p_of_megan2cbmz(156) = is_indole ; p_of_cbmz(156) = p_tol ; cbmz_per_megan(156) = 1.
    p_of_megan2cbmz(157) = is_jasmone ; p_of_cbmz(157) = p_iso ; cbmz_per_megan(157) = 2.
    p_of_megan2cbmz(158) = is_met_jasmonate ; p_of_cbmz(158) = p_iso ; cbmz_per_megan(158) = 3.
    p_of_megan2cbmz(159) = is_3met_3dctt ; p_of_cbmz(159) = p_iso ; cbmz_per_megan(159) = 3.
    p_of_megan2cbmz(160) = is_hexanal ; p_of_cbmz(160) = p_ald ; cbmz_per_megan(160) = 1.
    p_of_megan2cbmz(161) = is_hexanal ; p_of_cbmz(161) = p_par ; cbmz_per_megan(161) = 4.
    p_of_megan2cbmz(162) = is_hexanol_1 ; p_of_cbmz(162) = p_par ; cbmz_per_megan(162) = 6.
    p_of_megan2cbmz(163) = is_hexenal_c3 ; p_of_cbmz(163) = p_ald ; cbmz_per_megan(163) = 1.
    p_of_megan2cbmz(164) = is_hexenal_c3 ; p_of_cbmz(164) = p_par ; cbmz_per_megan(164) = 3.
    p_of_megan2cbmz(165) = is_hexenal_c3 ; p_of_cbmz(165) = p_oli ; cbmz_per_megan(165) = 1
    p_of_megan2cbmz(166) = is_hexenal_t2 ; p_of_cbmz(166) = p_ald ; cbmz_per_megan(166) = 1.
    p_of_megan2cbmz(167) = is_hexenal_t2 ; p_of_cbmz(167) = p_par ; cbmz_per_megan(167) = 6.
    p_of_megan2cbmz(168) = is_hexenal_t2 ; p_of_cbmz(168) = p_oli ; cbmz_per_megan(168) = 1.
    p_of_megan2cbmz(169) = is_hexenol_c3 ; p_of_cbmz(169) = p_par ; cbmz_per_megan(169) = 5.
    p_of_megan2cbmz(170) = is_hexenol_c3 ; p_of_cbmz(170) = p_oli ; cbmz_per_megan(170) = 1.
    p_of_megan2cbmz(171) = is_hexenyl_act_c3 ; p_of_cbmz(171) = p_ket ; cbmz_per_megan(171) = 1.
    p_of_megan2cbmz(172) = is_hexenyl_act_c3 ; p_of_cbmz(172) = p_par ; cbmz_per_megan(172) = 3.
    p_of_megan2cbmz(173) = is_hexenyl_act_c3 ; p_of_cbmz(173) = p_oli ; cbmz_per_megan(173) = 1.
  END SUBROUTINE get_megan2cbmz_table
  SUBROUTINE get_megan2radm2_table
    p_of_megan2radm2( 1) = is_isoprene ; p_of_radm2( 1) = p_iso ; radm2_per_megan( 1) = 1.000
    p_of_megan2radm2( 2) = is_myrcene ; p_of_radm2( 2) = p_olt ; radm2_per_megan( 2) = 0.5
    p_of_megan2radm2( 3) = is_myrcene ; p_of_radm2( 3) = p_oli ; radm2_per_megan( 3) = 0.5
    p_of_megan2radm2( 4) = is_sabinene ; p_of_radm2( 4) = p_olt ; radm2_per_megan( 4) = 1.000
    p_of_megan2radm2( 5) = is_limonene ; p_of_radm2( 5) = p_olt ; radm2_per_megan( 5) = 0.5
    p_of_megan2radm2( 6) = is_limonene ; p_of_radm2( 6) = p_oli ; radm2_per_megan( 6) = 0.5
    p_of_megan2radm2( 7) = is_carene_3 ; p_of_radm2( 7) = p_oli ; radm2_per_megan( 7) = 1.000
    p_of_megan2radm2( 8) = is_ocimene_t_b ; p_of_radm2( 8) = p_olt ; radm2_per_megan( 8) = 0.5
    p_of_megan2radm2( 9) = is_ocimene_t_b ; p_of_radm2( 9) = p_oli ; radm2_per_megan( 9) = 0.5
    p_of_megan2radm2( 10) = is_pinene_b ; p_of_radm2( 10) = p_olt ; radm2_per_megan( 10) = 1.000
    p_of_megan2radm2( 11) = is_pinene_a ; p_of_radm2( 11) = p_oli ; radm2_per_megan( 11) = 1.000
    p_of_megan2radm2( 12) = is_2met_styrene ; p_of_radm2( 12) = p_tol ; radm2_per_megan( 12) = 1.000
    p_of_megan2radm2( 13) = is_cymene_p ; p_of_radm2( 13) = p_tol ; radm2_per_megan( 13) = 1.000
    p_of_megan2radm2( 14) = is_cymene_o ; p_of_radm2( 14) = p_tol ; radm2_per_megan( 14) = 1.000
    p_of_megan2radm2( 15) = is_phellandrene_a ; p_of_radm2( 15) = p_oli ; radm2_per_megan( 15) = 1.000
    p_of_megan2radm2( 16) = is_thujene_a ; p_of_radm2( 16) = p_oli ; radm2_per_megan( 16) = 1.000
    p_of_megan2radm2( 17) = is_terpinene_a ; p_of_radm2( 17) = p_oli ; radm2_per_megan( 17) = 1.000
    p_of_megan2radm2( 18) = is_terpinene_g ; p_of_radm2( 18) = p_oli ; radm2_per_megan( 18) = 1.000
    p_of_megan2radm2( 19) = is_terpinolene ; p_of_radm2( 19) = p_oli ; radm2_per_megan( 19) = 1.000
    p_of_megan2radm2( 20) = is_phellandrene_b ; p_of_radm2( 20) = p_olt ; radm2_per_megan( 20) = 0.5
    p_of_megan2radm2( 21) = is_phellandrene_b ; p_of_radm2( 21) = p_oli ; radm2_per_megan( 21) = 0.5
    p_of_megan2radm2( 22) = is_camphene ; p_of_radm2( 22) = p_olt ; radm2_per_megan( 22) = 1.000
    p_of_megan2radm2( 23) = is_bornene ; p_of_radm2( 23) = p_oli ; radm2_per_megan( 23) = 1.000
    p_of_megan2radm2( 24) = is_fenchene_a ; p_of_radm2( 24) = p_olt ; radm2_per_megan( 24) = 1.000
    p_of_megan2radm2( 25) = is_ocimene_al ; p_of_radm2( 25) = p_oli ; radm2_per_megan( 25) = 1.000
    p_of_megan2radm2( 26) = is_ocimene_c_b ; p_of_radm2( 26) = p_olt ; radm2_per_megan( 26) = 0.5
    p_of_megan2radm2( 27) = is_ocimene_c_b ; p_of_radm2( 27) = p_oli ; radm2_per_megan( 27) = 0.5
    p_of_megan2radm2( 28) = is_tricyclene ; p_of_radm2( 28) = non_react ; radm2_per_megan( 28) = 1.000
    p_of_megan2radm2( 29) = is_estragole ; p_of_radm2( 29) = p_olt ; radm2_per_megan( 29) = 1.000
    p_of_megan2radm2( 30) = is_camphor ; p_of_radm2( 30) = p_hc8 ; radm2_per_megan( 30) = 0.388
    p_of_megan2radm2( 31) = is_fenchone ; p_of_radm2( 31) = non_react ; radm2_per_megan( 31) = 1.000
    p_of_megan2radm2( 32) = is_piperitone ; p_of_radm2( 32) = p_olt ; radm2_per_megan( 32) = 1.000
    p_of_megan2radm2( 33) = is_thujone_a ; p_of_radm2( 33) = non_react ; radm2_per_megan( 33) = 1.000
    p_of_megan2radm2( 34) = is_thujone_b ; p_of_radm2( 34) = non_react ; radm2_per_megan( 34) = 1.000
    p_of_megan2radm2( 35) = is_cineole_1_8 ; p_of_radm2( 35) = p_hc8 ; radm2_per_megan( 35) = 0.755
    p_of_megan2radm2( 36) = is_borneol ; p_of_radm2( 36) = non_react ; radm2_per_megan( 36) = 1.000
    p_of_megan2radm2( 37) = is_linalool ; p_of_radm2( 37) = p_olt ; radm2_per_megan( 37) = 0.5
    p_of_megan2radm2( 38) = is_linalool ; p_of_radm2( 38) = p_oli ; radm2_per_megan( 38) = 0.5
    p_of_megan2radm2( 39) = is_terpineol_4 ; p_of_radm2( 39) = p_oli ; radm2_per_megan( 39) = 1.000
    p_of_megan2radm2( 40) = is_terpineol_a ; p_of_radm2( 40) = p_oli ; radm2_per_megan( 40) = 1.000
    p_of_megan2radm2( 41) = is_linalool_oxd_c ; p_of_radm2( 41) = p_olt ; radm2_per_megan( 41) = 1.000
    p_of_megan2radm2( 42) = is_linalool_oxd_t ; p_of_radm2( 42) = p_olt ; radm2_per_megan( 42) = 1.000
    p_of_megan2radm2( 43) = is_ionone_b ; p_of_radm2( 43) = p_oli ; radm2_per_megan( 43) = 1.000
    p_of_megan2radm2( 44) = is_bornyl_act ; p_of_radm2( 44) = non_react ; radm2_per_megan( 44) = 1.000
    p_of_megan2radm2( 45) = is_farnescene_a ; p_of_radm2( 45) = p_olt ; radm2_per_megan( 45) = 0.5
    p_of_megan2radm2( 46) = is_farnescene_a ; p_of_radm2( 46) = p_oli ; radm2_per_megan( 46) = 0.5
    p_of_megan2radm2( 47) = is_caryophyllene_b ; p_of_radm2( 47) = p_olt ; radm2_per_megan( 47) = 0.5
    p_of_megan2radm2( 48) = is_caryophyllene_b ; p_of_radm2( 48) = p_oli ; radm2_per_megan( 48) = 0.5
    p_of_megan2radm2( 49) = is_acoradiene ; p_of_radm2( 49) = p_olt ; radm2_per_megan( 49) = 0.5
    p_of_megan2radm2( 50) = is_acoradiene ; p_of_radm2( 50) = p_oli ; radm2_per_megan( 50) = 0.5
    p_of_megan2radm2( 51) = is_aromadendrene ; p_of_radm2( 51) = p_olt ; radm2_per_megan( 51) = 1.000
    p_of_megan2radm2( 52) = is_bergamotene_a ; p_of_radm2( 52) = p_oli ; radm2_per_megan( 52) = 1.000
    p_of_megan2radm2( 53) = is_bergamotene_b ; p_of_radm2( 53) = p_olt ; radm2_per_megan( 53) = 0.5
    p_of_megan2radm2( 54) = is_bergamotene_b ; p_of_radm2( 54) = p_oli ; radm2_per_megan( 54) = 0.5
    p_of_megan2radm2( 55) = is_bisabolene_a ; p_of_radm2( 55) = p_oli ; radm2_per_megan( 55) = 1.000
    p_of_megan2radm2( 56) = is_bisabolene_b ; p_of_radm2( 56) = p_olt ; radm2_per_megan( 56) = 0.5
    p_of_megan2radm2( 57) = is_bisabolene_b ; p_of_radm2( 57) = p_oli ; radm2_per_megan( 57) = 0.5
    p_of_megan2radm2( 58) = is_bourbonene_b ; p_of_radm2( 58) = p_olt ; radm2_per_megan( 58) = 1.000
    p_of_megan2radm2( 59) = is_cadinene_d ; p_of_radm2( 59) = p_oli ; radm2_per_megan( 59) = 1.000
    p_of_megan2radm2( 60) = is_cadinene_g ; p_of_radm2( 60) = p_olt ; radm2_per_megan( 60) = 0.5
    p_of_megan2radm2( 61) = is_cadinene_g ; p_of_radm2( 61) = p_oli ; radm2_per_megan( 61) = 0.5
    p_of_megan2radm2( 62) = is_cedrene_a ; p_of_radm2( 62) = p_oli ; radm2_per_megan( 62) = 1.000
    p_of_megan2radm2( 63) = is_copaene_a ; p_of_radm2( 63) = p_oli ; radm2_per_megan( 63) = 1.000
    p_of_megan2radm2( 64) = is_cubebene_a ; p_of_radm2( 64) = p_oli ; radm2_per_megan( 64) = 1.000
    p_of_megan2radm2( 65) = is_cubebene_b ; p_of_radm2( 65) = p_olt ; radm2_per_megan( 65) = 1.000
    p_of_megan2radm2( 66) = is_elemene_b ; p_of_radm2( 66) = p_olt ; radm2_per_megan( 66) = 1.000
    p_of_megan2radm2( 67) = is_farnescene_b ; p_of_radm2( 67) = p_olt ; radm2_per_megan( 67) = 0.5
    p_of_megan2radm2( 68) = is_farnescene_b ; p_of_radm2( 68) = p_oli ; radm2_per_megan( 68) = 0.5
    p_of_megan2radm2( 69) = is_germacrene_b ; p_of_radm2( 69) = p_oli ; radm2_per_megan( 69) = 1.000
    p_of_megan2radm2( 70) = is_germacrene_d ; p_of_radm2( 70) = p_olt ; radm2_per_megan( 70) = 0.5
    p_of_megan2radm2( 71) = is_germacrene_d ; p_of_radm2( 71) = p_oli ; radm2_per_megan( 71) = 0.5
    p_of_megan2radm2( 72) = is_gurjunene_b ; p_of_radm2( 72) = p_oli ; radm2_per_megan( 72) = 1.000
    p_of_megan2radm2( 73) = is_humulene_a ; p_of_radm2( 73) = p_oli ; radm2_per_megan( 73) = 1.000
    p_of_megan2radm2( 74) = is_humulene_g ; p_of_radm2( 74) = p_oli ; radm2_per_megan( 74) = 1.000
    p_of_megan2radm2( 75) = is_isolongifolene ; p_of_radm2( 75) = p_oli ; radm2_per_megan( 75) = 1.000
    p_of_megan2radm2( 76) = is_longifolene ; p_of_radm2( 76) = p_olt ; radm2_per_megan( 76) = 1.000
    p_of_megan2radm2( 77) = is_longipinene ; p_of_radm2( 77) = p_oli ; radm2_per_megan( 77) = 1.000
    p_of_megan2radm2( 78) = is_muurolene_a ; p_of_radm2( 78) = p_oli ; radm2_per_megan( 78) = 1.000
    p_of_megan2radm2( 79) = is_muurolene_g ; p_of_radm2( 79) = p_olt ; radm2_per_megan( 79) = 0.5
    p_of_megan2radm2( 80) = is_muurolene_g ; p_of_radm2( 80) = p_oli ; radm2_per_megan( 80) = 0.5
    p_of_megan2radm2( 81) = is_selinene_b ; p_of_radm2( 81) = p_olt ; radm2_per_megan( 81) = 1.000
    p_of_megan2radm2( 82) = is_selinene_d ; p_of_radm2( 82) = p_oli ; radm2_per_megan( 82) = 1.000
    p_of_megan2radm2( 83) = is_nerolidol_c ; p_of_radm2( 83) = p_olt ; radm2_per_megan( 83) = 0.5
    p_of_megan2radm2( 84) = is_nerolidol_c ; p_of_radm2( 84) = p_oli ; radm2_per_megan( 84) = 0.5
    p_of_megan2radm2( 85) = is_nerolidol_t ; p_of_radm2( 85) = p_olt ; radm2_per_megan( 85) = 0.5
    p_of_megan2radm2( 86) = is_nerolidol_t ; p_of_radm2( 86) = p_oli ; radm2_per_megan( 86) = 0.5
    p_of_megan2radm2( 87) = is_cedrol ; p_of_radm2( 87) = non_react ; radm2_per_megan( 87) = 1.000
    p_of_megan2radm2( 88) = is_mbo_2m3e2ol ; p_of_radm2( 88) = p_iso ; radm2_per_megan( 88) = 1.000
    p_of_megan2radm2( 89) = is_methanol ; p_of_radm2( 89) = p_hc3 ; radm2_per_megan( 89) = 0.402
    p_of_megan2radm2( 90) = is_acetone ; p_of_radm2( 90) = p_ket ; radm2_per_megan( 90) = 0.253
    p_of_megan2radm2( 91) = is_methane ; p_of_radm2( 91) = p_ch4 ; radm2_per_megan( 91) = 1.000
    p_of_megan2radm2( 92) = is_ammonia ; p_of_radm2( 92) = p_nh3 ; radm2_per_megan( 92) = 1.000
    p_of_megan2radm2( 93) = is_nitrous_oxd ; p_of_radm2( 93) = p_no2 ; radm2_per_megan( 93) = 1.000
    p_of_megan2radm2( 94) = is_nitric_oxd ; p_of_radm2( 94) = p_no ; radm2_per_megan( 94) = 1.000
    p_of_megan2radm2( 95) = is_acetaldehyde ; p_of_radm2( 95) = p_ald ; radm2_per_megan( 95) = 1.000
    p_of_megan2radm2( 96) = is_ethanol ; p_of_radm2( 96) = p_hc3 ; radm2_per_megan( 96) = 1.198
    p_of_megan2radm2( 97) = is_formic_acid ; p_of_radm2( 97) = p_ora1 ; radm2_per_megan( 97) = 1.000
    p_of_megan2radm2( 98) = is_formaldehyde ; p_of_radm2( 98) = p_hcho ; radm2_per_megan( 98) = 1.000
    p_of_megan2radm2( 99) = is_acetic_acid ; p_of_radm2( 99) = p_ora2 ; radm2_per_megan( 99) = 1.000
    p_of_megan2radm2(100) = is_mbo_3m2e1ol ; p_of_radm2(100) = p_iso ; radm2_per_megan(100) = 1.000
    p_of_megan2radm2(101) = is_mbo_3m3e1ol ; p_of_radm2(101) = p_iso ; radm2_per_megan(101) = 1.000
    p_of_megan2radm2(102) = is_benzaldehyde ; p_of_radm2(102) = non_react ; radm2_per_megan(102) = 1.000
    p_of_megan2radm2(103) = is_butanone_2 ; p_of_radm2(103) = p_ket ; radm2_per_megan(103) = 1.000
    p_of_megan2radm2(104) = is_decanal ; p_of_radm2(104) = p_ald ; radm2_per_megan(104) = 1.000
    p_of_megan2radm2(105) = is_dodecene_1 ; p_of_radm2(105) = p_olt ; radm2_per_megan(105) = 1.000
    p_of_megan2radm2(106) = is_geranyl_acetone ; p_of_radm2(106) = p_oli ; radm2_per_megan(106) = 1.000
    p_of_megan2radm2(107) = is_heptanal ; p_of_radm2(107) = p_ald ; radm2_per_megan(107) = 1.000
    p_of_megan2radm2(108) = is_heptane ; p_of_radm2(108) = p_hc5 ; radm2_per_megan(108) = 1.226
    p_of_megan2radm2(109) = is_hexane ; p_of_radm2(109) = p_hc5 ; radm2_per_megan(109) = 1.049
    p_of_megan2radm2(110) = is_met_benzoate ; p_of_radm2(110) = p_hc8 ; radm2_per_megan(110) = 1.000
    p_of_megan2radm2(111) = is_met_heptenone ; p_of_radm2(111) = p_oli ; radm2_per_megan(111) = 1.000
    p_of_megan2radm2(112) = is_neryl_acetone ; p_of_radm2(112) = p_oli ; radm2_per_megan(112) = 1.000
    p_of_megan2radm2(113) = is_nonanal ; p_of_radm2(113) = p_ald ; radm2_per_megan(113) = 1.000
    p_of_megan2radm2(114) = is_nonenal ; p_of_radm2(114) = p_ald ; radm2_per_megan(114) = 1.000
    p_of_megan2radm2(115) = is_nonenal ; p_of_radm2(115) = p_hc8 ; radm2_per_megan(115) = 1.000
    p_of_megan2radm2(116) = is_octanal ; p_of_radm2(116) = p_ald ; radm2_per_megan(116) = 1.000
    p_of_megan2radm2(117) = is_octanol ; p_of_radm2(117) = p_hc8 ; radm2_per_megan(117) = 1.119
    p_of_megan2radm2(118) = is_octenol_1e3ol ; p_of_radm2(118) = p_olt ; radm2_per_megan(118) = 1.000
    p_of_megan2radm2(119) = is_oxopentanal ; p_of_radm2(119) = p_ald ; radm2_per_megan(119) = 1.000
    p_of_megan2radm2(120) = is_pentane ; p_of_radm2(120) = p_hc5 ; radm2_per_megan(120) = 0.847
    p_of_megan2radm2(121) = is_phenyl_cco ; p_of_radm2(121) = non_react ; radm2_per_megan(121) = 1.000
    p_of_megan2radm2(122) = is_pyruvic_acid ; p_of_radm2(122) = p_ora2 ; radm2_per_megan(122) = 1.000
    p_of_megan2radm2(123) = is_terpinyl_act_a ; p_of_radm2(123) = p_oli ; radm2_per_megan(123) = 1.000
    p_of_megan2radm2(124) = is_tetradecene_1 ; p_of_radm2(124) = p_olt ; radm2_per_megan(124) = 1.000
    p_of_megan2radm2(125) = is_toluene ; p_of_radm2(125) = p_tol ; radm2_per_megan(125) = 1.000
    p_of_megan2radm2(126) = is_carbon_monoxide ; p_of_radm2(126) = p_co ; radm2_per_megan(126) = 1.000
    p_of_megan2radm2(127) = is_butene ; p_of_radm2(127) = p_olt ; radm2_per_megan(127) = 1.000
    p_of_megan2radm2(128) = is_ethane ; p_of_radm2(128) = p_eth ; radm2_per_megan(128) = 1.000
    p_of_megan2radm2(129) = is_ethene ; p_of_radm2(129) = p_ol2 ; radm2_per_megan(129) = 1.000
    p_of_megan2radm2(130) = is_hydrogen_cyanide; p_of_radm2(130) = non_react ; radm2_per_megan(130) = 1.000
    p_of_megan2radm2(131) = is_propane ; p_of_radm2(131) = p_hc3 ; radm2_per_megan(131) = 0.519
    p_of_megan2radm2(132) = is_propene ; p_of_radm2(132) = p_olt ; radm2_per_megan(132) = 1.000
    p_of_megan2radm2(133) = is_carbon_2s ; p_of_radm2(133) = non_react ; radm2_per_megan(133) = 1.000
    p_of_megan2radm2(134) = is_carbonyl_s ; p_of_radm2(134) = non_react ; radm2_per_megan(134) = 1.000
    p_of_megan2radm2(135) = is_diallyl_2s ; p_of_radm2(135) = p_oli ; radm2_per_megan(135) = 1.000
    p_of_megan2radm2(136) = is_diallyl_2s ; p_of_radm2(136) = p_so2 ; radm2_per_megan(136) = 2.000
    p_of_megan2radm2(137) = is_2met_2s ; p_of_radm2(137) = p_eth ; radm2_per_megan(137) = 1.000
    p_of_megan2radm2(138) = is_2met_2s ; p_of_radm2(138) = p_so2 ; radm2_per_megan(138) = 2.000
    p_of_megan2radm2(139) = is_2met_s ; p_of_radm2(139) = p_eth ; radm2_per_megan(139) = 1.000
    p_of_megan2radm2(140) = is_2met_s ; p_of_radm2(140) = p_so2 ; radm2_per_megan(140) = 1.000
    p_of_megan2radm2(141) = is_met_chloride ; p_of_radm2(141) = non_react ; radm2_per_megan(141) = 1.000
    p_of_megan2radm2(142) = is_met_bromide ; p_of_radm2(142) = non_react ; radm2_per_megan(142) = 1.000
    p_of_megan2radm2(143) = is_met_iodide ; p_of_radm2(143) = non_react ; radm2_per_megan(143) = 1.000
    p_of_megan2radm2(144) = is_hydrogen_s ; p_of_radm2(144) = p_so2 ; radm2_per_megan(144) = 1.000
    p_of_megan2radm2(145) = is_met_mercaptan ; p_of_radm2(145) = p_ch4 ; radm2_per_megan(145) = 1.000
    p_of_megan2radm2(146) = is_met_mercaptan ; p_of_radm2(146) = p_so2 ; radm2_per_megan(146) = 1.000
    p_of_megan2radm2(147) = is_met_propenyl_2s ; p_of_radm2(147) = p_oli ; radm2_per_megan(147) = 1.000
    p_of_megan2radm2(148) = is_met_propenyl_2s ; p_of_radm2(148) = p_so2 ; radm2_per_megan(148) = 2.000
    p_of_megan2radm2(149) = is_pppp_2s ; p_of_radm2(149) = p_oli ; radm2_per_megan(149) = 1.000
    p_of_megan2radm2(150) = is_pppp_2s ; p_of_radm2(150) = p_so2 ; radm2_per_megan(150) = 2.000
    p_of_megan2radm2(151) = is_2met_nonatriene ; p_of_radm2(151) = p_olt ; radm2_per_megan(151) = 0.500
    p_of_megan2radm2(152) = is_2met_nonatriene ; p_of_radm2(152) = p_oli ; radm2_per_megan(152) = 0.500
    p_of_megan2radm2(153) = is_met_salicylate ; p_of_radm2(153) = p_hc8 ; radm2_per_megan(153) = 1.000
    p_of_megan2radm2(154) = is_indole ; p_of_radm2(154) = p_hc8 ; radm2_per_megan(154) = 1.238
    p_of_megan2radm2(155) = is_indole ; p_of_radm2(155) = p_hno3 ; radm2_per_megan(155) = 1.000
    p_of_megan2radm2(156) = is_jasmone ; p_of_radm2(156) = p_oli ; radm2_per_megan(156) = 1.000
    p_of_megan2radm2(157) = is_met_jasmonate ; p_of_radm2(157) = p_oli ; radm2_per_megan(157) = 1.000
    p_of_megan2radm2(158) = is_3met_3dctt ; p_of_radm2(158) = p_olt ; radm2_per_megan(158) = 0.500
    p_of_megan2radm2(159) = is_3met_3dctt ; p_of_radm2(159) = p_oli ; radm2_per_megan(159) = 0.500
    p_of_megan2radm2(160) = is_hexanal ; p_of_radm2(160) = p_ald ; radm2_per_megan(160) = 1.000
    p_of_megan2radm2(161) = is_hexanol_1 ; p_of_radm2(161) = p_hc5 ; radm2_per_megan(161) = 1.697
    p_of_megan2radm2(162) = is_hexenal_c3 ; p_of_radm2(162) = p_oli ; radm2_per_megan(162) = 1.000
    p_of_megan2radm2(163) = is_hexenal_t2 ; p_of_radm2(163) = p_oli ; radm2_per_megan(163) = 1.000
    p_of_megan2radm2(164) = is_hexenol_c3 ; p_of_radm2(164) = p_olt ; radm2_per_megan(164) = 1.000
    p_of_megan2radm2(165) = is_hexenyl_act_c3 ; p_of_radm2(165) = p_oli ; radm2_per_megan(165) = 1.000
  END SUBROUTINE get_megan2radm2_table
  SUBROUTINE get_megan2racm_table
    p_of_megan2racm( 1) = is_isoprene ; p_of_racm( 1) = p_iso ; racm_per_megan( 1) = 1.000
    p_of_megan2racm( 2) = is_myrcene ; p_of_racm( 2) = p_lim ; racm_per_megan( 2) = 1.000
    p_of_megan2racm( 3) = is_sabinene ; p_of_racm( 3) = p_api ; racm_per_megan( 3) = 1.000
    p_of_megan2racm( 4) = is_limonene ; p_of_racm( 4) = p_lim ; racm_per_megan( 4) = 1.000
    p_of_megan2racm( 5) = is_carene_3 ; p_of_racm( 5) = p_api ; racm_per_megan( 5) = 1.000
    p_of_megan2racm( 6) = is_ocimene_t_b ; p_of_racm( 6) = p_lim ; racm_per_megan( 6) = 1.000
    p_of_megan2racm( 7) = is_pinene_b ; p_of_racm( 7) = p_api ; racm_per_megan( 7) = 1.000
    p_of_megan2racm( 8) = is_pinene_a ; p_of_racm( 8) = p_api ; racm_per_megan( 8) = 1.000
    p_of_megan2racm( 9) = is_2met_styrene ; p_of_racm( 9) = p_lim ; racm_per_megan( 9) = 1.000
    p_of_megan2racm( 10) = is_cymene_p ; p_of_racm( 10) = p_lim ; racm_per_megan( 10) = 1.000
    p_of_megan2racm( 11) = is_cymene_o ; p_of_racm( 11) = p_lim ; racm_per_megan( 11) = 1.000
    p_of_megan2racm( 12) = is_phellandrene_a ; p_of_racm( 12) = p_lim ; racm_per_megan( 12) = 1.000
    p_of_megan2racm( 13) = is_thujene_a ; p_of_racm( 13) = p_api ; racm_per_megan( 13) = 1.000
    p_of_megan2racm( 14) = is_terpinene_a ; p_of_racm( 14) = p_lim ; racm_per_megan( 14) = 1.000
    p_of_megan2racm( 15) = is_terpinene_g ; p_of_racm( 15) = p_lim ; racm_per_megan( 15) = 1.000
    p_of_megan2racm( 16) = is_terpinolene ; p_of_racm( 16) = p_lim ; racm_per_megan( 16) = 1.000
    p_of_megan2racm( 17) = is_phellandrene_b ; p_of_racm( 17) = p_lim ; racm_per_megan( 17) = 1.000
    p_of_megan2racm( 18) = is_camphene ; p_of_racm( 18) = p_api ; racm_per_megan( 18) = 1.000
    p_of_megan2racm( 19) = is_bornene ; p_of_racm( 19) = p_api ; racm_per_megan( 19) = 1.000
    p_of_megan2racm( 20) = is_fenchene_a ; p_of_racm( 20) = p_api ; racm_per_megan( 20) = 1.000
    p_of_megan2racm( 21) = is_ocimene_al ; p_of_racm( 21) = p_lim ; racm_per_megan( 21) = 1.000
    p_of_megan2racm( 22) = is_ocimene_c_b ; p_of_racm( 22) = p_lim ; racm_per_megan( 22) = 1.000
    p_of_megan2racm( 23) = is_tricyclene ; p_of_racm( 23) = non_react ; racm_per_megan( 23) = 1.000
    p_of_megan2racm( 24) = is_estragole ; p_of_racm( 24) = p_lim ; racm_per_megan( 24) = 1.000
    p_of_megan2racm( 25) = is_camphor ; p_of_racm( 25) = p_hc8 ; racm_per_megan( 25) = 0.380
    p_of_megan2racm( 26) = is_fenchone ; p_of_racm( 26) = non_react ; racm_per_megan( 26) = 1.000
    p_of_megan2racm( 27) = is_piperitone ; p_of_racm( 27) = p_api ; racm_per_megan( 27) = 1.000
    p_of_megan2racm( 28) = is_thujone_a ; p_of_racm( 28) = non_react ; racm_per_megan( 28) = 1.000
    p_of_megan2racm( 29) = is_thujone_b ; p_of_racm( 29) = non_react ; racm_per_megan( 29) = 1.000
    p_of_megan2racm( 30) = is_cineole_1_8 ; p_of_racm( 30) = p_hc8 ; racm_per_megan( 30) = 0.738
    p_of_megan2racm( 31) = is_borneol ; p_of_racm( 31) = non_react ; racm_per_megan( 31) = 1.000
    p_of_megan2racm( 32) = is_linalool ; p_of_racm( 32) = p_lim ; racm_per_megan( 32) = 1.000
    p_of_megan2racm( 33) = is_terpineol_4 ; p_of_racm( 33) = p_api ; racm_per_megan( 33) = 1.000
    p_of_megan2racm( 34) = is_terpineol_a ; p_of_racm( 34) = p_api ; racm_per_megan( 34) = 1.000
    p_of_megan2racm( 35) = is_linalool_oxd_c ; p_of_racm( 35) = p_lim ; racm_per_megan( 35) = 1.000
    p_of_megan2racm( 36) = is_linalool_oxd_t ; p_of_racm( 36) = p_lim ; racm_per_megan( 36) = 1.000
    p_of_megan2racm( 37) = is_ionone_b ; p_of_racm( 37) = p_lim ; racm_per_megan( 37) = 1.000
    p_of_megan2racm( 38) = is_bornyl_act ; p_of_racm( 38) = non_react ; racm_per_megan( 38) = 1.000
    p_of_megan2racm( 39) = is_farnescene_a ; p_of_racm( 39) = p_lim ; racm_per_megan( 39) = 1.000
    p_of_megan2racm( 40) = is_caryophyllene_b ; p_of_racm( 40) = p_lim ; racm_per_megan( 40) = 1.000
    p_of_megan2racm( 41) = is_acoradiene ; p_of_racm( 41) = p_lim ; racm_per_megan( 41) = 1.000
    p_of_megan2racm( 42) = is_aromadendrene ; p_of_racm( 42) = p_api ; racm_per_megan( 42) = 1.000
    p_of_megan2racm( 43) = is_bergamotene_a ; p_of_racm( 43) = p_lim ; racm_per_megan( 43) = 1.000
    p_of_megan2racm( 44) = is_bergamotene_b ; p_of_racm( 44) = p_lim ; racm_per_megan( 44) = 1.000
    p_of_megan2racm( 45) = is_bisabolene_a ; p_of_racm( 45) = p_lim ; racm_per_megan( 45) = 1.000
    p_of_megan2racm( 46) = is_bisabolene_b ; p_of_racm( 46) = p_lim ; racm_per_megan( 46) = 1.000
    p_of_megan2racm( 47) = is_bourbonene_b ; p_of_racm( 47) = p_api ; racm_per_megan( 47) = 1.000
    p_of_megan2racm( 48) = is_cadinene_d ; p_of_racm( 48) = p_lim ; racm_per_megan( 48) = 1.000
    p_of_megan2racm( 49) = is_cadinene_g ; p_of_racm( 49) = p_lim ; racm_per_megan( 49) = 1.000
    p_of_megan2racm( 50) = is_cedrene_a ; p_of_racm( 50) = p_api ; racm_per_megan( 50) = 1.000
    p_of_megan2racm( 51) = is_copaene_a ; p_of_racm( 51) = p_api ; racm_per_megan( 51) = 1.000
    p_of_megan2racm( 52) = is_cubebene_a ; p_of_racm( 52) = p_api ; racm_per_megan( 52) = 1.000
    p_of_megan2racm( 53) = is_cubebene_b ; p_of_racm( 53) = p_api ; racm_per_megan( 53) = 1.000
    p_of_megan2racm( 54) = is_elemene_b ; p_of_racm( 54) = p_lim ; racm_per_megan( 54) = 1.000
    p_of_megan2racm( 55) = is_farnescene_b ; p_of_racm( 55) = p_lim ; racm_per_megan( 55) = 1.000
    p_of_megan2racm( 56) = is_germacrene_b ; p_of_racm( 56) = p_lim ; racm_per_megan( 56) = 1.000
    p_of_megan2racm( 57) = is_germacrene_d ; p_of_racm( 57) = p_lim ; racm_per_megan( 57) = 1.000
    p_of_megan2racm( 58) = is_gurjunene_b ; p_of_racm( 58) = p_api ; racm_per_megan( 58) = 1.000
    p_of_megan2racm( 59) = is_humulene_a ; p_of_racm( 59) = p_lim ; racm_per_megan( 59) = 1.000
    p_of_megan2racm( 60) = is_humulene_g ; p_of_racm( 60) = p_lim ; racm_per_megan( 60) = 1.000
    p_of_megan2racm( 61) = is_isolongifolene ; p_of_racm( 61) = p_api ; racm_per_megan( 61) = 1.000
    p_of_megan2racm( 62) = is_longifolene ; p_of_racm( 62) = p_api ; racm_per_megan( 62) = 1.000
    p_of_megan2racm( 63) = is_longipinene ; p_of_racm( 63) = p_api ; racm_per_megan( 63) = 1.000
    p_of_megan2racm( 64) = is_muurolene_a ; p_of_racm( 64) = p_lim ; racm_per_megan( 64) = 1.000
    p_of_megan2racm( 65) = is_muurolene_g ; p_of_racm( 65) = p_lim ; racm_per_megan( 65) = 1.000
    p_of_megan2racm( 66) = is_selinene_b ; p_of_racm( 66) = p_lim ; racm_per_megan( 66) = 1.000
    p_of_megan2racm( 67) = is_selinene_d ; p_of_racm( 67) = p_lim ; racm_per_megan( 67) = 1.000
    p_of_megan2racm( 68) = is_nerolidol_c ; p_of_racm( 68) = p_lim ; racm_per_megan( 68) = 1.000
    p_of_megan2racm( 69) = is_nerolidol_t ; p_of_racm( 69) = p_lim ; racm_per_megan( 69) = 1.000
    p_of_megan2racm( 70) = is_cedrol ; p_of_racm( 70) = non_react ; racm_per_megan( 70) = 1.000
    p_of_megan2racm( 71) = is_mbo_2m3e2ol ; p_of_racm( 71) = p_iso ; racm_per_megan( 71) = 1.000
    p_of_megan2racm( 72) = is_methanol ; p_of_racm( 72) = p_hc3 ; racm_per_megan( 72) = 0.490
    p_of_megan2racm( 73) = is_acetone ; p_of_racm( 73) = p_ket ; racm_per_megan( 73) = 0.330
    p_of_megan2racm( 74) = is_methane ; p_of_racm( 74) = p_ch4 ; racm_per_megan( 74) = 1.000
    p_of_megan2racm( 75) = is_ammonia ; p_of_racm( 75) = p_nh3 ; racm_per_megan( 75) = 1.000
    p_of_megan2racm( 76) = is_nitrous_oxd ; p_of_racm( 76) = p_no2 ; racm_per_megan( 76) = 1.000
    p_of_megan2racm( 77) = is_nitric_oxd ; p_of_racm( 77) = p_no ; racm_per_megan( 77) = 1.000
    p_of_megan2racm( 78) = is_acetaldehyde ; p_of_racm( 78) = p_ald ; racm_per_megan( 78) = 1.000
    p_of_megan2racm( 79) = is_ethanol ; p_of_racm( 79) = p_hc3 ; racm_per_megan( 79) = 1.370
    p_of_megan2racm( 80) = is_formic_acid ; p_of_racm( 80) = p_ora1 ; racm_per_megan( 80) = 1.000
    p_of_megan2racm( 81) = is_formaldehyde ; p_of_racm( 81) = p_hcho ; racm_per_megan( 81) = 1.000
    p_of_megan2racm( 82) = is_acetic_acid ; p_of_racm( 82) = p_ora2 ; racm_per_megan( 82) = 1.000
    p_of_megan2racm( 83) = is_mbo_3m2e1ol ; p_of_racm( 83) = p_iso ; racm_per_megan( 83) = 1.000
    p_of_megan2racm( 84) = is_mbo_3m3e1ol ; p_of_racm( 84) = p_iso ; racm_per_megan( 84) = 1.000
    p_of_megan2racm( 85) = is_benzaldehyde ; p_of_racm( 85) = non_react ; racm_per_megan( 85) = 1.000
    p_of_megan2racm( 86) = is_butanone_2 ; p_of_racm( 86) = p_ket ; racm_per_megan( 86) = 1.610
    p_of_megan2racm( 87) = is_decanal ; p_of_racm( 87) = p_ald ; racm_per_megan( 87) = 1.000
    p_of_megan2racm( 88) = is_dodecene_1 ; p_of_racm( 88) = p_olt ; racm_per_megan( 88) = 1.000
    p_of_megan2racm( 89) = is_geranyl_acetone ; p_of_racm( 89) = p_lim ; racm_per_megan( 89) = 1.000
    p_of_megan2racm( 90) = is_heptanal ; p_of_racm( 90) = p_ald ; racm_per_megan( 90) = 1.000
    p_of_megan2racm( 91) = is_heptane ; p_of_racm( 91) = p_hc5 ; racm_per_megan( 91) = 1.236
    p_of_megan2racm( 92) = is_hexane ; p_of_racm( 92) = p_hc5 ; racm_per_megan( 92) = 1.058
    p_of_megan2racm( 93) = is_met_benzoate ; p_of_racm( 93) = p_hc8 ; racm_per_megan( 93) = 1.000
    p_of_megan2racm( 94) = is_met_heptenone ; p_of_racm( 94) = p_oli ; racm_per_megan( 94) = 1.000
    p_of_megan2racm( 95) = is_neryl_acetone ; p_of_racm( 95) = p_oli ; racm_per_megan( 95) = 1.000
    p_of_megan2racm( 96) = is_nonanal ; p_of_racm( 96) = p_ald ; racm_per_megan( 96) = 1.000
    p_of_megan2racm( 97) = is_nonenal ; p_of_racm( 97) = p_hc8 ; racm_per_megan( 97) = 1.000
    p_of_megan2racm( 98) = is_nonenal ; p_of_racm( 98) = p_ald ; racm_per_megan( 98) = 1.000
    p_of_megan2racm( 99) = is_octanal ; p_of_racm( 99) = p_ald ; racm_per_megan( 99) = 1.000
    p_of_megan2racm(100) = is_octanol ; p_of_racm(100) = p_hc8 ; racm_per_megan(100) = 1.092
    p_of_megan2racm(101) = is_octenol_1e3ol ; p_of_racm(101) = p_olt ; racm_per_megan(101) = 1.000
    p_of_megan2racm(102) = is_oxopentanal ; p_of_racm(102) = p_ald ; racm_per_megan(102) = 1.000
    p_of_megan2racm(103) = is_pentane ; p_of_racm(103) = p_hc5 ; racm_per_megan(103) = 0.854
    p_of_megan2racm(104) = is_phenyl_cco ; p_of_racm(104) = non_react ; racm_per_megan(104) = 1.000
    p_of_megan2racm(105) = is_pyruvic_acid ; p_of_racm(105) = p_ora2 ; racm_per_megan(105) = 1.000
    p_of_megan2racm(106) = is_terpinyl_act_a ; p_of_racm(106) = p_api ; racm_per_megan(106) = 1.000
    p_of_megan2racm(107) = is_tetradecene_1 ; p_of_racm(107) = p_olt ; racm_per_megan(107) = 1.000
    p_of_megan2racm(108) = is_toluene ; p_of_racm(108) = p_tol ; racm_per_megan(108) = 1.000
    p_of_megan2racm(109) = is_carbon_monoxide ; p_of_racm(109) = p_co ; racm_per_megan(109) = 1.000
    p_of_megan2racm(110) = is_butene ; p_of_racm(110) = p_olt ; racm_per_megan(110) = 1.000
    p_of_megan2racm(111) = is_ethane ; p_of_racm(111) = p_eth ; racm_per_megan(111) = 1.000
    p_of_megan2racm(112) = is_ethene ; p_of_racm(112) = p_ete ; racm_per_megan(112) = 1.000
    p_of_megan2racm(113) = is_hydrogen_cyanide ; p_of_racm(113) = non_react ; racm_per_megan(113) = 1.000
    p_of_megan2racm(114) = is_propane ; p_of_racm(114) = p_hc3 ; racm_per_megan(114) = 0.570
    p_of_megan2racm(115) = is_propene ; p_of_racm(115) = p_olt ; racm_per_megan(115) = 1.000
    p_of_megan2racm(116) = is_carbon_2s ; p_of_racm(116) = non_react ; racm_per_megan(116) = 1.000
    p_of_megan2racm(117) = is_carbonyl_s ; p_of_racm(117) = non_react ; racm_per_megan(117) = 1.000
    p_of_megan2racm(118) = is_diallyl_2s ; p_of_racm(118) = p_oli ; racm_per_megan(118) = 1.000
    p_of_megan2racm(119) = is_diallyl_2s ; p_of_racm(119) = p_so2 ; racm_per_megan(119) = 2.000
    p_of_megan2racm(120) = is_2met_2s ; p_of_racm(120) = p_eth ; racm_per_megan(120) = 1.000
    p_of_megan2racm(121) = is_2met_2s ; p_of_racm(121) = p_so2 ; racm_per_megan(121) = 2.000
    p_of_megan2racm(122) = is_2met_s ; p_of_racm(122) = p_eth ; racm_per_megan(122) = 1.000
    p_of_megan2racm(123) = is_2met_s ; p_of_racm(123) = p_so2 ; racm_per_megan(123) = 1.000
    p_of_megan2racm(124) = is_met_chloride ; p_of_racm(124) = non_react ; racm_per_megan(124) = 1.000
    p_of_megan2racm(125) = is_met_bromide ; p_of_racm(125) = non_react ; racm_per_megan(125) = 1.000
    p_of_megan2racm(126) = is_met_iodide ; p_of_racm(126) = non_react ; racm_per_megan(126) = 1.000
    p_of_megan2racm(127) = is_hydrogen_s ; p_of_racm(127) = p_so2 ; racm_per_megan(127) = 1.000
    p_of_megan2racm(128) = is_met_mercaptan ; p_of_racm(128) = p_ch4 ; racm_per_megan(128) = 1.000
    p_of_megan2racm(129) = is_met_mercaptan ; p_of_racm(129) = p_so2 ; racm_per_megan(129) = 1.000
    p_of_megan2racm(130) = is_met_propenyl_2s ; p_of_racm(130) = p_oli ; racm_per_megan(130) = 1.000
    p_of_megan2racm(131) = is_met_propenyl_2s ; p_of_racm(131) = p_so2 ; racm_per_megan(131) = 2.000
    p_of_megan2racm(132) = is_pppp_2s ; p_of_racm(132) = p_oli ; racm_per_megan(132) = 1.000
    p_of_megan2racm(133) = is_pppp_2s ; p_of_racm(133) = p_so2 ; racm_per_megan(133) = 2.000
    p_of_megan2racm(134) = is_2met_nonatriene ; p_of_racm(134) = p_lim ; racm_per_megan(134) = 1.000
    p_of_megan2racm(135) = is_met_salicylate ; p_of_racm(135) = p_hc8 ; racm_per_megan(135) = 1.000
    p_of_megan2racm(136) = is_indole ; p_of_racm(136) = p_hc8 ; racm_per_megan(136) = 1.201
    p_of_megan2racm(137) = is_indole ; p_of_racm(137) = p_hno3 ; racm_per_megan(137) = 1.000
    p_of_megan2racm(138) = is_jasmone ; p_of_racm(138) = p_lim ; racm_per_megan(138) = 1.000
    p_of_megan2racm(139) = is_met_jasmonate ; p_of_racm(139) = p_lim ; racm_per_megan(139) = 1.000
    p_of_megan2racm(140) = is_3met_3dctt ; p_of_racm(140) = p_oli ; racm_per_megan(140) = 0.500
    p_of_megan2racm(141) = is_3met_3dctt ; p_of_racm(141) = p_olt ; racm_per_megan(141) = 0.500
    p_of_megan2racm(142) = is_hexanal ; p_of_racm(142) = p_ald ; racm_per_megan(142) = 1.000
    p_of_megan2racm(143) = is_hexanol_1 ; p_of_racm(143) = p_hc5 ; racm_per_megan(143) = 1.710
    p_of_megan2racm(144) = is_hexenal_c3 ; p_of_racm(144) = p_oli ; racm_per_megan(144) = 1.000
    p_of_megan2racm(145) = is_hexenal_t2 ; p_of_racm(145) = p_oli ; racm_per_megan(145) = 1.000
    p_of_megan2racm(146) = is_hexenol_c3 ; p_of_racm(146) = p_olt ; racm_per_megan(146) = 1.000
    p_of_megan2racm(147) = is_hexenyl_act_c3 ; p_of_racm(147) = p_oli ; racm_per_megan(147) = 1.000
  END SUBROUTINE get_megan2racm_table
  SUBROUTINE get_megan2racmSOA_table
    p_of_megan2racmSOA( 1) = is_isoprene ; p_of_racmSOA( 1) = p_iso ; racmSOA_per_megan( 1) = 1.000
    p_of_megan2racmSOA( 2) = is_myrcene ; p_of_racmSOA( 2) = p_lim ; racmSOA_per_megan( 2) = 1.000
    p_of_megan2racmSOA( 3) = is_sabinene ; p_of_racmSOA( 3) = p_api ; racmSOA_per_megan( 3) = 1.000
    p_of_megan2racmSOA( 4) = is_limonene ; p_of_racmSOA( 4) = p_lim ; racmSOA_per_megan( 4) = 1.000
    p_of_megan2racmSOA( 5) = is_carene_3 ; p_of_racmSOA( 5) = p_api ; racmSOA_per_megan( 5) = 1.000
    p_of_megan2racmSOA( 6) = is_ocimene_t_b ; p_of_racmSOA( 6) = p_lim ; racmSOA_per_megan( 6) = 1.000
    p_of_megan2racmSOA( 7) = is_pinene_b ; p_of_racmSOA( 7) = p_api ; racmSOA_per_megan( 7) = 1.000
    p_of_megan2racmSOA( 8) = is_pinene_a ; p_of_racmSOA( 8) = p_api ; racmSOA_per_megan( 8) = 1.000
    p_of_megan2racmSOA( 9) = is_2met_styrene ; p_of_racmSOA( 9) = p_lim ; racmSOA_per_megan( 9) = 1.000
    p_of_megan2racmSOA( 10) = is_cymene_p ; p_of_racmSOA( 10) = p_lim ; racmSOA_per_megan( 10) = 1.000
    p_of_megan2racmSOA( 11) = is_cymene_o ; p_of_racmSOA( 11) = p_lim ; racmSOA_per_megan( 11) = 1.000
    p_of_megan2racmSOA( 12) = is_phellandrene_a ; p_of_racmSOA( 12) = p_lim ; racmSOA_per_megan( 12) = 1.000
    p_of_megan2racmSOA( 13) = is_thujene_a ; p_of_racmSOA( 13) = p_api ; racmSOA_per_megan( 13) = 1.000
    p_of_megan2racmSOA( 14) = is_terpinene_a ; p_of_racmSOA( 14) = p_lim ; racmSOA_per_megan( 14) = 1.000
    p_of_megan2racmSOA( 15) = is_terpinene_g ; p_of_racmSOA( 15) = p_lim ; racmSOA_per_megan( 15) = 1.000
    p_of_megan2racmSOA( 16) = is_terpinolene ; p_of_racmSOA( 16) = p_lim ; racmSOA_per_megan( 16) = 1.000
    p_of_megan2racmSOA( 17) = is_phellandrene_b ; p_of_racmSOA( 17) = p_lim ; racmSOA_per_megan( 17) = 1.000
    p_of_megan2racmSOA( 18) = is_camphene ; p_of_racmSOA( 18) = p_api ; racmSOA_per_megan( 18) = 1.000
    p_of_megan2racmSOA( 19) = is_bornene ; p_of_racmSOA( 19) = p_api ; racmSOA_per_megan( 19) = 1.000
    p_of_megan2racmSOA( 20) = is_fenchene_a ; p_of_racmSOA( 20) = p_api ; racmSOA_per_megan( 20) = 1.000
    p_of_megan2racmSOA( 21) = is_ocimene_al ; p_of_racmSOA( 21) = p_lim ; racmSOA_per_megan( 21) = 1.000
    p_of_megan2racmSOA( 22) = is_ocimene_c_b ; p_of_racmSOA( 22) = p_lim ; racmSOA_per_megan( 22) = 1.000
    p_of_megan2racmSOA( 23) = is_tricyclene ; p_of_racmSOA( 23) = non_react ; racmSOA_per_megan( 23) = 1.000
    p_of_megan2racmSOA( 24) = is_estragole ; p_of_racmSOA( 24) = p_lim ; racmSOA_per_megan( 24) = 1.000
    p_of_megan2racmSOA( 25) = is_camphor ; p_of_racmSOA( 25) = p_hc8 ; racmSOA_per_megan( 25) = 0.380
    p_of_megan2racmSOA( 26) = is_fenchone ; p_of_racmSOA( 26) = non_react ; racmSOA_per_megan( 26) = 1.000
    p_of_megan2racmSOA( 27) = is_piperitone ; p_of_racmSOA( 27) = p_api ; racmSOA_per_megan( 27) = 1.000
    p_of_megan2racmSOA( 28) = is_thujone_a ; p_of_racmSOA( 28) = non_react ; racmSOA_per_megan( 28) = 1.000
    p_of_megan2racmSOA( 29) = is_thujone_b ; p_of_racmSOA( 29) = non_react ; racmSOA_per_megan( 29) = 1.000
    p_of_megan2racmSOA( 30) = is_cineole_1_8 ; p_of_racmSOA( 30) = p_hc8 ; racmSOA_per_megan( 30) = 0.738
    p_of_megan2racmSOA( 31) = is_borneol ; p_of_racmSOA( 31) = non_react ; racmSOA_per_megan( 31) = 1.000
    p_of_megan2racmSOA( 32) = is_linalool ; p_of_racmSOA( 32) = p_lim ; racmSOA_per_megan( 32) = 1.000
    p_of_megan2racmSOA( 33) = is_terpineol_4 ; p_of_racmSOA( 33) = p_api ; racmSOA_per_megan( 33) = 1.000
    p_of_megan2racmSOA( 34) = is_terpineol_a ; p_of_racmSOA( 34) = p_api ; racmSOA_per_megan( 34) = 1.000
    p_of_megan2racmSOA( 35) = is_linalool_oxd_c ; p_of_racmSOA( 35) = p_lim ; racmSOA_per_megan( 35) = 1.000
    p_of_megan2racmSOA( 36) = is_linalool_oxd_t ; p_of_racmSOA( 36) = p_lim ; racmSOA_per_megan( 36) = 1.000
    p_of_megan2racmSOA( 37) = is_ionone_b ; p_of_racmSOA( 37) = p_lim ; racmSOA_per_megan( 37) = 1.000
    p_of_megan2racmSOA( 38) = is_bornyl_act ; p_of_racmSOA( 38) = non_react ; racmSOA_per_megan( 38) = 1.000
    p_of_megan2racmSOA( 39) = is_farnescene_a ; p_of_racmSOA( 39) = p_lim ; racmSOA_per_megan( 39) = 1.000
    p_of_megan2racmSOA( 40) = is_caryophyllene_b ; p_of_racmSOA( 40) = p_lim ; racmSOA_per_megan( 40) = 1.000
    p_of_megan2racmSOA( 41) = is_acoradiene ; p_of_racmSOA( 41) = p_lim ; racmSOA_per_megan( 41) = 1.000
    p_of_megan2racmSOA( 42) = is_aromadendrene ; p_of_racmSOA( 42) = p_api ; racmSOA_per_megan( 42) = 1.000
    p_of_megan2racmSOA( 43) = is_bergamotene_a ; p_of_racmSOA( 43) = p_lim ; racmSOA_per_megan( 43) = 1.000
    p_of_megan2racmSOA( 44) = is_bergamotene_b ; p_of_racmSOA( 44) = p_lim ; racmSOA_per_megan( 44) = 1.000
    p_of_megan2racmSOA( 45) = is_bisabolene_a ; p_of_racmSOA( 45) = p_lim ; racmSOA_per_megan( 45) = 1.000
    p_of_megan2racmSOA( 46) = is_bisabolene_b ; p_of_racmSOA( 46) = p_lim ; racmSOA_per_megan( 46) = 1.000
    p_of_megan2racmSOA( 47) = is_bourbonene_b ; p_of_racmSOA( 47) = p_api ; racmSOA_per_megan( 47) = 1.000
    p_of_megan2racmSOA( 48) = is_cadinene_d ; p_of_racmSOA( 48) = p_lim ; racmSOA_per_megan( 48) = 1.000
    p_of_megan2racmSOA( 49) = is_cadinene_g ; p_of_racmSOA( 49) = p_lim ; racmSOA_per_megan( 49) = 1.000
    p_of_megan2racmSOA( 50) = is_cedrene_a ; p_of_racmSOA( 50) = p_api ; racmSOA_per_megan( 50) = 1.000
    p_of_megan2racmSOA( 51) = is_copaene_a ; p_of_racmSOA( 51) = p_api ; racmSOA_per_megan( 51) = 1.000
    p_of_megan2racmSOA( 52) = is_cubebene_a ; p_of_racmSOA( 52) = p_api ; racmSOA_per_megan( 52) = 1.000
    p_of_megan2racmSOA( 53) = is_cubebene_b ; p_of_racmSOA( 53) = p_api ; racmSOA_per_megan( 53) = 1.000
    p_of_megan2racmSOA( 54) = is_elemene_b ; p_of_racmSOA( 54) = p_lim ; racmSOA_per_megan( 54) = 1.000
    p_of_megan2racmSOA( 55) = is_farnescene_b ; p_of_racmSOA( 55) = p_lim ; racmSOA_per_megan( 55) = 1.000
    p_of_megan2racmSOA( 56) = is_germacrene_b ; p_of_racmSOA( 56) = p_lim ; racmSOA_per_megan( 56) = 1.000
    p_of_megan2racmSOA( 57) = is_germacrene_d ; p_of_racmSOA( 57) = p_lim ; racmSOA_per_megan( 57) = 1.000
    p_of_megan2racmSOA( 58) = is_gurjunene_b ; p_of_racmSOA( 58) = p_api ; racmSOA_per_megan( 58) = 1.000
    p_of_megan2racmSOA( 59) = is_humulene_a ; p_of_racmSOA( 59) = p_lim ; racmSOA_per_megan( 59) = 1.000
    p_of_megan2racmSOA( 60) = is_humulene_g ; p_of_racmSOA( 60) = p_lim ; racmSOA_per_megan( 60) = 1.000
    p_of_megan2racmSOA( 61) = is_isolongifolene ; p_of_racmSOA( 61) = p_api ; racmSOA_per_megan( 61) = 1.000
    p_of_megan2racmSOA( 62) = is_longifolene ; p_of_racmSOA( 62) = p_api ; racmSOA_per_megan( 62) = 1.000
    p_of_megan2racmSOA( 63) = is_longipinene ; p_of_racmSOA( 63) = p_api ; racmSOA_per_megan( 63) = 1.000
    p_of_megan2racmSOA( 64) = is_muurolene_a ; p_of_racmSOA( 64) = p_lim ; racmSOA_per_megan( 64) = 1.000
    p_of_megan2racmSOA( 65) = is_muurolene_g ; p_of_racmSOA( 65) = p_lim ; racmSOA_per_megan( 65) = 1.000
    p_of_megan2racmSOA( 66) = is_selinene_b ; p_of_racmSOA( 66) = p_lim ; racmSOA_per_megan( 66) = 1.000
    p_of_megan2racmSOA( 67) = is_selinene_d ; p_of_racmSOA( 67) = p_lim ; racmSOA_per_megan( 67) = 1.000
    p_of_megan2racmSOA( 68) = is_nerolidol_c ; p_of_racmSOA( 68) = p_lim ; racmSOA_per_megan( 68) = 1.000
    p_of_megan2racmSOA( 69) = is_nerolidol_t ; p_of_racmSOA( 69) = p_lim ; racmSOA_per_megan( 69) = 1.000
    p_of_megan2racmSOA( 70) = is_cedrol ; p_of_racmSOA( 70) = non_react ; racmSOA_per_megan( 70) = 1.000
    p_of_megan2racmSOA( 71) = is_mbo_2m3e2ol ; p_of_racmSOA( 71) = p_mbo ; racmSOA_per_megan( 71) = 1.000
    p_of_megan2racmSOA( 72) = is_methanol ; p_of_racmSOA( 72) = p_hc3 ; racmSOA_per_megan( 72) = 0.490
    p_of_megan2racmSOA( 73) = is_acetone ; p_of_racmSOA( 73) = p_ket ; racmSOA_per_megan( 73) = 0.330
    p_of_megan2racmSOA( 74) = is_methane ; p_of_racmSOA( 74) = p_ch4 ; racmSOA_per_megan( 74) = 1.000
    p_of_megan2racmSOA( 75) = is_ammonia ; p_of_racmSOA( 75) = p_nh3 ; racmSOA_per_megan( 75) = 1.000
    p_of_megan2racmSOA( 76) = is_nitrous_oxd ; p_of_racmSOA( 76) = p_no2 ; racmSOA_per_megan( 76) = 1.000
    p_of_megan2racmSOA( 77) = is_nitric_oxd ; p_of_racmSOA( 77) = p_no ; racmSOA_per_megan( 77) = 1.000
    p_of_megan2racmSOA( 78) = is_acetaldehyde ; p_of_racmSOA( 78) = p_ald ; racmSOA_per_megan( 78) = 1.000
    p_of_megan2racmSOA( 79) = is_ethanol ; p_of_racmSOA( 79) = p_hc3 ; racmSOA_per_megan( 79) = 1.370
    p_of_megan2racmSOA( 80) = is_formic_acid ; p_of_racmSOA( 80) = p_ora1 ; racmSOA_per_megan( 80) = 1.000
    p_of_megan2racmSOA( 81) = is_formaldehyde ; p_of_racmSOA( 81) = p_hcho ; racmSOA_per_megan( 81) = 1.000
    p_of_megan2racmSOA( 82) = is_acetic_acid ; p_of_racmSOA( 82) = p_ora2 ; racmSOA_per_megan( 82) = 1.000
    p_of_megan2racmSOA( 83) = is_mbo_3m2e1ol ; p_of_racmSOA( 83) = p_mbo ; racmSOA_per_megan( 83) = 1.000
    p_of_megan2racmSOA( 84) = is_mbo_3m3e1ol ; p_of_racmSOA( 84) = p_mbo ; racmSOA_per_megan( 84) = 1.000
    p_of_megan2racmSOA( 85) = is_benzaldehyde ; p_of_racmSOA( 85) = non_react ; racmSOA_per_megan( 85) = 1.000
    p_of_megan2racmSOA( 86) = is_butanone_2 ; p_of_racmSOA( 86) = p_ket ; racmSOA_per_megan( 86) = 1.610
    p_of_megan2racmSOA( 87) = is_decanal ; p_of_racmSOA( 87) = p_ald ; racmSOA_per_megan( 87) = 1.000
    p_of_megan2racmSOA( 88) = is_dodecene_1 ; p_of_racmSOA( 88) = p_olt ; racmSOA_per_megan( 88) = 1.000
    p_of_megan2racmSOA( 89) = is_geranyl_acetone ; p_of_racmSOA( 89) = p_lim ; racmSOA_per_megan( 89) = 1.000
    p_of_megan2racmSOA( 90) = is_heptanal ; p_of_racmSOA( 90) = p_ald ; racmSOA_per_megan( 90) = 1.000
    p_of_megan2racmSOA( 91) = is_heptane ; p_of_racmSOA( 91) = p_hc5 ; racmSOA_per_megan( 91) = 1.236
    p_of_megan2racmSOA( 92) = is_hexane ; p_of_racmSOA( 92) = p_hc5 ; racmSOA_per_megan( 92) = 1.058
    p_of_megan2racmSOA( 93) = is_met_benzoate ; p_of_racmSOA( 93) = p_hc8 ; racmSOA_per_megan( 93) = 1.000
    p_of_megan2racmSOA( 94) = is_met_heptenone ; p_of_racmSOA( 94) = p_oli ; racmSOA_per_megan( 94) = 1.000
    p_of_megan2racmSOA( 95) = is_neryl_acetone ; p_of_racmSOA( 95) = p_oli ; racmSOA_per_megan( 95) = 1.000
    p_of_megan2racmSOA( 96) = is_nonanal ; p_of_racmSOA( 96) = p_ald ; racmSOA_per_megan( 96) = 1.000
    p_of_megan2racmSOA( 97) = is_nonenal ; p_of_racmSOA( 97) = p_hc8 ; racmSOA_per_megan( 97) = 1.000
    p_of_megan2racmSOA( 98) = is_nonenal ; p_of_racmSOA( 98) = p_ald ; racmSOA_per_megan( 98) = 1.000
    p_of_megan2racmSOA( 99) = is_octanal ; p_of_racmSOA( 99) = p_ald ; racmSOA_per_megan( 99) = 1.000
    p_of_megan2racmSOA(100) = is_octanol ; p_of_racmSOA(100) = p_hc8 ; racmSOA_per_megan(100) = 1.092
    p_of_megan2racmSOA(101) = is_octenol_1e3ol ; p_of_racmSOA(101) = p_olt ; racmSOA_per_megan(101) = 1.000
    p_of_megan2racmSOA(102) = is_oxopentanal ; p_of_racmSOA(102) = p_ald ; racmSOA_per_megan(102) = 1.000
    p_of_megan2racmSOA(103) = is_pentane ; p_of_racmSOA(103) = p_hc5 ; racmSOA_per_megan(103) = 0.854
    p_of_megan2racmSOA(104) = is_phenyl_cco ; p_of_racmSOA(104) = non_react ; racmSOA_per_megan(104) = 1.000
    p_of_megan2racmSOA(105) = is_pyruvic_acid ; p_of_racmSOA(105) = p_ora2 ; racmSOA_per_megan(105) = 1.000
    p_of_megan2racmSOA(106) = is_terpinyl_act_a ; p_of_racmSOA(106) = p_api ; racmSOA_per_megan(106) = 1.000
    p_of_megan2racmSOA(107) = is_tetradecene_1 ; p_of_racmSOA(107) = p_olt ; racmSOA_per_megan(107) = 1.000
    p_of_megan2racmSOA(108) = is_toluene ; p_of_racmSOA(108) = p_tol ; racmSOA_per_megan(108) = 1.000
    p_of_megan2racmSOA(109) = is_carbon_monoxide ; p_of_racmSOA(109) = p_co ; racmSOA_per_megan(109) = 1.000
    p_of_megan2racmSOA(110) = is_butene ; p_of_racmSOA(110) = p_olt ; racmSOA_per_megan(110) = 1.000
    p_of_megan2racmSOA(111) = is_ethane ; p_of_racmSOA(111) = p_eth ; racmSOA_per_megan(111) = 1.000
    p_of_megan2racmSOA(112) = is_ethene ; p_of_racmSOA(112) = p_ete ; racmSOA_per_megan(112) = 1.000
    p_of_megan2racmSOA(113) = is_hydrogen_cyanide ; p_of_racmSOA(113) = non_react ; racmSOA_per_megan(113) = 1.000
    p_of_megan2racmSOA(114) = is_propane ; p_of_racmSOA(114) = p_hc3 ; racmSOA_per_megan(114) = 0.570
    p_of_megan2racmSOA(115) = is_propene ; p_of_racmSOA(115) = p_olt ; racmSOA_per_megan(115) = 1.000
    p_of_megan2racmSOA(116) = is_carbon_2s ; p_of_racmSOA(116) = non_react ; racmSOA_per_megan(116) = 1.000
    p_of_megan2racmSOA(117) = is_carbonyl_s ; p_of_racmSOA(117) = non_react ; racmSOA_per_megan(117) = 1.000
    p_of_megan2racmSOA(118) = is_diallyl_2s ; p_of_racmSOA(118) = p_oli ; racmSOA_per_megan(118) = 1.000
    p_of_megan2racmSOA(119) = is_diallyl_2s ; p_of_racmSOA(119) = p_so2 ; racmSOA_per_megan(119) = 2.000
    p_of_megan2racmSOA(120) = is_2met_2s ; p_of_racmSOA(120) = p_eth ; racmSOA_per_megan(120) = 1.000
    p_of_megan2racmSOA(121) = is_2met_2s ; p_of_racmSOA(121) = p_so2 ; racmSOA_per_megan(121) = 2.000
    p_of_megan2racmSOA(122) = is_2met_s ; p_of_racmSOA(122) = p_eth ; racmSOA_per_megan(122) = 1.000
    p_of_megan2racmSOA(123) = is_2met_s ; p_of_racmSOA(123) = p_so2 ; racmSOA_per_megan(123) = 1.000
    p_of_megan2racmSOA(124) = is_met_chloride ; p_of_racmSOA(124) = non_react ; racmSOA_per_megan(124) = 1.000
    p_of_megan2racmSOA(125) = is_met_bromide ; p_of_racmSOA(125) = non_react ; racmSOA_per_megan(125) = 1.000
    p_of_megan2racmSOA(126) = is_met_iodide ; p_of_racmSOA(126) = non_react ; racmSOA_per_megan(126) = 1.000
    p_of_megan2racmSOA(127) = is_hydrogen_s ; p_of_racmSOA(127) = p_so2 ; racmSOA_per_megan(127) = 1.000
    p_of_megan2racmSOA(128) = is_met_mercaptan ; p_of_racmSOA(128) = p_ch4 ; racmSOA_per_megan(128) = 1.000
    p_of_megan2racmSOA(129) = is_met_mercaptan ; p_of_racmSOA(129) = p_so2 ; racmSOA_per_megan(129) = 1.000
    p_of_megan2racmSOA(130) = is_met_propenyl_2s ; p_of_racmSOA(130) = p_oli ; racmSOA_per_megan(130) = 1.000
    p_of_megan2racmSOA(131) = is_met_propenyl_2s ; p_of_racmSOA(131) = p_so2 ; racmSOA_per_megan(131) = 2.000
    p_of_megan2racmSOA(132) = is_pppp_2s ; p_of_racmSOA(132) = p_oli ; racmSOA_per_megan(132) = 1.000
    p_of_megan2racmSOA(133) = is_pppp_2s ; p_of_racmSOA(133) = p_so2 ; racmSOA_per_megan(133) = 2.000
    p_of_megan2racmSOA(134) = is_2met_nonatriene ; p_of_racmSOA(134) = p_lim ; racmSOA_per_megan(134) = 1.000
    p_of_megan2racmSOA(135) = is_met_salicylate ; p_of_racmSOA(135) = p_hc8 ; racmSOA_per_megan(135) = 1.000
    p_of_megan2racmSOA(136) = is_indole ; p_of_racmSOA(136) = p_hc8 ; racmSOA_per_megan(136) = 1.201
    p_of_megan2racmSOA(137) = is_indole ; p_of_racmSOA(137) = p_hno3 ; racmSOA_per_megan(137) = 1.000
    p_of_megan2racmSOA(138) = is_jasmone ; p_of_racmSOA(138) = p_lim ; racmSOA_per_megan(138) = 1.000
    p_of_megan2racmSOA(139) = is_met_jasmonate ; p_of_racmSOA(139) = p_lim ; racmSOA_per_megan(139) = 1.000
    p_of_megan2racmSOA(140) = is_3met_3dctt ; p_of_racmSOA(140) = p_oli ; racmSOA_per_megan(140) = 0.500
    p_of_megan2racmSOA(141) = is_3met_3dctt ; p_of_racmSOA(141) = p_olt ; racmSOA_per_megan(141) = 0.500
    p_of_megan2racmSOA(142) = is_hexanal ; p_of_racmSOA(142) = p_ald ; racmSOA_per_megan(142) = 1.000
    p_of_megan2racmSOA(143) = is_hexanol_1 ; p_of_racmSOA(143) = p_hc5 ; racmSOA_per_megan(143) = 1.710
    p_of_megan2racmSOA(144) = is_hexenal_c3 ; p_of_racmSOA(144) = p_oli ; racmSOA_per_megan(144) = 1.000
    p_of_megan2racmSOA(145) = is_hexenal_t2 ; p_of_racmSOA(145) = p_oli ; racmSOA_per_megan(145) = 1.000
    p_of_megan2racmSOA(146) = is_hexenol_c3 ; p_of_racmSOA(146) = p_olt ; racmSOA_per_megan(146) = 1.000
    p_of_megan2racmSOA(147) = is_hexenyl_act_c3 ; p_of_racmSOA(147) = p_oli ; racmSOA_per_megan(147) = 1.000
    p_of_megan2racmSOA(148) = is_farnescene_a ; p_of_racmSOA(148) = p_sesq ; racmSOA_per_megan(148) = 1.000
    p_of_megan2racmSOA(149) = is_caryophyllene_b ; p_of_racmSOA(149) = p_sesq ; racmSOA_per_megan(149) = 1.000
    p_of_megan2racmSOA(150) = is_aromadendrene ; p_of_racmSOA(150) = p_sesq ; racmSOA_per_megan(150) = 1.000
    p_of_megan2racmSOA(151) = is_bergamotene_a ; p_of_racmSOA(151) = p_sesq ; racmSOA_per_megan(151) = 1.000
    p_of_megan2racmSOA(152) = is_bergamotene_b ; p_of_racmSOA(152) = p_sesq ; racmSOA_per_megan(152) = 1.000
    p_of_megan2racmSOA(153) = is_bisabolene_a ; p_of_racmSOA(153) = p_sesq ; racmSOA_per_megan(153) = 1.000
    p_of_megan2racmSOA(154) = is_bisabolene_b ; p_of_racmSOA(154) = p_sesq ; racmSOA_per_megan(154) = 1.000
    p_of_megan2racmSOA(155) = is_bourbonene_b ; p_of_racmSOA(155) = p_sesq ; racmSOA_per_megan(155) = 1.000
    p_of_megan2racmSOA(156) = is_cadinene_d ; p_of_racmSOA(156) = p_sesq ; racmSOA_per_megan(156) = 1.000
    p_of_megan2racmSOA(157) = is_cadinene_g ; p_of_racmSOA(157) = p_sesq ; racmSOA_per_megan(157) = 1.000
    p_of_megan2racmSOA(158) = is_cedrene_a ; p_of_racmSOA(158) = p_sesq ; racmSOA_per_megan(158) = 1.000
    p_of_megan2racmSOA(159) = is_copaene_a ; p_of_racmSOA(159) = p_sesq ; racmSOA_per_megan(159) = 1.000
    p_of_megan2racmSOA(160) = is_cubebene_a ; p_of_racmSOA(160) = p_sesq ; racmSOA_per_megan(160) = 1.000
    p_of_megan2racmSOA(161) = is_cubebene_b ; p_of_racmSOA(161) = p_sesq ; racmSOA_per_megan(161) = 1.000
    p_of_megan2racmSOA(162) = is_elemene_b ; p_of_racmSOA(162) = p_sesq ; racmSOA_per_megan(162) = 1.000
    p_of_megan2racmSOA(163) = is_farnescene_b ; p_of_racmSOA(163) = p_sesq ; racmSOA_per_megan(163) = 1.000
    p_of_megan2racmSOA(164) = is_germacrene_b ; p_of_racmSOA(164) = p_sesq ; racmSOA_per_megan(164) = 1.000
    p_of_megan2racmSOA(165) = is_germacrene_d ; p_of_racmSOA(165) = p_sesq ; racmSOA_per_megan(165) = 1.000
    p_of_megan2racmSOA(166) = is_gurjunene_b ; p_of_racmSOA(166) = p_sesq ; racmSOA_per_megan(166) = 1.000
    p_of_megan2racmSOA(167) = is_humulene_a ; p_of_racmSOA(167) = p_sesq ; racmSOA_per_megan(167) = 1.000
    p_of_megan2racmSOA(168) = is_humulene_g ; p_of_racmSOA(168) = p_sesq ; racmSOA_per_megan(168) = 1.000
    p_of_megan2racmSOA(169) = is_isolongifolene ; p_of_racmSOA(169) = p_sesq ; racmSOA_per_megan(169) = 1.000
    p_of_megan2racmSOA(170) = is_longifolene ; p_of_racmSOA(170) = p_sesq ; racmSOA_per_megan(170) = 1.000
    p_of_megan2racmSOA(171) = is_longipinene ; p_of_racmSOA(171) = p_sesq ; racmSOA_per_megan(171) = 1.000
    p_of_megan2racmSOA(172) = is_muurolene_a ; p_of_racmSOA(172) = p_sesq ; racmSOA_per_megan(172) = 1.000
    p_of_megan2racmSOA(173) = is_muurolene_g ; p_of_racmSOA(173) = p_sesq ; racmSOA_per_megan(173) = 1.000
    p_of_megan2racmSOA(174) = is_selinene_b ; p_of_racmSOA(174) = p_sesq ; racmSOA_per_megan(174) = 1.000
    p_of_megan2racmSOA(175) = is_selinene_d ; p_of_racmSOA(175) = p_sesq ; racmSOA_per_megan(175) = 1.000
    p_of_megan2racmSOA(176) = is_nerolidol_c ; p_of_racmSOA(176) = p_sesq ; racmSOA_per_megan(176) = 1.000
    p_of_megan2racmSOA(177) = is_nerolidol_t ; p_of_racmSOA(177) = p_sesq ; racmSOA_per_megan(177) = 1.000
    p_of_megan2racmSOA(178) = is_cedrol ; p_of_racmSOA(178) = p_sesq ; racmSOA_per_megan(178) = 1.000
  END SUBROUTINE get_megan2racmSOA_table
  SUBROUTINE get_megan2saprcnov_table
    p_of_megan2saprcnov( 1) = is_isoprene ; p_of_saprcnov( 1) = p_isoprene ; saprcnov_per_megan( 1) = 1.000
    p_of_megan2saprcnov( 2) = is_myrcene ; p_of_saprcnov( 2) = p_terp ; saprcnov_per_megan( 2) = 1.000
    p_of_megan2saprcnov( 3) = is_sabinene ; p_of_saprcnov( 3) = p_terp ; saprcnov_per_megan( 3) = 1.000
    p_of_megan2saprcnov( 4) = is_limonene ; p_of_saprcnov( 4) = p_terp ; saprcnov_per_megan( 4) = 1.000
    p_of_megan2saprcnov( 5) = is_carene_3 ; p_of_saprcnov( 5) = p_terp ; saprcnov_per_megan( 5) = 1.000
    p_of_megan2saprcnov( 6) = is_ocimene_t_b ; p_of_saprcnov( 6) = p_terp ; saprcnov_per_megan( 6) = 1.000
    p_of_megan2saprcnov( 7) = is_pinene_b ; p_of_saprcnov( 7) = p_terp ; saprcnov_per_megan( 7) = 1.000
    p_of_megan2saprcnov( 8) = is_pinene_a ; p_of_saprcnov( 8) = p_terp ; saprcnov_per_megan( 8) = 1.000
    p_of_megan2saprcnov( 9) = is_2met_styrene ; p_of_saprcnov( 9) = p_ole2 ; saprcnov_per_megan( 9) = 1.000
    p_of_megan2saprcnov( 10) = is_cymene_p ; p_of_saprcnov( 10) = p_aro2 ; saprcnov_per_megan( 10) = 1.000
    p_of_megan2saprcnov( 11) = is_cymene_o ; p_of_saprcnov( 11) = p_aro2 ; saprcnov_per_megan( 11) = 1.000
    p_of_megan2saprcnov( 12) = is_phellandrene_a ; p_of_saprcnov( 12) = p_terp ; saprcnov_per_megan( 12) = 1.000
    p_of_megan2saprcnov( 13) = is_thujene_a ; p_of_saprcnov( 13) = p_terp ; saprcnov_per_megan( 13) = 1.000
    p_of_megan2saprcnov( 14) = is_terpinene_a ; p_of_saprcnov( 14) = p_terp ; saprcnov_per_megan( 14) = 1.000
    p_of_megan2saprcnov( 15) = is_terpinene_g ; p_of_saprcnov( 15) = p_terp ; saprcnov_per_megan( 15) = 1.000
    p_of_megan2saprcnov( 16) = is_terpinolene ; p_of_saprcnov( 16) = p_terp ; saprcnov_per_megan( 16) = 1.000
    p_of_megan2saprcnov( 17) = is_phellandrene_b ; p_of_saprcnov( 17) = p_terp ; saprcnov_per_megan( 17) = 1.000
    p_of_megan2saprcnov( 18) = is_camphene ; p_of_saprcnov( 18) = p_terp ; saprcnov_per_megan( 18) = 1.000
    p_of_megan2saprcnov( 19) = is_bornene ; p_of_saprcnov( 19) = p_terp ; saprcnov_per_megan( 19) = 1.000
    p_of_megan2saprcnov( 20) = is_fenchene_a ; p_of_saprcnov( 20) = p_terp ; saprcnov_per_megan( 20) = 1.000
    p_of_megan2saprcnov( 21) = is_ocimene_al ; p_of_saprcnov( 21) = p_terp ; saprcnov_per_megan( 21) = 1.000
    p_of_megan2saprcnov( 22) = is_ocimene_c_b ; p_of_saprcnov( 22) = p_terp ; saprcnov_per_megan( 22) = 1.000
    p_of_megan2saprcnov( 23) = is_tricyclene ; p_of_saprcnov( 23) = p_alk5 ; saprcnov_per_megan( 23) = 1.000
    p_of_megan2saprcnov( 24) = is_estragole ; p_of_saprcnov( 24) = p_terp ; saprcnov_per_megan( 24) = 1.000
    p_of_megan2saprcnov( 25) = is_camphor ; p_of_saprcnov( 25) = p_terp ; saprcnov_per_megan( 25) = 1.000
    p_of_megan2saprcnov( 26) = is_fenchone ; p_of_saprcnov( 26) = p_alk5 ; saprcnov_per_megan( 26) = 1.000
    p_of_megan2saprcnov( 27) = is_piperitone ; p_of_saprcnov( 27) = p_terp ; saprcnov_per_megan( 27) = 1.000
    p_of_megan2saprcnov( 28) = is_thujone_a ; p_of_saprcnov( 28) = p_alk5 ; saprcnov_per_megan( 28) = 1.000
    p_of_megan2saprcnov( 29) = is_thujone_b ; p_of_saprcnov( 29) = p_alk5 ; saprcnov_per_megan( 29) = 1.000
    p_of_megan2saprcnov( 30) = is_cineole_1_8 ; p_of_saprcnov( 30) = p_alk5 ; saprcnov_per_megan( 30) = 1.000
    p_of_megan2saprcnov( 31) = is_borneol ; p_of_saprcnov( 31) = p_alk5 ; saprcnov_per_megan( 31) = 1.000
    p_of_megan2saprcnov( 32) = is_linalool ; p_of_saprcnov( 32) = p_terp ; saprcnov_per_megan( 32) = 1.000
    p_of_megan2saprcnov( 33) = is_terpineol_4 ; p_of_saprcnov( 33) = p_terp ; saprcnov_per_megan( 33) = 1.000
    p_of_megan2saprcnov( 34) = is_terpineol_a ; p_of_saprcnov( 34) = p_terp ; saprcnov_per_megan( 34) = 1.000
    p_of_megan2saprcnov( 35) = is_linalool_oxd_c ; p_of_saprcnov( 35) = p_terp ; saprcnov_per_megan( 35) = 1.000
    p_of_megan2saprcnov( 36) = is_linalool_oxd_t ; p_of_saprcnov( 36) = p_terp ; saprcnov_per_megan( 36) = 1.000
    p_of_megan2saprcnov( 37) = is_ionone_b ; p_of_saprcnov( 37) = p_terp ; saprcnov_per_megan( 37) = 1.000
    p_of_megan2saprcnov( 38) = is_bornyl_act ; p_of_saprcnov( 38) = p_alk5 ; saprcnov_per_megan( 38) = 1.000
    p_of_megan2saprcnov( 39) = is_farnescene_a ; p_of_saprcnov( 39) = p_sesq ; saprcnov_per_megan( 39) = 1.000
    p_of_megan2saprcnov( 40) = is_caryophyllene_b ; p_of_saprcnov( 40) = p_sesq ; saprcnov_per_megan( 40) = 1.000
    p_of_megan2saprcnov( 41) = is_acoradiene ; p_of_saprcnov( 41) = p_terp ; saprcnov_per_megan( 41) = 1.000
    p_of_megan2saprcnov( 42) = is_aromadendrene ; p_of_saprcnov( 42) = p_sesq ; saprcnov_per_megan( 42) = 1.000
    p_of_megan2saprcnov( 43) = is_bergamotene_a ; p_of_saprcnov( 43) = p_sesq ; saprcnov_per_megan( 43) = 1.000
    p_of_megan2saprcnov( 44) = is_bergamotene_b ; p_of_saprcnov( 44) = p_sesq ; saprcnov_per_megan( 44) = 1.000
    p_of_megan2saprcnov( 45) = is_bisabolene_a ; p_of_saprcnov( 45) = p_sesq ; saprcnov_per_megan( 45) = 1.000
    p_of_megan2saprcnov( 46) = is_bisabolene_b ; p_of_saprcnov( 46) = p_sesq ; saprcnov_per_megan( 46) = 1.000
    p_of_megan2saprcnov( 47) = is_bourbonene_b ; p_of_saprcnov( 47) = p_sesq ; saprcnov_per_megan( 47) = 1.000
    p_of_megan2saprcnov( 48) = is_cadinene_d ; p_of_saprcnov( 48) = p_sesq ; saprcnov_per_megan( 48) = 1.000
    p_of_megan2saprcnov( 49) = is_cadinene_g ; p_of_saprcnov( 49) = p_sesq ; saprcnov_per_megan( 49) = 1.000
    p_of_megan2saprcnov( 50) = is_cedrene_a ; p_of_saprcnov( 50) = p_sesq ; saprcnov_per_megan( 50) = 1.000
    p_of_megan2saprcnov( 51) = is_copaene_a ; p_of_saprcnov( 51) = p_sesq ; saprcnov_per_megan( 51) = 1.000
    p_of_megan2saprcnov( 52) = is_cubebene_a ; p_of_saprcnov( 52) = p_sesq ; saprcnov_per_megan( 52) = 1.000
    p_of_megan2saprcnov( 53) = is_cubebene_b ; p_of_saprcnov( 53) = p_sesq ; saprcnov_per_megan( 53) = 1.000
    p_of_megan2saprcnov( 54) = is_elemene_b ; p_of_saprcnov( 54) = p_sesq ; saprcnov_per_megan( 54) = 1.000
    p_of_megan2saprcnov( 55) = is_farnescene_b ; p_of_saprcnov( 55) = p_sesq ; saprcnov_per_megan( 55) = 1.000
    p_of_megan2saprcnov( 56) = is_germacrene_b ; p_of_saprcnov( 56) = p_sesq ; saprcnov_per_megan( 56) = 1.000
    p_of_megan2saprcnov( 57) = is_germacrene_d ; p_of_saprcnov( 57) = p_sesq ; saprcnov_per_megan( 57) = 1.000
    p_of_megan2saprcnov( 58) = is_gurjunene_b ; p_of_saprcnov( 58) = p_sesq ; saprcnov_per_megan( 58) = 1.000
    p_of_megan2saprcnov( 59) = is_humulene_a ; p_of_saprcnov( 59) = p_sesq ; saprcnov_per_megan( 59) = 1.000
    p_of_megan2saprcnov( 60) = is_humulene_g ; p_of_saprcnov( 60) = p_sesq ; saprcnov_per_megan( 60) = 1.000
    p_of_megan2saprcnov( 61) = is_isolongifolene ; p_of_saprcnov( 61) = p_sesq ; saprcnov_per_megan( 61) = 1.000
    p_of_megan2saprcnov( 62) = is_longifolene ; p_of_saprcnov( 62) = p_sesq ; saprcnov_per_megan( 62) = 1.000
    p_of_megan2saprcnov( 63) = is_longipinene ; p_of_saprcnov( 63) = p_sesq ; saprcnov_per_megan( 63) = 1.000
    p_of_megan2saprcnov( 64) = is_muurolene_a ; p_of_saprcnov( 64) = p_sesq ; saprcnov_per_megan( 64) = 1.000
    p_of_megan2saprcnov( 65) = is_muurolene_g ; p_of_saprcnov( 65) = p_sesq ; saprcnov_per_megan( 65) = 1.000
    p_of_megan2saprcnov( 66) = is_selinene_b ; p_of_saprcnov( 66) = p_sesq ; saprcnov_per_megan( 66) = 1.000
    p_of_megan2saprcnov( 67) = is_selinene_d ; p_of_saprcnov( 67) = p_sesq ; saprcnov_per_megan( 67) = 1.000
    p_of_megan2saprcnov( 68) = is_nerolidol_c ; p_of_saprcnov( 68) = p_sesq ; saprcnov_per_megan( 68) = 1.000
    p_of_megan2saprcnov( 69) = is_nerolidol_t ; p_of_saprcnov( 69) = p_sesq ; saprcnov_per_megan( 69) = 1.000
    p_of_megan2saprcnov( 70) = is_cedrol ; p_of_saprcnov( 70) = p_sesq ; saprcnov_per_megan( 70) = 1.000
    p_of_megan2saprcnov( 71) = is_mbo_2m3e2ol ; p_of_saprcnov( 71) = p_isoprene ; saprcnov_per_megan( 71) = 1.000
    p_of_megan2saprcnov( 72) = is_methanol ; p_of_saprcnov( 72) = p_meoh ; saprcnov_per_megan( 72) = 1.000
    p_of_megan2saprcnov( 73) = is_acetone ; p_of_saprcnov( 73) = p_acet ; saprcnov_per_megan( 73) = 1.000
    p_of_megan2saprcnov( 74) = is_methane ; p_of_saprcnov( 74) = p_ch4 ; saprcnov_per_megan( 74) = 1.000
    p_of_megan2saprcnov( 75) = is_ammonia ; p_of_saprcnov( 75) = non_react ; saprcnov_per_megan( 75) = 1.000
    p_of_megan2saprcnov( 76) = is_nitrous_oxd ; p_of_saprcnov( 76) = non_react ; saprcnov_per_megan( 76) = 1.000
    p_of_megan2saprcnov( 77) = is_nitric_oxd ; p_of_saprcnov( 77) = p_no ; saprcnov_per_megan( 77) = 1.000
    p_of_megan2saprcnov( 78) = is_acetaldehyde ; p_of_saprcnov( 78) = p_ccho ; saprcnov_per_megan( 78) = 1.000
    p_of_megan2saprcnov( 79) = is_ethanol ; p_of_saprcnov( 79) = p_alk3 ; saprcnov_per_megan( 79) = 1.000
    p_of_megan2saprcnov( 80) = is_formic_acid ; p_of_saprcnov( 80) = p_hcooh ; saprcnov_per_megan( 80) = 1.000
    p_of_megan2saprcnov( 81) = is_formaldehyde ; p_of_saprcnov( 81) = p_hcho ; saprcnov_per_megan( 81) = 1.000
    p_of_megan2saprcnov( 82) = is_acetic_acid ; p_of_saprcnov( 82) = p_cco_oh ; saprcnov_per_megan( 82) = 1.000
    p_of_megan2saprcnov( 83) = is_mbo_3m2e1ol ; p_of_saprcnov( 83) = p_isoprene ; saprcnov_per_megan( 83) = 1.000
    p_of_megan2saprcnov( 84) = is_mbo_3m3e1ol ; p_of_saprcnov( 84) = p_isoprene ; saprcnov_per_megan( 84) = 1.000
    p_of_megan2saprcnov( 85) = is_benzaldehyde ; p_of_saprcnov( 85) = p_bald ; saprcnov_per_megan( 85) = 1.000
    p_of_megan2saprcnov( 86) = is_butanone_2 ; p_of_saprcnov( 86) = p_mek ; saprcnov_per_megan( 86) = 1.000
    p_of_megan2saprcnov( 87) = is_decanal ; p_of_saprcnov( 87) = p_rcho ; saprcnov_per_megan( 87) = 1.000
    p_of_megan2saprcnov( 88) = is_dodecene_1 ; p_of_saprcnov( 88) = p_ole1 ; saprcnov_per_megan( 88) = 1.000
    p_of_megan2saprcnov( 89) = is_geranyl_acetone ; p_of_saprcnov( 89) = p_terp ; saprcnov_per_megan( 89) = 1.000
    p_of_megan2saprcnov( 90) = is_heptanal ; p_of_saprcnov( 90) = p_rcho ; saprcnov_per_megan( 90) = 1.000
    p_of_megan2saprcnov( 91) = is_heptane ; p_of_saprcnov( 91) = p_alk5 ; saprcnov_per_megan( 91) = 1.000
    p_of_megan2saprcnov( 92) = is_hexane ; p_of_saprcnov( 92) = p_alk4 ; saprcnov_per_megan( 92) = 1.000
    p_of_megan2saprcnov( 93) = is_met_benzoate ; p_of_saprcnov( 93) = p_aro1 ; saprcnov_per_megan( 93) = 1.000
    p_of_megan2saprcnov( 94) = is_met_heptenone ; p_of_saprcnov( 94) = p_ole2 ; saprcnov_per_megan( 94) = 1.000
    p_of_megan2saprcnov( 95) = is_neryl_acetone ; p_of_saprcnov( 95) = p_ole2 ; saprcnov_per_megan( 95) = 1.000
    p_of_megan2saprcnov( 96) = is_nonanal ; p_of_saprcnov( 96) = p_rcho ; saprcnov_per_megan( 96) = 1.000
    p_of_megan2saprcnov( 97) = is_nonenal ; p_of_saprcnov( 97) = p_ole1 ; saprcnov_per_megan( 97) = 1.000
    p_of_megan2saprcnov( 98) = is_octanal ; p_of_saprcnov( 98) = p_rcho ; saprcnov_per_megan( 98) = 1.000
    p_of_megan2saprcnov( 99) = is_octanol ; p_of_saprcnov( 99) = p_alk5 ; saprcnov_per_megan( 99) = 1.000
    p_of_megan2saprcnov(100) = is_octenol_1e3ol ; p_of_saprcnov(100) = p_ole1 ; saprcnov_per_megan(100) = 1.000
    p_of_megan2saprcnov(101) = is_oxopentanal ; p_of_saprcnov(101) = p_rcho ; saprcnov_per_megan(101) = 1.000
    p_of_megan2saprcnov(102) = is_pentane ; p_of_saprcnov(102) = p_alk4 ; saprcnov_per_megan(102) = 1.000
    p_of_megan2saprcnov(103) = is_phenyl_cco ; p_of_saprcnov(103) = p_aro1 ; saprcnov_per_megan(103) = 1.000
    p_of_megan2saprcnov(104) = is_pyruvic_acid ; p_of_saprcnov(104) = p_rco_oh ; saprcnov_per_megan(104) = 1.000
    p_of_megan2saprcnov(105) = is_terpinyl_act_a ; p_of_saprcnov(105) = p_terp ; saprcnov_per_megan(105) = 1.000
    p_of_megan2saprcnov(106) = is_tetradecene_1 ; p_of_saprcnov(106) = p_ole1 ; saprcnov_per_megan(106) = 1.000
    p_of_megan2saprcnov(107) = is_toluene ; p_of_saprcnov(107) = p_aro1 ; saprcnov_per_megan(107) = 1.000
    p_of_megan2saprcnov(108) = is_carbon_monoxide ; p_of_saprcnov(108) = p_co ; saprcnov_per_megan(108) = 1.000
    p_of_megan2saprcnov(109) = is_butene ; p_of_saprcnov(109) = p_ole1 ; saprcnov_per_megan(109) = 1.000
    p_of_megan2saprcnov(110) = is_ethane ; p_of_saprcnov(110) = p_c2h6 ; saprcnov_per_megan(110) = 1.000
    p_of_megan2saprcnov(111) = is_ethene ; p_of_saprcnov(111) = p_ethene ; saprcnov_per_megan(111) = 1.000
    p_of_megan2saprcnov(112) = is_hydrogen_cyanide ; p_of_saprcnov(112) = non_react ; saprcnov_per_megan(112) = 1.000
    p_of_megan2saprcnov(113) = is_propane ; p_of_saprcnov(113) = p_c3h8 ; saprcnov_per_megan(113) = 1.000
    p_of_megan2saprcnov(114) = is_propene ; p_of_saprcnov(114) = p_c3h6 ; saprcnov_per_megan(114) = 1.000
    p_of_megan2saprcnov(115) = is_carbon_2s ; p_of_saprcnov(115) = non_react ; saprcnov_per_megan(115) = 1.000
    p_of_megan2saprcnov(116) = is_carbonyl_s ; p_of_saprcnov(116) = non_react ; saprcnov_per_megan(116) = 1.000
    p_of_megan2saprcnov(117) = is_diallyl_2s ; p_of_saprcnov(117) = p_ole1 ; saprcnov_per_megan(117) = 1.000
    p_of_megan2saprcnov(118) = is_2met_2s ; p_of_saprcnov(118) = p_alk5 ; saprcnov_per_megan(118) = 1.000
    p_of_megan2saprcnov(119) = is_2met_s ; p_of_saprcnov(119) = p_alk4 ; saprcnov_per_megan(119) = 1.000
    p_of_megan2saprcnov(120) = is_met_chloride ; p_of_saprcnov(120) = non_react ; saprcnov_per_megan(120) = 1.000
    p_of_megan2saprcnov(121) = is_met_bromide ; p_of_saprcnov(121) = non_react ; saprcnov_per_megan(121) = 1.000
    p_of_megan2saprcnov(122) = is_met_iodide ; p_of_saprcnov(122) = non_react ; saprcnov_per_megan(122) = 1.000
    p_of_megan2saprcnov(123) = is_hydrogen_s ; p_of_saprcnov(123) = non_react ; saprcnov_per_megan(123) = 1.000
    p_of_megan2saprcnov(124) = is_met_mercaptan ; p_of_saprcnov(124) = p_alk5 ; saprcnov_per_megan(124) = 1.000
    p_of_megan2saprcnov(125) = is_met_propenyl_2s ; p_of_saprcnov(125) = p_ole1 ; saprcnov_per_megan(125) = 1.000
    p_of_megan2saprcnov(126) = is_pppp_2s ; p_of_saprcnov(126) = p_ole1 ; saprcnov_per_megan(126) = 1.000
    p_of_megan2saprcnov(127) = is_2met_nonatriene ; p_of_saprcnov(127) = p_terp ; saprcnov_per_megan(127) = 1.000
    p_of_megan2saprcnov(128) = is_met_salicylate ; p_of_saprcnov(128) = p_aro1 ; saprcnov_per_megan(128) = 1.000
    p_of_megan2saprcnov(129) = is_indole ; p_of_saprcnov(129) = p_aro2 ; saprcnov_per_megan(129) = 1.000
    p_of_megan2saprcnov(130) = is_jasmone ; p_of_saprcnov(130) = p_terp ; saprcnov_per_megan(130) = 1.000
    p_of_megan2saprcnov(131) = is_met_jasmonate ; p_of_saprcnov(131) = p_terp ; saprcnov_per_megan(131) = 1.000
    p_of_megan2saprcnov(132) = is_3met_3dctt ; p_of_saprcnov(132) = p_terp ; saprcnov_per_megan(132) = 1.000
    p_of_megan2saprcnov(133) = is_hexanal ; p_of_saprcnov(133) = p_rcho ; saprcnov_per_megan(133) = 1.000
    p_of_megan2saprcnov(134) = is_hexanol_1 ; p_of_saprcnov(134) = p_alk5 ; saprcnov_per_megan(134) = 1.000
    p_of_megan2saprcnov(135) = is_hexenal_c3 ; p_of_saprcnov(135) = p_ole2 ; saprcnov_per_megan(135) = 1.000
    p_of_megan2saprcnov(136) = is_hexenal_t2 ; p_of_saprcnov(136) = p_ole2 ; saprcnov_per_megan(136) = 1.000
    p_of_megan2saprcnov(137) = is_hexenol_c3 ; p_of_saprcnov(137) = p_ole2 ; saprcnov_per_megan(137) = 1.000
    p_of_megan2saprcnov(138) = is_hexenyl_act_c3 ; p_of_saprcnov(138) = p_ole2 ; saprcnov_per_megan(138) = 1.000
  END SUBROUTINE get_megan2saprcnov_table
  SUBROUTINE get_megan2crimech_table
    p_of_megan2crimech( 1) = is_isoprene ; p_of_crimech( 1) = p_c5h8 ; crimech_per_megan( 1) = 1.000
    p_of_megan2crimech( 2) = is_myrcene ; p_of_crimech( 2) = p_bpinene ; crimech_per_megan( 2) = 1.000
    p_of_megan2crimech( 3) = is_sabinene ; p_of_crimech( 3) = p_bpinene ; crimech_per_megan( 3) = 1.000
    p_of_megan2crimech( 4) = is_limonene ; p_of_crimech( 4) = p_apinene ; crimech_per_megan( 4) = 0.667
    p_of_megan2crimech( 5) = is_limonene ; p_of_crimech( 5) = p_bpinene ; crimech_per_megan( 5) = 0.333
    p_of_megan2crimech( 6) = is_carene_3 ; p_of_crimech( 6) = p_apinene ; crimech_per_megan( 6) = 1.000
    p_of_megan2crimech( 7) = is_ocimene_t_b ; p_of_crimech( 7) = p_bpinene ; crimech_per_megan( 7) = 1.000
    p_of_megan2crimech( 8) = is_pinene_b ; p_of_crimech( 8) = p_apinene ; crimech_per_megan( 8) = 1.000
    p_of_megan2crimech( 9) = is_pinene_a ; p_of_crimech( 9) = p_bpinene ; crimech_per_megan( 9) = 1.000
    p_of_megan2crimech( 10) = is_2met_styrene ; p_of_crimech( 10) = p_apinene ; crimech_per_megan( 10) = 1.000
    p_of_megan2crimech( 11) = is_2met_styrene ; p_of_crimech( 11) = p_bpinene ; crimech_per_megan( 11) = 1.000
    p_of_megan2crimech( 12) = is_cymene_p ; p_of_crimech( 12) = p_toluene ; crimech_per_megan( 12) = 1.500
    p_of_megan2crimech( 13) = is_cymene_o ; p_of_crimech( 13) = p_toluene ; crimech_per_megan( 13) = 1.500
    p_of_megan2crimech( 14) = is_phellandrene_a ; p_of_crimech( 14) = p_apinene ; crimech_per_megan( 14) = 1.000
    p_of_megan2crimech( 15) = is_thujene_a ; p_of_crimech( 15) = p_apinene ; crimech_per_megan( 15) = 1.000
    p_of_megan2crimech( 16) = is_terpinene_a ; p_of_crimech( 16) = p_apinene ; crimech_per_megan( 16) = 1.000
    p_of_megan2crimech( 17) = is_terpinene_g ; p_of_crimech( 17) = p_apinene ; crimech_per_megan( 17) = 1.000
    p_of_megan2crimech( 18) = is_terpinolene ; p_of_crimech( 18) = p_apinene ; crimech_per_megan( 18) = 0.667
    p_of_megan2crimech( 19) = is_terpinolene ; p_of_crimech( 19) = p_bpinene ; crimech_per_megan( 19) = 0.333
    p_of_megan2crimech( 20) = is_phellandrene_b ; p_of_crimech( 20) = p_bpinene ; crimech_per_megan( 20) = 1.000
    p_of_megan2crimech( 21) = is_camphene ; p_of_crimech( 21) = p_bpinene ; crimech_per_megan( 21) = 1.000
    p_of_megan2crimech( 22) = is_bornene ; p_of_crimech( 22) = p_apinene ; crimech_per_megan( 22) = 1.000
    p_of_megan2crimech( 23) = is_fenchene_a ; p_of_crimech( 23) = p_apinene ; crimech_per_megan( 23) = 1.000
    p_of_megan2crimech( 24) = is_ocimene_al ; p_of_crimech( 24) = p_apinene ; crimech_per_megan( 24) = 1.000
    p_of_megan2crimech( 25) = is_ocimene_c_b ; p_of_crimech( 25) = p_bpinene ; crimech_per_megan( 25) = 1.000
    p_of_megan2crimech( 26) = is_tricyclene ; p_of_crimech( 26) = p_apinene ; crimech_per_megan( 26) = 0.667
    p_of_megan2crimech( 27) = is_tricyclene ; p_of_crimech( 27) = p_bpinene ; crimech_per_megan( 27) = 0.333
    p_of_megan2crimech( 28) = is_estragole ; p_of_crimech( 28) = p_apinene ; crimech_per_megan( 28) = 0.667
    p_of_megan2crimech( 29) = is_estragole ; p_of_crimech( 29) = p_bpinene ; crimech_per_megan( 29) = 0.333
    p_of_megan2crimech( 30) = is_camphor ; p_of_crimech( 30) = p_c5h8 ; crimech_per_megan( 30) = 2.000
    p_of_megan2crimech( 31) = is_fenchone ; p_of_crimech( 31) = p_apinene ; crimech_per_megan( 31) = 0.667
    p_of_megan2crimech( 32) = is_fenchone ; p_of_crimech( 32) = p_bpinene ; crimech_per_megan( 32) = 0.333
    p_of_megan2crimech( 33) = is_piperitone ; p_of_crimech( 33) = p_apinene ; crimech_per_megan( 33) = 1.000
    p_of_megan2crimech( 34) = is_thujone_a ; p_of_crimech( 34) = p_apinene ; crimech_per_megan( 34) = 1.000
    p_of_megan2crimech( 35) = is_thujone_b ; p_of_crimech( 35) = p_bpinene ; crimech_per_megan( 35) = 1.000
    p_of_megan2crimech( 36) = is_cineole_1_8 ; p_of_crimech( 36) = p_c5h8 ; crimech_per_megan( 36) = 2.000
    p_of_megan2crimech( 37) = is_borneol ; p_of_crimech( 37) = p_c5h8 ; crimech_per_megan( 37) = 2.000
    p_of_megan2crimech( 38) = is_linalool ; p_of_crimech( 38) = p_apinene ; crimech_per_megan( 38) = 0.667
    p_of_megan2crimech( 39) = is_linalool ; p_of_crimech( 39) = p_bpinene ; crimech_per_megan( 39) = 0.333
    p_of_megan2crimech( 40) = is_terpineol_4 ; p_of_crimech( 40) = p_apinene ; crimech_per_megan( 40) = 1.000
    p_of_megan2crimech( 41) = is_terpineol_a ; p_of_crimech( 41) = p_apinene ; crimech_per_megan( 41) = 1.000
    p_of_megan2crimech( 42) = is_linalool_oxd_c ; p_of_crimech( 42) = p_apinene ; crimech_per_megan( 42) = 0.833
    p_of_megan2crimech( 43) = is_linalool_oxd_c ; p_of_crimech( 43) = p_bpinene ; crimech_per_megan( 43) = 0.416
    p_of_megan2crimech( 44) = is_linalool_oxd_t ; p_of_crimech( 44) = p_apinene ; crimech_per_megan( 44) = 0.833
    p_of_megan2crimech( 45) = is_linalool_oxd_t ; p_of_crimech( 45) = p_bpinene ; crimech_per_megan( 45) = 0.416
    p_of_megan2crimech( 46) = is_ionone_b ; p_of_crimech( 46) = p_apinene ; crimech_per_megan( 46) = 0.934
    p_of_megan2crimech( 47) = is_ionone_b ; p_of_crimech( 47) = p_bpinene ; crimech_per_megan( 47) = 0.462
    p_of_megan2crimech( 48) = is_bornyl_act ; p_of_crimech( 48) = p_apinene ; crimech_per_megan( 48) = 0.667
    p_of_megan2crimech( 49) = is_bornyl_act ; p_of_crimech( 49) = p_bpinene ; crimech_per_megan( 49) = 0.333
    p_of_megan2crimech( 50) = is_farnescene_a ; p_of_crimech( 50) = p_apinene ; crimech_per_megan( 50) = 1.500
    p_of_megan2crimech( 51) = is_caryophyllene_b ; p_of_crimech( 51) = p_apinene ; crimech_per_megan( 51) = 1.000
    p_of_megan2crimech( 52) = is_caryophyllene_b ; p_of_crimech( 52) = p_bpinene ; crimech_per_megan( 52) = 0.500
    p_of_megan2crimech( 53) = is_acoradiene ; p_of_crimech( 53) = p_apinene ; crimech_per_megan( 53) = 1.000
    p_of_megan2crimech( 54) = is_acoradiene ; p_of_crimech( 54) = p_bpinene ; crimech_per_megan( 54) = 0.500
    p_of_megan2crimech( 55) = is_aromadendrene ; p_of_crimech( 55) = p_apinene ; crimech_per_megan( 55) = 1.000
    p_of_megan2crimech( 56) = is_aromadendrene ; p_of_crimech( 56) = p_bpinene ; crimech_per_megan( 56) = 0.500
    p_of_megan2crimech( 57) = is_bergamotene_a ; p_of_crimech( 57) = p_apinene ; crimech_per_megan( 57) = 1.500
    p_of_megan2crimech( 58) = is_bergamotene_b ; p_of_crimech( 58) = p_bpinene ; crimech_per_megan( 58) = 1.500
    p_of_megan2crimech( 59) = is_bisabolene_a ; p_of_crimech( 59) = p_apinene ; crimech_per_megan( 59) = 1.500
    p_of_megan2crimech( 60) = is_bisabolene_b ; p_of_crimech( 60) = p_bpinene ; crimech_per_megan( 60) = 1.500
    p_of_megan2crimech( 61) = is_bourbonene_b ; p_of_crimech( 61) = p_apinene ; crimech_per_megan( 61) = 1.000
    p_of_megan2crimech( 62) = is_bourbonene_b ; p_of_crimech( 62) = p_bpinene ; crimech_per_megan( 62) = 0.500
    p_of_megan2crimech( 63) = is_cadinene_d ; p_of_crimech( 63) = p_apinene ; crimech_per_megan( 63) = 1.000
    p_of_megan2crimech( 64) = is_cadinene_d ; p_of_crimech( 64) = p_bpinene ; crimech_per_megan( 64) = 0.500
    p_of_megan2crimech( 65) = is_cadinene_g ; p_of_crimech( 65) = p_apinene ; crimech_per_megan( 65) = 1.000
    p_of_megan2crimech( 66) = is_cadinene_g ; p_of_crimech( 66) = p_bpinene ; crimech_per_megan( 66) = 1.500
    p_of_megan2crimech( 67) = is_cedrene_a ; p_of_crimech( 67) = p_apinene ; crimech_per_megan( 67) = 1.000
    p_of_megan2crimech( 68) = is_cedrene_a ; p_of_crimech( 68) = p_bpinene ; crimech_per_megan( 68) = 0.500
    p_of_megan2crimech( 69) = is_copaene_a ; p_of_crimech( 69) = p_apinene ; crimech_per_megan( 69) = 1.000
    p_of_megan2crimech( 70) = is_copaene_a ; p_of_crimech( 70) = p_bpinene ; crimech_per_megan( 70) = 0.500
    p_of_megan2crimech( 71) = is_cubebene_a ; p_of_crimech( 71) = p_apinene ; crimech_per_megan( 71) = 1.000
    p_of_megan2crimech( 72) = is_cubebene_a ; p_of_crimech( 72) = p_bpinene ; crimech_per_megan( 72) = 0.500
    p_of_megan2crimech( 73) = is_cubebene_b ; p_of_crimech( 73) = p_apinene ; crimech_per_megan( 73) = 1.000
    p_of_megan2crimech( 74) = is_cubebene_b ; p_of_crimech( 74) = p_bpinene ; crimech_per_megan( 74) = 0.500
    p_of_megan2crimech( 75) = is_elemene_b ; p_of_crimech( 75) = p_apinene ; crimech_per_megan( 75) = 1.000
    p_of_megan2crimech( 76) = is_elemene_b ; p_of_crimech( 76) = p_bpinene ; crimech_per_megan( 76) = 0.500
    p_of_megan2crimech( 77) = is_farnescene_b ; p_of_crimech( 77) = p_apinene ; crimech_per_megan( 77) = 1.000
    p_of_megan2crimech( 78) = is_farnescene_b ; p_of_crimech( 78) = p_bpinene ; crimech_per_megan( 78) = 0.500
    p_of_megan2crimech( 79) = is_germacrene_b ; p_of_crimech( 79) = p_apinene ; crimech_per_megan( 79) = 1.000
    p_of_megan2crimech( 80) = is_germacrene_b ; p_of_crimech( 80) = p_bpinene ; crimech_per_megan( 80) = 0.500
    p_of_megan2crimech( 81) = is_germacrene_d ; p_of_crimech( 81) = p_apinene ; crimech_per_megan( 81) = 1.000
    p_of_megan2crimech( 82) = is_germacrene_d ; p_of_crimech( 82) = p_bpinene ; crimech_per_megan( 82) = 0.500
    p_of_megan2crimech( 83) = is_gurjunene_b ; p_of_crimech( 83) = p_apinene ; crimech_per_megan( 83) = 1.000
    p_of_megan2crimech( 84) = is_gurjunene_b ; p_of_crimech( 84) = p_bpinene ; crimech_per_megan( 84) = 0.500
    p_of_megan2crimech( 85) = is_humulene_a ; p_of_crimech( 85) = p_apinene ; crimech_per_megan( 85) = 1.000
    p_of_megan2crimech( 86) = is_humulene_a ; p_of_crimech( 86) = p_bpinene ; crimech_per_megan( 86) = 0.500
    p_of_megan2crimech( 87) = is_humulene_g ; p_of_crimech( 87) = p_apinene ; crimech_per_megan( 87) = 1.000
    p_of_megan2crimech( 88) = is_humulene_g ; p_of_crimech( 88) = p_bpinene ; crimech_per_megan( 88) = 0.500
    p_of_megan2crimech( 89) = is_isolongifolene ; p_of_crimech( 89) = p_apinene ; crimech_per_megan( 89) = 1.000
    p_of_megan2crimech( 90) = is_isolongifolene ; p_of_crimech( 90) = p_bpinene ; crimech_per_megan( 90) = 0.500
    p_of_megan2crimech( 91) = is_longifolene ; p_of_crimech( 91) = p_apinene ; crimech_per_megan( 91) = 1.000
    p_of_megan2crimech( 92) = is_longifolene ; p_of_crimech( 92) = p_bpinene ; crimech_per_megan( 92) = 0.500
    p_of_megan2crimech( 93) = is_longipinene ; p_of_crimech( 93) = p_apinene ; crimech_per_megan( 93) = 1.000
    p_of_megan2crimech( 94) = is_longipinene ; p_of_crimech( 94) = p_bpinene ; crimech_per_megan( 94) = 0.500
    p_of_megan2crimech( 95) = is_muurolene_a ; p_of_crimech( 95) = p_apinene ; crimech_per_megan( 95) = 1.000
    p_of_megan2crimech( 96) = is_muurolene_a ; p_of_crimech( 96) = p_bpinene ; crimech_per_megan( 96) = 0.500
    p_of_megan2crimech( 97) = is_muurolene_g ; p_of_crimech( 97) = p_apinene ; crimech_per_megan( 97) = 1.000
    p_of_megan2crimech( 98) = is_muurolene_g ; p_of_crimech( 98) = p_bpinene ; crimech_per_megan( 98) = 0.500
    p_of_megan2crimech( 99) = is_selinene_b ; p_of_crimech( 99) = p_apinene ; crimech_per_megan( 99) = 1.000
    p_of_megan2crimech(100) = is_selinene_b ; p_of_crimech(100) = p_bpinene ; crimech_per_megan(100) = 0.500
    p_of_megan2crimech(101) = is_selinene_d ; p_of_crimech(101) = p_apinene ; crimech_per_megan(101) = 1.000
    p_of_megan2crimech(102) = is_selinene_d ; p_of_crimech(102) = p_bpinene ; crimech_per_megan(102) = 0.500
    p_of_megan2crimech(103) = is_nerolidol_c ; p_of_crimech(103) = p_apinene ; crimech_per_megan(103) = 1.000
    p_of_megan2crimech(104) = is_nerolidol_c ; p_of_crimech(104) = p_bpinene ; crimech_per_megan(104) = 0.500
    p_of_megan2crimech(105) = is_nerolidol_t ; p_of_crimech(105) = p_apinene ; crimech_per_megan(105) = 1.000
    p_of_megan2crimech(106) = is_nerolidol_t ; p_of_crimech(106) = p_bpinene ; crimech_per_megan(106) = 0.500
    p_of_megan2crimech(107) = is_cedrol ; p_of_crimech(107) = p_c5h8 ; crimech_per_megan(107) = 3.000
    p_of_megan2crimech(108) = is_mbo_2m3e2ol ; p_of_crimech(108) = p_c5h8 ; crimech_per_megan(108) = 2.400
    p_of_megan2crimech(109) = is_methanol ; p_of_crimech(109) = p_ch3oh ; crimech_per_megan(109) = 1.000
    p_of_megan2crimech(110) = is_acetone ; p_of_crimech(110) = p_ket ; crimech_per_megan(110) = 1.000
    p_of_megan2crimech(111) = is_methane ; p_of_crimech(111) = p_ch4 ; crimech_per_megan(111) = 1.000
    p_of_megan2crimech(112) = is_ammonia ; p_of_crimech(112) = p_nh3 ; crimech_per_megan(112) = 1.000
    p_of_megan2crimech(113) = is_nitrous_oxd ; p_of_crimech(113) = p_no2 ; crimech_per_megan(113) = 1.000
    p_of_megan2crimech(114) = is_nitric_oxd ; p_of_crimech(114) = p_no ; crimech_per_megan(114) = 1.000
    p_of_megan2crimech(115) = is_acetaldehyde ; p_of_crimech(115) = p_ch3cho ; crimech_per_megan(115) = 1.000
    p_of_megan2crimech(116) = is_ethanol ; p_of_crimech(116) = p_c2h5oh ; crimech_per_megan(116) = 1.000
    p_of_megan2crimech(117) = is_formic_acid ; p_of_crimech(117) = p_hcooh ; crimech_per_megan(117) = 1.000
    p_of_megan2crimech(118) = is_formaldehyde ; p_of_crimech(118) = p_hcho ; crimech_per_megan(118) = 1.000
    p_of_megan2crimech(119) = is_acetic_acid ; p_of_crimech(119) = p_ch3co2h ; crimech_per_megan(119) = 1.000
    p_of_megan2crimech(120) = is_mbo_3m2e1ol ; p_of_crimech(120) = p_c5h8 ; crimech_per_megan(120) = 1.260
    p_of_megan2crimech(121) = is_mbo_3m3e1ol ; p_of_crimech(121) = p_c5h8 ; crimech_per_megan(121) = 1.260
    p_of_megan2crimech(122) = is_benzaldehyde ; p_of_crimech(122) = p_benzene ; crimech_per_megan(122) = 0.883
    p_of_megan2crimech(123) = is_butanone_2 ; p_of_crimech(123) = p_mek ; crimech_per_megan(123) = 1.000
    p_of_megan2crimech(124) = is_decanal ; p_of_crimech(124) = p_c2h5cho ; crimech_per_megan(124) = 2.690
    p_of_megan2crimech(125) = is_dodecene_1 ; p_of_crimech(125) = p_tbut2ene; crimech_per_megan(125) = 3.000
    p_of_megan2crimech(126) = is_geranyl_acetone ; p_of_crimech(126) = p_apinene ; crimech_per_megan(126) = 1.400
    p_of_megan2crimech(127) = is_geranyl_acetone ; p_of_crimech(127) = p_bpinene ; crimech_per_megan(127) = 1.400
    p_of_megan2crimech(128) = is_heptanal ; p_of_crimech(128) = p_c2h5cho ; crimech_per_megan(128) = 1.966
    p_of_megan2crimech(129) = is_heptane ; p_of_crimech(129) = p_nc4h10 ; crimech_per_megan(129) = 1.724
    p_of_megan2crimech(130) = is_hexane ; p_of_crimech(130) = p_nc4h10 ; crimech_per_megan(130) = 1.483
    p_of_megan2crimech(131) = is_met_benzoate ; p_of_crimech(131) = p_toluene ; crimech_per_megan(131) = 1.478
    p_of_megan2crimech(132) = is_met_heptenone ; p_of_crimech(132) = p_tbut2ene; crimech_per_megan(132) = 3.000
    p_of_megan2crimech(133) = is_neryl_acetone ; p_of_crimech(133) = p_ket ; crimech_per_megan(133) = 3.379
    p_of_megan2crimech(134) = is_nonanal ; p_of_crimech(134) = p_c2h5cho ; crimech_per_megan(134) = 2.448
    p_of_megan2crimech(135) = is_nonenal ; p_of_crimech(135) = p_nc4h10 ; crimech_per_megan(135) = 2.413
    p_of_megan2crimech(136) = is_nonenal ; p_of_crimech(136) = p_c2h5cho ; crimech_per_megan(136) = 2.413
    p_of_megan2crimech(137) = is_octanal ; p_of_crimech(137) = p_c2h5cho ; crimech_per_megan(137) = 2.207
    p_of_megan2crimech(138) = is_octanol ; p_of_crimech(138) = p_hc8 ; crimech_per_megan(138) = 2.207
    p_of_megan2crimech(139) = is_octenol_1e3ol ; p_of_crimech(139) = p_tbut2ene; crimech_per_megan(139) = 2.286
    p_of_megan2crimech(140) = is_oxopentanal ; p_of_crimech(140) = p_mek ; crimech_per_megan(140) = 1.400
    p_of_megan2crimech(141) = is_pentane ; p_of_crimech(141) = p_nc4h10 ; crimech_per_megan(141) = 1.241
    p_of_megan2crimech(142) = is_phenyl_cco ; p_of_crimech(142) = p_toluene ; crimech_per_megan(142) = 1.300
    p_of_megan2crimech(143) = is_pyruvic_acid ; p_of_crimech(143) = p_ch3co2h ; crimech_per_megan(143) = 1.158
    p_of_megan2crimech(144) = is_terpinyl_act_a ; p_of_crimech(144) = p_apinene ; crimech_per_megan(144) = 0.934
    p_of_megan2crimech(145) = is_terpinyl_act_a ; p_of_crimech(145) = p_bpinene ; crimech_per_megan(145) = 0.466
    p_of_megan2crimech(146) = is_tetradecene_1 ; p_of_crimech(146) = p_tbut2ene; crimech_per_megan(146) = 3.500
    p_of_megan2crimech(147) = is_toluene ; p_of_crimech(147) = p_toluene ; crimech_per_megan(147) = 1.000
    p_of_megan2crimech(148) = is_carbon_monoxide ; p_of_crimech(148) = p_co ; crimech_per_megan(148) = 1.000
    p_of_megan2crimech(149) = is_butene ; p_of_crimech(149) = p_tbut2ene; crimech_per_megan(149) = 1.000
    p_of_megan2crimech(150) = is_ethane ; p_of_crimech(150) = p_c2h6 ; crimech_per_megan(150) = 1.000
    p_of_megan2crimech(151) = is_ethene ; p_of_crimech(151) = p_c2h4 ; crimech_per_megan(151) = 1.000
    p_of_megan2crimech(152) = is_hydrogen_cyanide ; p_of_crimech(152) = non_react ; crimech_per_megan(152) = 1.000
    p_of_megan2crimech(153) = is_propane ; p_of_crimech(153) = p_c3h8 ; crimech_per_megan(153) = 1.000
    p_of_megan2crimech(154) = is_propene ; p_of_crimech(154) = p_c3h6 ; crimech_per_megan(154) = 1.000
    p_of_megan2crimech(155) = is_carbon_2s ; p_of_crimech(155) = non_react ; crimech_per_megan(155) = 1.000
    p_of_megan2crimech(156) = is_carbonyl_s ; p_of_crimech(156) = non_react ; crimech_per_megan(156) = 1.000
    p_of_megan2crimech(157) = is_diallyl_2s ; p_of_crimech(157) = p_tbut2ene; crimech_per_megan(157) = 0.660
    p_of_megan2crimech(158) = is_diallyl_2s ; p_of_crimech(158) = p_so2 ; crimech_per_megan(158) = 2.000
    p_of_megan2crimech(159) = is_2met_2s ; p_of_crimech(159) = p_c2h6 ; crimech_per_megan(159) = 1.000
    p_of_megan2crimech(160) = is_2met_2s ; p_of_crimech(160) = p_so2 ; crimech_per_megan(160) = 2.000
    p_of_megan2crimech(161) = is_2met_s ; p_of_crimech(161) = p_c2h6 ; crimech_per_megan(161) = 1.000
    p_of_megan2crimech(162) = is_2met_s ; p_of_crimech(162) = p_so2 ; crimech_per_megan(162) = 1.000
    p_of_megan2crimech(163) = is_met_chloride ; p_of_crimech(163) = non_react ; crimech_per_megan(163) = 1.000
    p_of_megan2crimech(164) = is_met_bromide ; p_of_crimech(164) = non_react ; crimech_per_megan(164) = 1.000
    p_of_megan2crimech(165) = is_met_iodide ; p_of_crimech(165) = non_react ; crimech_per_megan(165) = 1.000
    p_of_megan2crimech(166) = is_hydrogen_s ; p_of_crimech(166) = p_so2 ; crimech_per_megan(166) = 1.000
    p_of_megan2crimech(167) = is_met_mercaptan ; p_of_crimech(167) = p_ch4 ; crimech_per_megan(167) = 1.000
    p_of_megan2crimech(168) = is_met_mercaptan ; p_of_crimech(168) = p_so2 ; crimech_per_megan(168) = 1.000
    p_of_megan2crimech(169) = is_met_propenyl_2s ; p_of_crimech(169) = p_c3h6 ; crimech_per_megan(169) = 2.800
    p_of_megan2crimech(170) = is_met_propenyl_2s ; p_of_crimech(170) = p_so2 ; crimech_per_megan(170) = 2.000
    p_of_megan2crimech(171) = is_pppp_2s ; p_of_crimech(171) = p_c3h6 ; crimech_per_megan(171) = 3.500
    p_of_megan2crimech(172) = is_pppp_2s ; p_of_crimech(172) = p_so2 ; crimech_per_megan(172) = 1.800
    p_of_megan2crimech(173) = is_2met_nonatriene ; p_of_crimech(173) = p_apinene ; crimech_per_megan(173) = 0.667
    p_of_megan2crimech(174) = is_2met_nonatriene ; p_of_crimech(174) = p_bpinene ; crimech_per_megan(174) = 0.333
    p_of_megan2crimech(175) = is_met_salicylate ; p_of_crimech(175) = p_toluene ; crimech_per_megan(175) = 1.652
    p_of_megan2crimech(176) = is_indole ; p_of_crimech(176) = p_toluene ; crimech_per_megan(176) = 1.200
    p_of_megan2crimech(177) = is_indole ; p_of_crimech(177) = p_hno3 ; crimech_per_megan(177) = 1.000
    p_of_megan2crimech(178) = is_jasmone ; p_of_crimech(178) = p_apinene ; crimech_per_megan(178) = 0.804
    p_of_megan2crimech(179) = is_jasmone ; p_of_crimech(179) = p_bpinene ; crimech_per_megan(179) = 0.402
    p_of_megan2crimech(180) = is_met_jasmonate ; p_of_crimech(180) = p_apinene ; crimech_per_megan(180) = 1.098
    p_of_megan2crimech(181) = is_met_jasmonate ; p_of_crimech(181) = p_bpinene ; crimech_per_megan(181) = 0.548
    p_of_megan2crimech(182) = is_3met_3dctt ; p_of_crimech(182) = p_tbut2ene; crimech_per_megan(182) = 3.893
    p_of_megan2crimech(183) = is_hexanal ; p_of_crimech(183) = p_c2h5cho ; crimech_per_megan(183) = 1.720
    p_of_megan2crimech(184) = is_hexanol_1 ; p_of_crimech(184) = p_nc4h10 ; crimech_per_megan(184) = 1.759
    p_of_megan2crimech(185) = is_hexenal_c3 ; p_of_crimech(185) = p_tbut2ene; crimech_per_megan(185) = 1.750
    p_of_megan2crimech(186) = is_hexenal_t2 ; p_of_crimech(186) = p_tbut2ene; crimech_per_megan(186) = 1.750
    p_of_megan2crimech(187) = is_hexenol_c3 ; p_of_crimech(187) = p_tbut2ene; crimech_per_megan(187) = 1.786
    p_of_megan2crimech(188) = is_hexenyl_act_c3 ; p_of_crimech(188) = p_tbut2ene; crimech_per_megan(188) = 2.536
  END SUBROUTINE get_megan2crimech_table
END MODULE module_data_mgn2mech
