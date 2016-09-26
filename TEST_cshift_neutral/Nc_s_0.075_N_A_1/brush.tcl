####################################################
#                                                  #
# (BRUSH) POLYMER SIMULATION                       #
# BY JACOB                                         #
# 26.9.2016                                        #
#                                                  #
####################################################
#TODO N LJ interakce polymeru
#TODO
####################################################
set cas   [clock seconds]
#################### SETUP #########################

set steps 1000

### PARAMETERS OF STRUCTURE
set N_P 1	;#number of polymers
set N_B 31 	;#monomer units of main chain
set N_A 1 	;#monomer units of side chain
set h	2  	;#space between brush chains
set N_M [expr $N_P*($N_B+($N_B+1)/$h*$N_A)]
set Mr          72.     ;#molar mass of monomer unit FOR BOX SIZE
set conc_pol    1.      ;#[g/L] concentration of polymer FOR BOX SIZE
set box_size	 [expr 2*$N_B]
#set box_size    [expr round(3.38339*($Mr*$N_M/$conc_pol)**(1./3.))];#size of cubic box
set lhalf	[expr $box_size/2]
set V	[expr $box_size**(3.)]
set r_m_count 0.0;#e.g. 0.2;ratio of multivalent counterions and counterions
set q_m_count 2 ;#charge of multivalent counterions
set N_S	0	;#number of salt units
#set N_S [expr $N_M*2];#number of salts c=N_S/V
### INTERACTIONS
set l_B 4  	;#CI
set p3m_acc 1e-4;#CI
set r_0 1.0	;#HB
set k_B 1000.0	;#HB
set K   7.0     ;#FENE
set dr  2.0     ;#FENE
set r0  0.0     ;#FENE
set fi_0 [expr 109.5/180*3.14];#HA
set k_A 10.0	;#HA
set HA  "yes"	;#use HA?	
set eps 1.0	;#LJ
set sgm 1.0	;#LJ
### GOOD SOLVENT CONDITION
set Gr_c [expr 2**(1/6)];#LJ
set Gc_s 0.25	;#LJ
### THETA SOLVENT CONDITION
set Tr_c 1.29	;#LJ
set Tc_s 0.169912;#LJ
### POOR SOLVENT CONDITION
set Pr_c 1.75	;#LJ
set Pc_s 0.0336033;#LJ
### TESTING SOLVENT
set Nc_s 0.075   ;# just for testing
set Nr_c [expr $sgm*((1+(1-4*$Nc_s)**0.5)/2/$Nc_s)**(1/6.)]


### OTHER PARAMETERS

set kbT		1.0;#fce teploty
set gamma	1.0;#LG
set friction    13.0    ;#LB
set agrid       1       ;#LB
set dens        0.8     ;#LB
set visc        2.8     ;#LB
set tau         0.0125  ;#LB
set bd_visc	2.0	;#BD
set bd_radius   1.0	;#BD
###
set test_mod    no    	;#yes/no        kratka simulace + vmd
set vmd         yes     ;#yes/no
set centering	no     	;#yes/no move molecule to center of box
set term        "LG"    ;#LG/LB/BD (Langevin/lattice Boltzmann/Brownian dynamics)
###
set unit        cpu     ;#cpu/gpu
set P3M         yes     ;#yes/no (electrostatics)
set pretuned    no      ;#yes/no (no for finding new parameters)
set mesh        12      ;
set cao         4       ;
set r_cut       0.275252;
set alpha       7.82625 ;
setmd box_l	$box_size $box_size $box_size
setmd time_step 0.0125
setmd skin	0.4
t_random seed [expr abs([clock clicks]%100000)]

### INTERACTIONS

inter 0 harmonic $k_B $r_0
inter 1 angle_harmonic $k_A $fi_0
inter 2 fene $K $dr $r0
inter 0 0 lennard-jones $eps $sgm $Nr_c $Nc_s
inter 0 1 lennard-jones $eps $sgm $Nr_c $Nc_s
inter 1 1 lennard-jones $eps $sgm $Nr_c $Nc_s
inter 0 2 lennard-jones $eps $sgm $Nr_c $Nc_s
inter 1 2 lennard-jones $eps $sgm $Nr_c $Nc_s
inter 2 2 lennard-jones $eps $sgm $Nr_c $Nc_s

inter 0 3 lennard-jones $eps $sgm $Gr_c $Gc_s
inter 1 3 lennard-jones $eps $sgm $Gr_c $Gc_s
inter 2 3 lennard-jones $eps $sgm $Gr_c $Gc_s
inter 3 3 lennard-jones $eps $sgm $Gr_c $Gc_s

### STRUCTURE
set t0 1;if {$P3M==yes && $h==1 && $N_A>0} {set t0 0}
set t1 1;if {$N_A<2} {set t1 0}
set t2 1;
set t3 1;if {$P3M==no} {set t3 0}
set t4 1;if {$N_S==0} {set t4 0}
if { $N_A==0 } { part 0 pos 0.0 0.0 0.0 type 2 q 1.0
} else { part 0 pos 0.0 0.0 0.0 type 0 }
for { set i 1 } { $i < $N_B } { incr i 1 } {
 if { $i%$h==0 && $N_A==0 } {
  part $i pos [expr $i*$r_0] 0.0 0.0 type 2 bond 0 [expr $i-1] q 1.0
 } else {
  part $i pos [expr $i*$r_0] 0.0 0.0 type 0 bond 0 [expr $i-1]
 }
 if { $i>1 && $HA == "yes" } { part [expr $i-1] bond 1 [expr $i-2] $i }
}
set NP $N_B
if { $N_A > 0 } {
 for { set i 0 } { $i < $N_B } { incr i $h } {
  if { $N_A > 1 } {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 type 1 bond 0 $i
   if { [expr $i-1] >= 0 && $HA == "yes"} {part [expr $i] bond 1 [expr $i-1] $NP}
   if { [expr $i+1] <= $N_B && $HA == "yes" } {part [expr $i] bond 1 [expr $i+1] $NP}
   incr NP 1
   for { set j 2 } { $j <= [expr $N_A-1] } { incr j 1 } {
    part $NP pos [expr $i*$r_0] [expr $j*$r_0] 0.0 type 1 bond 0 [expr $NP-1]
    if { $HA == "yes" } { if { $j>2 } {part [expr $NP-1] bond 1 [expr $NP-2] $NP 
     } else { part [expr $NP-1] bond 1 [expr $i] $NP } 
    }
    incr NP 1
   }
   part $NP pos [expr $i*$r_0] [expr $N_A*$r_0] 0.0 type 2 bond 0 [expr $NP-1] q 1.0
   if { $HA == "yes" } { if { $N_A>2 } { part [expr $NP-1] bond 1 [expr $NP-2] $NP 
    } else { part [expr $NP-1] bond 1 [expr $i-1] $NP }
   }
   incr NP 1
  } else {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 type 2 bond 0 [expr $i] q 1.0
   if { [expr $i-1] >= 0 && $HA == "yes"} {part [expr $i] bond 1 [expr $i-1] $NP}
   if { [expr $i+1] <= $N_B && $HA == "yes" } {part [expr $i] bond 1 [expr $i+1] $NP}
   incr NP 1
  }
 }
}
#COUNTERIONS
for { set i 1 } { $i <= $N_B } { incr i $h } {
 part $NP pos [expr $i*$r_0] [expr ($N_A+2)*$r_0] 0.0 type 3 q -1.0
 incr NP 1
}
#MOVE
set xc [lindex [analyze centermass 0] 0]
set yc [lindex [analyze centermass 0] 1]
set zc [lindex [analyze centermass 0] 2]
for { set i 0 } { $i < [expr $N_M+$N_S+($N_B+1)/$h*$r_m_count/$q_m_count+($N_B+1)/$h*(1-$r_m_count)] } { incr i } {
  set nx [expr [lindex [part $i print pos] 0]-$xc+$lhalf]
  set ny [expr [lindex [part $i print pos] 1]-$yc+$lhalf]
  set nz [expr [lindex [part $i print pos] 2]-$zc+$lhalf]
  part $i pos $nx $ny $nz
}
###########################################
###########################################
########### PRESIMULATION #################
###########################################
###########################################
if {$term=="LB"} {
        if {$unit=="cpu"} {
                lbfluid agrid $agrid dens $dens visc $visc tau $tau friction $friction
	} else {
                lbfluid $unit agrid $agrid dens $dens visc $visc tau $tau friction $friction
        }
}
if {$term=="LG"} {
 #TODO cellsystem nejde scitat po jednom
 if {$P3M=="no"} {cellsystem nsquare}
 thermostat langevin $kbT $gamma }
if {$term=="LB"} {
 cellsystem domain_decomposition -no_verlet_list
 thermostat lb $kbT }
if {$term=="BD"} {
 cellsystem domain_decomposition -no_verlet_list
 setmd sd_visc $bd_visc
 setmd sd_radius $bd_radius
 thermostat bd $kbT}
#USER OUTPUT
if {$test_mod=="yes"} {set vmd "yes"}
puts "box size: $box_size lb units"
puts "koncentrace polymeru: [expr 3.38339**3.*($Mr*$N_M/$box_size**3.)] g/L"
puts "koncentrace co-iontu: [expr 3.38339**3.*($N_S/$box_size**3.)] mol/L"
if {$vmd==yes} {
        set vtf_file    [open "STRUCTURES/vmdsimulation.vtf" "w"]
        writevsf $vtf_file}
#-----------------------------------
#|       l1        |      l2
#|  m1 |     |     |  m2  |      | integrate*time_step
#||||||||||||||||||||||||||||||||| time_step
#-----------------------------------
if {$test_mod==yes} {set m1 80  ;set l1 100      ;
                     set m2 80  ;set l2 $steps   ;
} else {             set m1 80  ;set l1 1000     ;
                     set m2 8000;set l2 $steps;}
####### WARMING UP
for { set krok 1 } { $krok <= 10 } { incr krok 1 } {
 puts "$krok/$l1 - warming up"
 inter forcecap [expr $krok*2]
 integrate $m1
 if { $centering == yes } {
  #MOVE
  set xc [lindex [analyze centermass 0] 0]
  set yc [lindex [analyze centermass 0] 1]
  set zc [lindex [analyze centermass 0] 2]
  for { set i 0 } { $i < [expr $N_M+$N_S+($N_B+1)/$h*$r_m_count/$q_m_count+($N_B+1)/$h*(1-$r_m_count)] } { incr i } {
    set nx [expr [lindex [part $i print pos] 0]-$xc+$lhalf]
    set ny [expr [lindex [part $i print pos] 1]-$yc+$lhalf]
    set nz [expr [lindex [part $i print pos] 2]-$zc+$lhalf]
    part $i pos $nx $ny $nz
  }
 }
 if { $vmd == "yes" } {writevcf $vtf_file folded}
}
inter forcecap 0
####### COULOMBIC INTERACTIONS
if {$P3M=="yes"} {
	puts "Coulombic interaction - ON"
	set casik [expr [clock seconds]-$cas]
	puts "CAS=$casik s (tj. [expr $casik/60] min)"
        if {$pretuned=="no"} {
                if {$unit=="cpu"} {
                        puts [inter coulomb $l_B p3m tunev2 accuracy $p3m_acc]
                } else {
                        puts [inter coulomb $l_B p3m $unit tunev2 accuracy $p3m_acc]
                }
        } else {
                if {$unit=="cpu"} {
                        inter coulomb $l_B p3m $r_cut $mesh $cao $alpha
                } else {
                        inter coulomb $l_B p3m $unit $r_cut $mesh $cao $alpha
                }
        }
}
set casik [expr [clock seconds]-$cas]
puts "CAS=$casik s (tj. [expr $casik/60] min)"
####### WARMING UP 2
for { set krok 11 } { $krok <= 20 } { incr krok 1 } {
 puts "$krok/$l1 - warming up 2"
 inter forcecap [expr $krok*2]
 integrate $m1
 if { $centering == yes } {
  #MOVE
  set xc [lindex [analyze centermass 0] 0]
  set yc [lindex [analyze centermass 0] 1]
  set zc [lindex [analyze centermass 0] 2]
  for { set i 0 } { $i < [expr $N_M+$N_S+($N_B+1)/$h*$r_m_count/$q_m_count+($N_B+1)/$h*(1-$r_m_count)] } { incr i } {
    set nx [expr [lindex [part $i print pos] 0]-$xc+$lhalf]
    set ny [expr [lindex [part $i print pos] 1]-$yc+$lhalf]
    set nz [expr [lindex [part $i print pos] 2]-$zc+$lhalf]
    part $i pos $nx $ny $nz
  }
 }
 if { $vmd == "yes" } {writevcf $vtf_file folded}
}
inter forcecap 0
####### EQULIBRATION
for { set krok 20 } { $krok <= $l1 } { incr krok 1 } {
 puts "$krok/$l1 - equlibrating"
 integrate [expr $m1*10]
 if { $centering == yes } {
  #MOVE
  set xc [lindex [analyze centermass 0] 0]
  set yc [lindex [analyze centermass 0] 1]
  set zc [lindex [analyze centermass 0] 2]
  for { set i 0 } { $i < [expr $N_M+$N_S+($N_B+1)/$h*$r_m_count/$q_m_count+($N_B+1)/$h*(1-$r_m_count)] } { incr i } {
    set nx [expr [lindex [part $i print pos] 0]-$xc+$lhalf]
    set ny [expr [lindex [part $i print pos] 1]-$yc+$lhalf]
    set nz [expr [lindex [part $i print pos] 2]-$zc+$lhalf]
    part $i pos $nx $ny $nz
  }
 }
 if { $vmd == "yes" } {writevcf $vtf_file folded}
}
#####################################################
#####                                           #####
#####################################################
############        SIMULATION       ################
#####################################################
#####                                           #####
#####################################################
set tot_time [expr $l2*$m2*0.0125]
set FILEchd [open "RESULTS/ch_distances" "w"]
set FILEa_backbone [open "RESULTS/anisotropy_backbone.txt" "w"]
set FILEa_all [open "RESULTS/anisotropy_all.txt" "w"]
#RDF
set FILErdf  [open "FUNCTIONS/rdf_2_3" "w"]
#STRF
if {$t0==yes} 				{set FILEstrf0   [open "FUNCTIONS/strf_0"   "w"]}
if {$t1==yes} 				{set FILEstrf1   [open "FUNCTIONS/strf_1"   "w"]}
if {$t2==yes} 				{set FILEstrf2   [open "FUNCTIONS/strf_2"   "w"]}
if {$t0==yes && $t2==yes} 		{set FILEstrf02  [open "FUNCTIONS/strf_02"  "w"]}
if {$t0==yes && $t1==yes && $t2==yes} 	{set FILEstrf012 [open "FUNCTIONS/strf_012" "w"]}
if {$t1==yes && $t2==yes} 		{set FILEstrf12  [open "FUNCTIONS/strf_12"  "w"]}
#CORR
for { set i 0 } { $i <= 3 } { incr i } { if { [set t$i] == 1 } {
 set FILEpos$i  [open "STRUCTURES/coordinates$i" "w"]
 set pos$i      [observable new particle_positions types $i]
 set center$i   [observable new com_position types $i]
 set c_pos$i    [correlation new obs1 [set pos$i] corr_operation square_distance_componentwise tau_lin 16 tau_max $tot_time dt 1]
 set c_center$i [correlation new obs1 [set center$i] corr_operation square_distance_componentwise tau_lin 16 tau_max $tot_time dt 1]
 correlation [set c_pos$i]    autoupdate start
 correlation [set c_center$i] autoupdate start
}}
if {$term=="LB"} {set inf2_file    [open "RESULTS/t_re_rg_rh_toten_kinen_coulen_flmass_flmomx_flmomy_flmomz.out" "w"]}
if {$term=="LG"} {set inf2_file    [open "RESULTS/t_re_rg_rh_toten_kinen_coulen.out" "w"]}
if {$term=="LB"} {puts $inf2_file "# 1_time 2_re 3_rg 4_rh 5_tot.en. 6_kin.en 7_coul.en 8_fl.mass 9_fl.mom_x_y_z ||$m1/$l1 $m2/$l2"}
if {$term=="LG"} {puts $inf2_file "# 1_time 2_re 3_rg 4_rh 5_tot.en. 6_kin.en 7_coul.en ||$m1/$l1 $m2/$l2"}
set mytime [clock seconds]
set casik [expr [clock seconds]-$cas]
puts "CAS=$casik s (tj. [expr $casik/60] min)"
puts "Simulation"

for { set krok 1 } { $krok <= $l2 } { incr krok 1 } {
 if { $term == "BD" } {integrate_sd $m2} else {integrate $m2}
 if { $centering == yes } {
  #MOVE
  set xc [lindex [analyze centermass 0] 0]
  set yc [lindex [analyze centermass 0] 1]
  set zc [lindex [analyze centermass 0] 2]
  for { set i 0 } { $i < [expr $N_M+$N_S+($N_B+1)/$h*$r_m_count/$q_m_count+($N_B+1)/$h*(1-$r_m_count)] } { incr i } {
    set nx [expr [lindex [part $i print pos] 0]-$xc+$lhalf]
    set ny [expr [lindex [part $i print pos] 1]-$yc+$lhalf]
    set nz [expr [lindex [part $i print pos] 2]-$zc+$lhalf]
    part $i pos $nx $ny $nz
  }
 }
 if { $vmd == "yes" } {writevcf $vtf_file folded}
 set doba [expr [clock seconds]-$mytime]
 set doba_kroku [expr $doba/($krok+0.000001)]
 set doba_celkem [expr $l2*$doba_kroku]
 puts "[format "%.4g" $krok]/[format "%.4g" $l2]\t| [format "%.3e" $doba]s/ [format "%.2g" $doba_celkem]s \t= [format "%.2g" [expr $doba_celkem/60]]min \t= [format "%.2g" [expr $doba_celkem/60./60.]]h \t|| [format "%.4g" [expr ($doba/($doba_celkem+0.000001))*100.]] % ||"

 ##############################
 ########## ANALYSE ###########
 ##############################
#POSITIONS
 for { set i 0 } { $i <= 3 } { incr i } { if { [set t$i] == 1 } {
  puts [set FILEpos$i] [observable [set pos$i] print] 
 }}
#RDF
 set rdf [lindex [analyze rdf 2 3 0.0 10 200] 1]
 if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {lappend list [lindex $rdf $i 0]};puts [set FILErdf] $list}
 set list "";for {set i 1} {$i<100} {incr i} {lappend list [lindex $rdf $i 1]};puts [set FILErdf] $list
#STRF
 if {$t0==yes} 				{set strf [analyze structurefactor 0 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf0] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf0] $list}
 if {$t1==yes} 				{set strf [analyze structurefactor 1 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf1] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf1] $list}
 if {$t2==yes} 				{set strf [analyze structurefactor 2 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf2] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf2] $list}
 if {$t0==yes && $t2==yes} 		{set strf [analyze structurefactor {0 2} 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf02] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf02] $list}
 if {$t0==yes && $t1==yes && $t2==yes} 	{set strf [analyze structurefactor {0 1 2} 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf012] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf012] $list}
 if {$t1==yes && $t2==yes} 		{set strf [analyze structurefactor {1 2} 100];if {$krok==1} {set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 0];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf12] $list}; set list "";for {set i 1} {$i<100} {incr i} {set lnew [lindex $strf $i 1];if {$lnew!=""} {lappend list $lnew}};puts [set FILEstrf02] $list}
##################
 set t [expr 0.0125*$m2*$krok]
 set xc [lindex [analyze centermass 0] 0]
 set yc [lindex [analyze centermass 0] 1]
 set zc [lindex [analyze centermass 0] 2]
## Rg ##
 if {$krok==1} {set FILErg [open "RESULTS/rg.txt" "w"];puts $FILErg "# t backbone all sidechains"}
 puts $FILErg "$t [lindex [analyze rg 0 1 $N_B] 0] [lindex [analyze rg $N_B 1 $N_M] 0] [lindex [analyze rg 0 1 [expr $N_M-$N_B]] 0]"
## Re ##
 if {$krok==1} {set FILEre [open "RESULTS/re.txt" "w"];puts $FILEre "# t backbone"}
 puts $FILEre "$t [lindex [analyze re 0 1 $N_B] 0]"
## gyration tensor backbone ##
 set A11 0;set A12 0;set A13 0;set A22 0;set A23 0;set A33 0;
 for { set i 0 } { $i <  $N_B } {incr i} {
  set xi [expr [lindex [part $i print pos] 0]-$xc]
  set yi [expr [lindex [part $i print pos] 1]-$yc]
  set zi [expr [lindex [part $i print pos] 2]-$zc]
  set A11 [expr $A11+$xi**2]
  set A12 [expr $A12+$xi*$yi]
  set A13 [expr $A13+$xi*$zi]
  set A22 [expr $A22+$yi**2]
  set A23 [expr $A23+$yi*$zi]
  set A33 [expr $A33+$zi**2]
 }
  set A11 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A11]
  set A12 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A12]
  set A13 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A13]
  set A22 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A22]
  set A23 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A23]
  set A33 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A33]
  set A21 $A12;set A31 $A13;set A32 $A23;
  puts $FILEa_backbone "$t $A11 $A21 $A31 $A12 $A22 $A32 $A13 $A23 $A33"
## gyration tensor all ##
 set A11 0;set A12 0;set A13 0;set A22 0;set A23 0;set A33 0;
 for { set i 0 } { $i < [expr $N_B+($N_B+1)/$h*$N_A] } {incr i} {
  set xi [expr [lindex [part $i print pos] 0]-$xc]
  set yi [expr [lindex [part $i print pos] 1]-$yc]
  set zi [expr [lindex [part $i print pos] 2]-$zc]
  set A11 [expr $A11+$xi**2]
  set A12 [expr $A12+$xi*$yi]
  set A13 [expr $A13+$xi*$zi]
  set A22 [expr $A22+$yi**2]
  set A23 [expr $A23+$yi*$zi]
  set A33 [expr $A33+$zi**2]
 }
  set A11 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A11]
  set A12 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A12]
  set A13 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A13]
  set A22 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A22]
  set A23 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A23]
  set A33 [expr 1./($N_B+($N_B+1)/$h*$N_A)*$A33]
  set A21 $A12;set A31 $A13;set A32 $A23;
  puts $FILEa_all "$t $A11 $A21 $A31 $A12 $A22 $A32 $A13 $A23 $A33"
## CHARGE DISTANCES ##
 set first -1;set second -1;set third -1;
 set ch1  0;set ch2  0;set ch3  0;
 set div1 0;set div2 0;set div3 0;
 for { set i 0 } { $i < [expr $N_B+($N_B+1)/$h*$N_A] } {incr i} {
  if { [part $i print q] == 1 } {
   if { $first == -1 } {
    set first $i
   } else { 
    if { $second == -1 } {
     set dx [expr [lindex [part $i print pos] 0]-[lindex [part $first print pos] 0]]
     set dy [expr [lindex [part $i print pos] 1]-[lindex [part $first print pos] 1]]
     set dz [expr [lindex [part $i print pos] 2]-[lindex [part $first print pos] 2]]
     set ch1 [expr $ch1+($dx**2+$dy**2+$dz**2)**0.5]
     incr div1
     set second $first
     set first  $i
    } else {
     if { $third == -1 } {
      set dx [expr [lindex [part $i print pos] 0]-[lindex [part $first print pos] 0]]
      set dy [expr [lindex [part $i print pos] 1]-[lindex [part $first print pos] 1]]
      set dz [expr [lindex [part $i print pos] 2]-[lindex [part $first print pos] 2]]
      set ch1 [expr $ch1+($dx**2+$dy**2+$dz**2)**0.5]
      set dx [expr [lindex [part $i print pos] 0]-[lindex [part $second print pos] 0]]
      set dy [expr [lindex [part $i print pos] 1]-[lindex [part $second print pos] 1]]
      set dz [expr [lindex [part $i print pos] 2]-[lindex [part $second print pos] 2]]
      set ch2 [expr $ch2+($dx**2+$dy**2+$dz**2)**0.5]
      incr div1;incr div2;
      set third $second;set second $first;set first  $i
     } else {
      set dx [expr [lindex [part $i print pos] 0]-[lindex [part $first print pos] 0]]
      set dy [expr [lindex [part $i print pos] 1]-[lindex [part $first print pos] 1]]
      set dz [expr [lindex [part $i print pos] 2]-[lindex [part $first print pos] 2]]
      set ch1 [expr $ch1+($dx**2+$dy**2+$dz**2)**0.5]
      set dx [expr [lindex [part $i print pos] 0]-[lindex [part $second print pos] 0]]
      set dy [expr [lindex [part $i print pos] 1]-[lindex [part $second print pos] 1]]
      set dz [expr [lindex [part $i print pos] 2]-[lindex [part $second print pos] 2]]
      set ch2 [expr $ch2+($dx**2+$dy**2+$dz**2)**0.5]
      set dx [expr [lindex [part $i print pos] 0]-[lindex [part $third print pos] 0]]
      set dy [expr [lindex [part $i print pos] 1]-[lindex [part $third print pos] 1]]
      set dz [expr [lindex [part $i print pos] 2]-[lindex [part $third print pos] 2]]
      set ch3 [expr $ch2+($dx**2+$dy**2+$dz**2)**0.5]
      incr div1;incr div2;incr div3;
      set third $second;set second $first;set first  $i 
     }
    }
   }
  }
 }
 if {$krok==1} {puts $FILEchd "#t ch_1 ch_2 ch_3"}
 puts $FILEchd "t [expr $ch1/$div1] [expr $ch2/$div2] [expr $ch3/$div3]"
}
close $FILEchd
close $FILErdf
close $FILEa_backbone
close $FILEa_all
if {$t0==yes} 				{close $FILEstrf0   }
if {$t1==yes} 				{close $FILEstrf1   }
if {$t2==yes} 				{close $FILEstrf2   }
if {$t0==yes && $t2==yes} 		{close $FILEstrf02  }
if {$t0==yes && $t1==yes && $t2==yes} 	{close $FILEstrf012 }
if {$t1==yes && $t2==yes} 		{close $FILEstrf12  }
for { set i 0 } { $i <= 3 } { incr i } { if { [set t$i] == 1 } {
 close [set FILEpos$i]
 correlation [set c_pos$i]    finalize
 correlation [set c_pos$i]    write_to_file "CORRELATIONS/t_sl_CORRpos$i.dat"
 correlation [set c_center$i] finalize
 correlation [set c_center$i] write_to_file "CORRELATIONS/t_sl_CORRcenter$i.dat"
}}

puts "FINISHED"
