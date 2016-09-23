# Brush simulation
# Jakub Kubecka
# 12.9.2016

set steps 1000

### PARAMETERS

set N_B 31 	;#monomer units of main chain
set N_A 2 	;#monomer units of side chain
set h	2  	;#space between brush chains
set l_B 4  	;#CI
set acc 1e-04	;#CI
set r_0 1.0	;#HB
set k_B 1000.0	;#HB
set fi0 [expr 109.5/180*3.14];#HA
set k_A 10.0	;#HA
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


### OTHER PARAMETERS

set kbT		1.0;#langevin
set gamma	1.0;#langevin
set box_size	[expr $N_B*2];
setmd box_l	$box_size $box_size $box_size
setmd time_step 0.0125
setmd skin	0.4

### INTERACTIONS

inter 0 harmonic $k_B $r_0
inter 1 angle_harmonic $k_A $fi0
inter 0 0 lennard-jones $eps $sgm $Pr_c $Pc_s
inter 0 1 lennard-jones $eps $sgm $Gr_c $Gc_s
inter 1 1 lennard-jones $eps $sgm $Gr_c $Gc_s

### STRUCTURE

part 0 pos 0.0 0.0 0.0 type 0
for { set i 1 } { $i < $N_B } { incr i 1 } {
 if { $i%2==1 && $N_A==0 } {
  part $i pos [expr $i*$r_0] 0.0 0.0 type 0 bond 0 [expr $i-1] q 1.0
 } else {
  part $i pos [expr $i*$r_0] 0.0 0.0 type 0 bond 0 [expr $i-1]
 }
 if { $i>2 } { part [expr $i-1] type 0 bond 1 [expr $i-2] $i }
}
set NP $N_B
if { $N_A > 0 } {
 for { set i 1 } { $i <= $N_B } { incr i $h } {
  if { $N_A > 1 } {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 type 0 bond 0 [expr $i-1]
   if {[expr $i-2]>=0} {part [expr $i-1] type 0 bond 1 [expr $i-2] $NP}
   if {[expr $i]<=$N_B} {part [expr $i-1] type 0 bond 1 [expr $i] $NP}
   incr NP 1
   for { set j 2 } { $j <= [expr $N_A-1] } { incr j 1 } {
    part $NP pos [expr $i*$r_0] [expr $j*$r_0] 0.0 type 0 bond 0 [expr $NP-1]
    if { $j>2 } {part [expr $NP-1] type 0 bond 1 [expr $NP-2] $NP 
    } else { part [expr $NP-1] type 0 bond 1 [expr $i-1] $NP }
    incr NP 1
   }
   part $NP pos [expr $i*$r_0] [expr $N_A*$r_0] 0.0 type 0 bond 0 [expr $NP-1] q 1.0
   if { $N_A>2 } {part [expr $NP-1] type 0 bond 1 [expr $NP-2] $NP 
    } else { part [expr $NP-1] type 0 bond 1 [expr $i-1] $NP }
   incr NP 1
  } else {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 type 0 bond 0 [expr $i-1] q 1.0
   incr NP 1
  }
 }
 for { set i 1 } { $i <= $N_B } { incr i $h } {
  part $NP pos [expr $i*$r_0] [expr ($N_A+2)*$r_0] 0.0 type 1 q -1.0
  incr NP 1
 }
}

inter coulomb $l_B p3m tunev2 accuracy $acc
### ANALYSIS

set ifile [open "results.txt" "w"] 
puts $ifile "# t rg re s3 s4 s6"
set afile [open "anisotropy.txt" "w"]

###########################################
###########################################
############ SIMULATION	###################
###########################################
###########################################

thermostat langevin $kbT $gamma
set vtf_file [open "vmdsimulation.vtf" "w"]
writevsf $vtf_file

####### WARMING UP
for { set krok 1 } { $krok <= 10 } { incr krok 1 } {
 puts "Warm $krok/10"
 inter forcecap [expr $krok*2]
 integrate 80
 writevcf $vtf_file folded
}
inter forcecap 0
####### EQULIBRATION
for { set krok 1 } { $krok <= 100 } { incr krok 1 } {
 puts "Equilibration $krok/100"
 integrate 800
 writevcf $vtf_file folded
}
####### SIMULATION
for { set krok 1 } { $krok <= $steps } { incr krok 1 } {
 puts "Simulation $krok/$steps"
 integrate 8000
 writevcf $vtf_file folded
 #############
 ## ANALYSE ##
 #############
 set t [expr 0.0125*8000*$krok]
 set xc [lindex [analyze centermass 0] 0]
 set yc [lindex [analyze centermass 0] 1]
 set zc [lindex [analyze centermass 0] 2]
 ## 1 ##
 set rg [lindex [analyze rg 0 1 $N_B] 0]
 ## 2 ##
 set re [lindex [analyze re 0 1 $N_B] 0]
 ## 3 ##
 set sum3 0
 set div3 0
 for { set i 0 } { $i < [expr $N_B+($N_B+1)/$h*$N_A] } {incr i} {
  if { [part $i print q] == 1 } {
   set dx [expr [lindex [part $i print pos] 0]-$xc]
   set dy [expr [lindex [part $i print pos] 1]-$yc]
   set dz [expr [lindex [part $i print pos] 2]-$zc]
   set sum3 [expr $sum3+($dx**2+$dy**2+$dz**2)**0.5]
   incr div3
  }
 }
 set sum3 [expr $sum3/$div3]
 ## 4 ##
 set sum4 0
 set div4 0
 for { set i 0 } { $i < [expr $N_B+($N_B+1)/$h*$N_A] } {incr i} {
  set dx [expr [lindex [part $i print pos] 0]-$xc]
  set dy [expr [lindex [part $i print pos] 1]-$yc]
  set dz [expr [lindex [part $i print pos] 2]-$zc]
  set sum4 [expr $sum4+($dx**2+$dy**2+$dz**2)**0.5]
  incr div4
 }
 set sum4 [expr $sum4/$div4]
 ## 5 ##
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
  puts $afile "$t $A11 $A21 $A31 $A12 $A22 $A32 $A13 $A23 $A33"
 ## 6 ##
 set first -1
 set sum6 0
 set div6 0
 for { set i 0 } { $i < [expr $N_B+($N_B+1)/$h*$N_A] } {incr i} {
  if { [part $i print q] == 1 } {
   if { $first != -1 } {
    set dx [expr [lindex [part $i print pos] 0]-[lindex [part $first print pos] 0]]
    set dy [expr [lindex [part $i print pos] 1]-[lindex [part $first print pos] 1]]
    set dz [expr [lindex [part $i print pos] 2]-[lindex [part $first print pos] 2]]
    set sum6 [expr $sum6+($dx**2+$dy**2+$dz**2)**0.5]
    incr div6
    set first $i
   } else {
    set first $i
   }
  }
 }
 set sum6 [expr $sum6/$div6]
 #######
 #puts [ analyze structurefactor 0 5 ]
 puts $ifile "$t $rg $re $sum3 $sum4 $sum6" 
}
exit
