# Brush simulation
# Jakub Kubecka
# 12.9.2016

set steps 10000

### PARAMETERS

set N_B	31 	;#monomer units of main chain
set N_A	0 	;#monomer units of side chain
set h	2  	;#space between brush chains
set l_B 1  	;#CI
set acc 1e-04	;#CI
set r_0 1.0	;#HB
set k_B 1000.0	;#HB
set eps 1.0	;#LJ
set sgm 1.0	;#LJ
set r_c [expr 2**(1/6)];#LJ
set c_s 0.25	;#LJ


### INTERACTIONS

inter 0 harmonic $k_B $r_0
inter 0 0 lennard-jones $eps $sgm $r_c $c_s

### OTHER PARAMETERS

set kbT		1.0;#langevin
set gamma	1.0;#langevin
set box_size	[expr $N_B+1];
setmd box_l	$box_size $box_size $box_size
setmd time_step 0.0125
setmd skin	0.4

### STRUCTURE

part 0 pos 0.0 0.0 0.0 
for { set i 1 } { $i < $N_B } { incr i 1 } {
	part $i pos [expr $i*$r_0] 0.0 0.0 bond 0 [expr $i-1]
}
set NP $N_B
if { $N_A > 0 } {
 for { set i 1 } { $i <= $N_B } { incr i $h } {
  if { $N_A > 1 } {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 bond 0 [expr $i-1]
   incr NP 1
   for { set j 2 } { $j <= [expr $N_A-1] } { incr n 1 } {
    part $NP pos [expr $i*$r_0] [expr $j*$r_0] 0.0 bond 0 [expr $NP-1]
    incr NP 1
   }
   part $NP pos [expr $i*$r_0] [expr $N_A*$r_0] 0.0 bond 0 [expr $NP-1] q 1.0
   incr NP 1
  } else {
   part $NP pos [expr $i*$r_0] [expr $r_0] 0.0 bond 0 [expr $i-1] q 1.0
   incr NP 1
  }
 }
 for { set i 1 } { $i <= $N_B } { incr i $h } {
  part $NP pos [expr $i*$r_0] [expr ($N_A+2)*$r_0] 0.0 bond 0 [expr $i-1] q -1.0
  incr NP 1
 }
 inter coulomb $l_B p3m tunev2 accuracy $acc
}

### ANALYSIS

set N [ expr $N_B+($N_B+1)/$h*$N_A ]
set ifile [open "results.txt" "w"] 
puts $ifile "# t rg re1 re2"

### SIMULATION	

thermostat langevin $kbT $gamma

for { set krok 1 } { $krok <= $steps } { incr krok 1 } {
 puts "step $krok"
 integrate 80
 set rg [lindex [analyze rg 0 1 $N] 0]
 set re1 [lindex [analyze re 0 1 $N_B] 0]
 set re2 [lindex [analyze re 0 1 $N] 0]
 #puts [ analyze structurefactor 0 5 ]
 puts $ifile "$krok $rg $re1 $re2" 
}
exit
