#VERLET
#SIMULACION II - 2017

#SE ELIGEN TODAS LAS AGUAS
set aguas [atomselect top "name OH2"]

set id [$aguas get resid]
set n [llength $id]
puts "moleculas: $n"

set ep 0.661
set sg 0.3345

#CONTADOR DE ITERACIONES
set k 1

#CONDICIONES INICIALES 
for { set s 0 } {$s < $n } { incr s } {

    set Fx($s) 0
    set Fy($s) 0
    set Fz($s) 0

    set Vx($s) 0
    set Vy($s) 0
    set Vz($s) 0
  
    set Ax_($s) 0
    set Ay_($s) 0
    set Az_($s) 0
    
    set Ax($s) 0
    set Ay($s) 0
    set Az($s) 0
}

#PARA CADA PAR DE MOLECULAS
for { set i 0 } {$i < $n } { incr i } {
  for { set j 0 } {$j < $n } { incr j } {
    set idx1 [lindex $id $i]
    set idx2 [lindex $id $j]
    
    if {$idx1!=$idx2} {
      
      set p1 [atomselect top "name OH2 and resid $idx1"]
      set p2 [atomselect top "name OH2 and resid $idx2"]
      
      if {$k!=1} {
       lappend particulas $p1
      }
      
      #CALCULO DISTANCIAS ENTRE PARTICULAS
      set x1 [scan [$p1 get {x}] %f]
      set y1 [scan [$p1 get {y}] %f]
      set z1 [scan [$p1 get {z}] %f]
      
      #######################################
      set x1_($i) $x1
      set y1_($i) $y1
      set z1_($i) $z1
      puts "** $i -> $x1_($i) $y1_($i) $z1_($i)"
      ########################################
      
      set x2 [scan [$p2 get {x}] %f]
      set y2 [scan [$p2 get {y}] %f]
      set z2 [scan [$p2 get {z}] %f]
    
      set dx [expr {$x1-$x2}]
      set dy [expr {$y1-$y2}]
      set dz [expr {$z1-$z2}]

      set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
      set theta [expr {atan($dy/$dx)}]
      
      puts "comparando: $k ->$i con $j ($idx1 con $idx2)"
      puts "x1:$x1 y1:$y1 z1:$z1 - x2:$x2 y2:$y2 z2:$z2"
      puts "dx:$dx dy:$dy dz:$dz - dist:$dist"

      #LENNARD JONES CON CAMBIO DE VARIABLES
      set A [expr {4*$ep*pow($sg,2)}]
      set B [expr {4*$ep*pow($sg,6)}]
      
      set r $dist
      
      #COMO F=-dVr SE TIENE QUE F=12A(r^-13)-6B(r^-7)
      set Vr [expr {($A/pow($r,12))-($B/pow($r,6))}]
      set Fr [expr {(12*$A/pow($r,13))-(6*$B/pow($r,7))}]
      
      #DESCOMPOSICION VECTORIAL DE LA FUERZA
      set Frx [expr {($Fr/$r)*$dx}] 
      set Fry [expr {($Fr/$r)*$dy}] 
      set Frz [expr {($Fr/$r)*$dz}] 
      
      puts "$i ->A:$A B:$B Vr:$Vr Fr:$Fr"
      puts "$i ->Frx:$Frx Fry:$Fry Frz:$Frz"
      puts "=================================="
      
      #ACTUALIZO LOS VECTORES DE FUERZA
      if {$k==0} {
        set Fx($i) $Frx
        set Fy($i) $Fry
        set Fz($i) $Frz
      } else {
        set Fx($i) [expr {$Fx($i)+$Frx}]
        set Fy($i) [expr {$Fy($i)+$Fry}]
        set Fz($i) [expr {$Fz($i)+$Frz}]
      }

      incr k

    }
 
  }

  #MUESTRO VECTOR FUERZA
  puts "&&&&&&&&&&&&&&&"
  puts "FUERZAS TOTALES"
  
  for { set s 0 } {$s < $n } { incr s } {
      puts "$s -> Fx:$Fx($s) Fy:$Fy($s) Fz:$Fz($s)"
      set Ax($s) 0
      set Ay($s) 0
      set Az($s) 0
  } 
  
  puts "&&&&&&&&&&&&&&&"
  puts "ACELERACIONES TOTALES"
  
  set m [$p1 get mass]
  for { set s 0 } {$s < $n } { incr s } {
    set Ax($s) [expr {$Fx($s)*$m}]
    set Ay($s) [expr {$Fy($s)*$m}]
    set Az($s) [expr {$Fz($s)*$m}]
    puts "$s -> Ax:$Ax($s) Ay:$Ay($s) Az:$Az($s)"
  }
  
}


  puts "&&&&&&&&&&&&&&&"
  puts "VERLET"
  
  #TIMESTEP
  set t 1
  
  puts "&&&&&&&&&&&&&&&"
  puts "VELOCIDADES TOTALES"
  
  for { set s 0 } {$s < $n } { incr s } {
    set Vx($s) [expr {$Vx($s)+($Ax($s)+$Ax($s))*$t}]
    set Vy($s) [expr {$Vy($s)+($Ay($s)+$Ay($s))*$t}]
    set Vz($s) [expr {$Vz($s)+($Az($s)+$Az($s))*$t}]
    puts "$s -> Vx:$Vx($s) Vy:$Vy($s) Vz:$Vz($s)"
  } 

  puts "&&&&&&&&&&&&&&&"
  puts "POSICIONES TOTALES"
  
  set id [$aguas get resid]
  
  for { set s 0 } {$s < $n } { incr s } {
    set Rx($s) [expr {$x1_($s)+$Vx($s)*$t+$Ax($s)*$t*$t/2}]
    set Ry($s) [expr {$y1_($s)+$Vy($s)*$t+$Ax($s)*$t*$t/2}]
    set Rz($s) [expr {$z1_($s)+$Vz($s)*$t+$Ax($s)*$t*$t/2}]
    puts "$s -> Rx:$Rx($s) Ry:$Ry($s) Rz:$Rz($s)"
    
    #ACTUALIZO POSICIONES
    set xf $Rx($s)
    set yf $Ry($s)
    set zf $Rz($s)
    
    set idx1 [lindex $id $s]
    set p1 [atomselect top "name OH2 and resid $idx1"]
    $p1 moveby [list $xf $yf $zf]
  }  





#CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
#ATOM      1  OH2 TIP3W  48      -2.250   3.844   4.629  1.00  0.00      WT1  O
#ATOM      2  H1  TIP3W  48      -2.890   3.290   4.998  1.00  0.00      WT1  H
#ATOM      3  H2  TIP3W  48      -2.836   4.245   3.962  1.00  0.00      WT1  H
#ATOM      4  OH2 TIP3W  87      -3.339   1.681  -3.057  1.00  0.00      WT1  O
#ATOM      5  H1  TIP3W  87      -3.057   0.739  -3.019  1.00  0.00      WT1  H
#ATOM      6  H2  TIP3W  87      -4.222   1.599  -3.517  1.00  0.00      WT1  H
#ATOM      7  OH2 TIP3W 485      -0.382   3.125   2.777  1.00  0.00      WT1  O
#ATOM      8  H1  TIP3W 485       0.405   2.919   3.314  1.00  0.00      WT1  H
#ATOM      9  H2  TIP3W 485      -0.998   3.458   3.422  1.00  0.00      WT1  H
#ATOM     10  OH2 TIP3W 486      -1.597  -0.313   4.815  1.00  0.00      WT1  O
#ATOM     11  H1  TIP3W 486      -1.891   0.272   4.097  1.00  0.00      WT1  H
#ATOM     12  H2  TIP3W 486      -1.183  -1.027   4.327  1.00  0.00      WT1  H
#ATOM     13  OH2 TIP3W 594      -2.290   1.185   2.540  1.00  0.00      WT1  O
#ATOM     14  H1  TIP3W 594      -2.533   1.236   1.561  1.00  0.00      WT1  H
#ATOM     15  H2  TIP3W 594      -1.590   1.855   2.580  1.00  0.00      WT1  H
#END