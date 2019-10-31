#VERLET
#SIMULACION II - 2017

#SE ELIGEN TODAS LAS PARTICULAS
set aguas [atomselect top all]

set id [$aguas get resid]
set n [llength $id]
puts "moleculas: $n"

set ep 1
set sg 1

#CONTADOR DE ITERACIONES
set k 1

#CONDICIONES INICIALES
# Ver en el tutorial la secciion
# Temperatura y el teorema de la Equipartición de la Energía
# 3.6. Implementación de un programa de Dinámica Molecular: 1 Inicializacion.

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


#VELOCIDADES INICIALES
set Vx(0) 0
set Vy(0) 0
set Vz(0) 0

set Vx(1) 0
set Vy(1) 1
set Vz(1) 0

set Vx(2) 1
set Vy(2) 0
set Vz(2) 0



#PARA CADA PAR DE MOLECULAS
;# por que no recorrer en dos ciclos uno en i > N, y el otro i+1 >N OK
;# luego recalcule el par de atomos con Fi = -Fj                    OK
;# el doble ciclo no se topa dos veces con el mismo par             OK


##############################################################esto esta ok
set sel [atomselect top all]
set coord [$sel get {x y z} ]
;#set m [$sel get mass ]

set m 39.95

puts "==================="
for {set i 0 } {$i < $n } { incr i } {
  set p1 [lindex $coord $i]

  puts "(P$i:$p1) (Vx$i:$Vx($i) Vy$i:$Vy($i) Vz$i:$Vz($i)) " 
}
puts "==================="

#for { set k 0 } {$k < 1 } { incr k } {

# puts "***PASO $k"

 for { set i 0 } {$i < $n } { incr i } {
  for { set j 0 } {$j < $i } { incr j } {
  
    puts ""
    puts "comparando: $i con $j"
  
    set p1 [lindex $coord $i]
    set p2 [lindex $coord $j]  
    
    puts "P1:$p1" 
    puts "P2:$p2" 
    
    #CALCULO DE DISTANCIAS
    set x1 [lindex $p1 0]
    set y1 [lindex $p1 1]
    set z1 [lindex $p1 2]   
    
    set x2 [lindex $p2 0]
    set y2 [lindex $p2 1]
    set z2 [lindex $p2 2]

    set dx [expr {$x1-$x2}]
    set dy [expr {$y1-$y2}]
    set dz [expr {$z1-$z2}]
    
    set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
    
    puts "(x1:$x1 y1:$y1 z1:$z1) - (x2:$x2 y2:$y2 z2:$z2)"
    puts "(dx:$dx dy:$dy dz:$dz) - dist:$dist"    
    
    
    #CALCULO DE FUERZAS A PARTIR DEL POTENCIAL DE LENNARD JONES
    set r2 [expr {$dx + $dy + $dz}]
    
    set r6i [expr {1.0/($r2*$r2*$r2)}]
    set f [expr {48*($r6i*$r6i-0.5*$r6i)}]
    
    set Fx($i) [expr {$Fx($i)+$dx*$f*double(1/$r2)}]
    set Fx($j) [expr {$Fx($j)-$dx*$f*double(1/$r2)}]
    
    set Fy($i) [expr {$Fy($i)+$dy*$f*double(1/$r2)}]
    set Fy($j) [expr {$Fy($j)-$dy*$f*double(1/$r2)}]
    
    set Fz($i) [expr {$Fz($i)+$dz*$f*double(1/$r2)}]
    set Fz($j) [expr {$Fz($j)-$dz*$f*double(1/$r2)}]    
    

    puts "r2:$r2 -> r6i:$r6i"
    puts "Fx($i):$Fx($i) Fy($i):$Fy($i) Fz($i):$Fz($i)"

  }
 }

 puts ""
 puts "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PASO 0"
 for { set i 0 } {$i < $n } { incr i } {
   puts "Fx($i):$Fx($i) Fy($i):$Fy($i) Fz($i):$Fz($i)"
 }
 
#}


#-----------------------------------------------------[CORTAR AQUI 3
#-----------------------------------------------------[CORTAR AQUI 3


#CRYST1   34.739   34.739   34.739  90.00  90.00  90.00 P 1           1
#ATOM      1  AR  AR      1     0.000   0.000   0.000  0.00  0.00            
#ATOM      2  AR  AR      2     4.000   0.000   0.000  0.00  0.00            
#ATOM      3  AR  AR      3     0.000   4.000   0.000  0.00  0.00             
#END


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