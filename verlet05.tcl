#VERLET
#SIMULACION II - NP FB 2017

#CONSTANTES
#set ep 0.661
#set sg 0.3345

set ep 145
set sg 3.8
set m 39.95
set L 10

#CONTADOR DE ITERACIONES
set k 1

#ENERGIA INICIAL
set e1 0
set e2 0

#TIMESTEP
set t 0.05


########################################################
#SE ELIGEN TODAS LAS PARTICULAS
set particulas [atomselect top all]
set id [$particulas get resid]
set n [llength $id]

set systemTime [clock seconds]

print "INICIANDO VERLET - NP,FB 2017: [clock format $systemTime]\n"

print "particulas: $n"


print "==================="
print "ALGORITMO DE VERLET"
print "ep=$ep sg=$sg"

print ""


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
    
    set Rx($s) 0
    set Ry($s) 0
    set Rz($s) 0    
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



print "==================="
for {set i 0 } {$i < $n } { incr i } {
  set p1 [lindex $coord $i]

  print "(P$i:$p1) (Vx$i:$Vx($i) Vy$i:$Vy($i) Vz$i:$Vz($i)) " 
}
print "==================="

for { set k 0 } {$k < 3 } { incr k } {

 print ""
 print "***PASO $k"

 
 generarXYZ $n
 generarXYZ $k

 for { set i 0 } {$i < $n } { incr i } {
  for { set j 0 } {$j < $i } { incr j } {
  
    print ""
    print "comparando: $i con $j"
  
    set p1 [lindex $coord $i]
    set p2 [lindex $coord $j]  
    
    print "P1:$p1" 
    print "P2:$p2" 
    
    #CALCULO DE DISTANCIAS
    set x1 [lindex $p1 0]
    set y1 [lindex $p1 1]
    set z1 [lindex $p1 2]   
    
    set Rx($i) $x1
    set Ry($i) $y1
    set Rz($i) $z1
    
    set x2 [lindex $p2 0]
    set y2 [lindex $p2 1]
    set z2 [lindex $p2 2]

    set dx [expr {$x1-$x2}]
    set dy [expr {$y1-$y2}]
    set dz [expr {$z1-$z2}]
    
    #CONDICIONES PERIODICAS DE BORDE
    
    set dx [expr {$dx - $L*round(double($dx)/$L)}]
    set dy [expr {$dy - $L*round(double($dy)/$L)}]
    set dz [expr {$dz - $L*round(double($dz)/$L)}]


    set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
    lappend d $dist
    
    print "(x1:$x1 y1:$y1 z1:$z1) - (x2:$x2 y2:$y2 z2:$z2)"
    print "(dx:$dx dy:$dy dz:$dz) - dist:$dist"    
    
    
    #CALCULO DE FUERZAS A PARTIR DEL POTENCIAL DE LENNARD JONES
    set r2 [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
    
    ####################################CALCULO DE ENERGIA Y GRADIENTE
    set r6 [expr {$r2*$r2*$r2}]
    set r12 [expr {$r6*$r6}]
    set r6i [expr {1.0/($r2*$r2*$r2)}]
    
    #PRIMERA FORMULA DE CALCULAR e
    set e1 [expr {$e1+4*(1.0/$r12-1.0/$r6)}]
    
    #SEGUNDA FORMA DE CALCULAR e
    set e2 [expr {$e2-4*$r6i*(1-$r6i)}]
    
    
    print "r6:$r6 r12:$r12 e1:$e1 e2:$e2"
    ####################################
    
    #Clj(r)
    set f [expr {(48/$r2)*($r6i*$r6i-0.5*$r6i)}]
    
    set Fx($i) [expr {$Fx($i)+$dx*$f*double(1/$r2)}]
    set Fx($j) [expr {$Fx($j)-$dx*$f*double(1/$r2)}]
    
    set Fy($i) [expr {$Fy($i)+$dy*$f*double(1/$r2)}]
    set Fy($j) [expr {$Fy($j)-$dy*$f*double(1/$r2)}]
    
    set Fz($i) [expr {$Fz($i)+$dz*$f*double(1/$r2)}]
    set Fz($j) [expr {$Fz($j)-$dz*$f*double(1/$r2)}]    
    

    print "r2:$r2 -> r6i:$r6i"
    print "Fx($i):$Fx($i) Fy($i):$Fy($i) Fz($i):$Fz($i)"

  }
 }

 print ""
 print "%%%%%%%%%%%%%%%%%%%%%%%%RESUMEN PASO $k"
 for { set i 0 } {$i < $n } { incr i } {
   print "Fx($i):$Fx($i) Fy($i):$Fy($i) Fz($i):$Fz($i)"
 }
 
 ###################ACELERACIONES
 print ""
 print "ACELERACIONES TOTALES"
  
 ;#set m [$p1 get mass]
 for { set s 0 } {$s < $n } { incr s } {
    set Ax($s) [expr {$Fx($s)*$m}]
    set Ay($s) [expr {$Fy($s)*$m}]
    set Az($s) [expr {$Fz($s)*$m}]
    print "$s -> Ax:$Ax($s) Ay:$Ay($s) Az:$Az($s)"
 }
 
 print "&&&&&&&&&&&&&&&"
 print "VELOCIDADES TOTALES"
  
 for { set s 0 } {$s < $n } { incr s } {
    set Vx($s) [expr {$Vx($s)+($Ax($s)+$Ax($s))*$t}]
    set Vy($s) [expr {$Vy($s)+($Ay($s)+$Ay($s))*$t}]
    set Vz($s) [expr {$Vz($s)+($Az($s)+$Az($s))*$t}]
    print "$s -> Vx:$Vx($s) Vy:$Vy($s) Vz:$Vz($s)"
    
    #CALCULO DE ENERGIA CINETICA
    set K [expr {$m*($Vx($s)*$Vx($s)+$Vy($s)*$Vy($s)+$Vy($s)*$Vy($s))}]
    lappend Ki $K
 } 
 
 print ""
 
 #SUMATORIA DE ENERGIAS CINETICAS
 set K_total [expr {[lindex $Ki 0]+[lindex $Ki 1]+[lindex $Ki 2]}]
 print "ENERGIA CINETICA K:$K_total"
 
 
 set k_ 1
 
 #CALCULO DE TEMPERATURA CINETICA INSTANTANEA
 set Temp [expr {(2*$K_total)/($n*$k_)}]
 print "TEMPERATURA CINETICA INSTANTANEA T:$Temp"
 print ""
 
 
 print "&&&&&&&&&&&&&&&"
 print "POSICIONES TOTALES"
  
 for { set s 0 } {$s < $n } { incr s } {
    set Rx($s) [expr {$Rx($s)+$Vx($s)*$t+$Ax($s)*$t*$t/2}]
    set Ry($s) [expr {$Ry($s)+$Vy($s)*$t+$Ax($s)*$t*$t/2}]
    set Rz($s) [expr {$Rz($s)+$Vz($s)*$t+$Ax($s)*$t*$t/2}]
    print "$s -> Rx:$Rx($s) Ry:$Ry($s) Rz:$Rz($s)"
    
    set aux "Ar $Rx($s)   $Rx($s)   $Rx($s)"
    generarXYZ $aux
    
    #ACTUALIZO POSICIONES
    ;#set xf $Rx($s)
    ;#set yf $Ry($s)
    ;#set zf $Rz($s)
    
    ;#set idx1 [lindex $id $s]
    ;#set p1 [atomselect top "name OH2 and resid $idx1"]
    ;#$p1 moveby [list $xf $yf $zf]
 }
 
}

GdR $d $n
unset d

puts ""
parray Rx

puts ""
parray Ry

puts ""
parray Rz

print "\n\n\n"

#FUNCION PARA GUARDAR EN TEXTO
proc grabar {dato} {
  set outfile [open output1.txt a]
  puts $outfile $dato
  close $outfile
}

proc print {dato} {
  puts $dato
  grabar "$dato"
}

proc GdR {d n} {
  set max_r 8
  set dr 0.2
  
  set vol [expr {4*3.14*$max_r*$max_r*$dr}]
  set den [expr {$n/$vol}]
  
  puts ""
  puts ">>>d"
  puts $n
  puts $d
  puts ""
  
  for { set j 0 } {$j < $max_r/$dr } { incr j } { 
    set a1 [expr {$dr*$j}]
    set a2 [expr {$dr*($j+1)}]

    set cont 0
    foreach r $d {
      #SI LA DISTANCIA ES MENOR AL RADIO MAXIMO DE LA ESFERA
      if {$r<$max_r} {
        #SI LA DISTANCIA ESTA DENTRO DEL RANGO
        if {$r>$a1 && $r<$a2} {
          incr cont
        }
      }
    }
    set cont [expr {$cont/double($n*$vol)}]
    puts "j:$a1 j+dr:$a2 cont:$cont"
    lappend xDist $a1
    lappend yProb $cont
    set cont 0
  }
  graficar $xDist $yProb
}

proc graficar {xDist yProb} {
  set plothandle [multiplot -x $xDist -y $yProb -title "GdR" -lines -linewidth 3 -marker point -plot]
  $plothandle configure -fillcolor yellow -radius 3
  $plothandle configure -vline {3 -width 2 -fill red -dash "."}
  $plothandle replot;
}

proc generarXYZ {dato} {
  set outfile2 [open trayectoria.xyz a]
  puts $outfile2 $dato
  close $outfile2
}

#-----------------------------------------------------[CORTAR AQUI 1
#3
#archivo XYZ argon

#Ar 0.0 0.0 0.0
#Ar 1.0 0.0 0.0
#Ar 0.0 1.0 0.0

#Ar 0.1 0.1 0.1
#Ar 1.2 0.2 0.2
#Ar 0.1 1.2 0.1

#-----------------------------------------------------[CORTAR AQUI 2
#CRYST1   34.739   34.739   34.739  90.00  90.00  90.00 P 1           1
#ATOM      1  AR  AR      1     0.000   0.000   0.000  0.00  0.00            
#ATOM      2  AR  AR      2     4.000   0.000   0.000  0.00  0.00            
#ATOM      3  AR  AR      3     0.000   4.000   0.000  0.00  0.00             
#END

#-----------------------------------------------------[CORTAR AQUI 3
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