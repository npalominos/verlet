#FUNCION DISTRIBUCION RADIAL
#SIMULACION II - 2017

#PARAMETROS INICIALES
set max_r 10
set dr 0.1

#OBTENGO NUMERO DE FRAMES
set nf [molinfo top get numframes]
puts "numframes:$nf"

#REDUZCO LA CANTIDAD DE FRAMES A ANALIZAR
set nf 1 

#RECORRO TODOS LOS FRAMES
for { set i 0 } {$i < $nf } { incr i } {
  
  #TOMO EL PRIMER FRAME COMO REFERENCIA
  set sel [atomselect top all frame $i]

  #VUELCO EL PRIMER FRAME A PDB
  #$sel writepdb frame0.pdb

  #SE ELIGEN TODAS LAS AGUAS DEL FRAME
  set aguas [atomselect top "name OH2"]
  
  #SE ALMACENAN SUS COORDENADAS 
  set coord [$aguas get {x y z}]
  set id [$aguas get resid]
  #puts "**$id"
  
  set k 0

  #MUESTRO COORDENADAS X,Y,Z DE CADA AGUA
  foreach e_id $id e_coord $coord {
  
    #EXTRAIGO COORDENADAS DE LA PRIMERA PARTICULA COMO REFERENCIA
    if {$k == 0} {
      #puts "primero!"
      set x_ref [scan [lindex $e_coord 0] %f]
      set y_ref [scan [lindex $e_coord 1] %f]
      set z_ref [scan [lindex $e_coord 2] %f]
      
      lappend v1 $x_ref $y_ref $z_ref
      puts "%%$v1"
    }

    #EXTRAIGO COORDENADAS DE LA PARTICULA ACTUAL
    set x [scan [lindex $e_coord 0] %f]
    set y [scan [lindex $e_coord 1] %f]
    set z [scan [lindex $e_coord 2] %f]
    
    set dx [expr {$x_ref-$x}]
    set dy [expr {$y_ref-$y}]
    set dz [expr {$z_ref-$z}]
    
    #CALCULO DISTANCIA
    set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
    
    lappend d $dist
    lappend k_ $k
    
    #set dist 0
    set texto "fr:$i id:$k resid:$e_id elemento:$e_coord d:$dist"
    set texto2 "fr:$i id:$k resid:$e_id x:$x x_ref:$x_ref  d:$dist"
    
    puts $texto2
    grabar $texto
    incr k
  }
  }
  
  #EN ESTE PUNTO, ESTA SETEADO UNA LISTA "D" CON LAS DISTANCIAS RESPECTO A LA REFERENCIA
  #ADEMAS SE DEFINIERON MAX_R Y DR AL INICIO
  #POR LO TANTO, DENTRO DE UN FOR SE HARAN CALZAR LAS DISTANCIAS HASTA EL R_MAX
  #EN INTERVALOS [r,r+dr] ->CASCARAS 
  
  set max_r 40
  set dr 0.5
  set num [molinfo top get numatoms]
  set vol [expr {4*3.14*$max_r*$max_r*$dr}]
  set den [expr {$num/$vol}]
  
  puts "max_r:$max_r"
  puts "dr:$dr"
  
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
    set cont [expr {$cont/($num*$vol)}]
    puts "j:$a1 j+dr:$a2 cont:$cont"
    lappend x2 $a1
    lappend y2 $cont
    set cont 0
  }
  puts "no desespere"

  #BORRO LA SELECCION ANTERIOR
  foreach atom [atomselect list] {
   #puts $atom
   $atom delete
  }
  #clear
}

proc graficar {}{
  set x_ {-2 -1 0 1 2 3 4 5 6 7 8 9 10}
  set y_ {-2  0 2 3 4 5 5 4 3 2 1 0 1}
  set plothandle [multiplot -x $x2 -y $y2 -title "Funcion Distribucion Radial de probabilidades g(r) MAX_R:$max_r DR:$dr" -lines -linewidth 3 -marker point -plot]
  $plothandle configure -fillcolor yellow -radius 6
  $plothandle configure -vline {3 -width 2 -fill red -dash "."}
  $plothandle replot;
}

#FUNCION PARA GUARDAR EN TEXTO
proc grabar {dato} {
  set outfile [open output1.txt a]
  puts $outfile $dato
  close $outfile
}


