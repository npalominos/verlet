set nf [molinfo top get numframes]
puts "numframes:$nf"

set k 0
#set nf 10

#vector distancias
#set d [veczero] 

for { set i 0 } {$i < $nf } { incr i } {
  set sel [atomselect top all frame $i]
  #$sel writepdb $i.pdb
  #clear
  set numatoms [molinfo top get numatoms]
  
  set ref [atomselect top "index 0"]
  set x_ref [$ref get {x}]
  set y_ref [$ref get {y}]
  set z_ref [$ref get {z}]  
  
  #set numatoms 10
  set outfile [open output.txt a]
  for { set j 0 } {$j < $numatoms } { incr j } { 
  #foreach atom [atomselect top "name OH2"] {
  
    set coordenadas [atomselect top "index $j"]
    $coordenadas get {x y z}
    
    #geom_center $coordenadas
    
    set x [$coordenadas get {x}]
    set y [$coordenadas get {y}]
    set z [$coordenadas get {z}]

    set dx [expr {$x_ref-$x}]
    set dy [expr {$y_ref-$y}]
    set dz [expr {$z_ref-$z}]
    
    set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
    
    puts $outfile "fr:$i at:$j >x:$x y:$y z:$z -> d:$dist"
    puts "fr:$i at:$j >x:$x y:$y z:$z -> d:$dist"

    $coordenadas delete
    
    unset x
    unset y
    unset z
    
    unset dx
    unset dy
    unset dz
    
    unset dist
    
  }
  $sel delete
  unset x_ref
  unset y_ref
  unset z_ref
  unset ref

} 
close $outfile

puts "OK"

foreach atom [atomselect list] {
   puts $atom
   $atom delete
}

atomselect list
#puts info vars

proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}