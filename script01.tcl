set nf [molinfo top get numframes]

puts "numframes:$nf"

set k 0
set outfile [open output.txt a]
for { set i 0 } {$i < $nf } { incr i } {
  set sel [atomselect top all frame $i]
  #$sel writepdb $i.pdb
  clear
  set numatoms [molinfo top get numatoms]
  for { set j 0 } {$j < $numatoms } { incr j } { 
    set coordenadas [atomselect top "index $j"]
    $coordenadas get {x y z}
    set ref [atomselect top "index 0"]
    
    set x_ref [$ref get {x}]
    set y_ref [$ref get {y}]
    set z_ref [$ref get {z}]
    
    set x [$coordenadas get {x}]
    set y [$coordenadas get {y}]
    set z [$coordenadas get {z}]

    set dx [expr {$x_ref-$x}]
    set dy [expr {$y_ref-$y}]
    set dz [expr {$z_ref-$z}]
    
    set dist [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
 
    
    puts $outfile "atom:$j >x:$x y:$y z:$z -> dist:$dist"

    puts "atom:$j >x:$x y:$y z:$z -> dist:$dist"
    
    #if($k>10){
      #clear
      #set k 0
      #puts "=)"
    #}
    #else{
      #incr $k
      #puts "&"
    #}


  }
} 
close $outfile


#for {set x 0} {$x < 10} {incr x} {
#  for {set y 0} {$y < 10} {incr y} {
#    puts "x:$x y:$y"
#  }
#}

#set outfile [open output.dat w]
#puts $outfile [measure bond {143 387} first 0 last 42]
#close $outfile