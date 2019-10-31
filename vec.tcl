#MANEJO DE VECTORES TCL

set d [veczero] 

set v1 {1 2 3}
set v2 {4 5 6}
puts "v1: {$v1} v2: {$v2}"


set v3 [vecadd $v1 $v2]
puts "vecadd: {$v1} + {$v2} = {$v3}"

set v4 [vecsub $v1 $v2]
puts "vecsub: {$v1} - {$v2} = {$v4}"

set dot [vecdot $v1 $v2] 
puts "vecdot: {$v1} . {$v2} = {$dot}"

set cross [veccross $v1 $v2] 
puts "veccross: {$v1} x {$v2} = {$cross}"

set length [veclength $v1] 
puts "veclength: L{$v1} = {$length}"

set dist [vecdist $v1 $v2] 
puts "vecdist: D{$v1 $v2} = {$dist}"

set norm [vecnorm $v1] 
puts "vecnorm: N{$v1} = {$norm}"

set vinvert [vecinvert $v1] 
puts "vecinvert: I{$v1} = {$vinvert}"
