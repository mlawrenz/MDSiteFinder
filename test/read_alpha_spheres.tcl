
proc read_alpha {file} {

set xcoors {}
set ycoors {}
set zcoors {}
set radii {}

set inStream [open $file r]
puts "opened"
foreach line [split [read $inStream] \n] {
    set atom [ string range $line 0 3 ]
    set atomnum [ string range $line 6 10 ]
    set atomname [ string range $line 12 15 ]
    set resname [ string range $line 16 20 ]
    set chain [ string range $line 21 21 ]
    set resnum [ string range $line 22 25 ]
    set xcoor [ string range $line 30 37 ]
    set ycoor [ string range $line 38 45 ]
    set zcoor [ string range $line 46 53 ]
    set occupancy [ string range $line 58 61 ]
    #set beta [ string range $line 67 70 ]
    #set occupancy [string range $line 54 59]
    set beta [string range $line 60 65]
    set dist [ llength $line ]
    if { [llength $line ] > 1 } { 
    lappend xcoors $xcoor
    lappend ycoors $ycoor
    lappend zcoors $zcoor
    lappend radii $beta
    }
}
close $inStream
puts "getting spheres"
set total [ llength $radii ]
for {set x 0 } { $x < $total } { incr x } {
set xcoor [ lindex $xcoors $x ]
set ycoor [ lindex $ycoors $x ]
set zcoor [ lindex $zcoors $x ]
set coor "$xcoor $ycoor $zcoor"
set rad [ lindex $radii $x ]
draw sphere $coor radius $rad
}
}

