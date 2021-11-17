load "myStyleEps.gnu"

set style data image

set colorbox origin 0.75,0.6
set cbtics offset -0.5,0

set label 7 at graph 0.95,0.9 right front

#
set xtics -0.5,0.1
set ytics -0.5,0.1

set size ratio -1

set cbrange [0:1]

do for [topo in "TB BT"] {
 array tb[2]
 tb[1] = (topo eq "TB" ? "T" : "B")
 tb[2] = (topo eq "BT" ? "T" : "B")

 set label 7 sprintf("%s",topo)

 set xlabel sprintf("p_{1%s,y} [GeV]",tb[1])
 set ylabel sprintf("p_{2%s,y} [GeV]",tb[2])

 if(topo eq "TB") { xlo = 0.1; xhi = 0.5 } else { xlo = -0.5; xhi = -0.1 }
 if(topo eq "BT") { ylo = 0.1; yhi = 0.5 } else { ylo = -0.5; yhi = -0.1 }

 set xtics -0.5,0.1
 set ytics -0.5,0.1

 if(topo eq "TB") { set label 7 at graph 0.05,0.9 left front }
 if(topo eq "BT") { set label 7 at graph 0.05,0.1 left front }

 set parametric
 set trange [-0.5:0.5]
 set xrange [xlo:xhi]
 set yrange [ylo:yhi]

 w = 0.07
 d = 0.01

 if(topo eq "TB") {
   do for [i=2:5] {
     set object i+1 rectangle \
         from  0.015-d+i*w,-0.025+d-i*w rto +w+2*d,-w-2*d front \
         fs empty border @green lw 3 dt 3
   }
 } else {
   do for [i=2:5] {
     set object i+1 rectangle \
         from -0.025+d-i*w, 0.015-d+i*w rto -w-2*d,+w+2*d front \
         fs empty border @green lw 3 dt 3
   }
 }

 set output sprintf("../eps/rp/veto_eff_py_%s.eps",topo)
 plot \
   sprintf("../out/rp/2part/veto_py.his") u 1:2:3
#   -0.19,t ls 2, t,-0.19 ls 2, \
#    0.65,t ls 2, t, 0.65 ls 2
}

