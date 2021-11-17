load "myStyleEps.gnu"

set style data his

#
set xlabel "HLT-tracking-PID efficiency"
set ylabel "Density" offset 1.5,0
set ytics 0,1e-2

do for [p in "pi ka pr"] {

#
set output sprintf("../eps/track/comb_efficiency_%s.eps",p)

print " distrts ",p

plot [0-1e-3:1+1e-3][0:4e-2] \
 sprintf("<zcat ../out/track/cmEff/%s_2_2_TB.his.gz | \
           histogram 5 0. 1. 100 -",p) u 1:($2/$3) t "TB" ls 1, \
 sprintf("<zcat ../out/track/cmEff/%s_2_2_BT.his.gz | \
           histogram 5 0. 1. 100 -",p) u 1:($2/$3) t "BT" ls 2, \
 sprintf("<zcat ../out/track/cmEff/%s_2_2_TT.his.gz | \
           histogram 5 0. 1. 100 -",p) u 1:($2/$3) t "TT" ls 3, \
 sprintf("<zcat ../out/track/cmEff/%s_2_2_BB.his.gz | \
           histogram 5 0. 1. 100 -",p) u 1:($2/$3) t "BB" ls 4
}
set ytics auto

#
set colorbox vertical user origin 0.2,0.2 size 0.03,0.3 front

set xlabel "{/Symbol f}"
set xrange [0:pi]

set ylabel "m [GeV]"

set log cb ; set format cb "10^{%T}"
#set cbrange [1e-2:1] # FIXME VERY
set cbrange [1e-3:1] # FIXME VERY

set label 12 at graph 0.1,0.9 left font "Helvetica,40" front

do for [p in "pi ka pr"] {

print " min/max ",p

 if(p eq "pi") { set yrange [0.0:2.5] ; set label 12 "{/Symbol p}"}
 if(p eq "ka") { set yrange [0.5:3.0] ; set label 12 "K"}
 if(p eq "pr") { set yrange [1.0:3.5] ; set label 12 "p"}

# do for [b in "1_1 2_4 4_7"] {
#   if(b eq "1_1") { ext = "a" ; i=1; j=1; topo="TB" }
#   if(b eq "2_4") { ext = "b" ; i=2; j=4; topo="BT" }
#   if(b eq "4_7") { ext = "c" ; i=4; j=7; topo="TT" }
  do for [topo in "TB BT TT BB"] {
# FIXME VERY
#  i = 4; j = 6; b = "4_6"
   i = 2; j = 2; b = "2_2"
   if(topo eq "TB") { ext = "a" }
   if(topo eq "BT") { ext = "b" }
   if(topo eq "TT") { ext = "c" }
   if(topo eq "BB") { ext = "d" }

   #
   set output sprintf("../eps/track/comb_efficiency_min_%s_%s.eps",p,ext)
   plot sprintf("<zcat ../out/track/cmEff/%s_%s_%s.his.gz | \
         ./cm_combined.awk",p,b,topo) u 1:2:3 w image

   if(p eq "pr") {
    set label 13 sprintf("%.2f<p_{1,T}<%.2f GeV\n%.2f<p_{1,T}<%.2f GeV", \
        0.2+i*0.05,0.2+(i+1)*0.05, \
        0.2+j*0.05,0.2+(j+1)*0.05) at graph 0.9,0.15 right front
    set label 14 topo at graph 0.9,0.9 right front
   }

   #
   set output sprintf("../eps/track/comb_efficiency_max_%s_%s.eps",p,ext)
   plot sprintf("<zcat ../out/track/cmEff/%s_%s_%s.his.gz | \
         ./cm_combined.awk",p,b,topo) u 1:2:4 w image

   if(p eq "pr") {
    unset label 13
    unset label 14
   }
  }
}

