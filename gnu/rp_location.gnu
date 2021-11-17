load "myStyleEps.gnu"

unset colorbox

lab(a,s,r) = sprintf("%d%s%s",a+1,s==0?"n":"f",r==0?"T":"B")

set xtics -2,1
set ytics -40,20

set tics scale 1

#
set output "../eps/rp/rp_location.eps"
set style data image

set term post eps enh color dashed dl 2 "Helvetica" 15 size 0.9*5in,2*4.5in
set multiplot layout 8,4 margins 0.08,0.98,0.08,0.98 spacing 0.008

set size square

set object 5 rectangle from screen 0.525,screen 0 rto 0.5,1 \
  fillc rgb "yellow" fs solid 0.1 noborder behind

do for [arm=0:1] {
do for [sta=0:1] {
do for [rpt=0:1] {
 if(arm==1 && sta==1 && rpt==1) { \
  set xlabel "x [mm]" offset 0,0.5; set xtics } \
  else { unset xlabel; unset xtics }

 set label 7 lab(arm,sta,rpt) at graph 0.05,0.55 left front

 do for [n in "2 0"] {
 do for [topo in "TB BT TT BB"] {

  if(n==2 && topo eq "TB") { set ylabel "y [mm]" offset 1,0; \
                             set ytics nomirror }

  t = substr(topo,arm+1,arm+1)
  ok = (t eq "T" && rpt==0) || (t eq "B" && rpt==1)

  if(ok) {
    set label 6 topo at graph 0.95,0.9 right front

    y_low =  7
    y_hig = 24

    if(sta == 1) {
      y_low = y_low * 1.14
      y_hig = y_hig * 1.14
    }

    if(rpt == 1) {
      y_low = -y_low
      y_hig = -y_hig
    }

    plot \
      sprintf("../out/rp/%spart/hits/loc_%d%d%d_%s.his",n, arm,sta,rpt, topo), \
      0 lt 0, \
      y_low lt 2 dt 3, \
      y_hig lt 2 dt 3
    

    unset ylabel; unset ytics
    unset object 5
    unset label 7
  }

 }}
}}}

unset multiplot

###############################################################################
f(x) = a*exp(-((x-m)/s)**2/2)
w(x) = (x > 1 ? sqrt(x) : 1)

do for [dir in "x y"] {
set output sprintf("../eps/rp/rp_location_%s.eps",dir)
set style data his

if(dir eq "x") {
 set xrange [-1:1]
} else {
 set auto x
 set log y
}

set multiplot layout 8,4 margins 0.08,0.98,0.08,0.98 spacing 0.008

set size square

set object 5 rectangle from screen 0.525,screen 0 rto 0.5,1 \
  fillc rgb "yellow" fs solid 0.1 noborder behind

do for [arm=0:1] {
do for [sta=0:1] {
do for [rpt=0:1] {
 if(arm==1 && sta==1 && rpt==1) { \
   set xlabel sprintf("%s [mm]",dir) offset 0,0.5;
   if(dir eq "x") { set xtics auto -1,0.5 } else { set xtics -40,20 } } \
 else { unset xlabel; unset xtics }

 set label 7 lab(arm,sta,rpt) at graph -0.20,0.5 left front

 do for [n in "2 0"] {

# if(n == 2) { set yrange [1e+3:1e+5] } else { set yrange [1e+2:1e+5] }

 if(n == 0) {
  set object 6 rectangle from graph 0,graph 0 rto 1,1 \
   fillc rgb "white" fs solid noborder behind
 }

 do for [topo in "TB BT TT BB"] {
  
  t = substr(topo,arm+1,arm+1)
  ok = (t eq "T" && rpt==0) || (t eq "B" && rpt==1)

  diag = (topo eq "TB") || (topo eq "BT")
  para = (topo eq "TT") || (topo eq "BB")

  oppo = (topo eq "TB" ? "BT" : \
         (topo eq "BT" ? "TB" : \
         (topo eq "TT" ? "BB" : "TT")))

  if(ok) {
   set label 6 topo at graph 0.95,0.9 right front

   if(dir eq "x") {
     a = 1e+4; m = 1e-3; s = 0.2
     mx = 0.4

     fit [-mx:mx] f(x) sprintf("../out/rp/%spart/hits/loc_%s_%d%d%d_%s.his",\
                      n,dir, arm,sta,rpt, topo) u 1:2:(w($2)) zerror via a
     fit [-mx:mx] f(x) "" u 1:2:(w($2)) zerror via a,m,s

     set label 8 at graph 0.05,0.7 left front
     if(n == 0 && diag) {
      set label 8 sprintf("%.0f {/Symbol m}m",m*1e+3)   tc rgb "black"
     } else {
      set label 8 sprintf("(%.0f {/Symbol m}m)",m*1e+3) tc rgb "dark-gray"
     }

     #
     plot "" u 1:2 ls 1, (abs(x) < mx ? f(x) : 1/0) ls 3 dt 1 lw 1
     unset label 8
   } else {
     y_low =  7
     y_hig = 24

     if(sta == 1) {
       y_low = y_low * 1.14
       y_hig = y_hig * 1.14
     }

     if(rpt == 1) {
       y_low = -y_low
       y_hig = -y_hig
     }

     set arrow from y_low,graph 0 rto 0,graph 1 lt 2 dt 3 nohead
     set arrow from y_hig,graph 0 rto 0,graph 1 lt 2 dt 3 nohead

     plot \
       sprintf("../out/rp/%spart/hits/loc_%s_%d%d%d_%s.his", \
               n,dir, arm,sta,  rpt, topo) u 1:2 ls 1
#      sprintf("../out/rp/%spart/hits/loc_%s_%d%d%d_%s.his", \
#              n,dir, arm,sta,1-rpt, oppo) u 1:2 ls 3

     unset arrow
   }

   unset ylabel; unset ytics
   unset object 5
   unset label 7
  }

 }}
 unset object 6
}}}

unset multiplot
}
