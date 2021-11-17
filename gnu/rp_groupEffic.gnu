load "myStyleEps.gnu"

lab(a,s,r,u) = \
 sprintf("%d%s%s%s",a+1,(s==0?"n":"f"),(r==0?"T":"B"),(u==0?"u":"v"))

###############################################################################
unset xlabel; unset xtics
unset ylabel; unset ytics

set output "../eps/rp/groupEffic/dist.eps"

set term post eps enh color dashed dl 2 "Helvetica" 25 size 2*5in,2*4.5in

set multiplot layout 5,6 margins 0.08,0.98,0.08,0.98 spacing 0.005

do for [run in "319104 319124 319125 319159 319160 \
                319174 319175 319176 319177 319190 \
                319222 319223 319254 319255 319256 \
                319260 319262 319263 319264 319265 \
                319266 319267 319268 319270 319300 319311"] {

lc(a,s,r,u) = \
 (a==0 && s==0 && r==1 && u==0 ? 1 : \
 (a==0 && s==0 && r==1 && u==1 ? 2 : \
 (a==1 && s==0 && r==1 && u==0 ? 3 : \
 (a==1 && s==0 && r==1 && u==1 ? 4 : -1))))

set xrange [-0.01:1.01]
set yrange [0.5:*]

set log y; set format y "10^{%T}"
set key at graph 0.6,graph 0.8 samplen 2

set label 6 sprintf("run #%s",run) at graph 0.05,0.9 left

if(run >= 319266) {
  set xlabel "tracklet efficiency"; set xtics auto
}

if(run==319104 || run==319175 || run==319254 || run==319264 || run==319300) {
  set ylabel "counts"; set ytics auto
}

plot for [a=0:1] \
     for [s=0:1] \
     for [r=0:1] \
     for [u=0:1] \
 sprintf("<awk '{if($1>=10 && $1<512-10 && $2>-0.3 && $2<0.3) print}' \
  ../out/rp/2part/groupEffic/%s_%d%d%d_%d.his | histogram 3 -0.01 1.01 102 -", \
  run, a,s,r, u) t (lc(a,s,r,u) != -1 ? lab(a,s,r,u) : "") w his ls lc(a,s,r,u)

unset xtics
unset ytics

unset xlabel
unset ylabel
}

unset multiplot
unset log y; set format y

unset label 6

set auto x
set auto y

###############################################################################
set key top left

do for [run in "319104 319124 319125 319159 319160 \
                319174 319175 319176 319177 319190 \
                319222 319223 319254 319255 319256 \
                319260 319262 319263 319264 319265 \
                319266 319267 319268 319270 319300 319311"] {

set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,3*4.5in
set output sprintf("../eps/rp/groupEffic/%s.eps",run)

set label 6 sprintf("run #%s",run) at graph 0.98,0.8 right front

set multiplot layout 16,1 margins 0.09,0.98,0.05,0.99 spacing 0.005

#
set label 5 at graph 0.015,0.8 left front

set xrange [0:512]

set yrange [-0.3:0.3]; set ytics -0.4,0.2
set ylabel "slope" offset 1,0

#
do for [x=0+128:512-128:128] {
 set arrow from x,graph 0 rto 0,graph 1 lt 2 dt 2 nohead front
}

set colorbox origin screen 0.9, screen 0.06 \
               size screen 0.02,screen 0.09 front
set cbrange [0:1]

set style data image

#
do for [a=0:1] {
do for [s=0:1] {
do for [r=0:1] {
do for [u=0:1] {
 unset xlabel
 unset xtics

 if(a==1 && s==1 && r==1 && u==1) {
  set xlabel "tracklet location at RP center [pitch units]"; set xtics 0,128
 } 

 set label 5 lab(a,s,r,u)

 plot sprintf("../out/rp/2part/groupEffic/%s_%d%d%d_%d.his",run, a,s,r, u) 

 unset label 6
}}}}

unset multiplot

unset label 5

}
