load "myStyleEps.gnu"

set xtics

dat(p) = \
 sprintf("<grep '%d|%d|%d|u%d|p%d|r%d' .r%d.dat", a,s,r,u, p,x, x)

set xrange [-1:26]
#set xlabel "run number"

set yrange [0:1]
# set ylabel "strip efficiency"

set style data e

#
set term post eps enh color dashed dl 2 "Helvetica" 30 size 2*5in,3*4.5in

# FIXME 2part
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r200" > .r200.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r225" > .r225.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r250" > .r250.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r275" > .r275.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r300" > .r300.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r325" > .r325.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r350" > .r350.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r375" > .r375.dat
! zcat ../out/rp/2part/stripEffic/319*.dat.gz | grep "r400" > .r400.dat

set label 8 at graph 0.9,graph 0.1 right front

do for [x=200:400:25] {

set key at graph 0.7,graph 0.375

unset xtics
unset ytics

set label 7 sprintf("strip #%d",x) at screen 0.06,screen 0.975 left front

set output sprintf("../eps/rp/stripEffic/%d_vs_run.eps",x)
set multiplot layout 4,4 margins 0.05,0.99,0.05,0.99 spacing 0.01

lab(a,s,r,u) = sprintf("%d%s%s%s",a+1,s==0?"n":"f",r==0?"T":"B",u==0?"u":"v")

do for [a=0:1] {
do for [s=0:1] {

if(a==1 && s==1) {
set xtics ( \
 "#319104"  0, "#319124"  1, "#319125"  2, "#319159"  3, "#319160"  4, \
 "#319174"  5, "#319175"  6, "#319176"  7, "#319177"  8, "#319190"  9, \
 "#319222" 10, "#319223" 11, "#319254" 12, "#319255" 13, "#319256" 14, \
 "#319260" 15, "#319262" 16, "#319263" 17, "#319264" 18, "#319265" 19, \
 "#319266" 20, "#319267" 21, "#319268" 22, "#319270" 23, "#319300" 24, \
 "#319311" 25) font "Helvetica,20" rotate by 90 offset 0,-2.25
}

do for [r=0:1] {
do for [u=0:1] {

unset ytics
unset mytics

if(r==0 && u==0) {
 set ytics
} else {
 set ytics ("" 0, "" 0.2, "" 0.4, "" 0.6, "" 0.8, "" 1)
}

 set label 8 lab(a,s,r,u)
# for [p=0:4] dat(0) u 0:2:3 \
#     w errorlines t sprintf("pla #%d",p+1) lt (p+1)
 plot \
  for [p=0:4] dat(p) u 0:2:3 \
      w e t sprintf("plane #%d",p+1) lt (p+1), \
  for [p=0:4] dat(p) u 0:2:($2>1e-2 ? 1e-1 : 1e-6) sm acs w l lt (p+1)

 unset label 7
 unset key

}}}}

unset multiplot

}

! rm .r*.dat
