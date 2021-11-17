load "myStyleEps.gnu"

set lmargin
set rmargin
set tmargin
set bmargin

set tics scale 1,0.5 front

set term post eps enh color dashed dl 2 "Helvetica" 20 size 6in,2*4.5in

set style data xe

dat(var,fac,npars) = \
 sprintf("<./tunes.awk s=\"%s\" npars=%d ../tune/tunes_%s/results.txt",var,npars,fac)

dim(var,nsoft) = \
 sprintf("<awk '{if($1==\"%s\") print $%d}' soft.dat",var,nsoft+1)


unset ytics

set output "../eps/tunes.eps"

set multiplot layout 17,1

set lmargin 9
set tmargin 0
set bmargin 1.75
set rmargin 3+10

set xtics offset 0,0.4

set datafile missing "0"

set label 9 at graph -0.17,first 0 left
unset key

set border 1
set xtics nomirror

do for [var in "bjac bexp aoo ao bo sigo ep asp gd2 pp0(1) bb0(1) bb0(2) cc0(1) cc0(2) bex(1) bex(2) cfsi"] {

  #
  name = var

  if(name eq "cfsi") { name = "c_{sec}" }

  if(name eq "bexp") { name = "b_{exp}" }
  if(name eq "aoo" ) { name = "b_{pow}" }
  if(name eq "ao"  ) { name = "a_{or}" }
  if(name eq "bo"  ) { name = "b_{or}" }

# if(name eq "pp0(1)") name = "$2 |a_1|^2$"

  if(name eq "bb0(1)") { name = "c_1 [GeV^2]" }
  if(name eq "bb0(2)") { name = "c_2 [GeV^2]" }

  if(name eq "cc0(1)") { name = "d_1" }
  if(name eq "cc0(2)") { name = "d_2" }

  if(name eq "bex(1)") { name = "b_1 [GeV^2]" }
  if(name eq "bex(2)") { name = "b_2 [GeV^2]" }

  if(name eq "sigo") { name = "{/Symbol s}_P [mb]" }
  if(name eq "asp")  { name = "{/Symbol a}@_P' [GeV^{-2}]" }
  if(name eq "ep")   { name = "{/Symbol a}_P" }

  if(name eq "bjac") { name = "B_P [GeV^{-2}]" }

  set label 9 name

  #
  set xtics auto

  if(var eq "bjac") { set xrange [5.6:6.4] }
  if(var eq "ao")   { set xrange [0.8:1.0] ; set xtics 0.8,0.05 }
  if(var eq "cfsi") { set xtics 0.2,0.2 }
  if(var eq "ep")   { set xtics 0.08,0.02 }

  if(strstrt(var,"bb0(")>0) { set xrange [0.0:0.5] }
  if(strstrt(var,"cc0(")>0) { set xrange [0.4:0.6] ; set xtics 0.4,0.05 }
  if(strstrt(var,"bex(")>0) { set xrange [2:10] }

  if(var eq "cfsi") {
    set key at screen 0.99,0.5 horizontal maxcols 1
  }

  plot [][-3:3] \
    dim(var,1) u 1:(0):(1.5) w ye lt 1 lc 9 pt 0 ps 0 lw 10, \
    dim(var,2) u 1:(0):(1.5) w ye lt 1 lc 9 pt 0 ps 0 lw 10, \
    dim(var,3) u 1:(0):(1.5) w ye lt 1 lc 9 pt 0 ps 0 lw 10, \
    dim(var,4) u 1:(0):(1.5) w ye lt 1 lc 9 pt 0 ps 0 lw 10, \
    dat(var,"exp"  ,14) u 2:( 1):3 t "exp" ls 3 ps 2, \
    dat(var,"power",14) u 2:( 0):3 t "power" ls 1 ps 2, \
    dat(var,"orear",15) u 2:(-1):3 t "Orear" ls 2 ps 2

  unset key

  set auto x
}

unset multiplot

