load "myStyleEps.gnu"

g(x) = a*erf(x/sqrt(2)/s)
f(x) = g(x-m+0.01)-g(x-m-0.01)

dat(mas,b,p) = sprintf("<awk 'function abs(x) { return (x>0 ? x : -x) } \
 {if($1==%.2f) print}' ../out/track/massResol/%s_%s.dat",mas,b,p)

set label 12 at graph 0.1,0.9 center

do for [p in "pi ka pr"] {

 if(p eq "pi") { mlo =  29 ; m_ = mpi; set label 12 "{/Symbol p}"}
 if(p eq "ka") { mlo =  99 ; m_ = mka; set label 12 "K" }
 if(p eq "pr") { mlo = 191 ; m_ = mpr; set label 12 "p" }

 do for [b in "1_1 2_4 4_7"] {

  if(b eq "1_1") { ext = "a" ; i=1; j=1 }
  if(b eq "2_4") { ext = "b" ; i=2; j=4 }
  if(b eq "4_7") { ext = "c" ; i=4; j=7 }

  a = 1; m = -2e-3; s = 10e-3

  set print "a.out"
  do for [im=mlo:250:2] {
#  a = 1; m = 0.001; s = 0.01
   mas = im*1e-2
   fit f(x) dat(mas,b,p) u 2:3 via a
   fit f(x) ""           u 2:3 via a,m,s

   print mas, m,m_err, abs(s),s_err
  }

  set print

  #
  v(E) = sqrt((E/2)**2 - m_**2)

  set xlabel "p* [GeV]"
  set style data e

  set output sprintf("../eps/track/massResol_del_%s_%s.eps",p,ext) 
  set ylabel  "{/Symbol D}m [MeV]"
  set yrange [-10:10]
  plot [0:1.5] "a.out" u (v($1)):($2*1e+3):($3*1e+3) ls 1, 0 lt 0
  set auto y

  if(p eq "pr") {
    set label 13 sprintf("%.2f<p_{1,T}<%.2f GeV\n%.2f<p_{1,T}<%.2f GeV", \
        0.2+i*0.05,0.2+(i+1)*0.05, \
        0.2+j*0.05,0.2+(j+1)*0.05) at graph 0.9,0.15 right front
  }

  set output sprintf("../eps/track/massResol_sig_%s_%s.eps",p,ext) 
  set ylabel "{/Symbol s}_m [MeV]"
  set yrange [0:35]
  plot [0:1.5]      "" u (v($1)):($4*1e+3):($5*1e+3) ls 3 dt 1
  set auto y

  unset label 13
 }
}

! rm a.out

