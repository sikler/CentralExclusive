load "myStyleEps.gnu"

parts = "pi ka pr"
names = "'{/Symbol p}' K p"

set style data image

set xlabel "p_T [GeV]"
set ylabel "{/Symbol s}_{ln({/Symbol e}/[MeV/cm])}" offset 2,0

set cbtics offset -1,0

set label 5 at graph 0.5,0.95 center front

do for [i=1:words(parts)] {

 p = word(parts,i)
 t = word(names,i)

 if(p eq "pr") { set cbtics 0,0.1 }

 set label 1 t at graph 0.9,0.9 front right font "Helvetica,40"

 n = 30
 set label 5 "0.0 < {/Symbol h} < 0.1"
 set output sprintf("../eps/track/relsig%s_a.eps",p)
 plot sprintf("../out/track/relSigma_%s.his",p) \
         i n u 2:3:($4*1e-2)

 n = 35
 set label 5 "0.5 < {/Symbol h} < 0.6"
 set output sprintf("../eps/track/relsig%s_b.eps",p)
 plot "" i n u 2:3:($4*1e-2)

 n = 40
 set label 5 "1.0 < {/Symbol h} < 1.1"
 set output sprintf("../eps/track/relsig%s_c.eps",p)
 plot "" i n u 2:3:($4*1e-2)

 n = 45
 set label 5 "1.5 < {/Symbol h} < 1.6"
 set output sprintf("../eps/track/relsig%s_d.eps",p)
 plot "" i n u 2:3:($4*1e-2)

 n = 50
 set label 5 "2.0 < {/Symbol h} < 2.1"
 set output sprintf("../eps/track/relsig%s_e.eps",p)
 plot "" i n u 2:3:($4*1e-2)

 n = 55
 set label 5 "2.5 < {/Symbol h} < 2.6"
 set output sprintf("../eps/track/relsig%s_f.eps",p)
 plot "" i n u 2:3:($4*1e-2)

 set cbtics auto
}

