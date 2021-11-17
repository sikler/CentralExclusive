set fit quiet
set fit logfile "/dev/null"

_exp_(x) = (x < -100 ? 0 : exp(x))

f(x) = a * _exp_(-((x-m)/s)**2/2)

#
set print sprintf("../out/track/ptResol_%s.dat",p)

do for [ieta=-295:300:10] {
 do for [ipt=75:1500:50] {
  eta = ieta/100.
  pt  = ipt/1000.

  stats sprintf("<awk 'BEGIN {print 0} {if($1==%.2f && $2==%.3f && $4>0) \
                  for(i=0; i<$4; i++) print $3}' \
                  ../out/track/ptResol_%s.his", eta,pt, p) nooutput

  m = STATS_mean
  s = STATS_adev * sqrt(pi/2)

  print eta,pt,m,abs(s)
 }
 print ""

 eta = eta + 0.1
}

