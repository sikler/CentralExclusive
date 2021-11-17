#!/usr/bin/awk -f

BEGIN {
 o_phim=-1; o_phi=-1; min=1e+99; max=-1e+99
} {
 if(NF > 0)
 {
  phim = $1" "$2
  phi  = $1

  if(o_phim != -1 && o_phim != phim)
  {
   print o_phim, min, max
   min=1e+99; max = -1e+99

   if(ophi != -1 && o_phi != phi)
   {
     print ""
   }
  }

  if($5>max) max=$5
  if($5<min) min=$5

  o_phim = phim
  o_phi  = phi
 }
} END {
   print o_phim, min, max
}

