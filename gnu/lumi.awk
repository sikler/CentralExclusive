#!/usr/bin/awk -f

function bad(pred,meas)
{
  return (meas < pred - 3*sqrt(pred) || meas <= 1)
} BEGIN {
  bprev = 0
} {
  if(NR==1 || NFprev==0) 
  {
    print sprintf("set label \" %s\" at \"%s %s\", \
          graph 0.8 right front rotate by 90; ",$10,$12,$13);
  }
  NFprev = NF

  if(NF>0)
  {
    if(topo == "TB") b = bad($1,$2)
    if(topo == "BT") b = bad($3,$4)
    if(topo == "TT") b = bad($5,$6)
    if(topo == "BB") b = bad($7,$8)
  }
  else b = 0

  if(b)
  { 
    if(!bprev) { fromls=$11; fromdate=$12; fromtime=$13 }
                 lastls=$11
  }
  else
  {  
    if(bprev)
    {
      if(fromls == lastls)
        print sprintf("set label \"%d\" at \"%s %s\", \
              graph 0.1 rotate by 90 front font \"Helvetica,10\"; ",\
              fromls,       fromdate,fromtime)
      else
        print sprintf("set label \"%d-%d\" at \"%s %s\", \
              graph 0.1 rotate by 90 front font \"Helvetica,10\"; ",\
              fromls,lastls,fromdate,fromtime)
    }
  }

  bprev = b
}

