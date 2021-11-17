#!/usr/bin/awk -f

{
 if(NF==1) d=$1;
 else
 {
  n=$6

  for(i=1; i<=5; i++)
  if($i != -99)
  {
   if($i != int($i))
     a[d][i]+=n

   b[d][i]+=n
  }
 }
} END {
 for(x in a) 
 {
  printf(" %s",x)

  for(i=1; i<=5; i++)
  {
    num = a[x][i]+0 
    sum = b[x][i]+0 

    if(sum > 0) printf(" %f",num/sum);
           else printf(" 0");
  }

  printf("\n")
 }
}

