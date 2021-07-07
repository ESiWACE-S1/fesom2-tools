#!/usr/bin/env bash

set -eu -o pipefail

command -v m2gmetis # Converts a mesh into a graph that is compatible with METIS.
command -v mpmetis # Partitions a mesh into a specified number of parts
command -v ndmetis # fill-reducing ordering of the vertices of the graph using nested dissection

test -e elem2d.out.bak || cp elem2d.out elem2d.out.bak
test -e nod2d.out.bak  || cp nod2d.out nod2d.out.bak

### Prepare metis files and reorder permuttation for fill-reduction
m2gmetis              elem2d.out.bak elem2d.e.graph
ndmetis elem2d.e.graph   # -> elem2d.e.graph.iperm
#m2gmetis -gtype=nodal elem2d.out.bak elem2d.n.graph
#ndmetis elem2d.n.graph   # -> elem2d.n.graph.iperm

#gpmetis elem2d.e.graph 2 # -> elem2d.e.graph.part.2

npoints=$(head -n1 nod2d.out | awk '{print $1}')
ncells=$(head -n1 elem2d.out | awk '{print $1}')

# https://www.datafix.com.au/BASHing/2020-11-11.html
# awk -F"\t" 'FNR==NR {a[$1]=$2; next} FNR>1 {$3=a[$3]} 1' OFS="\t" lookup fileA

## renumber elements in elem2d.out
echo "$ncells" > elem2d.out
awk '
FNR == NR {
mapping[FNR]=$1
next
}
FNR > 1 {
print mapping[FNR-1]" "$0
}' elem2d.e.graph.iperm elem2d.out.bak | sort -n -k 1 | awk '{ print $2" "$3" "$4 }' >> elem2d.out

## renumber levels in elvls.out
test -e elvls.orig || cp elvls.out elvls.orig
awk '
FNR == NR {
mapping[FNR]=$1
next
}
FNR >= 1 {
print mapping[FNR]" "$0
}' elem2d.e.graph.iperm elvls.orig | sort -n -k 1 | awk '{ print $2" "$3" "$4 }' > elvls.out

## construct elem2d.n.iperm
awk '
BEGIN {
  for(i=1; i<='$npoints'; i++) {
	  mapping[i]=-1
  }
  points_id=0
}
NR > 1 {
  if(mapping[$1] == -1){
	  points_id++
	  mapping[$1]=points_id
  }
  if(mapping[$2] == -1){
	  points_id++
	  mapping[$2]=points_id
  }
  if(mapping[$3] == -1){
	  points_id++
	  mapping[$3]=points_id
  }
}
END {
  for(i=1; i<='$npoints'; i++) {
	  print mapping[i]
  }
}' elem2d.out > elem2d.n.iperm

# renumber nodes in nod2d.out
echo "$npoints" > nod2d.out
awk '
FNR == NR {
mapping[FNR]=$1
next
}
FNR > 1 {
print mapping[$1]" "$2" "$3" "$4
}' elem2d.n.iperm nod2d.out.bak | sort -n -k 1 >> nod2d.out

## renumber levels in nlvls.out
test -e nlvls.orig || cp nlvls.out nlvls.orig
awk '
FNR == NR {
mapping[FNR]=$1
next
}
FNR >= 1 {
print mapping[FNR]" "$0
}' elem2d.n.iperm nlvls.orig | sort -n -k 1 | awk '{ print $2" "$3" "$4 }' > nlvls.out

# renumber nodes in elements
cp elem2d.out elem2d.out.old
echo "$ncells" > elem2d.out
awk '
FNR == NR {
mapping[FNR]=$1
next
}
FNR > 1 {
print mapping[$1]" "mapping[$2]" "mapping[$3]
}' elem2d.n.iperm elem2d.out.old >> elem2d.out


