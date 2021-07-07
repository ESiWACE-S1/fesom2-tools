#!/usr/bin/env bash

set -eu -o pipefail

npoints=$(head -n1 nod2d.out)
ncells=$(head -n1 elem2d.out)
totalpointsincells=$(( $ncells * (1+3) ))
VTK_TRIANGLE=5

# https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
# flat lat/lon
#$(cat nod2d.out | tail -n+2 | sed "s/ [0-9]\+$//" | awk '{ $1=$1-1; print $2" "$3" 0" }')
cat > mesh.vtk <<EOF
# vtk DataFile Version 2.0
Really cool data
ASCII
DATASET UNSTRUCTURED_GRID
POINTS $npoints float
$(cat nod2d.out | tail -n+2 | sed "s/ [0-9]\+$//" | awk '
BEGIN {
	R=6371;
	pi = atan2(0, -1);
	rad = pi/180
}
{
	$1=$1-1;
	lon=$2*rad;
	lat=$3*rad;
	x = R*cos(lat)*cos(lon);
	y = R*cos(lat)*sin(lon);
	z = R*sin(lat);
	print x" "y" "z
}')
CELLS $ncells $totalpointsincells
$(cat elem2d.out | tail -n+2 | awk '{$1=$1-1; $2=$2-1; $3=$3-1; print "3 "$0}')
CELL_TYPES $ncells
$(for i in $(seq 1 $ncells); do echo $VTK_TRIANGLE; done)
EOF

### POINT data
(
echo "POINT_DATA $npoints"
echo "SCALARS id int"
echo "LOOKUP_TABLE default"
for i in $(seq 1 $npoints); do echo $i; done
if test -e elem2d.n.iperm; then
echo "SCALARS i_perm int"
echo "LOOKUP_TABLE default"
cat elem2d.n.iperm
fi
if test -e nlvls.out; then
echo "SCALARS nlvls int"
echo "LOOKUP_TABLE default"
cat nlvls.out
fi
# put partition data if they exist
for dir in $(ls --color=no -d dist_* || true); do
echo "SCALARS $dir int"
echo "LOOKUP_TABLE default"
awk '
BEGIN {
  npart=0
  part_id=0
  partitions[0]=0
  partitions_read=0
}
NR == 1 {
  npart=$1
}
partitions_read == 1 {
  entry_id=$1-1
  part_id=npart
  while(entry_id < partitions[part_id] && part_id >= 0) {
	  part_id--
  }
  print part_id
}
(NR > 1) && (partitions_read == 0) {
  for ( i=1; i <= NF; i++) {
    part_id++
    partitions[part_id]=$i
  }
  if (part_id >= npart) {
	  partitions_read=1
	  part_id = 0;
	  for ( i=1; i <= npart ; i++) {
	    partitions[i]=partitions[i-1]+partitions[i]
	  }
  }
  
}
END {
  #print "npart="npart
  #for ( i=0; i <= npart; i++) {
  #  print partitions[i]
  #}
}
' $dir/rpart.out
done
) >> mesh.vtk

### CELL data
(
echo "CELL_DATA $ncells"
echo "SCALARS id int"
echo "LOOKUP_TABLE default"
for i in $(seq 1 $ncells); do echo $i; done
if test -e elem2d.e.graph.iperm; then
echo "SCALARS i_perm int"
echo "LOOKUP_TABLE default"
cat elem2d.e.graph.iperm
fi

if test -e elvls.out; then
echo "SCALARS elvls int"
echo "LOOKUP_TABLE default"
cat elvls.out
fi
) >> mesh.vtk

