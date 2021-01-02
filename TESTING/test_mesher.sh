#!/bin/bash
./copytemplates.sh $buildstyle
set -e
cd MESHER
mkdir -p Diags
make -sj 4
FILES=$(find Models/ -type f -name '*.bm')
status='normal'
for file in $FILES; do 
  # Escape slashes in file name
  file=$(echo "$file" | sed 's/\//\\\//g'); 
  echo ""
  echo "****************************************************************"
  echo "Testing external model in $file"
  echo "****************************************************************"
  echo ""
  # Modify inparam_mesh to use this file
  sed -e "s/EXT_TEMPLATE/'$file'/" ../TESTING/inparam_mesh_external_template > inparam_mesh
  set +e
  ./xmesh &> OUTPUT 
  if [ $? == 0 ]; then
    tail -n 1 OUTPUT 
  else
    tail -n 20 OUTPUT 
    status='error'
  fi
  set -e
  echo "****************************************************************"
done
if [ $status == 'error' ]; then
  exit 1
else
  exit 0
fi
