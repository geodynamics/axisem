#!/bin/bash
set -e
./copytemplates.sh
cd MESHER
mkdir Diags
make -sj 4
FILES=$(find Models/ -type f -name '*.bm*')
for file in $FILES; do 
  # Escape slashes in file name
  file=$(echo "$file" | sed 's/\//\\\//g'); 
  echo ""
  echo "****************************************************************"
  echo "Testing external model in $file"
  echo "****************************************************************"
  echo ""
  # Modify inparam_mesh to use this file
  sed -e "s/EXT_TEMPLATE/'$file'/" ../TRAVIS/inparam_mesh_external_template > inparam_mesh
  ./xmesh > OUTPUT 
  tail -n 15 OUTPUT
  echo "****************************************************************"
done
