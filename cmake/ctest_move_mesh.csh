#!/bin/csh -f

mkdir -p $1/Diags
mkdir -p $1/Mesh
mkdir -p $1/Info
mkdir -p $1/Data
mv $1/meshdb* $1/Mesh

