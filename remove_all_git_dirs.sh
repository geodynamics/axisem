#!/bin/bash
find . | grep .git | xargs rm -rf
rm .gitignore
rm .gitmodules

