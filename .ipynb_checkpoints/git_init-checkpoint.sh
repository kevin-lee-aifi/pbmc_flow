#!/bin/bash

git init
git add .
git commit -m "first commit"
git branch -M master
git remote add origin git@github.com:kevin-lee-aifi/pbmc_flow.git
git push -u origin master