#!/usr/bin/env bash
touch demo/data/$1/log.txt
./LDRT demo/data/$1/$1.xml 2>&1 | tee demo/data/$1/log.txt
