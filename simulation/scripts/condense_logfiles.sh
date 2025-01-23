#!/bin/bash

for log in $(find out -name simulation.log); do
    lg="${log//\//_}"
    lg="${lg/out_/}"
    cp $log out/logfiles/$lg
done
