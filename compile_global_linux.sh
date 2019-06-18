#! /bin/bash

mpif90 -O3 FlowFill.f90 -o /usr/local/bin/flowfill -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -ffast-math -march=native -mtune=native
