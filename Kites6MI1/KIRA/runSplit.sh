#!/bin/bash
python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_0.m' -o './splitres'
if [ -e ./results/hxb/numerics_2147483647_1.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_1.m' -o './splitres'
fi