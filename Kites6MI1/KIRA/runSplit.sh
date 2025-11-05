#!/bin/bash
python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_0.m' -o './splitres'
if [ -e ./results/hxb/numerics_2147483647_1.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_1.m' -o './splitres'
fi
if [ -e ./results/hxb/numerics_2147483647_2.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_2.m' -o './splitres'
fi
if [ -e ./results/hxb/numerics_2147483647_3.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_3.m' -o './splitres'
fi
if [ -e ./results/hxb/numerics_2147483647_4.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_4.m' -o './splitres'
fi
if [ -e ./results/hxb/numerics_2147483647_5.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_5.m' -o './splitres'
fi
if [ -e ./results/hxb/numerics_2147483647_6.m ]
then
    python3 SplitNumericalOutput.py -i './results/hxb/numerics_2147483647_6.m' -o './splitres'
fi