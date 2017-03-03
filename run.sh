#!/bin/bash -ex

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/eos.out $@
python $script_dir/plot_timeseries.py

