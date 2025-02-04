#!/usr/bin/bash

dir=./figures
ffmpeg -i ${dir}/%3d.png ${dir}/output.gif
convert ${dir}/output.gif -rotate 180 ${dir}/rotated_output.gif

