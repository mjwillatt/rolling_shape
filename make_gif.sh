#!/usr/bin/bash

dir=./figures

#n = number of images files
n=$(ls -l ${dir}/* | awk '/[0-9]+.png/' | wc -l)
#n2 = n/2
n2=$(echo ${n} | awk '{printf "%.4i", $0/2}')
#generate palette from the middle image
ffmpeg -i ${dir}/${n2}.png -vf palettegen ${dir}/palette.png
#create the gif using the palette
ffmpeg -i ${dir}/%4d.png -i ${dir}/palette.png -r 50 -lavfi paletteuse ${dir}/output.gif 
#rotate the gif 180 degrees
convert ${dir}/output.gif -rotate 180 ${dir}/rotated_output.gif

