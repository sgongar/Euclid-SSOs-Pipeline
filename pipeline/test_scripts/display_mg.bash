#! /bin/bash


# $1 -> Display number
# $2 -> input image
# $3 -> output image


declare -i output_display=$1
# declare -s input_image=$2
# declare -s output_image=$3

echo $output_display


export DISPLAY=:$output_display

Xvfb :$output_display -ac -screen 0 1024x768x16 &
XVFB_PID=$!

echo "blabla"

echo $2
echo $3

# ds9 -fits $2 -saveimage $3 -exit

bash test.sh

# ds9 -version -exit

# ds9 $2 -zoom to fit -histequ -saveimage jpeg $3 -exit

kill $XVFB_PID

echo "out"

# echo out