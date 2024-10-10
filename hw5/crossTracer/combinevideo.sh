casename1="pbl_half_PU_1"
casename2="pbl_half_PU_2"
dmname="Y64pm16"
INPUT1="$casename1-$dmname.mp4"
INPUT2="$casename2-$dmname.mp4"

ffmpeg -i $INPUT1 -i $INPUT2 -filter_complex "hstack=inputs=2" -c:v libx264 -pix_fmt yuv420p -y "crossChem-$dmname-$casename1-$casename2.mp4"

