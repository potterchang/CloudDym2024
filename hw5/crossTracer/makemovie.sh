#crossTracer-pbl_half_PU_1-000272.png

casename="pbl_half_PU_1"
dmname="Y64pm16"

ffmpeg -framerate 2.5 -i "$casename/$dmname/crossTracer-${casename}-%06d.png" \
-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p -y "${casename}-${dmname}.mp4"

