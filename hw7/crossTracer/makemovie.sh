#crossTracer-pbl_half_PU_1-000272.png

casename="pbl_half_PU_uarea_1"
dmname="tr01"

ffmpeg -framerate 5 -i "$casename/$dmname/crossTracer-${casename}-%06d.png" \
-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p -y "${casename}-${dmname}.mp4"

