ffmpeg -framerate 2.5 -pattern_type glob -i 'CPcrossB5_tpe20190722cln_2mn-*.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p -y CPcrossB5_tpe20190722cln_2mn.mp4

