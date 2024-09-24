'reinit '
'set grads off'
casename=pbl_evergreen_qc
* open control file
'open ../VVM/DATA/'casename'/gs_ctl_files/radiation.ctl'
'open ../VVM/DATA/'casename'/gs_ctl_files/surface.ctl'
* set region and time
'set x 1'
'set y 1'
'set t 1 421'
* subplot (left)
'set parea 1 10 1 7.5'
* plot line
'set gxout line'
* TOA downward short wave
'set ccolor 2'
'set cstyle 1'
'd fdswtoa.1'
* Surface flux of th (convert K kg m-2 s-1 to W m-2)
'define wthwm=wth.2*1004'
'set ccolor 12'
'set cstyle 1'
'd wthwm'
* Surface flux of qv (convert kg m-2 s-1 to W m-2)
'define wqvwm=wqv.2*2.5e6'
'set ccolor 4'
'set cstyle 1'
'd wqvwm'
'cbar_line2 -x 8 -y 7 -c 2 12 4 -m 2 2 2 -l 1 1 1 -t "fdswtoa" "wth" "wqv"'
* add title
'draw title Time Series of Energy'
* set label
'draw ylab Energy [W m-2]'
'draw xlab Local Time'
* save figure
'set grads off'
'printim energySeries-'casename'.png x1100 y850 white'
* close file
* 'close 1'

