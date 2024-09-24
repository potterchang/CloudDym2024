'reinit '
'set grads off'
casename=pbl_evergreen_qc
* open control file
'open ../VVM/DATA/'casename'/gs_ctl_files/topo.ctl'
* set region and time
'set x 1 128'
'set y 1 128'
* subplot (left)
'set vpage 0 5.5 0 8.5'
'set parea 1 5 1 7.5'
* set ticks
'set xlint 0.05'
'set ylint 0.05'
* plot filled contours
'set gxout shaded'
'set cmin 0'
'set cmax 50'
'set cint 1'
'd topo'
'cbar'
* add title
'draw title TOPO'
* subplot (left)
'set vpage 5.5 11 0 8.5'
'set parea 1 5 1 7.5'
* set ticks
'set xlint 0.05'
'set ylint 0.05'
* plot filled contours
'set gxout shaded'
'set cmin 0'
'set cmax 15'
'set cint 1'
'd lu'
'cbar'
* add title
'draw title LU'
* save figure
'set grads off'
'set display white'
'printim TOPO-'casename'.png'
* close file
* 'close 1'



