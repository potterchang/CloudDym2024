'reinit '
'set grads off'
datasetdir='../VVM/DATA/'
casename=pbl_evergreen_qc
* open control file
'open 'datasetdir casename'/gs_ctl_files/dynamic.ctl'
'open 'datasetdir casename'/gs_ctl_files/thermodynamic.ctl'
* set region and time
'set x 1'
'set y 1'
'set lev 0 2000'

* subplot (1): th
'set vpage 0.5 3.25 1 7.5'
'set parea 1 3.25 0.5 7.5'
'set vrange 292 310'
'set xlint 4' 
* plot line
'set gxout line'
'set cstyle 1'
'set t 1'
'set ccolor 1'
'd th.2'
'set t 31'
'set ccolor 3'
'd th.2'
'set t 211'
'set ccolor 2'
'd th.2'
'set t 391'
'set ccolor 8'
'd th.2'
* 'set t 571'
* 'set ccolor 4'
* 'd th.2'
* add title
'draw title theta (K)'
'draw ylab height [m]'
* legend
* 'cbar_line2 -x 1.25 -y 7.25 -c 1 3 2 8 4 -m 2 3 4 5 1 -l 1 1 1 1 1 -t "Initial" "06LT" "12LT" "18LT" "24LT"'
'cbar_line2 -x 1.25 -y 7.25 -c 1 3 2 8 -m 2 3 4 5 -l 1 1 1 1 -t "Initial" "06LT" "12LT" "18LT"'

* subplot (2): qv
'set vpage 2.75 5.5 1 7.5'
'set parea 1 3.25 0.5 7.5'
'set vrange 0 20'
'set xlint 4'
* plot line
'set gxout line'
'set cstyle 1'
'set t 1'
'set ccolor 1'
'd qv.2*1000'
'set t 31'
'set ccolor 3'
'd qv.2*1000'
'set t 211'
'set ccolor 2'
'd qv.2*1000'
'set t 391'
'set ccolor 8'
'd qv.2*1000'
* 'set t 571'
* 'set ccolor 4'
* 'd qv.2*1000*1000'
* add title
'draw title qv (g/kg)'

* subplot (3): u
'set vpage 5.0 7.75 1 7.5'
'set parea 1 3.25 0.5 7.5'
'set vrange -1.5 1.5'
'set xlint 0.5'
* plot barb (BUT STILL PLOT LINE)
* 'set gxout barb'
* 'set t 211'
* 'd u.1;v.1'
* plot line
'set gxout barb'
'set cstyle 1'
'set t 1'
'set ccolor 1'
'd u.1'
'set t 31'
'set ccolor 3'
'd u.1'
'set t 211'
'set ccolor 2'
'd u.1'
'set t 391'
'set ccolor 8'
'd u.1'
* 'set t 571'
* 'set ccolor 4'
* 'd u.1'
* add title
'draw title u (m/s)'

* subplot (4): v
'set vpage 7.25 10 1 7.5'
'set parea 1 3.25 0.5 7.5'
'set vrange -1.5 1.5'
'set xlint 0.5'
* plot barb (BUT PLOT LINE INSTEAD)
* 'set gxout barb'
* 'set t 151'
* 'd u.1;v.1'
* plot line
'set gxout barb'
'set cstyle 1'
'set t 1'
'set ccolor 1'
'd v.1'
'set t 31'
'set ccolor 3'
'd v.1'
'set t 211'
'set ccolor 2'
'd v.1'
'set t 391'
'set ccolor 8'
'd v.1'
* 'set t 571'
* 'set ccolor 4'
* 'd v.1'
* add title
'draw title v (m/s)'

* save figure
'set grads off'
'printim profile-'casename'.png x1650 y1275 white'
* close file
* 'close 1'

