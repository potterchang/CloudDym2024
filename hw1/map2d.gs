'reinit '
'set grads off'
datasetdir='../VVM/DATA/'
casename=pbl_ctl
* datasetdir='/data/poyen/taiwanvvm/'
* casename=tpe20190722cln

* open control file
'open 'datasetdir casename'/gs_ctl_files/dynamic.ctl'
'open 'datasetdir casename'/gs_ctl_files/thermodynamic.ctl'
'open 'datasetdir casename'/gs_ctl_files/bar.ctl'

* set region and time
'set x 1 128'
'set y 1 128'
'set lev 1000'
'set t 1'
* 'define pibarlev=0.966089'
'define pibarlev=pibar.3'
'set t 20'

* subplot (left)
'set vpage 0 5.5 0 8.5'
'set parea 1 5 1 7.5'
'set xlint 0.05'
'set ylint 0.05'
* thermodynamic: th, qv
'set gxout shaded'
'set ccolor rainbow'
* 'set cmin 0'
* 'set cmax 12'
'set cint 0.0005'
'd qv.2*1000'
'cbar'
'set gxout contour'
'define t2=th.2*pibarlev'
'set cmin 285'
'set cmax 300'
'set cint 0.01'
'd t2'
* add title
'draw title qv (shaded); th (contour)'
* subplot (right)
'set vpage 5.5 11 0 8.5'
'set parea 1 5 1 7.5'
'set xlint 0.05'
'set ylint 0.05'
* dynamic: w
'set gxout shaded'
'set ccolor rainbow'
'set cmin -0.5'
'set cmax 0.5'
'set cint 0.01'
'd w.1'
'cbar'
* dynamic: u, v
'set ccolor 1'
'set cthick 4'
'd skip(u.1,8);skip(v.1,8)'
* add title
'draw title Horizontal Wind (vector)\Vertical Motion (shaded)'
* save figure
'set grads off'
'set display white'
'printim test.png'
* close file
* 'close 1'

