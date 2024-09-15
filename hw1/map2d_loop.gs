'reinit'
'set grads off'
datasetdir='../VVM/DATA/'
casename=pbl_ctl

* 打開資料檔
'open 'datasetdir casename'/gs_ctl_files/dynamic.ctl'
'open 'datasetdir casename'/gs_ctl_files/thermodynamic.ctl'
'open 'datasetdir casename'/gs_ctl_files/bar.ctl'

* 設定區域與固定變數
'set x 1 128'
'set y 1 128'
'set lev 1000'
'set t 1'
'define pibarlev = pibar.3'
'q defval pibarlev 1 1'
pibarval=subwrd(result, 3)

* 固定的圖例設定（colorbar 與 vector scale）
'set gxout shaded'
'set ccolor rainbow'
'set cmin -0.5'
'set cmax 0.5'
'set cint 0.01'

* 開始迴圈處理時間步驟
tmax = 721
i = 1
while (i <= tmax)
  'c'
  'set t 'i

* subplot (左圖)
  'set vpage 0 5.5 0 8.5'
  'set parea 1 5 1 7.5'
  'set xlint 0.05'
  'set ylint 0.05'
* annotate time
  'q time'
   datetime=subwrd(result, 3)

* 繪製 th (shaded)
  'set gxout shaded'
  'set ccolor rainbow'
  'set cmin 285'
  'set cmax 290'
  'set cint 0.2'
  'define t2 = th.2*const(th.2,'pibarval')'
  'd t2'
  'cbar'

* 繪製 qv (contour)
  'set gxout contour'
  'set ccolor 1'
*  'set cmin 0.007'
*  'set cmax 0.012'
  'set cint 0.002'
  'd qv.2*1000'
  
* 添加標題
  'draw title 'datetime'\qv (contour); th (shaded)'

* subplot (右圖)
  'set vpage 5.5 11 0 8.5'
  'set parea 1 5 1 7.5'
  'set xlint 0.05'
  'set ylint 0.05'

* 繪製 w (shaded)
  'set gxout shaded'
  'set clevs -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5'
  'd w.1'
  'cbar'

* 繪製 u, v (vector)
  'set ccolor 1'
  'set cthick 4'
  'd skip(u.1,8); skip(v.1,8)'

* 添加標題
  'draw title Horizontal Wind (vector)\Vertical Motion (shaded)'

* 圖片儲存
  if i < 100
    fname = 'graphs/map2d-'casename'-0000'i'.png'
  else
    fname = 'graphs/map2d-'casename'-000'i'.png'
  endif
  if i < 10
    fname = 'graphs/map2d-'casename'-00000'i'.png'
  endif
  'printim 'fname' white'

* 設定迴圈下一個時間點
  i = i + 1
endwhile

