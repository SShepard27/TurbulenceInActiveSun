;
function double_gaussian, x, c
; coeffs as follows:
; c0 - continuum (background) level 
; c1 - intensity of first gaussian
; c2 - centroid of first gaussian
; c3 - standard deviation of first gaussian
; c4 - intensity of second gaussian
; c5 - centroid of second gaussian
; c6 - standard deviation of second gaussian
gauss_1 = c(1)*exp(-(x-c(2))^2./(2.*c(3)^2.))
gauss_2 = c(4)*exp(-(x-c(5))^2./(2.*c(6)^2.))
model=c(0)+gauss_1+gauss_2
return,model
end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                 ;
;    MAIN PROGRAM                                 ;
;                                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                              
hard = 1 ; 0 for screen plot, 1 for hard copy  ;
;                                              
; name of folder where data FITS file is located, and name of the file
;
folder_name='c:\users\agr75681\Documents\IDL\IRIS\data_files\'
;
file_name='data_for_gordon_SiIV1394_oct242015_Region1.fits'
;file_name='data_for_gordon_OIV1402_oct242015_Region2.fits'
;
; name of the folder in which to place results
;
results_folder_name='c:\users\agr75681\Documents\IDL\IRIS\results\'
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; number of time intervals to average (to smooth the plots)
;
averaging_time_interval=4
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; get some colors
;
white=FSC_COLOR("White", !D.Table_Size-1)
black=FSC_COLOR("Black", !D.Table_Size-2)
!p.charsize=1
!p.charthick=1
!p.thick=1
!x.thick=1
!y.thick=1
!z.thick=1
backcolor=black
axiscolor=white
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                ;
; get atomic species and flare region associated with data file  ;
;                                                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
species=strmid(file_name,16,8)
region=strmid(file_name,strlen(file_name)-12,7)
event_date=strmid(file_name,strlen(file_name)-22,9)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                  ;
; properties of the spectral line under observation (e.g., Si IV)  ;
;                                                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
if species ne 'SiIV1394' and species ne 'OIV1402_' then begin
  print,'Impossible to Determine Ion Species from Filename'
  stop 
endif
;
if species eq 'SiIV1394' then begin
  t_formation=7.e4 ;formation temperature for Si IV (K)
  ion_mass=28.*1.67e-24 ; mass of Si ion in grams
  rest_wavelength=1393.75 ; rest wavelength in Angstroms
endif
;
if species eq 'OIV1402_' then begin
  t_formation=7.e4 ;formation temperature for Si IV (K)
  ion_mass=16.*1.67e-24 ; mass of O ion in grams
  rest_wavelength=1402.8 ; rest wavelength in Angstroms
  endif
;
v_thermal = sqrt(2.*1.38e-16*t_formation/ion_mass) ; v_{th} = \sqrt(2 k T/m}
;
inst_width = 0.01214 ; instrumental line width (1 standard deviation) in Angstroms
;
if region ne 'Region1' and region ne 'Region2' then begin
  print,'Impossible to Determine Flare Region from Filename'
  stop
endif
;
if region eq 'Region1' then begin
  i_time_min=180                               
  i_time_max=400
endif
;
if region eq 'Region2' then begin
  i_time_min=50                               ;
  i_time_max=100
endif
;
;;;;;;;;;;;;;;;;;;;;
;                  ;
; read in data     ;
;                  ; 
;;;;;;;;;;;;;;;;;;;;
;
wavelength_arr=readfits(folder_name+file_name,header,exten_no=0)
intensity=readfits(folder_name+file_name,header,exten_no=1)
time_arr=readfits(folder_name+file_name,header,exten_no=2)
;
; determine size of arrays
;
n_lambda=n_elements(wavelength_arr)
n_time=n_elements(time_arr)
n_pixel=n_elements(intensity(0,0,*))
intensity_pixel=dblarr(n_time,n_pixel) ; lightcurves for each pixel
intensity_pixel_plot=dblarr(n_time,n_pixel)
intensity_overall=dblarr(n_time) ; overall lightcurve
intensity_err=dblarr(n_lambda,n_time,n_pixel) ; uncertainty in intensity (Poissonian)
warning_flag=intarr(n_time,n_elements(intensity(0,0,*))) ; set to 1 if count in any channel exceeds 4000
fluence_pixel=dblarr(n_pixel) ; event-intgrated fluence for each pixel
;
;  set warning flag for saturated pixels/times
;  
for i_time=0,n_time-1 do begin
  for slit_pos=0,n_pixel-1 do begin
    for i_lambda=0,n_lambda-1 do begin
      if intensity(i_lambda,i_time,slit_pos) gt 4000. then warning_flag(i_time,slit_pos)=1
    endfor
  endfor
endfor
;
; determine wavelength-integrated intensity for each pixel and time
;
for slit_pos=0,n_pixel-1 do begin
  for i_time=0,n_time-1 do begin
    intensity_pixel(i_time,slit_pos)=total(intensity(*,i_time,slit_pos))
  endfor
endfor
;
; determine lightcurves for each pixel, excluding positions with warning flags
;
for i_time=0,n_time-1 do begin
  for slit_pos=0,n_pixel-1 do begin
    if warning_flag(i_time,slit_pos) eq 0 then $
      intensity_pixel_plot(i_time,slit_pos)=intensity_pixel(i_time,slit_pos) else $
      intensity_pixel_plot(i_time,slit_pos)=-999.
  endfor
endfor
;
; determine time-integrated fluence for each pixel in the image during intensity peak times
;
for slit_pos=0,n_pixel-1 do fluence_pixel(slit_pos)=total(intensity_pixel(*,slit_pos))
;
; determine lightcurve for the whole flare (including saturated 4000 values)
;
for i_time=0,n_time-2 do begin
  intensity_overall(i_time)=0.
  for slit_pos=0,n_pixel-1 do begin
    for i_lambda=0,n_lambda-2 do begin
      intensity_overall(i_time) = intensity_overall(i_time)+$
        intensity(i_lambda,i_time,slit_pos)
    endfor
  endfor
endfor
;
; find time of maximum overall intensity and brightest three pixels (with no warning flag) *at that time*
;
brightest_time=where(intensity_overall eq max(intensity_overall))
brightest_pixels_sorted=reverse(sort(intensity_pixel_plot(brightest_time,*)))
slit_pos_array_unsorted=[brightest_pixels_sorted(0),brightest_pixels_sorted(1),brightest_pixels_sorted(2)]
slit_pos_array=slit_pos_array_unsorted(sort(slit_pos_array_unsorted))
n_pos=n_elements(slit_pos_array)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                              ;
; plot pixel fluences, lightcurves for event and for pertinent slit positions  ;
;                                                                              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; set time ranges
time_plot_low=100.*floor(time_arr(i_time_min)/100.)
time_plot_high=100.*(floor(time_arr(i_time_max)/100.)+1)
time_range=[time_plot_low,time_plot_high]
;
if (hard eq 0) then window,0 else begin
  !p.charsize=1.5
  !p.charthick=3
  !x.thick=3
  !y.thick=3
  !z.thick=3
  !p.thick=3
  axiscolor=black ; axis color for hard copy
  backcolor=white ; background color for hard copy
  set_plot,'ps'
  device,/color,filename=results_folder_name+event_date+'_'+region+'_lightcurves.eps'
endelse
;
plot,fluence_pixel,xtitle='Pixel Number',position=[0.15,0.6,0.45,0.9],$
  title='Fluence (cts)',color=axiscolor,background=backcolor
plot,time_arr,intensity_overall,xtitle='Time (s)',position=[0.6,0.6,0.9,0.9],$
  xrange=time_range,xticks=4,xstyle=1,$
  title='LightCurve (c/s)',color=axiscolor,background=backcolor,/noerase
;
for i_pos=0,n_pos-1 do begin
  if i_pos eq 0 then y_label=replicate('',8) else y_label=replicate(' ',8)
  plot,time_arr,intensity_pixel_plot(*,slit_pos_array(i_pos)),$
    yrange=[0,max(intensity_pixel_plot)],$
    xrange=time_range,xticks=4,xtickformat='(a1)',$
    ytickname=y_label,$
    position=[0.1+0.3*i_pos,0.15,0.33+0.3*i_pos,0.45],/noerase
  xyouts,time_arr(i_time_min+5),0.8*max(intensity_pixel_plot),$
    'Pos = '+string(format='(i2)',slit_pos_array(i_pos))
endfor
;
if (hard eq 1) then begin
  !p.charsize=1
  !p.charthick=1
  !p.thick=1
  !x.thick=1
  !y.thick=1
  !z.thick=1
  axiscolor=white
  backcolor=black
  device,/close_file
  set_plot,'win'
endif ; restore plot parameters to defaults
;
;stop
;
; arrays for Gaussian fit coefficients and fitted line profile
;
; intensity_error
;
for i_lambda=0,n_lambda-1 do begin
  for i_time=0,n_time-1 do begin
    for slit_pos=0,n_pixel-1 do begin
      if intensity(i_lambda,i_time,slit_pos) gt 0 then $
        intensity_err(i_lambda,i_time,i_pos)=sqrt(intensity(i_lambda,i_time,slit_pos)) else $
        intensity_err(i_lambda,i_time,i_pos)=sqrt(-intensity(i_lambda,i_time,slit_pos))
    endfor
  endfor
endfor
;
coeff_fit=dblarr(n_time,n_pixel,7)
double_gauss_fit=dblarr(n_lambda,n_time,n_pixel)
gauss1_fit=dblarr(n_lambda,n_time,n_pixel)
gauss2_fit=dblarr(n_lambda,n_time,n_pixel)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                            ;
; main (slit position) loop  ;
;                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
for slit_pos=0,n_pixel-1 do begin
  ;
  print,'Slit Position = ',slit_pos ; just to keep track of what's going on
  ;
  ; do double Gaussian fits to each line profile
  ;
  coeff_info=replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 7)
  ;
  ; shift of principal component must be between -0.2 and +0.2 Angstroms
  ;
  coeff_info[2].limited[0] = 1
  coeff_info[2].limits[0] = rest_wavelength-0.2
  coeff_info[2].limited[1] = 1
  coeff_info[2].limits[1] = rest_wavelength+0.2
  ;
  ; width of principal component must be positive and less than 0.2 Angstroms
  ;
  coeff_info[3].limited[0] = 1
  coeff_info[3].limits[0] = 0.001
  coeff_info[3].limited[1] = 1
  coeff_info[3].limits[1] = 0.2
  ;
  ; intensity of shifted component must be positive 
  ;
  coeff_info[4].limited[0] = 1
  coeff_info[4].limits[0] = 1.
  ;
  ; shift of secondary component must be between -0.1 and +0.3 Angstroms
  ;
  coeff_info[5].limited[0] = 1
  coeff_info[5].limits[0] = rest_wavelength-0.1
  coeff_info[5].limited[1] = 1
  coeff_info[5].limits[1] = rest_wavelength+0.3
  ;
  ; width of shifted component must be positive and less than 0.2 Angstroms
  ;
  coeff_info[6].limited[0] = 1
  coeff_info[6].limits[0] = 0.001
  coeff_info[6].limited[1] = 1
  coeff_info[6].limits[1] = 0.2
  ;
  for i_time=0,n_time-1 do begin
    if warning_flag(i_time,slit_pos) ne 1 then begin
      start=[0.,max(intensity(*,i_time,slit_pos)),rest_wavelength,0.1,$
        max(intensity(*,i_time,slit_pos))/3.,rest_wavelength+0.2,0.1]
      result=mpfitfun('double_gaussian',wavelength_arr,intensity(*,i_time,slit_pos),$
        intensity_err(*,i_time,slit_pos),weights=replicate(1,n_lambda),start,/quiet,parinfo=coeff_info)
      for j=0,6 do coeff_fit(i_time,slit_pos,j)=result(j)
    endif 
  endfor
  ;
  for i_time=0,n_time-1 do begin
    for i_lambda=0,n_lambda-1 do begin
      gauss1_fit(i_lambda,i_time,slit_pos)=coeff_fit(i_time,slit_pos,0)+$
        coeff_fit(i_time,slit_pos,1)*exp(-(wavelength_arr(i_lambda)-coeff_fit(i_time,slit_pos,2))^2./$
        (2.*coeff_fit(i_time,slit_pos,3)^2.))
      gauss2_fit(i_lambda,i_time,slit_pos)=$
        coeff_fit(i_time,slit_pos,4)*exp(-(wavelength_arr(i_lambda)-coeff_fit(i_time,slit_pos,5))^2./(2.*coeff_fit(i_time,slit_pos,6)^2.))
      double_gauss_fit(i_lambda,i_time,slit_pos)=gauss1_fit(i_lambda,i_time,slit_pos)+gauss2_fit(i_lambda,i_time,slit_pos)
    endfor
  endfor
  ;
endfor
;
; plot sample line and superimposed Gaussians
;
for i_pos_plot=0,2 do begin
  if (hard eq 0) then begin
    window,i_pos_plot+1
    axiscolor=white ; axis color for screen print
    backcolor=black ; background color for screen print
  endif else begin
    !p.charsize=1.5
    !p.charthick=3
    !x.thick=3
    !y.thick=3
    !z.thick=3
    !p.thick=3
    axiscolor=black ; axis color for hard copy
    backcolor=white ; background color for hard copy
    set_plot,'ps'
    if brightest_time ge 100 then begin
      if slit_pos_array(i_pos_plot) lt 10 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i1)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i3)',brightest_time)+'.eps'
      if slit_pos_array(i_pos_plot) ge 10 and slit_pos_array(i_pos_plot) lt 100 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i2)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i3)',brightest_time)+'.eps'
      if slit_pos_array(i_pos_plot) ge 100 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i3)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i3)',brightest_time)+'.eps'
    endif
    if brightest_time lt 100 then begin
      if slit_pos_array(i_pos_plot) lt 10 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i1)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i2)',brightest_time)+'.eps'
      if slit_pos_array(i_pos_plot) ge 10 and slit_pos_array(i_pos_plot) lt 100 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i2)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i2)',brightest_time)+'.eps'
      if slit_pos_array(i_pos_plot) ge 100 then device,/color,filename=$
        results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
        string(format='(i3)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i2)',brightest_time)+'.eps'
    endif
  endelse
  ;
  plot,wavelength_arr,intensity(*,brightest_time,slit_pos_array(i_pos_plot)),$
    xtitle='Wavelength (!6!sA!r!u!9%!6!n )',ytitle='Intensity',$
    yrange=[0,max(intensity(*,brightest_time,slit_pos_array(i_pos_plot)))],$
    min_value=0.1,xticks=3
  oplot,wavelength_arr,gauss1_fit(*,brightest_time,slit_pos_array(i_pos_plot)),linestyle=1
  oplot,wavelength_arr,gauss2_fit(*,brightest_time,slit_pos_array(i_pos_plot)),linestyle=2
  ; oplot,wavelength,gauss1_fit(*,i_plot,i_pos)+gauss2_fit(*,i_plot,i_pos),linestyle=3
  xyouts,rest_wavelength-0.5,0.8*max(intensity(*,brightest_time,slit_pos_array(i_pos_plot))),$
    'Slit Pos = '+string(format='(i2)',slit_pos_array(i_pos_plot))
  xyouts,rest_wavelength-0.5,0.6*max(intensity(*,brightest_time,slit_pos_array(i_pos_plot))),$
    'Time = '+string(format='(i4)',round(time_arr(brightest_time))) 
  ;
  if (hard eq 1) then begin
    !p.charsize=1
    !p.charthick=1
    !p.thick=1
    !x.thick=1
    !y.thick=1
    !z.thick=1
    axiscolor=white
    backcolor=black
    device,/close_file
    set_plot,'win'
  endif ; restore plot parameters to defaults
;
endfor
;
if 2*(averaging_time_interval/2) ne averaging_time_interval then $
  averaging_time_interval=averaging_time_interval+1 ; make sure that # of times to average is +/- integer
;
coeff_fit_avg=dblarr(n_time,n_pixel,7)
;
for slit_pos=0,n_pixel-1 do begin
  for j=0,6 do begin
    for i_time=averaging_time_interval/2,n_time-1-averaging_time_interval/2 do begin
      values_to_be_averaged=coeff_fit([i_time-averaging_time_interval/2:i_time+averaging_time_interval/2],slit_pos,j)
      coeff_fit_avg(i_time,slit_pos,j)=mean(values_to_be_averaged [where(values_to_be_averaged ne 0)])
    endfor
  endfor
endfor
;
int_average_main=dblarr(n_pos)
int_average_shifted=dblarr(n_pos)
v_average_main=dblarr(n_pos)
v_average_shifted=dblarr(n_pos)
width_average_main=dblarr(n_pos)
width_average_shifted=dblarr(n_pos)
int_stdev_main=dblarr(n_pos)
int_stdev_shifted=dblarr(n_pos)
v_stdev_main=dblarr(n_pos)
v_stdev_shifted=dblarr(n_pos)
width_stdev_main=dblarr(n_pos)
width_stdev_shifted=dblarr(n_pos)
;
; set tickmarks for time (x) axis
;
time_plot_low=100*floor(time_arr(i_time_min-1+averaging_time_interval/2)/100)
time_plot_high=100*floor(time_arr(i_time_max-1-averaging_time_interval/2)/100)
time_range=[time_plot_low,time_plot_high]
;
for i_pos=0,n_pos-1 do begin
  if (hard eq 0) then begin
    window,i_pos+n_pos+1
    axiscolor=white ; axis color for screen print
    backcolor=black ; background color for screen print
  endif else begin
    !p.charsize=1.5
    !p.charthick=3
    !x.thick=3
    !y.thick=3
    !z.thick=3
    !p.thick=3
    axiscolor=black ; axis color for hard copy
    backcolor=white ; background color for hard copy
    set_plot,'ps'
    if slit_pos_array(i_pos) lt 10 then device,/color,filename=results_folder_name+event_date+'_'+region+'_fit_parameters_'+$
      string(format='(i1)',slit_pos_array(i_pos))+'.eps'
    if slit_pos_array(i_pos) ge 10 and slit_pos_array(i_pos) lt 100 then $
      device,/color,filename=results_folder_name+event_date+'_'+region+'_fit_parameters_'+$
      string(format='(i2)',slit_pos_array(i_pos))+'.eps'
    if slit_pos_array(i_pos) ge 100 then device,/color,filename=results_folder_name+event_date+'_'+region+'_fit_parameters_'+$
      string(format='(i3)',slit_pos_array(i_pos))+'.eps'
  endelse
  ;
  ; intensities
  ;
  plot,time_arr,coeff_fit_avg(*,slit_pos_array(i_pos),1),position=[0.1,0.15,0.35,0.9],$
    min_value=0.01,title='Intensity',xrange=time_range,xstyle=1,xticks=3,$
    yrange=[0,max(coeff_fit_avg(*,slit_pos_array(i_pos),1))],xtitle='t (s)'
  oplot,time_arr,coeff_fit_avg(*,slit_pos_array(i_pos),4),linestyle=2,min_value=0.01
  xyouts,time_arr(i_time_min+1),0.96*max(coeff_fit_avg(*,slit_pos_array(i_pos),1)),$
    'Slit Pos='+string(format='(i2)',slit_pos_array(i_pos))
  ;
  ; average of non-zero values
  ; 
  xx=coeff_fit_avg(*,slit_pos_array(i_pos),1)
  int_average_main(i_pos)=mean(xx(where(xx ne 0)))
  int_stdev_main(i_pos)=stdev(xx(where(xx ne 0)))
  xxx=coeff_fit_avg(*,slit_pos_array(i_pos),4)
  int_average_shifted(i_pos)=mean(xxx(where(xxx ne 0)))
  int_stdev_shifted(i_pos)=stdev(xxx(where(xxx ne 0)))
  ;
  ; mean wavelengths
  ; 
  plot,time_arr,(coeff_fit_avg(*,slit_pos_array(i_pos),2)-rest_wavelength)*$
    (3.e10/v_thermal)/rest_wavelength,position=[0.4,0.15,0.65,0.9],$
    xrange=time_range,min_value=0.001,xstyle=1,xticks=3,$
    title='Centroid Velocity',/noerase,yrange=[-2,8],xtitle='t (s)'
  oplot,time_arr,(coeff_fit_avg(*,slit_pos_array(i_pos),2)-rest_wavelength)*$
    (3.e10/v_thermal)/rest_wavelength,max_value=-0.001
  oplot,time_arr,(coeff_fit_avg(*,slit_pos_array(i_pos),5)-rest_wavelength)*$
    (3.e10/v_thermal)/rest_wavelength,min_value=0.001,linestyle=2
  oplot,time_arr,(coeff_fit_avg(*,slit_pos_array(i_pos),5)-rest_wavelength)*$
    (3.e10/v_thermal)/rest_wavelength,max_value=-0.001,linestyle=2
  ;
  ; average of non-zero values
  ;
  yy=coeff_fit_avg(*,slit_pos_array(i_pos),2)
  v_average_main(i_pos)=(mean(yy(where(yy ne 0)))-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength
  v_stdev_main(i_pos)=stdev(yy(where(yy ne 0)))*(3.e10/v_thermal)/rest_wavelength
  yyy=coeff_fit_avg(*,slit_pos_array(i_pos),5)
  v_average_shifted(i_pos)=(mean(yyy(where(yyy ne 0)))-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength
  v_stdev_shifted(i_pos)=stdev(yyy(where(yyy ne 0)))*(3.e10/v_thermal)/rest_wavelength
  ;
  ; widths
  ;
  plot,time_arr,sqrt(2.)*coeff_fit_avg(*,slit_pos_array(i_pos),3)*$
    (3.e10/v_thermal)/rest_wavelength,position=[0.7,0.15,0.95,0.9],$
    xrange=time_range,xstyle=1,xticks=3,$
    title='Width',/noerase,yrange=[0,6],xtitle='t (s)',min_value=0.01
  oplot,time_arr,sqrt(2.)*coeff_fit_avg(*,slit_pos_array(i_pos),6)*$
    (3.e10/v_thermal)/rest_wavelength,min_value=0.01,linestyle=2
;
  zz=coeff_fit_avg(*,slit_pos_array(i_pos),3)
  width_average_main(i_pos)=sqrt(2.)*mean(zz(where(zz ne 0)))*(3.e10/v_thermal)/rest_wavelength
  width_stdev_main(i_pos)=sqrt(2.)*stdev(zz(where(zz ne 0)))*(3.e10/v_thermal)/rest_wavelength
  zzz=coeff_fit_avg(*,slit_pos_array(i_pos),6)
  width_average_shifted(i_pos)=sqrt(2.)*mean(zzz(where(zzz ne 0)))*(3.e10/v_thermal)/rest_wavelength
  width_stdev_shifted(i_pos)=sqrt(2.)*stdev(zzz(where(zzz ne 0)))*(3.e10/v_thermal)/rest_wavelength
;
  if (hard eq 1) then begin
    !p.charsize=1
    !p.charthick=1
    !p.thick=1
    !x.thick=1
    !y.thick=1
    !z.thick=1
    axiscolor=white
    backcolor=black
    device,/close_file
    set_plot,'win'
  endif ; restore plot parameters to defaults
;
endfor
;
; print event-averaged values (and standard devations) of intensity, centroid velocity, and width
;    of both Gaussian components
;
print,'Main Component'
for i_pos=0,n_pos-1 do begin
  print,'Slit Position',slit_pos_array(i_pos)
  print,'Intensity',int_average_main(i_pos),'  p/m',int_stdev_main(i_pos)
  print,'Centroid Velocity',v_average_main(i_pos),'  p/m',v_stdev_main(i_pos)
  print,'Width',width_average_main(i_pos),'  p/m',width_stdev_main(i_pos) 
endfor
print,'Shifted Component'
for i_pos=0,n_pos-1 do begin
  print,'Slit Position',slit_pos_array(i_pos)
  print,'Intensity',int_average_shifted(i_pos),'  p/m',int_stdev_shifted(i_pos)
  print,'Centroid Velocity',v_average_shifted(i_pos),'  p/m',v_stdev_shifted(i_pos)
  print,'Width',width_average_shifted(i_pos),'  p/m',width_stdev_shifted(i_pos)
endfor
;
end