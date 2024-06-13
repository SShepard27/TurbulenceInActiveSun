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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                              ;
hard = 1 ; 0 for screen plot, 1 for hard copy  ;
;                                              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; name of folder where data FITS file is located, and name of the file
;
folder_name='c:\users\agr75681\Documents\IDL\IRIS\data_files\'
;
;file_name='data_for_gordon_SiIV1394_oct242015_Region1.fits'
file_name='data_for_gordon_OIV1402_oct242015_Region2.fits'
;
; name of the folder in which to place results
;
results_folder_name='c:\users\agr75681\Documents\IDL\IRIS\results\'
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
  i_time_min=200                               
  i_time_max=300
  i_time_plot=240 ; time value for illustrative double gaussiam fit plot
  slit_pos_array=[16,20,24]
  ;slit_pos_array=[24]
endif
;
if region eq 'Region2' then begin
  i_time_min=0                               ;
  i_time_max=50
  i_time_plot=25  ; time value for illustrative double gaussiam fit plot
  slit_pos_array=[16,17,18]
  ;slit_pos_array=[18]
endif
;
;;;;;;;;;;;;;;;;;;;;
;                  ;
; read in data     ;
;                  ; 
;;;;;;;;;;;;;;;;;;;;
;
wavelength=readfits(folder_name+file_name,header,exten_no=0)
intensity=readfits(folder_name+file_name,header,exten_no=1)
time=readfits(folder_name+file_name,header,exten_no=2)
;
;stop
;
; determine size of arrays
;
n_lambda=n_elements(wavelength)
n_time=n_elements(time)
n_pos=n_elements(slit_pos_array)
intensity_pixel=dblarr(n_time,n_pos) ; lightcurves for each pixel
intensity_overall=dblarr(n_time) ; for overall lightcurve
intensity_err=dblarr(n_lambda,n_time,n_pos)
warning_flag=intarr(n_time,n_elements(intensity(0,0,*)))
;
;  need to filter out saturated pixels/times
;  
for i_time=0,n_time-1 do begin
  for slit_pos=0,n_elements(intensity(0,0,*))-1 do begin
    for i_lambda=0,n_lambda-1 do begin
      if intensity(i_lambda,i_time,slit_pos) gt 4000. then begin
        ;intensity(i_lambda,i_time,slit_pos) = 4000.
        warning_flag(i_time,slit_pos)=1
      endif
    endfor
  endfor
endfor
;  
; determine lightcurves for each pixel and for the overall flare
;
for i_time=0,n_time-1 do begin
  for i_pos=0,n_pos-1 do begin
    intensity_pixel(i_time,i_pos)=total(intensity(*,i_time,slit_pos_array(i_pos)-1))
  endfor
  intensity_overall(i_time)=0.
  for slit_pos=0,n_elements(intensity(0,0,*))-1 do begin
    intensity_overall(i_time) = intensity_overall(i_time) + total(intensity(*,i_time,slit_pos))
  endfor
endfor
;
;stop
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                              ;
; plot lightcurves for event and for Pertinent Slit Positions  ;
;                                                              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
plot,intensity_overall,xtitle='Time (s)',ytitle='Intensity (Counts)',position=[0.2,0.6,0.8,0.9],$
  title='Overall',color=axiscolor,background=backcolor
for i_slit_pos=0,n_pos-1 do begin
  if i_slit_pos eq 0 then y_label=replicate('',6) else y_label=replicate(' ',6)
  plot,intensity_pixel(*,i_slit_pos),$
    xtitle='Time (s)',xrange=[i_time_min,i_time_max],yrange=[0,1.e5],$
    xticks=2,xtickv=[i_time_min,(i_time_min+i_time_max)/2,i_time_max],ytickname=y_label,$
    position=[0.1+0.3*i_slit_pos,0.15,0.33+0.3*i_slit_pos,0.45],/noerase
  xyouts,i_time_min+0.2*(i_time_max-i_time_min),8e4,$
    'Pos = '+string(format='(i2)',slit_pos_array(i_slit_pos))
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
    for i_pos=0,n_pos-1 do begin
      if intensity(i_lambda,i_time,slit_pos_array(i_pos)) gt 0 then $
        intensity_err(i_lambda,i_time,i_pos)=sqrt(intensity(i_lambda,i_time,slit_pos_array(i_pos))) else $
        intensity_err(i_lambda,i_time,i_pos)=sqrt(-intensity(i_lambda,i_time,slit_pos_array(i_pos)))
    endfor
  endfor
endfor
;
coeff_fit=dblarr(n_time,n_pos,7)
double_gauss_fit=dblarr(n_lambda,n_time,n_pos)
gauss1_fit=dblarr(n_lambda,n_time,n_pos)
gauss2_fit=dblarr(n_lambda,n_time,n_pos)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                            ;
; main (slit position) loop  ;
;                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
for i_pos=0,n_pos-1 do begin
  ;
  print,'Slit Position = ',slit_pos_array(i_pos) ; just to keep track of what's going on
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
    if warning_flag(i_time,i_pos) ne 1 then begin
      start=[0.,max(intensity(*,i_time,slit_pos_array(i_pos))),rest_wavelength,0.1,$
        max(intensity(*,i_time,slit_pos_array(i_pos)))/3.,rest_wavelength+0.2,0.1]
      result=mpfitfun('double_gaussian',wavelength,intensity(*,i_time,slit_pos_array(i_pos)),$
        intensity_err(*,i_time,i_pos),weights=replicate(1,n_lambda),start,/quiet,parinfo=coeff_info)
      for j=0,6 do coeff_fit(i_time,i_pos,j)=result(j)
    endif 
  endfor
  ;
  ; check for which Gaussian is the larger of the two and reverse labels if necessary
  ;
  for i_time=0,n_time-1 do begin
    if coeff_fit(i_time,i_pos,1) lt coeff_fit(i_time,i_pos,4) then begin
      temp_big=coeff_fit(i_time,i_pos,4)
      coeff_fit(i_time,i_pos,4)=coeff_fit(i_time,i_pos,1)
      coeff_fit(i_time,i_pos,1)=temp_big
      shift_comp=coeff_fit(i_time,i_pos,2)
      coeff_fit(i_time,i_pos,2)=coeff_fit(i_time,i_pos,5)
      coeff_fit(i_time,i_pos,5)=shift_comp
      width_shifted=coeff_fit(i_time,i_pos,3)
      coeff_fit(i_time,i_pos,3)=coeff_fit(i_time,i_pos,6)
      coeff_fit(i_time,i_pos,6)=width_shifted
    endif
  endfor
  ; 
  for i_time=0,n_time-1 do begin
    for i_lambda=0,n_lambda-1 do begin
      gauss1_fit(i_lambda,i_time,i_pos)=coeff_fit(i_time,i_pos,0)+$
        coeff_fit(i_time,i_pos,1)*exp(-(wavelength(i_lambda)-coeff_fit(i_time,i_pos,2))^2./(2.*coeff_fit(i_time,i_pos,3)^2.))
      gauss2_fit(i_lambda,i_time,i_pos)=$
        coeff_fit(i_time,i_pos,4)*exp(-(wavelength(i_lambda)-coeff_fit(i_time,i_pos,5))^2./(2.*coeff_fit(i_time,i_pos,6)^2.))
      double_gauss_fit(i_lambda,i_time,i_pos)=gauss1_fit(i_lambda,i_time,i_pos)+gauss2_fit(i_lambda,i_time,i_pos)
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
    if i_time_plot ge 100 then device,/color,filename=results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
      string(format='(i2)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i3)',i_time_plot)+'.eps'
    if i_time_plot lt 100 then device,/color,filename=results_folder_name+event_date+'_'+region+'_double_gaussian_fit_'+$
      string(format='(i2)',slit_pos_array(i_pos_plot))+'_time_'+string(format='(i2)',i_time_plot)+'.eps'
  endelse
  ;
  plot,wavelength,intensity(*,i_time_plot,slit_pos_array(i_pos_plot)),$
    xtitle='Wavelength',ytitle='Intensity',$
    yrange=[0,max(intensity(*,i_time_plot,slit_pos_array(i_pos_plot)))],$
    min_value=0.1,color=axiscolor,background=backcolor,xticks=3
  oplot,wavelength,gauss1_fit(*,i_time_plot,i_pos_plot),linestyle=1
  oplot,wavelength,gauss2_fit(*,i_time_plot,i_pos_plot),linestyle=2
  ; oplot,wavelength,gauss1_fit(*,i_plot,i_pos)+gauss2_fit(*,i_plot,i_pos),linestyle=3
  xyouts,rest_wavelength,0.8*max(intensity(*,i_time_plot,slit_pos_array(i_pos_plot))),$
    'Slit Pos = '+string(format='(i2)',slit_pos_array(i_pos_plot))
  xyouts,rest_wavelength,0.6*max(intensity(*,i_time_plot,slit_pos_array(i_pos_plot))),$
    'Time = '+string(format='(i3)',i_time_plot) 
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
; plot intensities, mean wavelengths (in thermal velocity units), and standard deviations (in thermal velocity units) 
; 
;loadct,74
;xx=findgen(256)
;cc=bytscl(xx,top=!D.n_colors-1)
;main_component_color=cc(30)
;shifted_component_color=cc(200)
;
averaging_time_interval=4
;
if 2*(averaging_time_interval/2) ne averaging_time_interval then $
  averaging_time_interval=averaging_time_interval+1 ; make sure that # of times to average is +/- integer
;
coeff_fit_avg=dblarr(n_time,n_pos,7)
;
for i_pos=0,n_pos-1 do begin
  for j=0,6 do begin
    for i_time=averaging_time_interval/2,n_time-1-averaging_time_interval/2 do begin
      values_to_be_averaged=coeff_fit([i_time-averaging_time_interval/2:i_time+averaging_time_interval/2],i_pos,j)
      coeff_fit_avg(i_time,i_pos,j)=mean(values_to_be_averaged [where(values_to_be_averaged ne 0)])
    endfor
  endfor
endfor
;
time_range=[i_time_min-1+averaging_time_interval/2,i_time_max-1-averaging_time_interval/2]
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
    device,/color,filename=results_folder_name+event_date+'_'+region+'_fit_parameters_'+$
      string(format='(i2)',slit_pos_array(i_pos))+'.eps'
  endelse
  ;
  ; intensities
  ;
  plot,coeff_fit_avg(*,i_pos,1),position=[0.1,0.05,0.35,0.95],min_value=0.01,title='Intensity',$
    xrange=time_range,xticks=2,xtickv=[i_time_min,(i_time_min+i_time_max)/2,i_time_max],$
    yrange=[0,max(coeff_fit_avg(*,i_pos,1))],xtitle='t (s)'
  ;print,'i_pos =',i_pos
  oplot,coeff_fit_avg(*,i_pos,4),linestyle=2,min_value=0.01
  xyouts,i_time_min+10,0.96*max(coeff_fit_avg(*,i_pos,1)),'Slit Pos='+string(format='(i2)',slit_pos_array(i_pos))
  ;
  ;stop
  ;
  ; mean wavelengths
  ; 
  plot,(coeff_fit_avg(*,i_pos,2)-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength,$
    position=[0.4,0.05,0.65,0.95],xrange=time_range,min_value=0.001,$
    title='Centroid Velocity',/noerase,yrange=[-2,8],$
    xticks=2,xtickv=[i_time_min,(i_time_min+i_time_max)/2,i_time_max],$
    xtitle='t (s)'
  oplot,(coeff_fit_avg(*,i_pos,2)-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength,max_value=-0.001
  oplot,(coeff_fit_avg(*,i_pos,5)-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength,min_value=0.001,$
    linestyle=2
  oplot,(coeff_fit_avg(*,i_pos,5)-rest_wavelength)*(3.e10/v_thermal)/rest_wavelength,max_value=-0.001,$
    linestyle=2
  ;
  ; widths
  ;
  plot,sqrt(2.)*coeff_fit_avg(*,i_pos,3)*(3.e10/v_thermal)/rest_wavelength,position=[0.7,0.05,0.95,0.95],$
    title='Width',/noerase,yrange=[0,6],$;color=main_component_color,background=black
    xticks=2,xtickv=[i_time_min,(i_time_min+i_time_max)/2,i_time_max],xrange=time_range,$
    xtitle='t (s)',min_value=0.01
  oplot,sqrt(2.)*coeff_fit_avg(*,i_pos,6)*(3.e10/v_thermal)/rest_wavelength,min_value=0.01,$
    linestyle=2;color=shifted_component_color
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
end