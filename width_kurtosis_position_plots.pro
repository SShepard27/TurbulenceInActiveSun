;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
hard = 1 ; 0 for screen plot, 1 for hard copy
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
folder_name='c:\users\agr75681\Documents\IDL\IRIS\'
file_name='data_for_gordon_SiIV1394_oct242015.fits'
wavelength=readfits(folder_name+file_name,header,exten_no=0)
intensity=readfits(folder_name+file_name,header,exten_no=1)
time=readfits(folder_name+file_name,header,exten_no=2)
n_lambda=n_elements(wavelength)
n_time=n_elements(time)
n_pos=30
coeff_fit=fltarr(n_time,n_pos,4) ;set of gaussian fit coefficients
int_gauss=fltarr(n_lambda,n_time,n_pos)
;
for i_pos=0,n_pos-1,3 do begin
;
if (hard eq 0) then window,i_pos/3 else begin
  !p.charsize=1.5
  !p.charthick=3
  !x.thick=3
  !y.thick=3
  !z.thick=3
  !p.thick=3
  set_plot,'ps'
  device,filename=folder_name+'width_and_kurtosis_position_'+string(format='(i2)',i_pos)+'.eps'
endelse
;
for i_time=0,n_time-1 do begin
    gfit=gaussfit(wavelength,intensity(*,i_time,i_pos),coeff,nterms=4)
    for k = 0,3 do coeff_fit(i_time,i_pos,k)=coeff(k)
    for i_lambda=0,n_lambda-1 do int_gauss(i_lambda,i_time,i_pos)=coeff(3) + coeff(0) * exp(-(wavelength(i_lambda)-coeff(1))^2./(2.*coeff(2)^2.))
endfor
;
; normalized moments
;
delta_lambda=fltarr(n_lambda,n_time,n_pos)
raw_moment=fltarr(n_time,n_pos,5)
norm_moment=fltarr(n_time,n_pos,5)
sdev=fltarr(n_time,n_pos)
v_excess=fltarr(n_time,n_pos)
red_kurt=fltarr(n_time,n_pos)
;
t_formation=7.e4 ;formation temperature for Si IV (K)
ion_mass=28.*1.67e-24 ; mass of Si ion in grams
rest_wavelength=1394. ; rest wavelength in Angstroms
;
v_thermal = sqrt(1.38e-16*t_formation/ion_mass)
;
; range of wavelengths for moment calculations
;
;i_lambda_low=0
;i_lambda_high=n_lambda-1
;
i_lambda_low=5
i_lambda_high=30
;
for i_time=0,ntime-1 do begin
    for i_lambda=i_lambda_low,i_lambda_high do begin
      delta_lambda(i_lambda,i_time,i_pos)=wavelength(i_lambda)-coeff_fit(i_time,i_pos,1) ; deviation from centroid position
    endfor
    for n_mom=0,4 do begin
      for i_lambda=i_lambda_low,i_lambda_high do begin
        raw_moment(i_time,i_pos,n_mom)=raw_moment(i_time,i_pos,n_mom)+(delta_lambda(i_lambda,i_time,i_pos))^n_mom * intensity(i_lambda,i_time,i_pos) 
      endfor
    endfor
    ;
    ; normalized moments
    ;
    for n_mom=0,4 do begin
      norm_moment(i_time,i_pos,n_mom)=raw_moment(i_time,i_pos,n_mom)/raw_moment(i_time,i_pos,0)
    endfor
    ;
    ; width and (reduced) kurtosis
    ;
    sdev(i_time,i_pos)=sqrt(norm_moment(i_time,i_pos,2)) ; sdev in Angstrom units
    v_excess(i_time,i_pos)=sqrt(((sdev(i_time,i_pos)/rest_wavelength)*3.e10/v_thermal)^2. - 1.) ; excess velocity width in units of v_thermal
    red_kurt(i_time,i_pos)=norm_moment(i_time,i_pos,4)/sdev(i_time,i_pos)^4. - 3.
endfor
;
plot,v_excess(*,i_pos),red_kurt(*,i_pos),xtitle='Excess Width (in units of v!dth!n)',ytitle='Reduced Kurtosis',xrange=[0,10],yrange=[-10,10]
xyouts,6,-5,'Slit Position = '+string(format='(i2)',i_pos)
xyouts,6,-7,'Times as indicated'
for i_time=200,300,10 do xyouts,v_excess(i_time,i_pos),red_kurt(i_time,i_pos),string(format='(i3)',i_time)
;
if (hard eq 1) then begin
  !p.charsize=1
  !p.charthick=1
  !p.thick=1
  !x.thick=1
  !y.thick=1
  !z.thick=1
  device,/close_file
  set_plot,'win'
endif ; restore plot parameters to defaults
;
endfor
;
end