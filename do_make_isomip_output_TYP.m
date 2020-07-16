%function do_make_isomip_output(hisname,grdname,outname) 
% read in file
if 0 %ocean0
hisname = 'TYP/ocean_his_ocean0_TYP.nc';
grdname = 'Ocean1/isomip_plus_ocean1.nc';
outname = 'TYP/TMP.Ocean0_TYP_ROMSUTAS.nc'
elseif 0 %ocean1
hisname = 'TYP/ocean_his_ocean1_TYP.nc';
grdname = 'Ocean1/isomip_plus_ocean1.nc';
outname = 'TYP/TMP.Ocean1_TYP_ROMSUTAS.nc'
elseif 1 %ocean2
hisname = 'TYP/ocean_his_ocean2_TYP.nc';
grdname = 'Ocean2/isomip_plus_ocean2.nc';
outname = 'TYP/TMP.Ocean2_TYP_ROMSUTAS.nc'
end
Vtransform=2;
Vstretching=1;

if 1%if 1 also redo analysis
% Get vertical coords
h=flipud(ncread(grdname,'h'));
zice=flipud(ncread(grdname,'zice'));
lat_rho=ncread(grdname,'lat_rho');
y_rho=ncread(grdname,'y_rho');
x_rho=ncread(grdname,'x_rho');
lon_rho=ncread(grdname,'lon_rho');
lat_u=ncread(grdname,'lat_u');
lon_u=ncread(grdname,'lon_u');
lat_v=ncread(grdname,'lat_v');
lon_v=ncread(grdname,'lon_v');
mask_rho=flipud(ncread(grdname,'mask_rho'));
mask_zice=flipud(ncread(grdname,'mask_zice'));
pm=flipud(ncread(grdname,'pm'));
pn=flipud(ncread(grdname,'pn'));
Cs_w=ncread(hisname,'Cs_w');
Cs_r=ncread(hisname,'Cs_r');

Z_interp=-717.5:5:-2.5; %m
N = length(Cs_r);
dx = 1./pm;
dy = 1./pn;
mask_closed=ones(size(mask_zice)); mask_closed([1,end],:)=NaN; mask_closed(:,[1,end])=NaN;
mask_zice(mask_zice==0)=NaN; mask_zice=mask_zice;
mask_rho(mask_rho==0)=NaN; mask_rho=mask_rho;
Area_ice = dx.*dy.*mask_zice.*mask_closed;
Area_water=dx.*dy.*mask_rho.*mask_closed;
rho_i=918;%kgm^-3
Cd=2.5e-3; %dimensionless
u_t=0.01;%m/s
addpath('/ds/projects/iomp/matlab_scripts/ROMS_MATLAB/utility/')
Zw=nan(size(lat_rho,1),size(lat_rho,2),N+1);
Z=nan(size(lat_rho,1),size(lat_rho,2),N);
for jj=1:size(lat_rho,1)
[z,~,~]=scoord(h,zice,lon_rho,lat_rho,Vtransform,Vstretching,4,0.9,20,N,0,1,jj,0); %z is 3d of depths of each cell.
[zw,~,~]=scoord(h,zice,lon_rho,lat_rho,Vtransform,Vstretching,4,0.9,20,N,1,1,jj,0); %z is 3d of depths of each cell.
Z(jj,:,:)=z;
Zw(jj,:,:)=zw;
end
rmpath('/ds/projects/iomp/matlab_scripts/ROMS_MATLAB/utility/')
%Z = bsxfun(@plus,bsxfun(@times,h+zice,shiftdim(Cs_r,-2)),zice); %manually make Z-levels with loaded Cs_r
%Zw = bsxfun(@plus,bsxfun(@times,h+zice,shiftdim(Cs_w,-2)),zice);
dz = Zw(:,:,2:end)-Zw(:,:,1:end-1);
dxdydz=bsxfun(@times,bsxfun(@times,dz,dx),dy);
dxdydz_masked=bsxfun(@times,dxdydz,mask_rho.*mask_closed);

% load data from output file
m = flipud(double(ncread(hisname,'m')));
zeta = flipud(ncread(hisname,'zeta'));
temp = flipud(ncread(hisname,'temp'));
salt = flipud(ncread(hisname,'salt'));
u = flipud(ncread(hisname,'u'));
v = flipud(ncread(hisname,'v'));
Tb = flipud(double(ncread(hisname,'Tb')));
Sb = flipud(double(ncread(hisname,'Sb')));
ubar=flipud(ncread(hisname,'ubar'));
vbar=flipud(ncread(hisname,'vbar'));
ocean_time=ncread(hisname,'ocean_time');
% make requested data

meanMeltRate = squeeze(nansum(nansum(bsxfun(@times,m,Area_ice),2),1)) / squeeze(nansum(nansum(Area_ice,2),1));
meanMeltFlux = squeeze(nansum(nansum(bsxfun(@times,m,Area_ice*rho_i),2),1));
sumOfVolume = squeeze(nansum(nansum(bsxfun(@rdivide,bsxfun(@plus,(h+zice).*mask_rho.*mask_closed,zeta),(pm.*pn)),2),1));
meanTemperature = squeeze(nansum(nansum(nansum(bsxfun(@times,temp,bsxfun(@times,dxdydz,mask_rho.*mask_closed)),3),2),1)./nansum(nansum(nansum(bsxfun(@times,dxdydz,mask_rho.*mask_closed),3),2),1));
meanSalinity = squeeze(nansum(nansum(nansum(bsxfun(@times,salt,bsxfun(@times,dxdydz,mask_rho.*mask_closed)),3),2),1)./nansum(nansum(nansum(bsxfun(@times,dxdydz,mask_rho.*mask_closed),3),2),1));
meltRate=bsxfun(@times,m,mask_zice.*mask_closed);
%u*
u_mod = (u(1:end-1,2:end-1,:,:)+u(2:end,2:end-1,:,:))/2;
v_mod = (v(2:end-1,1:end-1,:,:)+v(2:end-1,2:end,:,:))/2;
u_mod(end+1,:,:,:)=NaN; u_mod(:,end+1,:,:)=NaN; 
v_mod(end+1,:,:,:)=NaN;v_mod(:,end+1,:,:)=NaN;
u_mod(2:end+1,:,:,:)=u_mod;u_mod(1,:,:,:)=NaN;u_mod(:,2:end+1,:,:)=u_mod;u_mod(:,1,:,:)=NaN;
v_mod(2:end+1,:,:,:)=v_mod;v_mod(1,:,:,:)=NaN;v_mod(:,2:end+1,:,:)=v_mod;v_mod(:,1,:,:)=NaN;
Ustar = sqrt(Cd)*sqrt(sqrt(squeeze(u_mod(:,:,N,:)).^2 + squeeze(v_mod(:,:,N,:)).^2).^2+u_t^2); disp('ustar computed post-results from sqrt(cd)*sqrt(u^2+u_t^2)') 
Ustar = bsxfun(@times,Ustar,mask_zice.*mask_closed);
%Tb is in situ -> convert to pot
addpath(genpath('/ds/projects/iomp/matlab_scripts/GSW'))
p_top = gsw_p_from_z(squeeze(Zw(:,:,N)),lat_rho);
disp('Changing basal temperature to potential temperature')
for tt = 1:size(temp,4)
TbPot(:,:,tt) = gsw_pt0_from_t(gsw_SA_from_SP(squeeze(Sb(:,:,tt)),p_top,lon_rho,lat_rho),squeeze(Tb(:,:,tt)),p_top);
if ~rem(tt,12)
disp([num2str(tt/size(temp,4)*100) '%'])
end
end
rmpath(genpath('/ds/projects/iomp/matlab_scripts/GSW'))
Tstar = bsxfun(@times,squeeze(temp(:,:,N,:))-TbPot,mask_zice.*mask_closed); 
Sstar = bsxfun(@times,squeeze(salt(:,:,N,:))-Sb,mask_zice.*mask_closed);
u_top = bsxfun(@times,squeeze(u_mod(:,:,N,:)),mask_zice.*mask_closed);
v_top = bsxfun(@times,squeeze(v_mod(:,:,N,:)),mask_zice.*mask_closed);
%psi = -cumsum(ubarMod2.*(h_Mod+zice_Mod).*dy(2:end-1,2:end-1),2);
disp('calculate barotropic SF')
if 0 % method to move u/v to rho points for summation
ubar(:,end,:)=0;
ubar_mod = 1/2*(ubar(1:end-1,2:end-1,:)+ubar(2:end,2:end-1,:)); ubar_mod(isnan(ubar_mod))=0;
vbar_mod = 1/2*(vbar(1:end-1,2:end-1,:)+vbar(2:end,2:end-1,:)); vbar_mod(isnan(vbar_mod))=0;
baroSF = -cumsum(bsxfun(@times,ubar_mod,(h(2:end-1,2:end-1)+zice(2:end-1,2:end-1)).*dy(2:end-1,2:end-1)),2);
baroSF(end+1,:,:,:)=NaN; baroSF(:,end+1,:,:)=NaN;
baroSF(2:end+1,:,:,:)=baroSF;baroSF(1,:,:,:)=NaN;baroSF(:,2:end+1,:,:)=baroSF;baroSF(:,1,:,:)=NaN;
elseif 1 % method to move h/zice/pn/pm to u,v points
et_s=1;
et_e=240;
et_c=240;
xi_s=1;
xi_e=40;
xi_c=40;
 % water column thickness at velocity points
hv = (h(2:end-1,1:end-1)+h(2:end-1,2:end))/2 + (zice(2:end-1,1:end-1)+zice(2:end-1,2:end))/2;
hu = (h(1:end-1,2:end-1)+h(2:end,2:end-1))/2 + (zice(1:end-1,2:end-1)+zice(1:end-1,2:end-1))/2;
 % grid size at velocity points
dx = 2./(pm(:,1:end-1)+pm(:,2:end));
de = 2./(pn(1:end-1,:)+pn(2:end,:));
 % initialise psi
psi = zeros(size(h));
 % make nan ubars 0
vbar(isnan(vbar))=0;  ubar(isnan(ubar))=0;
 % compute transports at density points
tu = (ubar(1:end-1,et_s+1:et_e-1,:).*de(1:end-1,et_s+1:et_e-1).*hu(1:end-1,:)+ubar(2:end,et_s+1:et_e-1,:).*de(2:end,et_s+1:et_e-1).*hu(2:end,:))./2;
tv = (vbar(xi_s+1:xi_e-1,1:end-1,:).*dx(xi_s+1:xi_e-1,1:end-1).*hv(:,1:end-1)+vbar(xi_s+1:xi_e-1,2:end,:).*dx(xi_s+1:xi_e-1,2:end).*hv(:,2:end))./2;
 % calculate baro SF
psi=-cumsum(tu,2);
baroSF = nan([size(m)]);
baroSF(2:end-1,2:end-1,:)=psi;
end
%
if 0 %calculate OTSF
% dx*dz*v integrate into page, cumsum bottom-top
disp('interpolating v-velocity to z-levels')
[X_isurf,Y_isurf,Z_isurf]=meshgrid(squeeze(x_rho(1,:)),squeeze(y_rho(:,1)),Z_interp);
v_mod_i = nan(size(X_isurf,1),size(X_isurf,2),size(X_isurf,3),size(v_mod,4));
X_RHO = repmat(x_rho,[1 1 N]);
Y_RHO = repmat(y_rho,[1 1 N]);
for tt=1:size(v_mod,4)
v_mod_i(3:end-2,:,:,tt) = griddata(squeeze(X_RHO(3:end-2,:,:)),squeeze(Y_RHO(3:end-2,:,:)),squeeze(Z(3:end-2,:,:)),squeeze(v_mod(3:end-2,:,:,tt)),squeeze(X_isurf(3:end-2,:,:)),squeeze(Y_isurf(3:end-2,:,:)),squeeze(Z_isurf(3:end-2,:,:)));
v_mod_i(2,:,:,tt) = griddata(squeeze(X_RHO(2,:,:)),squeeze(Z(2,:,:)),squeeze(v_mod(2,:,:,tt)),squeeze(X_isurf(end-1,:,:)),squeeze(Z_isurf(end-1,:,:)));
v_mod_i(end-1,:,:,tt) = griddata(squeeze(X_RHO(end-1,:,:)),squeeze(Z(end-1,:,:)),squeeze(v_mod(end-1,:,:,tt)),squeeze(X_isurf(end-1,:,:)),squeeze(Z_isurf(end-1,:,:)));
if ~rem(tt,12)
disp([num2str(tt/size(v_mod,4)*100) '%'])
end
end
disp('masking interpolation within convex hull, but not within data')
%remove bad interp within convex hull
IN.maskICE=nan(size(v_mod_i,1),size(v_mod_i,2),size(v_mod_i,3));
disp('making ice mask for otsf')
for jj=1:size(v_mod,1)
IN.maskICE(jj,:,:) = double(inpolygon(squeeze(X_isurf(jj,:,:)),squeeze(Z_isurf(jj,:,:)),...
[squeeze(X_RHO(jj,:,N)),squeeze(X_RHO(jj,end,:))',fliplr(squeeze(X_RHO(jj,:,1))),fliplr(squeeze(X_RHO(jj,1,:))')],...
[squeeze(Z(jj,:,N)),squeeze(Z(jj,end,:))',fliplr(squeeze(Z(jj,:,1))),fliplr(squeeze(Z(jj,1,:))')]));
if ~rem(jj,5)
disp([num2str(jj/size(v_mod,1)*100) '%'])
end
end
IN.maskICE(IN.maskICE==0)=NaN;
v_mod_I = bsxfun(@times,v_mod_i,IN.maskICE); v_mod_I(isnan(v_mod_I))=0;
vbarx2= squeeze(nansum(bsxfun(@times,bsxfun(@times,v_mod_I,dx),(Z_interp(2)-Z_interp(1))*ones(size(dx))),1));
otSF = cumsum(vbarx2,2);
else
otSF = zeros(size(x_rho,2),length(Z_interp),length(ocean_time));
end

T_bot = bsxfun(@times,squeeze(temp(:,:,1,:)),mask_rho.*mask_closed);
S_bot = bsxfun(@times,squeeze(salt(:,:,1,:)),mask_rho.*mask_closed);

%interpolation of the form:
X_trans1=repmat(squeeze(x_rho(20,:))',[1 N]);
Z_trans1=squeeze(Z(20,:,:));
X_trans2=repmat(squeeze(x_rho(21,:))',[1 N]);
Z_trans2=squeeze(Z(21,:,:));
Y_trans1=repmat(squeeze(y_rho(:,100)),[1 N]);
Z_trans3=squeeze(Z(:,100,:));
Y_trans2=repmat(squeeze(y_rho(:,101)),[1 N]);
Z_trans4=squeeze(Z(:,101,:));
[X_isurf,Z_isurf]=meshgrid(squeeze(x_rho(20,:)),Z_interp);
tempXZ1=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
tempXZ2=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
saltXZ1=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
saltXZ2=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
[Y_isurf,Z2_isurf]=meshgrid(squeeze(y_rho(:,100)),Z_interp);
tempYZ1=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
tempYZ2=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
saltYZ1=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
saltYZ2=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
disp('interpolating transect data')
for tt=1:size(temp,4)
tempXZ1(:,:,tt) = griddata(X_trans1,Z_trans1,squeeze(temp(20,:,:,tt)),X_isurf,Z_isurf);
tempXZ2(:,:,tt) = griddata(X_trans2,Z_trans2,squeeze(temp(21,:,:,tt)),X_isurf,Z_isurf);
saltXZ1(:,:,tt) = griddata(X_trans1,Z_trans1,squeeze(salt(20,:,:,tt)),X_isurf,Z_isurf);
saltXZ2(:,:,tt) = griddata(X_trans2,Z_trans2,squeeze(salt(21,:,:,tt)),X_isurf,Z_isurf);
tempYZ1(:,:,tt) = griddata(Y_trans1,Z_trans3,squeeze(temp(:,100,:,tt)),Y_isurf,Z2_isurf);
tempYZ2(:,:,tt) = griddata(Y_trans2,Z_trans4,squeeze(temp(:,101,:,tt)),Y_isurf,Z2_isurf);
saltYZ1(:,:,tt) = griddata(Y_trans1,Z_trans3,squeeze(salt(:,100,:,tt)),Y_isurf,Z2_isurf);
saltYZ2(:,:,tt) = griddata(Y_trans2,Z_trans4,squeeze(salt(:,101,:,tt)),Y_isurf,Z2_isurf);
if ~rem(tt,12)
disp([num2str(tt/size(temp,4)*100) '%'])
end
end

%remove interp within convex hull
IN.maskXZ1 = double(inpolygon(X_isurf,Z_isurf,[X_trans1(:,N);X_trans1(end,:)';flipud(X_trans1(:,1));flipud(X_trans1(1,:)')],[Z_trans1(:,N);Z_trans1(end,:)';flipud(Z_trans1(:,1));flipud(Z_trans1(1,:)')])); IN.maskXZ1(IN.maskXZ1==0)=NaN;
IN.maskXZ2 = double(inpolygon(X_isurf,Z_isurf,[X_trans2(:,N);X_trans2(end,:)';flipud(X_trans2(:,1));flipud(X_trans2(1,:)')],[Z_trans2(:,N);Z_trans2(end,:)';flipud(Z_trans2(:,1));flipud(Z_trans2(1,:)')])); IN.maskXZ2(IN.maskXZ2==0)=NaN;
IN.maskYZ1 = double(inpolygon(Y_isurf,Z2_isurf,[Y_trans1(:,N);Y_trans1(end,:)';flipud(Y_trans1(:,1));flipud(Y_trans1(1,:)')],[Z_trans3(:,N);Z_trans3(end,:)';flipud(Z_trans3(:,1));flipud(Z_trans3(1,:)')])); IN.maskYZ1(IN.maskYZ1==0)=NaN;
IN.maskYZ2 = double(inpolygon(Y_isurf,Z2_isurf,[Y_trans2(:,N);Y_trans2(end,:)';flipud(Y_trans2(:,1));flipud(Y_trans2(1,:)')],[Z_trans4(:,N);Z_trans4(end,:)';flipud(Z_trans4(:,1));flipud(Z_trans4(1,:)')])); IN.maskYZ2(IN.maskYZ2==0)=NaN;

tempXZ=( bsxfun(@times,tempXZ1,IN.maskXZ1) + bsxfun(@times,tempXZ2,IN.maskXZ2) )/2;
saltXZ=( bsxfun(@times,saltXZ1,IN.maskXZ1) + bsxfun(@times,saltXZ2,IN.maskXZ2) )/2;
tempYZ=( bsxfun(@times,tempYZ1,IN.maskYZ1) + bsxfun(@times,tempYZ2,IN.maskYZ2) )/2;
saltYZ=( bsxfun(@times,saltYZ1,IN.maskYZ1) + bsxfun(@times,saltYZ2,IN.maskYZ2) )/2;

%end  %breakpoint
end %if just make netcdf

disp(' ')
disp([' Creating the file : ',outname])
disp(' ')

id = netcdf.create(outname, 'clobber');

% define dims
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;
nx_dim = netcdf.defDim(id, 'nx', Mp);
ny_dim = netcdf.defDim(id, 'ny', Lp);
nz_dim = netcdf.defDim(id, 'nz', length(Z_interp));
nTime_dim=netcdf.defDim(id, 'nTime', length(ocean_time));

%define vars
x_id = netcdf.defVar(id, 'x', 'double', nx_dim);
netcdf.putAtt(id, x_id, 'long_name', 'cell centre in x-direction');
netcdf.putAtt(id, x_id, 'units', 'meter');
y_id = netcdf.defVar(id, 'y', 'double', ny_dim);
netcdf.putAtt(id, y_id, 'long_name', 'cell centre in y-direction');
netcdf.putAtt(id, y_id, 'units', 'meter');
z_id = netcdf.defVar(id, 'z', 'double', nz_dim);
netcdf.putAtt(id, z_id, 'long_name', 'cell centre in z-direction');
netcdf.putAtt(id, z_id, 'units', 'meter');
time_id = netcdf.defVar(id, 'time', 'double', nTime_dim);
netcdf.putAtt(id, time_id, 'long_name', 'time from the start of the simulation');
netcdf.putAtt(id, time_id, 'units', 'second');


meanMeltRate_id = netcdf.defVar(id, 'meanMeltRate', 'double', nTime_dim);
netcdf.putAtt(id, meanMeltRate_id, 'long_name', 'water equivalent melt rate, positive for melting and negative for freezing, averaged over the ice-shelf base');
netcdf.putAtt(id, meanMeltRate_id, 'units', 'm second^-1');
disp('!!!! check for water equiv')

totalMeltFlux_id = netcdf.defVar(id, 'totalMeltFlux', 'double', nTime_dim);
netcdf.putAtt(id, totalMeltFlux_id, 'long_name', 'the total mass flux of freshwater across the ice-ocean interface, positive for melting and negative for freezing');
netcdf.putAtt(id, totalMeltFlux_id, 'units', 'kg second^-1');


totalOceanVolume_id = netcdf.defVar(id, 'totalOceanVolume', 'double', nTime_dim);
netcdf.putAtt(id, totalOceanVolume_id, 'long_name', 'the total volume of the ocean');
netcdf.putAtt(id, totalOceanVolume_id, 'units', 'm^3');

meanTemperature_id = netcdf.defVar(id, 'meanTemperature', 'double', nTime_dim);
netcdf.putAtt(id, meanTemperature_id, 'long_name', 'the potential temperature averaged over the ocean volume');
netcdf.putAtt(id, meanTemperature_id, 'units', 'Celsius');

meanSalinity_id = netcdf.defVar(id, 'meanSalinity', 'double', nTime_dim);
netcdf.putAtt(id, meanSalinity_id, 'long_name', 'the salinity averaged over the ocean volume');
netcdf.putAtt(id, meanSalinity_id, 'units', 'PSU');

iceDraft_id = netcdf.defVar(id, 'iceDraft', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, iceDraft_id, 'long_name', 'the elevation of the ice-ocean interface');
netcdf.putAtt(id, iceDraft_id, 'units', 'm');
netcdf.putAtt(id, iceDraft_id, 'missing_value', -1e34);
netcdf.putAtt(id, iceDraft_id, '_FillValue', -1e34);

bathymetry_id = netcdf.defVar(id, 'bathymetry', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bathymetry_id, 'long_name', 'the elevation of the bathymetry');
netcdf.putAtt(id, bathymetry_id, 'units', 'm');
netcdf.putAtt(id, bathymetry_id, 'missing_value', -1e34);
netcdf.putAtt(id, bathymetry_id, '_FillValue', -1e34);

meltRate_id = netcdf.defVar(id, 'meltRate', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, meltRate_id, 'long_name', 'the melt rate, positive for melting and negative for freezing');
netcdf.putAtt(id, meltRate_id, 'units', 'm s^-1');
netcdf.putAtt(id, meltRate_id, 'missing_value', -1e34);
netcdf.putAtt(id, meltRate_id, '_FillValue', -1e34);

frictionVelocity_id = netcdf.defVar(id, 'frictionVelocity', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, frictionVelocity_id, 'long_name', 'the friction velocity, used in melt calculations');
netcdf.putAtt(id, frictionVelocity_id, 'units', 'm s^-1');
netcdf.putAtt(id, frictionVelocity_id, 'missing_value', -1e34);
netcdf.putAtt(id, frictionVelocity_id, '_FillValue', -1e34);

thermalDriving_id = netcdf.defVar(id, 'thermalDriving', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, thermalDriving_id, 'long_name', 'the thermal driving used in melt calculations, calculated as the difference between theta in the boundary layer and freezing theta at the interface.');
netcdf.putAtt(id, thermalDriving_id, 'units', 'Celsius');
netcdf.putAtt(id, thermalDriving_id, 'missing_value', -1e34);
netcdf.putAtt(id, thermalDriving_id, '_FillValue', -1e34);

halineDriving_id = netcdf.defVar(id, 'halineDriving', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, halineDriving_id, 'long_name', 'the haline driving used in melt calculations, calculated as the difference between salinity in the boundary layer and salinity at the interface.');
netcdf.putAtt(id, halineDriving_id, 'units', 'PSU');
netcdf.putAtt(id, halineDriving_id, 'missing_value', -1e34);
netcdf.putAtt(id, halineDriving_id, '_FillValue', -1e34);


uBoundaryLayer_id = netcdf.defVar(id, 'uBoundaryLayer', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, uBoundaryLayer_id, 'long_name', 'the u-components of velocity in the boundary layer that were used to compute friction velocity');
netcdf.putAtt(id, uBoundaryLayer_id, 'units', 'm s^-1');
netcdf.putAtt(id, uBoundaryLayer_id, 'missing_value', -1e34);
netcdf.putAtt(id, uBoundaryLayer_id, '_FillValue', -1e34);

vBoundaryLayer_id = netcdf.defVar(id, 'vBoundaryLayer', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, vBoundaryLayer_id, 'long_name', 'the v-components of velocity in the boundary layer that were used to compute friction velocity');
netcdf.putAtt(id, vBoundaryLayer_id, 'units', 'm s^-1');
netcdf.putAtt(id, vBoundaryLayer_id, 'missing_value', -1e34);
netcdf.putAtt(id, vBoundaryLayer_id, '_FillValue', -1e34);

barotropicStreamfunction_id = netcdf.defVar(id, 'barotropicStreamfunction', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, barotropicStreamfunction_id, 'long_name', 'the barotropic streamfunction');
netcdf.putAtt(id, barotropicStreamfunction_id, 'units', 'm^3 s^-1');
netcdf.putAtt(id, barotropicStreamfunction_id, 'missing_value', -1e34);
netcdf.putAtt(id, barotropicStreamfunction_id, '_FillValue', -1e34);

overturningStreamfunction_id = netcdf.defVar(id, 'overturningStreamfunction', 'double', [nx_dim nz_dim nTime_dim]); 
netcdf.putAtt(id, overturningStreamfunction_id, 'long_name', 'the overturning streamfunction');
netcdf.putAtt(id, overturningStreamfunction_id, 'units', 'm^3 s^-1');
netcdf.putAtt(id, overturningStreamfunction_id, 'missing_value', -1e34);
netcdf.putAtt(id, overturningStreamfunction_id, '_FillValue', -1e34);

bottomTemperature_id = netcdf.defVar(id, 'bottomTemperature', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bottomTemperature_id, 'long_name', 'the potential temperature in the bottom-most cell in each ocean column');
netcdf.putAtt(id, bottomTemperature_id, 'units', 'Celsius');
netcdf.putAtt(id, bottomTemperature_id, 'missing_value', -1e34);
netcdf.putAtt(id, bottomTemperature_id, '_FillValue', -1e34);

bottomSalinity_id = netcdf.defVar(id, 'bottomSalinity', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bottomSalinity_id, 'long_name', 'the salinity in the bottom-most cell in each ocean column');
netcdf.putAtt(id, bottomSalinity_id, 'units', 'PSU');
netcdf.putAtt(id, bottomSalinity_id, 'missing_value', -1e34);
netcdf.putAtt(id, bottomSalinity_id, '_FillValue', -1e34);

temperatureXZ_id = netcdf.defVar(id, 'temperatureXZ', 'double', [nx_dim nz_dim nTime_dim]);
netcdf.putAtt(id, temperatureXZ_id, 'long_name', 'potential temperature transect in the x-z plane through the centre of the domain');
netcdf.putAtt(id, temperatureXZ_id, 'units', 'Celsius');
netcdf.putAtt(id, temperatureXZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, temperatureXZ_id, '_FillValue', -1e34);

salinityXZ_id = netcdf.defVar(id, 'salinityXZ', 'double', [nx_dim nz_dim nTime_dim]);
netcdf.putAtt(id, salinityXZ_id, 'long_name', 'salinity transect in the x-z plane through the centre of the domain, y=40km');
netcdf.putAtt(id, salinityXZ_id, 'units', 'PSU');
netcdf.putAtt(id, salinityXZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, salinityXZ_id, '_FillValue', -1e34);

temperatureYZ_id = netcdf.defVar(id, 'temperatureYZ', 'double', [ny_dim nz_dim nTime_dim]);
netcdf.putAtt(id, temperatureYZ_id, 'long_name', 'potential temperature transect in the y-z plane through the centre of the domain, x=520km');
netcdf.putAtt(id, temperatureYZ_id, 'units', 'Celsius');
netcdf.putAtt(id, temperatureYZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, temperatureYZ_id, '_FillValue', -1e34);

salinityYZ_id = netcdf.defVar(id, 'salinityYZ', 'double', [ny_dim nz_dim nTime_dim]);
netcdf.putAtt(id, salinityYZ_id, 'long_name', 'salinity transect in the y-z plane through the centre of the domain, x=520km');
netcdf.putAtt(id, salinityYZ_id, 'units', 'PSU');
netcdf.putAtt(id, salinityYZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, salinityYZ_id, '_FillValue', -1e34);

%add global data
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ISOMIP results: Regional Ocean Modelling System');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', outname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'RESULTS file');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ROMS');
netcdf.endDef(id);

%fill in variables
netcdf.putVar(id, x_id, x_rho(1,:));
netcdf.putVar(id, y_id, y_rho(:,1));
netcdf.putVar(id, z_id, flip(Z_interp)); %flip dim to start at -2.5 and go through to -717.5
netcdf.putVar(id, time_id, ocean_time);
netcdf.putVar(id, meanMeltRate_id, meanMeltRate);
netcdf.putVar(id, totalMeltFlux_id, meanMeltFlux);
netcdf.putVar(id, totalOceanVolume_id, sumOfVolume);
netcdf.putVar(id, meanTemperature_id, meanTemperature);
netcdf.putVar(id, meanSalinity_id, meanSalinity);
netcdf.putVar(id, iceDraft_id, permute(repmat(zice,[1 1 length(ocean_time)]),[2 1 3]));
netcdf.putVar(id, bathymetry_id, permute(repmat(h,[1 1 length(ocean_time)]),[2 1 3]));
netcdf.putVar(id, meltRate_id, permute(meltRate,[2 1 3]));
netcdf.putVar(id, frictionVelocity_id, permute(Ustar,[2 1 3]));
netcdf.putVar(id, thermalDriving_id, permute(Tstar,[2 1 3]));
netcdf.putVar(id, halineDriving_id, permute(Sstar,[2 1 3]));
netcdf.putVar(id, uBoundaryLayer_id, permute(v_top,[2 1 3]));
netcdf.putVar(id, vBoundaryLayer_id, permute(-u_top,[2 1 3]));
netcdf.putVar(id, barotropicStreamfunction_id, permute(baroSF,[2 1 3]));
netcdf.putVar(id, overturningStreamfunction_id, otSF);
netcdf.putVar(id, bottomTemperature_id, permute(T_bot,[2 1 3]));
netcdf.putVar(id, bottomSalinity_id, permute(S_bot,[2 1 3]));
netcdf.putVar(id, temperatureXZ_id, flip(permute(tempXZ,[2 1 3]),2));
netcdf.putVar(id, salinityXZ_id, flip(permute(saltXZ,[2 1 3]),2));
netcdf.putVar(id, temperatureYZ_id, flip(permute(tempYZ,[2 1 3]),2));
netcdf.putVar(id, salinityYZ_id, flip(permute(saltYZ,[2 1 3]),2));

netcdf.close(id);


