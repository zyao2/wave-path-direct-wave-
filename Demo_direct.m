% Wavepath 
% Theory:
% WP   =   P   •   Q
% P is the forward propagation field, and Q is the backward propagation field.

%Define a homogenerous velocity model, and the model parameters;
clear all
nz=81;nx=201;nt=2001;
vel=zeros(nz,nx)+1000; 
dx=5;dt=0.0005;
x = (0:nx-1)*dx; z = (0:nz-1)*dx;
gx=(0:2:(nx-1))*dx; gz=zeros(size(gx)); g=1:numel(gx);t=(0:nt-1)*dt;
%Define the source and receiver geometry;
rec_indx=80;
s_indx=50;
sx=s_indx*dx;sz = 200; 
%Setup FD parameters and source wavelet;

nbc=40; 
freq=20; 
s=ricker(freq,dt);
%Generate the synthetic data;
%seis=forward(vel,nbc,dx,nt,dt,s,sx,sz,gx,gz);
%trace=seis(:,rec_indx);
%save trace trace
load trace;

%Calculate and plot the wavepath for direct wave;
tic; 
wp=wavepath(trace,vel,nbc,dx,nt,dt,s,sx,sz,gx(rec_indx),gz(rec_indx)); 
toc;

imagesc(x,z,wp);colormap(gray);caxis([-10 10]);
xlabel('X (m)'); ylabel('Z (m)'); title('wave path for direct wave');

