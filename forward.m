function seis=forward(v,nbc,dx,nt,dt,s,sx,sz,gx,gz)
%  purpose:  2DTDFD solution to acoustic wave equation
%            use the absorbing boundary condition
%
%  IN   v(:,:) -- velocity,      nbc         -- grid number of boundary
%       dx     -- grid intervel, nt          -- number of sample
%       dt     -- time interval, s(:)        -- wavelet
%       sx,sz  -- src position,  gx(:),gz(:) -- rec position
%  OUT  seis(:,:)          -- Output seismogram
seis=zeros(nt,numel(gx));
ng=numel(gx); 
C=getC;
v=padvel(v,nbc);
abc=Coef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*C(1)*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);
[isx,isz,igx,igz]=index(sx,sz,gx,gz,dx,nbc);
p1=zeros(size(v)); p0=zeros(size(v));

% Time Looping
for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (C(2)*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +C(3)*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +C(4)*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +C(5)*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    for ig=1:ng
        seis(it,ig)=p(igz(ig),igx(ig));
    end
    p0=p1;
    p1=p;
end
