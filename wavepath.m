function wp=wavepath(seis,v,nbc,dx,nt,dt,s,sx,sz,gx,gz)
%  2DTDFD solver to calculate wavepath
%  IN   seis(:,:) -- seismogram,    v(:,:) -- velocity
%       nbc       -- grid number of boundary
%       dx        -- grid intervel, nt          -- number of sample
%       dt        -- time interval, s(:)        -- wavelet
%       sx,sz     -- src position,  gx(:),gz(:) -- rec position
%  OUT  wp(:,:)           -- Output wave path

[nz,nx]=size(v); ng=numel(gx); wp=zeros(nz,nx);

C=getC;
% setup ABC and temperary variables
v=padvel(v,nbc);
abc=Coef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*C(1)*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);

[isx,isz,igx,igz]=index(sx,sz,gx,gz,dx,nbc);

p1=zeros(size(v)); p0=zeros(size(v));

bc_top=zeros(5,nx,nt);
bc_bottom=zeros(5,nx,nt);
bc_left=zeros(nz,5,nt);
bc_right=zeros(nz,5,nt);

% Time Looping
for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (C(2)*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +C(3)*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +C(4)*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +C(5)*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    [bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it)]=save_boundary(p,nz,nx,nbc);
    p0=p1;
    p1=p;
end
% save final wavefield
bc_p_nt_1=p0;
bc_p_nt=p1;

bp1=bc_p_nt_1;
bp0=bc_p_nt;
q0=zeros(size(v)); q1=zeros(size(v));

% Time Loop
for it=nt-2:-1:1
    % subtrace source
    bp0(isz,isx)=bp0(isz,isx)-s(it+2)*beta_dt(isz,isx);
    bp=temp1.*bp1-temp2.*bp0+alpha.*...
        (C(2)*(circshift(bp1,[0,1,0])+circshift(bp1,[0,-1,0])+circshift(bp1,[1,0,0])+circshift(bp1,[-1,0,0]))...
        +C(3)*(circshift(bp1,[0,2,0])+circshift(bp1,[0,-2,0])+circshift(bp1,[2,0,0])+circshift(bp1,[-2,0,0]))...
        +C(4)*(circshift(bp1,[0,3,0])+circshift(bp1,[0,-3,0])+circshift(bp1,[3,0,0])+circshift(bp1,[-3,0,0]))...
        +C(5)*(circshift(bp1,[0,4,0])+circshift(bp1,[0,-4,0])+circshift(bp1,[4,0,0])+circshift(bp1,[-4,0,0])));
    bp=load_boundary(bp,bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it),nz,nx,nbc); 
    q=temp1.*q1-temp2.*q0+alpha.*...
        (C(2)*(circshift(q1,[0,1,0])+circshift(q1,[0,-1,0])+circshift(q1,[1,0,0])+circshift(q1,[-1,0,0]))...
        +C(3)*(circshift(q1,[0,2,0])+circshift(q1,[0,-2,0])+circshift(q1,[2,0,0])+circshift(q1,[-2,0,0]))...
        +C(4)*(circshift(q1,[0,3,0])+circshift(q1,[0,-3,0])+circshift(q1,[3,0,0])+circshift(q1,[-3,0,0]))...
        +C(5)*(circshift(q1,[0,4,0])+circshift(q1,[0,-4,0])+circshift(q1,[4,0,0])+circshift(q1,[-4,0,0])));
    % Add seismogram
    for ig=1:ng
        q(igz(ig),igx(ig))=q(igz(ig),igx(ig))+beta_dt(igz(ig),igx(ig))*seis(it,ig);
    end
    wp=image_condition(wp,bp1,q0,nz,nx,nbc);  
    % wf refresh
    bp0=bp1; bp1=bp;
    q0=q1; q1=q;
end
