function [ FBsPCA_data, z,zcost,info, Timing, imager, image0, ZA]=hMRA_uniform(data_raw, data, sigma)
%%
% TODO: write docsnois
% Based on hMRA_uniform from
%       https://github.com/chaom1026/2DhMRA/blob/master/hMRA_uniform.m
% 
% Depends on ASPIRE and Manopt
% 

%hMRA_uniform: heterogeneous MRA with uniform distribution. The optimization problem does not consider the distribution.
%
%inputs:
%  Nclus: number of classes
%  Ncopy: number of observations in each class
%  SNR: signal-to-noise ratio
%  numcoef: number of coefficients taken by sPCA
%  shift_range: the maximum length of shifts 
%  file_name: address of the data file
%outputs:
%  FBsPCA_data: results of sPCA
%  z: solution of the optimization problem
%  zcost: objective value of the optimization problem
%  info: information from the optimization solver
%  imager: recovered images
%  image0: groundtruth images after sPCA
%  ZA: permutation between imager and image0


Nclus = 1;

%%============================Steerable PCA==========================
T=cputime;
L = size(data, 1);
r_max = floor(L/2);
P=size(data, 3);
N=floor(L/2);
  
[x, y]=meshgrid(-N:N, -N:N);    % generate square grids
r=sqrt(x.^2+y.^2);              % compute the radius at each grid point
test=reshape(data, L^2, P);
test=test(r>r_max, :);
noise_variance=var(test(:));    % estimate the noise variance from corner region.

% generate basis
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(data, r_max);

% determine sPCA basis
[ UU, Freqs, rad_Freqs, W ] = FBSPCA_MP_rankEst(P, U, D, freqs, rad_freqs, noise_variance);

% compute sPCA coefficients
[ Coeff ] = WF_FBSPCA( data, Mean, r_max, [0, 0], UU, Freqs, W, 0);
[ Coeff_raw ] = WF_FBSPCA( data_raw, Mean, r_max, [0, 0], UU, Freqs, W, 0);

Timing.t_spca = cputime-T;

% save sPCA results
FBsPCA_data.r_max = r_max;     
FBsPCA_data.UU = UU;                 % eigenimages by sPCA
FBsPCA_data.Freqs = Freqs;           % angular frequencies
FBsPCA_data.rad_Freqs = rad_freqs;
FBsPCA_data.Coeff = Coeff;           % sPCA coefficients
FBsPCA_data.Coeff_raw = Coeff_raw;   % sPCA coefficients of groundtruth images
FBsPCA_data.Mean = Mean;             % mean image
FBsPCA_data.W = W;                   % linear filter weight

FBsPCA_data.data_raw=data_raw(:,:,1:Nclus); % groundtruth images
FBsPCA_data.data_sample=data(:,:,1:10);     % samples of noisy observations

clear U D freqs rad_freqs Mean UU
'SPCA finished'
size(Coeff,1)


%%=============================Mean Values=================================
T=cputime;

M=zeros(100,1);
count=1;
id=find(Freqs==0);
l=length(id);
for i=1:l
    M(count,:)=id(i);
    count=count+1;
end
M=M(1:count-1,:);
Mean=Coeff(M,:);
mixed_m=mean(Mean,2);


%%===========================Power Spectrums===============================
P=zeros(50000,3);
count=1;
max_freq=max(Freqs);
for i=0:max_freq
    id=find(Freqs==i);
    l=length(id);
    for j=1:l
        for k=1:j
            P(count,:)=[id(j),id(k),k==j];
            count=count+1;
        end
    end
end
P=P(1:count-1,:);
Powerspec=Coeff(P(:,1),:).*conj(Coeff(P(:,2),:));
mixed_p=mean(Powerspec,2)-sigma^2*P(:,3);


%%===========================Bispectrums===================================
L=length(Freqs);

B=zeros(500000,4);
cm=size(mixed_m,1);
count=1;
max_freq=max(Freqs);
for i1=0:max_freq
    for i2=0:i1
        i3=i1+i2;
        id1=find(Freqs==i1);
        id2=find(Freqs==i2);
        id3=find(Freqs==i3);
        if ~isempty(id3)
            for k1=1:length(id1)
                for k2=1:length(id2)
                    for k3=1:length(id3)
                        try 
                            B(count,:)=[id1(k1),id2(k2),id3(k3),(i1==0)*(k2==k3)*mixed_m(min(k1,cm))+(i2==0)*(k1==k3)*mixed_m(min(k2,cm))+(i3==0)*(k1==k2)*conj(mixed_m(min(k3,cm)))];
                        catch
                            B(count,:)=[id1(k1),id2(k2),id3(k3),0];
                        end
                        count=count+1;
                    end
                end
            end
        end
    end
end
B=B(1:count-1,:);
Bispec=Coeff(B(:,1),:).*Coeff(B(:,2),:).*conj(Coeff(B(:,3),:));
mixed_b=mean(Bispec,2)-sigma^2*B(:,4);


% estimate mixed invariants
All=[mixed_b;mixed_p;mixed_m];
lb=size(B,1);
lp=size(P,1);
lm=size(M,1);
Timing.t_invariants = cputime-T;
clear i j k K N r x x0 y y0
delete(gcp('nocreate'));

%%=======================formulate optimization problem====================
lf=sum(Freqs==0);
mani.a=euclideanfactory(lf,Nclus);           % coefficients with frequency 0
mani.b=euclideancomplexfactory(L-lf,Nclus);  % coefficients with other frequencies
manifold=productmanifold(mani);
problem.M=manifold;

problem.cost=@costf;
problem.egrad=@gradf;

f=@(z)([z(B(:,1),:).*z(B(:,2),:).*conj(z(B(:,3),:));z(P(:,1),:).*conj(z(P(:,2),:));z(M,:)]);
g=@(z)(Nclus*All-sum(f([z.a;z.b]),2));

figure;
checkgradient(problem);

T=cputime;
option.maxinner=25;
[z, zcost, info] = conjugategradient(problem);
Timing.t_opt = cputime-T;

% align recovered results with groundtruth images
ZA=zeros(1,Nclus);
z=[z.a;z.b];
for I=1:Nclus
    zz=z(:,I);
    zze=zeros(size(zz,1),3600);
    for II=0:3599
        zze(:,II+1)=zz.*exp(1i*(II*2*pi/3600)*Freqs);
    end
    Align_dis=zeros(2,2*Nclus);
    for II=1:2*Nclus
        aln_dis=sum((abs(zze-Coeff_raw(:,II)*ones(1,3600))).^2,1);
        [dis,ind]=min(aln_dis);
        Align_dis(:,II)=[dis;ind];
    end
    [~,J]=min(Align_dis(1,:));
    ZA(I)=J;
    ang=Align_dis(2,J);
    zz=zz.*exp(1i*((ang-1)*2*pi/3600)*Freqs);
    z(:,I)=zz;
end

% recovered images and groundtruth images after sPCA
LL=size(data, 1);
N=floor(LL/2);
[x, y] = meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
imager=zeros(LL^2,Nclus);
UU=FBsPCA_data.UU;
tmp=real(UU(:,Freqs==0)*z(Freqs==0,:))+2*real(UU(:,Freqs~=0)*z(Freqs~=0,:));
imager(r<=r_max,:)=imager(r<=r_max,:)+tmp;
FMean=zeros(LL,LL,Nclus);
for I=1:Nclus
    FMean(:,:,I)=FBsPCA_data.Mean;
end
imager=reshape(imager,LL,LL,Nclus)+FMean;
for I=1:Nclus
    if ZA(I)>Nclus
        imager(:,:,I)=flipud(imager(:,:,I));
        ZA(I)=ZA(I)-Nclus;
    end
end

image0=zeros(LL^2,Nclus);
tmp=real(UU(:,Freqs==0)*Coeff_raw(Freqs==0,1:Nclus))+2*real(UU(:,Freqs~=0)*Coeff_raw(Freqs~=0,1:Nclus));
image0(r<=r_max,:)=image0(r<=r_max,:)+tmp;
image0=reshape(image0,LL,LL,Nclus)+FMean;


%%=============================sub-functions=================================

%     % cost function
%     function y =costf(z)
%         temp=g(z);
%         temp(1:lb)=temp(1:lb)/sqrt(lb);
%         temp(lb+1:lb+lp)=temp(lb+1:lb+lp)/sqrt(lp);
%         temp(lb+lp+1:lb+lp+lm)=temp(lb+lp+1:lb+lp+lm)/sqrt(lm);
%         y=temp'*temp;
%     end
% 
%     % gradient
%     function GD = gradf(z)
%         gg=g(z);
%         z=[z.a;z.b];
%         ggb=gg(1:lb);
%         ggp=gg(lb+1:lb+lp);
%         ggm=gg(lb+lp+1:lb+lp+lm);
%         gd=zeros(L,Nclus);
%         for ii=1:lb
%             a=B(ii,:);
%             if a(1)==a(2) && a(2)~=a(3)
%                 gd([a(1),a(3)],:)= gd([a(1),a(3)],:)+[-4*ggb(ii)*conj(z(a(1),:)).*z(a(3),:);-2*conj(ggb(ii))*z(a(1),:).^2]/lb;
% 	        elseif a(1)==a(3) && a(2)~=a(3)
% 		        gd([a(1),a(2)],:)= gd([a(1),a(2)],:)+[-2*conj(ggb(ii))*z(a(1),:).*z(a(2),:)-2*ggb(ii)*conj(z(a(2),:)).*z(a(1),:);-2*ggb(ii)*z(a(1),:).*conj(z(a(1),:))]/lb;
% 	        elseif a(2)==a(3) && a(1)~=a(3)
% 		        gd([a(1),a(2)],:)= gd([a(1),a(2)],:)+[-2*ggb(ii)*z(a(2),:).*conj(z(a(2),:));-2*conj(ggb(ii))*z(a(2),:).*z(a(1),:)-2*ggb(ii)*conj(z(a(1),:)).*z(a(2),:)]/lb;
% 	        elseif a(1)==a(2) && a(2)==a(3)
% 		        gd(a(1),:)= gd(a(1),:)+(-2)*[2*ggb(ii)*z(a(1),:).*conj(z(a(1),:))+conj(ggb(ii))*z(a(1),:).*z(a(1),:)]/lb;
%                 %[-6*ggb(ii)*z(a(1),:).*z(a(1),:)];
%             else
%                 gd([a(1),a(2),a(3)],:)= gd([a(1),a(2),a(3)],:)+...
%                     (-2)*[ggb(ii)*conj(z(a(2),:)).*z(a(3),:);ggb(ii)*conj(z(a(1),:)).*z(a(3),:);conj(ggb(ii))*z(a(1),:).*z(a(2),:)]/lb;
%             end
%         end
%         for ii=1:lp
%             a=P(ii,:);
%             if a(1)==a(2)
%                 gd(a(1),:)=gd(a(1),:)-4*real(ggp(ii))*z(a(1),:)/lp;
%             else
%                 gd([a(1),a(2)],:)=gd([a(1),a(2)],:)-2*[ggp(ii)*z(a(2),:);conj(ggp(ii))*z(a(1),:)]/lp;
%             end
%         end
%         for ii=1:lm
%             a=M(ii);
%             gd(a,:)=gd(a,:)-2*ggm(a)/lm;
%         end
%         GD.a=gd(1:lf,:);
%         GD.b=gd(lf+1:L,:);
%     end


 % cost function
    function y =costf(z)
        temp=g(z);
        temp(1:lb)=temp(1:lb)/sqrt((1+sigma^2+sigma^4));
        temp(lb+1:lb+lp)=temp(lb+1:lb+lp)/sqrt((1+sigma^2));
        temp(lb+lp+1:lb+lp+lm)=temp(lb+lp+1:lb+lp+lm)/sqrt(1);
        y=temp'*temp;
    end

    % gradient
    function GD = gradf(z)
        gg=g(z);
        z=[z.a;z.b];
        ggb=gg(1:lb);
        ggp=gg(lb+1:lb+lp);
        ggm=gg(lb+lp+1:lb+lp+lm);
        gd=zeros(L,Nclus);
        for ii=1:lb
            a=B(ii,:);
            if a(1)==a(2) && a(2)~=a(3)
                gd([a(1),a(3)],:)= gd([a(1),a(3)],:)+[-4*ggb(ii)*conj(z(a(1),:)).*z(a(3),:);-2*conj(ggb(ii))*z(a(1),:).^2]/(1+sigma^2+sigma^4);
	        elseif a(1)==a(3) && a(2)~=a(3)
		        gd([a(1),a(2)],:)= gd([a(1),a(2)],:)+[-2*conj(ggb(ii))*z(a(1),:).*z(a(2),:)-2*ggb(ii)*conj(z(a(2),:)).*z(a(1),:);-2*ggb(ii)*z(a(1),:).*conj(z(a(1),:))]/(1+sigma^2+sigma^4);
	        elseif a(2)==a(3) && a(1)~=a(3)
		        gd([a(1),a(2)],:)= gd([a(1),a(2)],:)+[-2*ggb(ii)*z(a(2),:).*conj(z(a(2),:));-2*conj(ggb(ii))*z(a(2),:).*z(a(1),:)-2*ggb(ii)*conj(z(a(1),:)).*z(a(2),:)]/(1+sigma^2+sigma^4);
	        elseif a(1)==a(2) && a(2)==a(3)
		        gd(a(1),:)= gd(a(1),:)+(-2)*[2*ggb(ii)*z(a(1),:).*conj(z(a(1),:))+conj(ggb(ii))*z(a(1),:).*z(a(1),:)]/(1+sigma^2+sigma^4);
                %[-6*ggb(ii)*z(a(1),:).*z(a(1),:)];
            else
                gd([a(1),a(2),a(3)],:)= gd([a(1),a(2),a(3)],:)+...
                    (-2)*[ggb(ii)*conj(z(a(2),:)).*z(a(3),:);ggb(ii)*conj(z(a(1),:)).*z(a(3),:);conj(ggb(ii))*z(a(1),:).*z(a(2),:)]/(1+sigma^2+sigma^4);
            end
        end
        for ii=1:lp
            a=P(ii,:);
            if a(1)==a(2)
                gd(a(1),:)=gd(a(1),:)-4*real(ggp(ii))*z(a(1),:)/(1+sigma^2);
            else
                gd([a(1),a(2)],:)=gd([a(1),a(2)],:)-2*[ggp(ii)*z(a(2),:);conj(ggp(ii))*z(a(1),:)]/(1+sigma^2);
            end
        end
        for ii=1:lm
            a=M(ii);
            gd(a,:)=gd(a,:)-2*ggm(a)/1;
        end
        GD.a=gd(1:lf,:);
        GD.b=gd(lf+1:L,:);
    end


end

