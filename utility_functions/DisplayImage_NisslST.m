% function DisplayImage(im,EigInfo,ConvInfo,para)
%
% This function displays the different informations returned by function
% coherence_orientation. 
% INPUT : 
% - im : 2D image
% - EigInfo : 
% - ConvInfo : 
% - para : 
%    para.Step : stepsize for oreintation field display
%    para.scl : length of orientation segments
%
% OUTPUT : 
% - Figures will pop up 1) ground truth image, 2) convexhull, 3)
% orientation field, 4) anisotropoy map.
%
% Developer:  Wenxing Zhang, wenxing84@gmail.com, July 2014

function DisplayImage(im,EigInfo,ConvInfo,para,varargin)
 
if length(varargin)>0 % only plot the tensors close to cells in the image
    mask_near_cells_flag = varargin{1};
else
    mask_near_cells_flag = false;
end
%%%%%%%%% plots parameters
if isfield(para,'Step');    Step   =para.Step;    else Step   =20;   end
if isfield(para,'scl');     scl    =para.scl;     else scl  =10;    end
               
%%%%%%%%%%% Image Info and convex hull Info.
[nx,ny]    = size(im);
% im         = im/max(im(:)); %%%%% 
imconv     = ConvInfo.imconv;  

%%%%%%%%%%%%%% Distance function and Tensor info.
nu1 = EigInfo.nu1; w1 = EigInfo.w1;   
nu2 = EigInfo.nu2; w2 = EigInfo.w2;
J_rho = EigInfo.J_rho;

%%%%%%%%%%%%% cell center, direction and radii info.
A =zeros(nx,ny); A(1:Step:nx,1:Step:ny)=1;
% % % ROEYCOMMENTED ALL OF THESE 
[i,j] = find(imconv.*A>0);
Synt.Centre = [j,i]';
        
Centre= Synt.Centre; %%% Synthetic and real data with ground-truth center
if median(Centre(:))<1;  Centre=Centre*max(nx,ny);  end
Cent  = fix(Centre(1,:))*nx+fix(Centre(2,:))+1; 
Centrex = Centre(1,:);  Centrey = Centre(2,:); 

if mask_near_cells_flag
    bins_centers = [0:254]+0.5;
    counts = hist(im(:),bins_centers);
    [thresh,~] = otsuthresh(counts);
    im_bw = 1-im2bw(im,thresh);
    im_bw = imdilate(im_bw,strel('disk', round(23*(0.645/0.46))));%round(23*(0.645/0.97)))); % Using the same value as rho makes sense. &&& Change for different datasets
    mask_near_cells = double(im_bw); % a maked of enlarged cells
else 
    mask_near_cells = double(ones(size(im)));
end

clear EigInfo  ConvInfo  DistInfo  Synt

%%%%%%%%%%%%%%%%%%%%%%%%%%% displays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % figure; imshow(im,[]); colormap gray; caxis([0.1 .5]) ; title(['truth image',num2str(nx),'*',num2str(ny)]);
% % % ROEYCOMMENTED THIS figure; imshow(imconv,[]); colormap gray; title('convHull');hold on        

%% %%%%%%%%% comparisons the resutls
figure;   %%% orientations
a1=w1(:,:,2); a2=w1(:,:,1);  
%imshow(im,[]); % Roey: I commented this because I don't want to rescale the pixels intensities
imshow(im);
colormap gray;
% caxis([0.1 .5]);  % Roey: I commented this becase I don't want to rescale the color map
hold on;
P1=kron([1;1],Centrex)+scl*kron([1;-1],a1(Cent).*mask_near_cells(Cent));
P2=kron([1;1],Centrey)+scl*kron([1;-1],a2(Cent).*mask_near_cells(Cent));
P1short=kron([1;1],Centrex)+0.9*scl*kron([1;-1],a1(Cent).*mask_near_cells(Cent));
P2short=kron([1;1],Centrey)+0.9*scl*kron([1;-1],a2(Cent).*mask_near_cells(Cent));

   % plot(P1,P2,'b','linewidth',3);% title('tensor orientation');hold off
   % Plot it with colors
   for sI = 1:size(P1,2)
    theta(sI) = atand(-diff(P2(:,sI))/diff(P1(:,sI)));
   end
    theta(theta<0) = theta(theta<0)+180;
    vals = theta;
    vals = [vals, 0, 180];
    rgbTmp = vals2colormap(vals, hsv, []);
    rgbTmp(end-1:end,:) = [];
   
   for sI = 1:size(P1,2)
       plot(P1(:,sI),P2(:,sI),'k','linewidth',6*0.645/0.46); %Or for brainmaps: (6*0.65/0.46));% title('tensor orientation');hold off
       h = plot(P1short(:,sI),P2short(:,sI),'linewidth',4*0.645/0.46);% title('tensor orientation');hold off
       h.Color = rgbTmp(sI,:);
   end
   
   end
% % % % % % % figure; %%%% anisotropy
% % % % % % % aniso=sqrt(nu1./nu2); %ROEYCOMMENTED aniso(imconv<0.5)=1; 
% % % % % % % imshow(aniso,[]); 
% % % % % % % colorbar('location','eastoutside'); title('anisotropy')
% % % % % % % colormap jet; caxis([0.58 1.02])   
    





