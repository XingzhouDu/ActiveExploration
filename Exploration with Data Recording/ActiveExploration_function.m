function [PosSwm, xyhv_co_mat, xyROI, hvROI, Branch, Img_total_arr, ReconsData, all_end, pix_inc, frame_ct_n, mf_dir, dist_hvroi] = ActiveExploration_function(imgroi2_top, img_grey_side, img_grey_side_ini, ...
    len_x,len_y, ROI_len,ROI_hei, xyROI_p, hvROI_p, xyhv_cpre_mat, Branch_p, PosSwm_p, Region, Img_total_arr_p, ReconsData_p, pix_inc_p, num_frame, frame_ct,mf_dir_p)

% Copyright © 2023 DU Xingzhou.
% Licence: Apache license 2.0.
% 
% Input:
%  imgroi2_top - binarized top-view image with a resolution of 1400*900;
%  img_grey_side - greyscale side-view image with a resolution of 1400*900;
%  img_grey_side_ini - greyscale background side-view image with a resolution of 1400*900;
%  len_x, len_y - consts showing the dimension (in millimeter) of top view in x and y directions, respectively;
%  ROI_len - length of the dynamic ROI on top view in pixels;
%  ROI_hei - height of the dynamic ROI on side view in pixels;
%  xyROI_p, hvROI_p - dynamic ROI on top view and side view in last cycle;
%  xyhv_cpre_mat - position of the swarm in previous cycles;
%  Branch_p - branch map from last cycle;
%  PosSwm_p - position of the swarm in last cycle;
%  Region - region of exploration;
%  Img_total_arr_p - exploration map from last cycle;
%  ReconsData_p - data for reconstruction from last cycle;
%  pix_inc_p - incremental region of the swarm in last cycle;
%  num_frame - current number of image frames;
%  frame_ct - counter for pause;
%  mf_dir_p - yaw angle of the rotating magnetic field in last cycle. 
% 
% Output:
%  PosSwm - position of the swarm;
%  xyhv_co_mat - positions of the swarm till current frame;
%  xyROI, hvROI - dynamic ROI on top view and side view, respectively;
%  Branch - updated branch map;
%  Img_total_arr - updated exploration map;
%  ReconsData - updated data for reconstruction;
%  all_end - marker to complete exploration;
%  pix_inc - incremental region of the swarm;
%  frame_ct - updated counter for pause;
%  mf_dir_p - updataed yaw angle of the rotating magnetic field; 
%  dist_hvroi - moving direction of the ROI on side view.

% clearvars -except imgroi2_top img_grey_side img_grey_side_ini len_x len_y ROI_len ROI_hei xyROI_p hvROI_p xyhv_cpre_mat Branch_p PosSwm_p Region Img_total_arr_p ReconsData_p pix_inc_p num_frame frame_ct mf_dir_p
xy_cpre_mat = xyhv_cpre_mat(:,1:2);
hv_cpre_mat = xyhv_cpre_mat(:,3:4);
xy_cpre = xy_cpre_mat(end,:);
hv_cpre = hv_cpre_mat(end,:);
xyROI_p = round(xyROI_p);
all_end = 0;
PointType = 0;
POI = [];
if isempty(Img_total_arr_p)
    Img_total_arr_p = zeros(900,1400);
end
Img_total_arr = Img_total_arr_p;
ReconsData = ReconsData_p;
pix_inc = pix_inc_p;
mf_dir = mf_dir_p;

Img_total = Img_total_arr;
Img_total(Img_total~=0) = 1;
Img_total_p = Img_total_arr_p;
Img_total_p(Img_total_p~=0) = 1;

% TOP VIEW
img_top_p_o = imresize(Img_total_arr_p(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3)), [256, 256]);
mark_top_p = unique(img_top_p_o);
mark_top_p = sort(mark_top_p,'descend');
img_top_p = zeros(256,256);
for ii_mk = 1:length(mark_top_p)-1
    img_top_p_tmp = img_top_p_o;
    img_top_p_tmp(img_top_p_tmp~=mark_top_p(ii_mk))=0;
    img_top_p_tmp(img_top_p_tmp==mark_top_p(ii_mk))=1;
    if any(any(img_top_p_tmp&imgroi2_top)) && sum(sum(img_top_p_tmp))>50
        img_top_p = img_top_p_tmp;
        break;
    end
end
imlabel_top = bwlabel(imgroi2_top);                       %% find the largest area in the image, eliminate noises
stats_top = regionprops(logical(imlabel_top),'Area');
pos_top = regionprops(logical(imlabel_top),'Centroid');
area_top = cat(1,stats_top.Area);
[area_top_max_a,ii_t] = sort(area_top,'descend');
if ~isempty(ii_t)
    pos_top_max = pos_top(ii_t(1)).Centroid;
    img_area_top = zeros(256,256);
    area_top_max = area_top_max_a(1);
    for iii_t = 1:length(ii_t)
        pos_top_tmp = pos_top(ii_t(iii_t)).Centroid;
        pos_top_tmp = pos_top_tmp./256.*ROI_len;
        pos_top_c_tmp = [(xyROI_p(1)+xyROI_p(3))/2 + (pos_top_tmp(1)-ROI_len/2),...
            (xyROI_p(2)+xyROI_p(4))/2 + (pos_top_tmp(2)-ROI_len/2)];
        img_area_tmp = ismember(imlabel_top,ii_t(iii_t));
        if ~(pos_top_c_tmp(1)<=Region(1)-50 || pos_top_c_tmp(1)>=Region(3)+50 || ...
                pos_top_c_tmp(2)<=Region(2)-50 || pos_top_c_tmp(2)>=Region(4)+50) && ...
            (any(any(img_top_p & img_area_tmp))||length(ii_t)==1)
            pos_top_max = pos_top(ii_t(iii_t)).Centroid;
            img_area_top = ismember(imlabel_top,ii_t(iii_t));
            area_top_max = area_top_max_a(iii_t);
            break;
        end
    end
else
    pos_top_max = [128,128];
    img_area_top = imlabel_top;
    area_top_max = 0;
end
% BWDISTSC
bw = ~img_area_top;
aspect=[1 1 1];
shape=size(bw);
if(length(shape)==2) shape=[shape,1]; end
if(length(aspect)==2) aspect=[aspect,1]; end
bwXY=bw;
DXY=zeros(shape(1:2));
D1=zeros(shape(1:2));
xlower=repmat(Inf,shape(1:2));
xlower(1,find(bwXY(1,:)))=1;
for i=2:shape(1)
    xlower(i,:)=xlower(i-1,:);
    xlower(i,find(bwXY(i,:)))=i;
end
xupper=repmat(Inf,shape(1:2));
xupper(end,find(bwXY(end,:)))=shape(1);
for i=shape(1)-1:-1:1
    xupper(i,:)=xupper(i+1,:);
    xupper(i,find(bwXY(i,:)))=i;
end
idx=find(~bwXY); [sbx,~]=ind2sub(shape(1:2),idx);
DXY(idx)=aspect(1)^2*min((sbx-xlower(idx)).^2,(sbx-xupper(idx)).^2);
D1=repmat(Inf,shape(1:2));
p=shape(2);
for i=1:shape(2)
    d0=DXY(:,i);
    dtmp=d0+aspect(2)^2*(p-i)^2;
    L=D1(:,p)>dtmp;
    idx=find(L);
    D1(idx,p)=dtmp(L);
    map_lower=L;
    idx_lower=idx;
    for ii=p-1:-1:1
        dtmp=d0(idx_lower)+aspect(2)^2*(ii-i)^2;
        L=D1(idx_lower,ii)>dtmp;
        D1(idx_lower(L),ii)=dtmp(L);
        map_lower(idx_lower)=L;
        idx_lower=idx_lower(L);
        if(isempty(idx_lower)) break; end
    end
end
edtImg_top=sqrt(D1);
% BWDISTSC END
edtImg_top = 2 .* edtImg_top;
xy_c_roi = pos_top_max;
% SKELETONIZE IMAGE
threshold_skl = round(max(max(edtImg_top)));
skeletonizedImage_top = bwmorph(img_area_top, 'thin', inf);
filter = [1 1 1;1 0 1;1 1 1];
Img_disconnect = skeletonizedImage_top & ~(skeletonizedImage_top & conv2(double(skeletonizedImage_top), filter, 'same')>2);
cc = bwconncomp(Img_disconnect);
numPixels = cellfun(@numel,cc.PixelIdxList);
[sorted_px, ind] = sort(numPixels);
%Remove components shorter than threshold
for iii=ind(sorted_px<threshold_skl)
    cur_comp = cc.PixelIdxList{iii};
    skeletonizedImage_top(cur_comp) = 0;
    full_cc = bwconncomp(skeletonizedImage_top);
    if full_cc.NumObjects>1
        skeletonizedImage_top(cur_comp) = 1;
    end
end
skeletonizedImage_top = bwmorph(skeletonizedImage_top, 'spur');
% SKELETONIZE IMAGE END
skel_ends = bwmorph(skeletonizedImage_top,'endpoints');
[EndPoint(:,2),EndPoint(:,1)] = find(skel_ends);
skel_bran = bwmorph(skeletonizedImage_top,'branchpoints');
if (max(max(skel_bran))~=0 && size(EndPoint,1)>=3)
    [BranPoint(:,2),BranPoint(:,1)] = find(skel_bran);
    ind_ep = ones(size(EndPoint,1),1);
    for i_br = 1:size(BranPoint,1)
        ind_ep = (vecnorm((BranPoint(i_br,:)-EndPoint)')>10)' & ind_ep;
    end
    EndPoint = EndPoint(ind_ep,:);
end
diamImg_top = edtImg_top .* double(skeletonizedImage_top);
[SkelPoint_top(:,2),SkelPoint_top(:,1)] = find(diamImg_top>0);
peri_skel_top = regionprops(logical(diamImg_top),'Perimeter');
if ~isempty(peri_skel_top)
    perimeter_skel_top = peri_skel_top.Perimeter;
    perimeter_skel_top = perimeter_skel_top./2;
else
    perimeter_skel_top = 0;
end
% TOP VIEW END

% SIDE VIEW
AW_side = imcrop(img_grey_side,[hvROI_p(1), hvROI_p(2), hvROI_p(3)-hvROI_p(1), hvROI_p(4)-hvROI_p(2)]);
AW_side = imresize(AW_side, [400, 128]);
imgroi2_side = im2double(AW_side);
AW_side_ini = imcrop(img_grey_side_ini,[hvROI_p(1), hvROI_p(2), hvROI_p(3)-hvROI_p(1), hvROI_p(4)-hvROI_p(2)]);
AW_side_ini = imresize(AW_side_ini, [400, 128]);
imgroi2_side_ini = im2double(AW_side_ini);
imgroi2_side_dev = abs(imgroi2_side-imgroi2_side_ini);
thresh_side = 5;
imgroi2_side_dev(imgroi2_side_dev<thresh_side)=0;
imgroi2_side_dev(imgroi2_side_dev>=thresh_side)=255;
imgroi2_side_dev = imfill(imgroi2_side_dev,'holes');
imgroi2_side = imgroi2_side_dev;

imlabel_side = bwlabel(imgroi2_side);                       %% find the largest area in the image, eliminate noises
stats_side = regionprops(logical(imlabel_side),'Area');
pos_side = regionprops(logical(imlabel_side),'Centroid');
area_side = cat(1,stats_side.Area);
[area_side_max,ii_s] = sort(area_side,'descend');
if ~isempty(ii_s)
    pos_side_max = pos_side(ii_s(1)).Centroid;
    img_area_side = ismember(imlabel_side,ii_s(1));
    dst_pos_side_min = 999;
    for iii_s = 1:min(length(ii_s),3)
        pos_side_tmp = pos_side(ii_s(iii_s)).Centroid;
        dst_pos_side = abs( ((hvROI_p(2)+hvROI_p(4))/2 + (pos_side_tmp(2)/400*ROI_hei-ROI_hei/2))-hv_cpre(2) );
        if hv_cpre(2)~=9999 &&  dst_pos_side<dst_pos_side_min
            dst_pos_side_min = dst_pos_side;
            pos_side_max = pos_side(ii_s(iii_s)).Centroid;
            img_area_side = ismember(imlabel_side,ii_s(iii_s));
        end
    end % REVISED HERE MAYBE GOOD FOR EXP
else
    pos_side_max = [0,0];
    img_area_side = imlabel_side;
    area_side_max = 0;
end
hv_c_roi = pos_side_max;
hv_area = max(area_side_max);
% SIDE VIEW END

if (max(area_top_max)>300 && max(area_top_max)<256*256*0.75 && ~frame_ct)
    p_x = len_x/1400;
    p_y = len_y/900;
    barr = 500;
    
    vroi_c = 900;
    
    if (xy_cpre(1)==9999 || xy_cpre(2)==9999)
        img_inc_top = zeros(256,256);
    else
        img_inc_top = max( (img_area_top - imresize(Img_total(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3)), [256, 256])), zeros(256,256) );
    end
    img_inc_top = imopen(img_inc_top,strel('disk',1));
    imlabel_inc_top = bwlabel(img_inc_top);
    stats_inc_top = regionprops(logical(imlabel_inc_top),'Area');
    pos_inc_top = regionprops(logical(imlabel_inc_top),'Centroid');
    area_inc_top = cat(1,stats_inc_top.Area);
    [max_area_inc,ii_inc] = sort(area_inc_top,'descend');
    ii_inc = ii_inc(max_area_inc>10);
    if ~isempty(max_area_inc)
        pix_inc = [pix_inc,max_area_inc(1)];
    end
    pos_inc_top = cat(1,pos_inc_top.Centroid);
    pos_inc_top = pos_inc_top(ii_inc,:);
    pos_inc_top = pos_inc_top./256.*ROI_len;
    if ~isempty(pos_inc_top)
        pos_inc_top = [(xyROI_p(1)+xyROI_p(3))/2 + (pos_inc_top(:,1)-ROI_len/2), (xyROI_p(2)+xyROI_p(4))/2 + (pos_inc_top(:,2)-ROI_len/2)];
    else
        pos_inc_top = [(xyROI_p(1)+xyROI_p(3))/2, (xyROI_p(2)+xyROI_p(4))/2];
    end
    
    xy_c_roi = xy_c_roi./256.*ROI_len;
    
    x_c = (xyROI_p(1)+xyROI_p(3))/2 + (xy_c_roi(1)-ROI_len/2);
    y_c = (xyROI_p(2)+xyROI_p(4))/2 + (xy_c_roi(2)-ROI_len/2);
    
    if ( (abs(x_c-xy_cpre(1))<=barr && abs(y_c-xy_cpre(2))<=barr) || xy_cpre(1)==9999 || xy_cpre(2)==9999 ) 
        xy_co = [x_c, y_c];
    else
        xy_co = xy_cpre;
    end
    
    Img_total(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3)) = imresize(img_area_top, [ROI_len+1, ROI_len+1])|Img_total(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3));
    Img_total_arr(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3)) = (Img_total(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3))-Img_total_p(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3)))...
        .*num_frame + Img_total_arr(xyROI_p(2):xyROI_p(4),xyROI_p(1):xyROI_p(3));
    
    % MAPPING
    xy_cpre_onpath = xy_cpre;
    for ii_pre = max(1,size(xy_cpre_mat,1)-10):-1:max(1,size(xy_cpre_mat,1)-20)
        xy_cpre_tmp = xy_cpre_mat(ii_pre,:);
        if (norm(xy_cpre_tmp-xy_co)>ROI_len/5 && norm(xy_cpre_tmp-xy_co)<ROI_len && xy_cpre_tmp(1)~=9999)
            xy_cpre_onpath = xy_cpre_tmp;
            break;
        end
    end
    EndPoint_c = EndPoint./256.*ROI_len;
    EndPoint_c(:,1) = (xyROI_p(1)+xyROI_p(3))/2 + (EndPoint_c(:,1)-ROI_len/2);
    EndPoint_c(:,2) = (xyROI_p(2)+xyROI_p(4))/2 + (EndPoint_c(:,2)-ROI_len/2);
    Branch = Branch_p;
    if (max(max(skel_bran))~=0 && size(EndPoint,1)>=3)
        % DISCOVER BRANCH AND POTENTIAL DIRECTIONS
        BranPoint_c = BranPoint./256.*ROI_len;
        BranPoint_c = [(xyROI_p(1)+xyROI_p(3))/2 + (BranPoint_c(:,1)-ROI_len/2), ...
            (xyROI_p(2)+xyROI_p(4))/2 + (BranPoint_c(:,2)-ROI_len/2)]; 
        dist_bc = 9999; num_bc = 1;
        for i_brc = 1:size(BranPoint_c,1)
            if norm(BranPoint_c(i_brc,:)-xy_cpre_onpath)<dist_bc
                dist_bc = norm(BranPoint_c(i_brc,:)-xy_cpre_onpath); num_bc = i_brc;
            end
        end
        BranPoint_cn = BranPoint_c(num_bc,:);
        dist_ep = 9999; num_ep = 1;
        dist_ep_tmp = 9999; num_ep_tmp = 1;
        for i_ep = 1:size(EndPoint_c,1)
            dist_ep_i = norm(EndPoint_c(i_ep,:)-xy_cpre_onpath);
            if norm(EndPoint_c(i_ep,:)-BranPoint_cn)>(256/2-30)/256*ROI_len
                dist_ep = dist_ep_i; num_ep = i_ep;
            end
            if dist_ep_i<dist_ep_tmp
                if norm(EndPoint_c(i_ep,:)-BranPoint_cn)>(256/2-30)/256*ROI_len
                    dist_ep = dist_ep_i; num_ep = i_ep;
                end
                dist_ep_tmp = dist_ep_i; num_ep_tmp = i_ep;
            end
        end
        if dist_ep == 9999
            num_ep = num_ep_tmp;
        end
        PotDir = EndPoint_c-BranPoint_cn;
        % DISCOVER BRANCH AND POTENTIAL DIRECTIONS END
        % REFRESH BRANCHING POINT AND DIRECTIONS
        if ~isempty(Branch_p)
            ind_branpnt_recd = find(~Branch_p(:,3));
            flag_cover = 0;
            for i_brp = 1:length(ind_branpnt_recd)
                if norm(Branch_p(ind_branpnt_recd(i_brp),1:2)-BranPoint_cn)<ROI_len/4
                    Branch(ind_branpnt_recd(i_brp),1:2) = (Branch_p(ind_branpnt_recd(i_brp),1:2)+BranPoint_cn)./2;
                    Branch(ind_branpnt_recd(i_brp),5) = Branch_p(ind_branpnt_recd(i_brp),5)+1;
                    ind_bran_cover = ind_branpnt_recd(i_brp);
                    flag_cover = 1;
                    break;
                end
            end
            if flag_cover
                ind_potdir_onbran = find(Branch_p(:,4)>(Branch_p(ind_bran_cover,4)*10) & Branch_p(:,4)<((Branch_p(ind_bran_cover,4)+1)*10) ...
                    & Branch_p(:,3)==1);
                for i_potdir = 1:size(PotDir,1)
                    ang_potdirs = 9999;
                    for i_potdir_recd = 1:length(ind_potdir_onbran)
                        ang_tmp = acos( dot(Branch_p(ind_potdir_onbran(i_potdir_recd),1:2),PotDir(i_potdir,:))...
                            /(norm(Branch_p(ind_potdir_onbran(i_potdir_recd),1:2))*norm(PotDir(i_potdir,:))) );
                        if ang_tmp < ang_potdirs
                            ang_potdirs = ang_tmp; ind_closedir = ind_potdir_onbran(i_potdir_recd);
                        end
                    end
                    if ang_potdirs<=(30/180*pi)
                        Branch(ind_closedir,1:2) = (Branch_p(ind_closedir,1:2)+PotDir(i_potdir,:))./2;
                        Branch(ind_closedir,6) = ( norm(Branch_p(ind_closedir,1:2)-PotDir(i_potdir,1:2))...
                            + Branch_p(ind_closedir,6)*Branch(ind_closedir,5) )/(Branch(ind_closedir,5)+1); % RECORD THE SPEED
                        Branch(ind_closedir,5) = Branch(ind_closedir,5)+1;
                    else
                        Branch_addnewdir = [PotDir(i_potdir,1:2),1,max(Branch_p(ind_potdir_onbran,4))+1,1,0,0 ]; 
                        Branch = [Branch_p(1:ind_potdir_onbran(end),:); Branch_addnewdir; Branch_p(ind_potdir_onbran(end)+1:end,:)];
                    end
                end
            else
                % ADD NEW BRANCH TO TREE
                curr_bran = find(Branch_p(:,7)==1 & Branch_p(:,3)==1);
                [~,ind_dir] = max(Branch_p(curr_bran,4));
                curr_bran = curr_bran(ind_dir);
                if isempty(curr_bran)
                    Branch = [BranPoint_cn,0,0,1,0,0];
                    Branch_tmp = zeros(size(PotDir,1),7);
                    for i_ep = 1:size(PotDir,1)
                        if (i_ep~= num_ep)
                            Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep, 1,0,0];
                        else
                            Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep, 1,0,2];
                        end
                    end
                    Branch = [Branch;Branch_tmp];
                    Branch_p = Branch;
                else
                    if ~any(Branch_p(:,4)==Branch_p(curr_bran,4) & Branch_p(:,3)==0)
                        Branch = [Branch_p; BranPoint_cn,0,Branch_p(curr_bran,4),1,0,1];
                        Branch_tmp = zeros(size(PotDir,1),7);
                        for i_ep = 1:size(PotDir,1)
                            if (i_ep~= num_ep)
                                Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep+Branch_p(curr_bran,4)*10, 1,0,0];
                            else
                                Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep+Branch_p(curr_bran,4)*10, 1,0,2];
                            end
                        end
                        Branch = [Branch;Branch_tmp];
                    elseif Branch_p(Branch_p(:,4)==Branch_p(curr_bran,4)&Branch_p(:,3)==0,5)<=3
                        Branch = Branch_p;
                        Branch(Branch(:,4)==Branch(curr_bran,4)&Branch(:,3)==0,:) = [BranPoint_cn,0,Branch_p(curr_bran,4),1,0,1];
                        ind_rmv = find(Branch(:,4)>Branch(curr_bran,4)*10 & Branch(:,4)<(Branch(curr_bran,4)+1)*10 & Branch(:,3)==1);
                        Branch(ind_rmv,:) = [];
                        Branch_tmp = zeros(size(PotDir,1),7);
                        for i_ep = 1:size(PotDir,1)
                            if (i_ep~= num_ep)
                                Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep+Branch_p(curr_bran,4)*10, 1,0,0];
                            else
                                Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep+Branch_p(curr_bran,4)*10, 1,0,2];
                            end
                        end
                        Branch = [Branch(1:min(ind_rmv)-1,:);Branch_tmp;Branch(min(ind_rmv):end,:)];
                    else
                        % DO NOTHING
                    end
                end
                % ADD NEW BRANCH TO TREE END
            end
        else
            Branch = [BranPoint_cn,0,0,1,0,0];
            Branch_tmp = zeros(size(PotDir,1),7);
            for i_ep = 1:size(PotDir,1)
                if (i_ep~= num_ep)
                    Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep, 1,0,0]; 
                else
                    Branch_tmp(i_ep,:) = [PotDir(i_ep,:),1, i_ep, 1,0,2];
                end
            end
            Branch = [Branch;Branch_tmp];
        end
        % REFRESH BRANCHING POINT AND DIRECTIONS END
    elseif (xy_cpre(1)==9999 || xy_cpre(2)==9999) && max(max(skel_bran))==0
        PotDir = [0,0];
    else
        PotDir = xy_co-xy_cpre;
        if abs( mod(atan2(-PotDir(2),PotDir(1))/pi*180+360,360)-mod(mf_dir_p+360,360) )>=160
            PotDir = [0,0];
        end
        if (xy_cpre_onpath~=xy_cpre)
            [~,ind_EPfrn] = max(vecnorm(EndPoint_c'-xy_cpre_onpath'));
        else
            [~,ind_EPfrn] = min(vecnorm(EndPoint'-[128;128]));
        end
        if all((EndPoint(ind_EPfrn,1)<40|EndPoint(ind_EPfrn,1)>216) | (EndPoint(ind_EPfrn,2)<40|EndPoint(ind_EPfrn,2)>216))...
                && perimeter_skel_top>40 && perimeter_skel_top>50
            PotDir_EPdev = ([EndPoint(ind_EPfrn,1),EndPoint(ind_EPfrn,2)]-[128,128])./256*ROI_len .*0.6;
        else
            PotDir_EPdev = [0,0];
        end
    end
    % MAPPING END
    
    % DECISION MAKING
    if (xy_cpre(1)==9999 || xy_cpre(2)==9999)
        xyROI = xyROI_p;
    elseif (xy_co(1)<=Region(1) || xy_co(1)>=Region(3) || xy_co(2)<=Region(2) || xy_co(2)>=Region(4)) ...
            || (length(pix_inc)>20&&all(pix_inc(max(length(pix_inc)-20,1):end)<10)&&size(EndPoint,1)<3&&max(max(skel_bran))==0) % OR SPEEED IS ZERO
        pix_inc = [];
        ind_currdir = find(Branch(:,7)==1 & Branch(:,3)==1);
        [~,ind_dir] = max(Branch(ind_currdir,4));
        ind_currdir = ind_currdir(ind_dir);
        ind_bran_rmv = find((Branch(:,3)==0 & Branch(:,5)<=3) | (Branch(:,3)==0 & Branch(:,4)==ind_currdir));
        ind_rmv = zeros(size(Branch,1),1);
        for ii_rmv = 1:length(ind_bran_rmv)
            ind_rmv_tmp = (Branch(:,3)==0 & Branch(:,4)==Branch(ind_bran_rmv(ii_rmv),4)) ...
                | (Branch(:,3)==1 & Branch(:,4)>Branch(ind_bran_rmv(ii_rmv),4)*10 & Branch(:,4)<(Branch(ind_bran_rmv(ii_rmv),4)+1)*10);
            ind_rmv = ind_rmv|ind_rmv_tmp;
        end
        if any(ind_rmv)
            Branch(ind_rmv,:) = [];
        end
        ind_currdir = find(Branch(:,7)==1 & Branch(:,3)==1);
        [~,ind_dir] = max(Branch(ind_currdir,4));
        ind_currdir = ind_currdir(ind_dir);
        ind_currdir_rec = ind_currdir;
        Branch(ind_currdir,7) = 2;
        while ~isempty(ind_currdir)
            curr_bran = find(Branch(:,3)==0 & Branch(:,4)==floor(Branch(ind_currdir,4)/10));
            ind_potdir_onbran = find(Branch(:,4)>(Branch(curr_bran,4)*10) & Branch(:,4)<((Branch(curr_bran,4)+1)*10) ...
                & Branch(:,3)==1 & Branch(:,7)~=2);
            if ~isempty(ind_potdir_onbran)
                ind_rmv = Branch(ind_potdir_onbran,5)<4 & Branch(ind_potdir_onbran,7)~=2;
                Branch(ind_potdir_onbran(ind_rmv),:) = [];
                ind_potdir_onbran = find(Branch(:,4)>(Branch(curr_bran,4)*10) & Branch(:,4)<((Branch(curr_bran,4)+1)*10) ...
                    & Branch(:,3)==1 & Branch(:,7)~=2);
            end
            if ~isempty(ind_potdir_onbran)
                [~,ind_dir] = max(Branch(ind_potdir_onbran,6));
                Branch(ind_potdir_onbran(ind_dir),7) = 1;
                POI = Branch(curr_bran,1:2)+Branch(ind_potdir_onbran(ind_dir),1:2);
                mf_dir = atan2( -Branch(ind_potdir_onbran(ind_dir),2), Branch(ind_potdir_onbran(ind_dir),1) )/pi*180;
                frame_ct = 3;
                break;
            else
                Branch(curr_bran,7) = 2;
                ind_currdir = find(Branch(:,4)==Branch(curr_bran,4) & Branch(:,3)==1 & Branch(:,7)==1);
                Branch(ind_currdir,7) = 2;
            end
        end
        if isempty(ind_currdir)
            all_end = 1;
        end
        PointType = 2;
    elseif (max(max(skel_bran))~=0 && size(EndPoint,1)>=3)
        curr_bran = find(Branch(:,7)==1 & ~Branch(:,3));
        if isempty(curr_bran)
            curr_bran = 1;
        end
        curr_bran = max(curr_bran);
        ind_potdir_onbran = find(Branch(:,4)>(Branch(curr_bran,4)*10) & Branch(:,4)<((Branch(curr_bran,4)+1)*10) ...
            & Branch(:,3)==1 & Branch(:,7)~=2);
        if Branch(curr_bran,5)>8 && ~any(Branch(ind_potdir_onbran,7)==1)
            [~,ind_dir] = max(Branch(ind_potdir_onbran,6));  
            Branch(ind_potdir_onbran(ind_dir),7) = 1;
            Branch(curr_bran,7) = 1;
            for ii = 1:size(pos_inc_top,1)
                dir_tmp = pos_inc_top(ii,:)-Branch(curr_bran,1:2);
                if acos( dot(dir_tmp, Branch(ind_potdir_onbran(ind_dir),1:2))/...
                        (norm(dir_tmp)*norm(Branch(ind_potdir_onbran(ind_dir),1:2))) )/pi*180 < 10
                    POI = pos_inc_top(ii,:);
                    mf_dir = atan2( -dir_tmp(2), dir_tmp(1) )/pi*180;
                    break;
                end
            end
        elseif any(Branch(ind_potdir_onbran,7)==1)
            for ii = 1:size(pos_inc_top,1)
                dir_tmp = pos_inc_top(ii,:)-Branch(curr_bran,1:2);
                if acos( dot(dir_tmp, Branch(ind_potdir_onbran(Branch(ind_potdir_onbran,7)==1),1:2))/...
                        (norm(dir_tmp)*norm(Branch(ind_potdir_onbran(Branch(ind_potdir_onbran,7)==1),1:2))) )/pi*180 < 10
                    POI = pos_inc_top(ii,:);
                    mf_dir = atan2( -dir_tmp(2), dir_tmp(1) )/pi*180;
                    break;
                end
            end
        else
            POI = Branch(curr_bran,1:2);
            PointType = 1;
        end
    else
        if ~isempty(max_area_inc)
            POI = [(xyROI_p(1)+xyROI_p(3))/2, (xyROI_p(2)+xyROI_p(4))/2] + min(min(max_area_inc(1)/(threshold_skl*10),perimeter_skel_top/80),1).*PotDir ... 
                + PotDir_EPdev; 
%             POI = [(xyROI_p(1)+xyROI_p(3))/2, (xyROI_p(2)+xyROI_p(4))/2] + PotDir; 
            mf_dir = atan2( -(POI(2)-(xyROI_p(2)+xyROI_p(4))/2-PotDir_EPdev(2)), POI(1)-(xyROI_p(1)+xyROI_p(3))/2-PotDir_EPdev(1) )/pi*180;
        end
    end
    
    if ~isempty(POI)
        xyROI = [POI(1)-ROI_len/2, POI(2)-ROI_len/2, POI(1)+ROI_len/2, POI(2)+ROI_len/2]; 
    else
        xyROI = xyROI_p;
    end
    % DECISION MAKING END
    
    x = (xy_co(1)-1400/2)*p_x;
    y = -(xy_co(2)-900/2)*p_y;
    
    hroi_c = 700+1704.348*( x/(376+y) );
    hvROI = [hroi_c-ROI_len/2, vroi_c-ROI_hei, hroi_c+ROI_len/2, vroi_c];
    
    h_c_roi = hv_c_roi(1)/128*ROI_len;
    v_c_roi = hv_c_roi(2)/400*ROI_hei;
    h_c = (hvROI_p(1)+hvROI_p(3))/2 + (h_c_roi-ROI_len/2);
    v_c = (hvROI_p(2)+hvROI_p(4))/2 + (v_c_roi-ROI_hei/2);
    if ( abs(v_c-hv_cpre(2))<=barr || hv_cpre(2)==9999 || hv_area>200)
        hv_co = [h_c, v_c];
    else
        hv_co = hv_cpre;
    end
    z = 124 - ( (16/26325) * (hv_co(2)-450) * (376+y) );
    
    PosSwm = [x,y,z,PointType];
    
    if (size(xy_cpre_mat,1)==1 && xy_cpre_mat(1,1)==9999)
        xy_co_mat = xy_co;
    else
        xy_co_mat = [xy_cpre_mat;xy_co];
    end
    
    if (size(hv_cpre_mat,1)==1 && hv_cpre_mat(1,1)==9999)
        hv_co_mat = hv_co;
    else
        hv_co_mat = [hv_cpre_mat;hv_co];
    end
    
    % DATA RECORDING
    % GET THE LOCATIONS AND RADIUS
    if (max(max(skel_bran))~=0 && size(EndPoint,1)>=3) && all_end~=1
        [SkelPoint_side_tmp,~] = find(img_area_side(:,max(1,round(BranPoint(num_bc,1)/2)-5):min(round(BranPoint(num_bc,1)/2)+5,128))); 
        if ~isempty(SkelPoint_side_tmp)
            SkelPoint_side_tmp = mean(SkelPoint_side_tmp);
            SampPoint_pix = [BranPoint(num_bc,1:2), SkelPoint_side_tmp, diamImg_top(BranPoint(num_bc,2),BranPoint(num_bc,1)),0, 0]; 
        else
            SampPoint_pix = [BranPoint(num_bc,1:2), 0, diamImg_top(BranPoint(num_bc,2),BranPoint(num_bc,1)),0, 0];
        end
        SampPoint_tmp = zeros(size(EndPoint,1),6);
        for ii_tmp = 1:size(EndPoint,1)
            [SkelPoint_side_tmp,~] = find(img_area_side(:,max(1,round(EndPoint(ii_tmp,1)/2)-5):min(round(EndPoint(ii_tmp,1)/2)+5,128))); 
            if ~isempty(SkelPoint_side_tmp)
                SkelPoint_side_tmp = mean(SkelPoint_side_tmp);
                SampPoint_tmp(ii_tmp,:) = [EndPoint(ii_tmp,:), SkelPoint_side_tmp, ...
                    max(max(diamImg_top(max(1,EndPoint(ii_tmp,2)-10):min(EndPoint(ii_tmp,2)+10,256),max(1,EndPoint(ii_tmp,1)-10):min(EndPoint(ii_tmp,1)+10,256)))),1, 0]; 
            else
                SampPoint_tmp(ii_tmp,:) = [EndPoint(ii_tmp,:), 0, ...
                    max(max(diamImg_top(max(1,EndPoint(ii_tmp,2)-10):min(EndPoint(ii_tmp,2)+10,256),max(1,EndPoint(ii_tmp,1)-10):min(EndPoint(ii_tmp,1)+10,256)))),1, 0];
            end
        end
        SampPoint_pix = [SampPoint_pix;SampPoint_tmp];
        SampPoint(:,1) = ( (xyROI_p(1)+xyROI_p(3))/2 + (SampPoint_pix(:,1)./256.*ROI_len-ROI_len/2) - 1400/2 ).*p_x;
    	SampPoint(:,2) = -( (xyROI_p(2)+xyROI_p(4))/2 + (SampPoint_pix(:,2)./256.*ROI_len-ROI_len/2) - 900/2 ).*p_y;
        SampPoint(:,3) = 124 - ( (16/26325) * ((hvROI_p(2)+hvROI_p(4))/2 + (SampPoint_pix(:,3)./400.*ROI_hei-ROI_hei/2) - 450) .* (376+SampPoint(:,2)) );
        SampPoint(:,4) = SampPoint_pix(:,4)./256.*ROI_len.*p_x;
        SampPoint(:,5:6) = SampPoint_pix(:,5:6);
        if ~isempty(SampPoint(:,3)==0) || any(isnan(SampPoint(:,3))) || any(abs(SampPoint(:,3)-z)>5)
                SampPoint(SampPoint(:,3)==0|isnan(SampPoint(:,3))|abs(SampPoint(:,3)-z)>5,3) = z;
        end
        if ~isempty(ReconsData_p)
            if ~isempty(Branch_p) && flag_cover 
                SampPoint(1,6) = Branch_p(ind_bran_cover,4);
                ReconsData(ReconsData_p(:,6)==SampPoint(1,6) & ReconsData_p(:,5)==0, 1:4) = ...
                    [(ReconsData_p(ReconsData_p(:,6)==SampPoint(1,6) & ReconsData_p(:,5)==0, 1:3).*2+SampPoint(1,1:3))./3, ...
                    ReconsData_p(ReconsData_p(:,6)==SampPoint(1,6) & ReconsData_p(:,5)==0, 4).*0.1+SampPoint(1,4)*0.9]; 
                ind_recdata_onbran = find(ReconsData(:,6)>(Branch_p(ind_bran_cover,4)*10) & ReconsData(:,6)<((Branch_p(ind_bran_cover,4)+1)*10) ...
                    & ReconsData(:,5)==1);
                str_recdata_onbran = unique(ReconsData(ind_recdata_onbran,6));
                for ii_samp = 2:size(SampPoint,1)
                    ang_rec = 9999; ind_ang = 1;
                    for ii_tmp = 1:length(str_recdata_onbran)
                        ind_recdata_onthisbran = find(ReconsData(:,6)==str_recdata_onbran(ii_tmp));
                        [~,ind_merge] = min(vecnorm(SampPoint(ii_samp,1:3)'-ReconsData(ind_recdata_onthisbran,1:3)'));
                        ang_tmp = acos( dot(SampPoint(ii_samp,1:3), ReconsData(ind_recdata_onthisbran(ind_merge),1:3))/...
                            (norm(SampPoint(ii_samp,1:3))*norm(ReconsData(ind_recdata_onthisbran(ind_merge),1:3))) );
                        if ang_tmp<ang_rec
                            ang_rec = ang_tmp; ind_ang = ii_tmp;
                        end
                    end
                    if ang_rec<=(30/180*pi)
                        SampPoint(ii_samp,6) = ReconsData(ind_recdata_onbran(ind_ang),6);
                        if Branch(Branch(:,4)==SampPoint(ii_samp,6),7)~=2
                            ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                        end
                    else
                        ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                    end
                end
            elseif ~isempty(Branch_p) && ~flag_cover 
                ind_currdir_rec = find(Branch(:,7)==1 & Branch(:,3)==1);
                [~,ind] = max(Branch(ind_currdir_rec,4));
                ind_currdir_rec = ind_currdir_rec(ind);
                if isempty(ind_currdir_rec)
                    ReconsData(ReconsData(:,6)~=0) = 0;
                    ReconsData(ReconsData(:,5)==0) = 1;
                    ReconsData = [ReconsData; SampPoint(1,:)];
                    num_bran = 1;
                    for ii_samp = 2:size(SampPoint,1)
                        SampPoint(ii_samp,6) = num_bran;
                        ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                        num_bran = num_bran+1;
                    end
                else
                    SampPoint(1,6) = Branch_p(ind_currdir_rec,4);
                    if ~any(ReconsData_p(:,6)==SampPoint(1,6) & ReconsData_p(:,5)==0)
                        ReconsData = [ReconsData; SampPoint(1,:)];
                        num_bran = 1;
                        for ii_samp = 2:size(SampPoint,1)
                            SampPoint(ii_samp,6) = num_bran+Branch_p(ind_currdir_rec,4)*10;
                            ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                            num_bran = num_bran+1;
                        end
                    end
                end
            else
                ReconsData = [ReconsData; SampPoint(1,:)];
                num_bran = 1;
                for ii_samp = 2:size(SampPoint,1)
                    SampPoint(ii_samp,6) = num_bran;
                    ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                    num_bran = num_bran+1;
                end
            end
        else
            ReconsData = SampPoint(1,:);
            num_bran = 1;
            for ii_samp = 2:size(SampPoint,1)
                SampPoint(ii_samp,6) = num_bran;
                ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                num_bran = num_bran+1;
            end
        end
    elseif size(EndPoint,1)>=1 && all_end~=1
        [SkelPoint_side_tmp,~] = find(img_area_side(:,max(1,round(EndPoint(1,1)/2)-5):min(round(EndPoint(1,1)/2)+5,128))); 
        if ~isempty(SkelPoint_side_tmp)
            SkelPoint_side_tmp = mean(SkelPoint_side_tmp);
            SampPoint_pix = [EndPoint(1,1:2), SkelPoint_side_tmp, diamImg_top(EndPoint(1,2),EndPoint(1,1)),1, 0]; 
        else
            SampPoint_pix = [EndPoint(1,1:2), 0, diamImg_top(EndPoint(1,2),EndPoint(1,1)),1, 0];
        end
        num_samp = max(floor(perimeter_skel_top./30),2);
        dist_samp = perimeter_skel_top/(num_samp-1); num_samped = 1;
        while dist_samp < (perimeter_skel_top-perimeter_skel_top/(num_samp-1)) && num_samped<num_samp
            SkelPoint_top_tmp = SkelPoint_top(vecnorm(EndPoint(1,:)'-SkelPoint_top')>dist_samp,:);
            [~,ind_skel_tmp] = min(vecnorm(EndPoint(1,:)'-SkelPoint_top_tmp'));
            SkelPoint_top_tmp = SkelPoint_top_tmp(ind_skel_tmp,:);
            SampPoint_tmp = zeros(size(SkelPoint_top_tmp,1),6);
            for ii_tmp = 1:size(SkelPoint_top_tmp,1)
                [SkelPoint_side_tmp,~] = find(img_area_side(:,max(1,round(SkelPoint_top_tmp(ii_tmp,1)/2)-5):min(round(SkelPoint_top_tmp(ii_tmp,1)/2)+5,128))); 
                if ~isempty(SkelPoint_side_tmp)
                    SkelPoint_side_tmp = mean(SkelPoint_side_tmp);
                    SampPoint_tmp(ii_tmp,:) = [SkelPoint_top_tmp(ii_tmp,:), SkelPoint_side_tmp, ...
                        max(max(diamImg_top(max(1,EndPoint(ii_tmp,2)-10):min(EndPoint(ii_tmp,2)+10,256),max(1,EndPoint(ii_tmp,1)-10):min(EndPoint(ii_tmp,1)+10,256)))),1, 0]; 
                else
                    SampPoint_tmp(ii_tmp,:) = [SkelPoint_top_tmp(ii_tmp,:), 0, ...
                        max(max(diamImg_top(max(1,EndPoint(ii_tmp,2)-10):min(EndPoint(ii_tmp,2)+10,256),max(1,EndPoint(ii_tmp,1)-10):min(EndPoint(ii_tmp,1)+10,256)))),1, 0];
                end
            end
            SampPoint_pix = [SampPoint_pix;SampPoint_tmp];
            num_samped = num_samped+size(SkelPoint_top_tmp,1);
            dist_samp = dist_samp+perimeter_skel_top/(num_samp-1);
        end
        % CALCULATE THE COORDINATES
        SampPoint(:,1) = ( (xyROI_p(1)+xyROI_p(3))/2 + (SampPoint_pix(:,1)./256.*ROI_len-ROI_len/2) - 1400/2 ).*p_x;
    	SampPoint(:,2) = -( (xyROI_p(2)+xyROI_p(4))/2 + (SampPoint_pix(:,2)./256.*ROI_len-ROI_len/2) - 900/2 ).*p_y;
        SampPoint(:,3) = 124 - ( (16/26325) * ((hvROI_p(2)+hvROI_p(4))/2 + (SampPoint_pix(:,3)./400.*ROI_hei-ROI_hei/2) - 450) .* (376+SampPoint(:,2)) );
        SampPoint(:,4) = SampPoint_pix(:,4)./256.*ROI_len.*p_x;
        SampPoint(:,5:6) = SampPoint_pix(:,5:6);
        % DETERMINE WHETHER TO MERGE OR NOT
        for ii_samp = 1:size(SampPoint,1)
            if SampPoint(ii_samp,3)==0 || isnan(SampPoint(ii_samp,3)) || abs(SampPoint(ii_samp,3)-z)>5
                SampPoint(ii_samp,3) = z;
            end
            if ~isempty(ReconsData)
                if PointType==2
                    ind_currdir = find(Branch(:,7)==1 & Branch(:,3)==1);
                    [~,ind_dir] = max(Branch(ind_currdir,4));
                    ind_currdir = ind_currdir(ind_dir);
                    ind_currdir_rec = ind_currdir;
                    SampPoint(ii_samp,6) = Branch(ind_currdir_rec,4);
                    ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                elseif ~isempty(Branch) && any(Branch(:,7)==1 & Branch(:,3)==1) 
                    ind_currdir_rec = find(Branch(:,7)==1 & Branch(:,3)==1);
                    [~,ind] = max(Branch(ind_currdir_rec,4));
                    ind_currdir_rec = ind_currdir_rec(ind);
                    SampPoint(ii_samp,6) = Branch(ind_currdir_rec,4);
                    ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                else
                    ReconsData = [ReconsData; SampPoint(ii_samp,:)];
                end
            else
                ReconsData = SampPoint(ii_samp,:);
            end
        end
    else
        % DO NOTHING
    end
    % DATA RECORDING END
    
else
    PosSwm = PosSwm_p;
    xyROI = xyROI_p;
    xy_co_mat = xy_cpre_mat;
    hvROI = hvROI_p;
    hv_co_mat = hv_cpre_mat;
    Branch = Branch_p;
end

dist_hvroi = hvROI(1)-hvROI_p(1);
frame_ct_n = max(frame_ct-1,0);
xyhv_co_mat = [xy_co_mat,hv_co_mat];

end

