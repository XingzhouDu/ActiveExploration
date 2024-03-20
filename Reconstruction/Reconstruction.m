% Copyright © 2023 DU Xingzhou.
% Licence: Apache license 2.0.

clear
clc
close all

CurrPath = pwd;
cd('..');
Branch_Exp = readmatrix('Demo for Active Exploration\sample result\Branch.xlsx');
Recons_Exp = readmatrix('Demo for Active Exploration\sample result\ReconstructionData.xlsx');
cd(CurrPath);

ReconsData = Recons_Exp;

f9 = figure('Position',[800 200 400 240]);
figure(f9)
plot3(ReconsData(:,1),ReconsData(:,2),ReconsData(:,3),'*');
hold on;
str_bran_o = unique(ReconsData(:,6));
str_bran_o = sort(str_bran_o);
xlabel('x (mm)','interpreter','latex');
ylabel('y (mm)','interpreter','latex');
zlabel('z (mm)','interpreter','latex');
grid on;
axis equal;
set(gca,'ZTick',70:20:99)
view([20,50])

NumFrame = max(Branch_Exp(Branch_Exp(:,2)==999,1));
ind_bra = find(Branch_Exp(:,1)==(NumFrame-2) & Branch_Exp(:,2)==999);
ind_bra_next = find(Branch_Exp(:,1)==(NumFrame-1) & Branch_Exp(:,2)==999);
if ind_bra ~= ind_bra_next-1
    Branch_fin = Branch_Exp((ind_bra+1):(ind_bra_next-1),:);
else
    Branch_fin = [];
end

str_bran = unique(Branch_fin(:,4));
str_bran = sort(str_bran);
flag_remain = ones(size(ReconsData,1),1);
num_branplot = 1; flag_haveend = 0;
ReconsData_rev = [];
length_uni = [];
for ii_plot = 1:length(str_bran) 
    if length(find(ReconsData(:,5)==1 & ReconsData(:,6)==str_bran(ii_plot)))<=1
        continue;
    end
    ReconsData_braend = ReconsData(ReconsData(:,5)==0 & ReconsData(:,6)==str_bran(ii_plot),:);
    ind_endbran = find(ReconsData(:,5)==1 & ReconsData(:,6)>str_bran(ii_plot)*10 & ReconsData(:,6)<=(str_bran(ii_plot)+1)*10);
    if length(ind_endbran) <= length(unique(ReconsData(ind_endbran,6)))
        ReconsData_braend = [];
    end
    if str_bran(ii_plot)==0
        ind_onbran = find(ReconsData(:,6)==str_bran(ii_plot));
        ReconsData_rev_tmp = ReconsData(ind_onbran(1:2),:);
        dir_pre = ReconsData_rev_tmp(2,1:2)-ReconsData_rev_tmp(1,1:2);
        ind_remain = ind_onbran(3:end);
        flag_remain(ind_onbran) = 0;
    else
        ReconsData_rev_tmp = ReconsData(ReconsData(:,5)==0 & ReconsData(:,6)==floor(str_bran(ii_plot)/10),:);
        ReconsData_rev_tmp = ReconsData_rev_tmp(end,:);
        if ~isempty(ReconsData_braend)
            dir_pre = ReconsData_braend(end,1:2)-ReconsData_rev_tmp(1:2);
        else
            dir_pre = [Branch_fin(Branch_fin(:,3)==1 & Branch_fin(:,4)==str_bran(ii_plot), 1),...
                -Branch_fin(Branch_fin(:,3)==1 & Branch_fin(:,4)==str_bran(ii_plot), 2)];
        end
        ind_remain = find(ReconsData(:,5)==1 & ReconsData(:,6)==str_bran(ii_plot));
    end
    
    ReconsData_remain = ReconsData(ind_remain,:);
    stop_branch = 0;
    while stop_branch~=1
        dst_min = 999; connscore_min = 999999; ang_min = 180;
        for ii_rev = 1:size(ReconsData_remain,1)
            dst_tmp = vecnorm(ReconsData_remain(ii_rev,1:2)'-ReconsData_rev_tmp(end,1:2)');
            if ~isempty(ReconsData_braend)
                dir_tmp = ReconsData_remain(ii_rev,1:2)-ReconsData_rev_tmp(1,1:2);
                ang_tmp = real(acos( dot(dir_tmp,dir_pre)/(norm(dir_tmp)*norm(dir_pre)) ))/pi*180 /90*35;
            else
                dir_tmp = ReconsData_remain(ii_rev,1:2)-ReconsData_rev_tmp(end,1:2);
                ang_tmp = real(acos( dot(dir_tmp,dir_pre)/(norm(dir_tmp)*norm(dir_pre)) ))/pi*180;
            end
            connscore = dst_tmp^3*sqrt(1-cos(ang_tmp/180*pi))+10^3*max(0,ang_tmp-35)^2+10^3*max(0,dst_tmp-15)^2;
            if connscore<connscore_min 
                connscore_min = connscore; 
                dst_min = dst_tmp;
                ang_min = ang_tmp;
                ind_conn = ii_rev;
            end
        end

        if dst_min<30 && ang_min<135
            ReconsData_rev_tmp = [ReconsData_rev_tmp;ReconsData_remain(ind_conn,:)];
            if ~isempty(ReconsData_braend)
                dir_pre = ReconsData_braend(end,1:2)-ReconsData_rev_tmp(end,1:2);
            else
                dir_pre = ReconsData_rev_tmp(end,1:2)-ReconsData_rev_tmp(end-1,1:2);
            end
            if ~isempty(ReconsData_braend) && vecnorm(ReconsData_remain(ind_conn,1:2)'-ReconsData_braend(end,1:2)')<10
                ReconsData_rev_tmp = [ReconsData_rev_tmp;ReconsData_braend(end,:)];
                stop_branch = 1;
            end
            ind_remain = [ind_remain(1:ind_conn-1);ind_remain(ind_conn+1:end)];
            ReconsData_remain = ReconsData(ind_remain,:);
        elseif ~isempty(ReconsData_braend)
            ReconsData_rev_tmp = [ReconsData_rev_tmp;ReconsData_braend(end,:)];
            stop_branch = 1;
        else
            stop_branch = 1;
        end
    end

    for ii_zcorr = 2:size(ReconsData_rev_tmp,1)-1
        if abs(ReconsData_rev_tmp(ii_zcorr,3)-ReconsData_rev_tmp(ii_zcorr-1,3)) > ...
                tan(50/180*pi)*norm(ReconsData_rev_tmp(ii_zcorr,1:2)-ReconsData_rev_tmp(ii_zcorr-1,1:2))
            ReconsData_rev_tmp(ii_zcorr,3) = ReconsData_rev_tmp(ii_zcorr-1,3);
        end
        if ReconsData_rev_tmp(ii_zcorr,4) ==0 && ReconsData_rev_tmp(ii_zcorr-1,4) ~=0
            ReconsData_rev_tmp(ii_zcorr,4) = ReconsData_rev_tmp(ii_zcorr-1,4);
        elseif ReconsData_rev_tmp(ii_zcorr,4) ==0
            ReconsData_rev_tmp(ii_zcorr,4) = mean(ReconsData_rev_tmp(ReconsData_rev_tmp(:,4)~=0,4));
        else
        end
    end
       
    ReconsData_rev_tmp = [ReconsData_rev_tmp,str_bran(ii_plot)*ones(size(ReconsData_rev_tmp,1),1)];
    ReconsData_rev = [ReconsData_rev; ReconsData_rev_tmp];
end

flag_uni = zeros(length(str_bran),1);
ReconsData_uni = [];
for ii_uni = 1:length(str_bran)
    if flag_uni(ii_uni)==0
        ReconsData_uni_tmp = ReconsData_rev(ReconsData_rev(:,7)==str_bran(ii_uni),:);
        flag_merge = 0;
        for jj_uni = ii_uni+1:length(str_bran)
            ReconsData_uni_tmp_tmp = ReconsData_rev(ReconsData_rev(:,7)==str_bran(jj_uni),:);
            if flag_uni(jj_uni)==0 && ~isempty(ReconsData_uni_tmp_tmp)...
                    && ((all(ReconsData_uni_tmp_tmp(1,1:2)==ReconsData_uni_tmp(1,1:2)) && all(ReconsData_uni_tmp_tmp(end,1:2)==ReconsData_uni_tmp(end,1:2))) ...
                    || (all(ReconsData_uni_tmp_tmp(1,1:2)==ReconsData_uni_tmp(end,1:2)) && all(ReconsData_uni_tmp_tmp(end,1:2)==ReconsData_uni_tmp(1,1:2))))
                if ReconsData_uni_tmp(1,5)==0
                    ReconsData_uni_mer = ReconsData_uni_tmp;
                else
                    ReconsData_uni_mer = ReconsData_uni_tmp_tmp;
                    ReconsData_uni_tmp_tmp = ReconsData_uni_tmp;
                end
                for jjj_mer = 2:size(ReconsData_uni_tmp_tmp,1)-1
                    for iii_mer = 2:size(ReconsData_uni_mer,1)
                        if norm(ReconsData_uni_mer(end,1:2)-ReconsData_uni_mer(iii_mer-1,1:2)) ...
                                > norm(ReconsData_uni_mer(end,1:2)-ReconsData_uni_tmp_tmp(jjj_mer,1:2)) ...
                                && norm(ReconsData_uni_mer(end,1:2)-ReconsData_uni_mer(iii_mer,1:2)) ...
                                < norm(ReconsData_uni_mer(end,1:2)-ReconsData_uni_tmp_tmp(jjj_mer,1:2))...
                            ReconsData_uni_tmp_tmp(jjj_mer,end) = ReconsData_uni_mer(iii_mer-1,end);
                            ReconsData_uni_mer = [ReconsData_uni_mer(1:iii_mer-1,:);ReconsData_uni_tmp_tmp(jjj_mer,:);...
                                ReconsData_uni_mer(iii_mer:end,:)];
                        end
                    end
                end
                ReconsData_uni_mer(:,1) = smooth(ReconsData_uni_mer(:,1)); ReconsData_uni_mer(:,2) = smooth(ReconsData_uni_mer(:,2));
                ReconsData_uni_mer(:,3) = smooth(ReconsData_uni_mer(:,3)); ReconsData_uni_mer(:,4) = smooth(ReconsData_uni_mer(:,4));
                ReconsData_uni = [ReconsData_uni;ReconsData_uni_mer];
                flag_uni(jj_uni) = 1;
                flag_merge = 1;
                length_uni_tmp = 0;
                for iii_length = 2:size(ReconsData_uni_mer,1)
                    length_uni_tmp = length_uni_tmp+norm(ReconsData_uni_mer(iii_length,1:3)-ReconsData_uni_mer(iii_length-1,1:3));
                end
                length_uni = [length_uni;length_uni_tmp];
                disp(str_bran(ii_uni))
            elseif isempty(ReconsData_uni_tmp_tmp)
                flag_uni(jj_uni) = 1;
            end
        end
        if flag_merge==0
            ReconsData_uni_tmp(:,1) = smooth(ReconsData_uni_tmp(:,1)); ReconsData_uni_tmp(:,2) = smooth(ReconsData_uni_tmp(:,2));
            ReconsData_uni_tmp(:,3) = smooth(ReconsData_uni_tmp(:,3)); ReconsData_uni_tmp(:,4) = smooth(ReconsData_uni_tmp(:,4));
            ReconsData_uni = [ReconsData_uni;ReconsData_uni_tmp];
            length_uni_tmp = 0;
            for iii_length = 2:size(ReconsData_uni_tmp,1)
                length_uni_tmp = length_uni_tmp+norm(ReconsData_uni_tmp(iii_length,1:3)-ReconsData_uni_tmp(iii_length-1,1:3));
            end
            length_uni = [length_uni;length_uni_tmp];
        end
        flag_uni(ii_uni) = 1;
    end
end

ReconsData_plot = ReconsData_uni;
str_bran = unique(ReconsData_plot(:,7));
str_bran = sort(str_bran);
f0 = figure('Position',[500 200 400 200]);
figure(f0);
hold on;
num_branplot = 1;
length_uni_act_lym = 0;
for iii_plot = 1:length(str_bran)
    length_uni_act_lym = length_uni_act_lym+length_uni(iii_plot); num_color = num_branplot;
    switch mod(num_color,7)
        case 1
            color_tube = [128 158 194]/255;
        case 2
            color_tube = [243 164 71]/255;
        case 3
            color_tube = [165 181 146]/255;
        case 4
            color_tube = [208 146 167]/255;
        case 5
            color_tube = [128 158 194]/255;
        case 6
            color_tube = [231 188 41]/255;
        otherwise
            color_tube = [156 133 192]/255;
    end
    plot3(ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),1),ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),2),...
        ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),3),'DisplayName',string(str_bran(iii_plot)),'linewidth',1,'color',color_tube);
%     [x_tb, y_tb, z_tb] = tubeplot([ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),1)';ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),2)';...
%         ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),3)'],ReconsData_plot(ReconsData_plot(:,7)==str_bran(iii_plot),4)'./2); % demands tubeplot function: https://www.mathworks.com/matlabcentral/fileexchange/5562-tubeplot
%     surf(x_tb, y_tb, z_tb,'EdgeColor','None','FaceColor',color_tube);
    num_branplot = num_branplot+1;
end
xlabel('x (mm)','interpreter','latex');
ylabel('y (mm)','interpreter','latex');
zlabel('z (mm)','interpreter','latex');
grid on;
axis equal;
axis([-100,50,-70,70,60,90])
set(gca,'ZTick',60:15:90)
view([30,30])

