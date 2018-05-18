% function DOS_GGR
%% DOS calculation using generalized GR (GGR) method
% The program is for DOS calculation using GGR method, referring
% to article "Generalized Gilat-Raubenheimer Method for Density-of-States 
% Calculation in Photonic Crystals". For more information, please refer to our website:
% https://github.com/boyuanliuoptics/DOS-calculation/edit/master/DOS_GGR.m
% The first edition is finished in Nov. 20th, 2017.
%% Important notice for initial parameters!!!
% Necessary parameters: 
% 0. three files include band frequencies on high symmetry points,
% band frequencies in the whole Brillouin zone, group velocities in the
% whole Brillouin zone; 
% 1. the reciprocal vectors; 2. the number of k points; 3. number of bands.
% Optional parameters: 4. maximum and minimum of band frequency (w_max, w_min); 
% 5. resolution about the frequency  (N_w); 6. inter number of points between 
% two high symmetry k points; 7. parameters about plot like color, fontsize, etc;

draw_band=1;    % draw the band structure with DOS
% draw_band=0;    % only calculate DOS and output the DOS data

file_bandline='band.txt';
file_bandmap='frequency_GGR.txt';
file_velocity='velocity_GGR.txt';

file_DOSdata='output.txt';    % save file for Density of states data

% 2D examle: trirods.ctl
% reciprocalvector1=[sqrt(3)/2 1/2 0];
% reciprocalvector2=[sqrt(3)/2 -1/2 0];
% reciprocalvector3=[0 0 0];
% 3D examle: DG.ctl
reciprocalvector1=[0 1 1];
reciprocalvector2=[1 0 1];
reciprocalvector3=[1 1 0];

% num_kpoints(i) is the number of k points along the bi axis. 
% The calculation at present requires the sampling to be uniform in BZ,
% which is num_kpoints(1)=num_kpoints(2)=num_kpoints(3) in 3D structures.
% If you use half of the BZ (T symmetry), num_kpoints(1) should be half of
% num_kpoints(2).
% if it is 2D structure, num_kpoints(3) should be 1.
num_kpoints=[6,12,12];

N_band=20;      % the total number of frequency bands

w_max_custom=-1;   % the range of frequency, '-1' denotes default settings
w_min_custom=-1;
w_max_dsp_coefficient=0.9;  % the displaying maximum frequency is the product of w_max_dsp_coefficient and w_max
N_w_custom=20000;       % denotes the resolution of frequency : dw = (w_max - w_min) / N_w

kinter = 30;       % the inter quantity of k points between two high symmetry points
maxDOS_custom=100;        % the parameters about plot, '-1' denotes default settings
sequence_points={'H','\Gamma','N','P','\Gamma'};    % sequence of high symmetry points in 3D example
DOSratio=1;        % the plot ratio of DOS in x axis 
%% Initialization and import data
% the reciprocal vectors initialization
vectorb1=reciprocalvector1;
vectorb2=reciprocalvector2;
vectorb3=reciprocalvector3;
vectorsb=[vectorb1;vectorb2;vectorb3];

% Nx,Ny,Nz is the number of k points along the x,y,z axis. N_kpoints is the
% total number of k points in Brillouin zone
n_kpoints=prod(num_kpoints);

N_kpoints=N_band*n_kpoints;

% import data
% the two importing txt files are arranged as matrix of N*1 and N*3
dataw=importdata(file_bandmap);
datav_original=importdata(file_velocity);   % the real group velocities
datav=(vectorsb*datav_original')';     % the transformed group velocities

if w_max_custom==-1
    w_max=1.05*max(dataw); % the maximum of frequency should be larger than max(dataw) a little
else
    w_max=w_max_custom;
end

if w_min_custom==-1
    w_min=0;
else
    w_min=w_min_custom;
end

itmd_v=sort(abs(datav),2,'descend');   % intermediate velocity: v1 >= v2 >= v3

% other parmeters

% N_w=20*N_kpoints;       % divide the frequency region into N_w part
N_w=N_w_custom;
step_w=(w_max-w_min)/N_w;      % the resolution of frequency
hside=1/num_kpoints(2)/2;    % half of side length of one transfromed cube
DOSarray=zeros(N_w+1,1);       % initialze the density of states array

w1=hside*abs(itmd_v(:,1)-itmd_v(:,2)-itmd_v(:,3));
w2=hside*(itmd_v(:,1)-itmd_v(:,2)+itmd_v(:,3));
w3=hside*(itmd_v(:,1)+itmd_v(:,2)-itmd_v(:,3));
w4=hside*(itmd_v(:,1)+itmd_v(:,2)+itmd_v(:,3));

%% DOS calculation
% principle of calculation process can be found in our article 
% "Generalized Gilat-Raubenheimer Method for Density-of-States Calculation 
% in Photonic Crystals"
for num_k=1:N_kpoints
    n_w_kcenter=round((dataw(num_k)-w_min)/step_w);
    v=norm(datav(num_k,:));
    v1=itmd_v(num_k,1);
    v2=itmd_v(num_k,2);
    v3=itmd_v(num_k,3);
    
    flag_delta_n_w=0;       % first time compute delta_n_w = 1
    for vdirection=0:1      % two velocity directions denote w-w_k0 > 0 and <0
        for delta_n_w=1:N_w
            n_tmpt=n_w_kcenter+(-1)^vdirection*(delta_n_w-1);
            delta_w=abs(dataw(num_k)-(n_tmpt*step_w+w_min));
            if delta_w<=w1(num_k)
                if v1>=v2+v3
                    DOScontribution=4*hside^2/v1;
                else
                    DOScontribution=(2*hside^2*(v1*v2+v2*v3+v3*v1)-...
                        (delta_w^2+(hside*v)^2))/v1/v2/v3;
                end
            elseif delta_w<w2(num_k)
                DOScontribution=(hside^2*(v1*v2+3*v2*v3+v3*v1)-...
                    hside*delta_w*(-v1+v2+v3)-(delta_w^2+hside^2*v^2)/2)/v1/v2/v3;
            elseif delta_w<w3(num_k)
                DOScontribution=2*(hside^2*(v1+v2)-hside*delta_w)/v1/v2;
            elseif delta_w<w4(num_k)
                DOScontribution=(hside*(v1+v2+v3)-delta_w)^2/v1/v2/v3/2;
            else
                break;
            end
            if DOScontribution>8*hside^3/step_w
                DOScontribution=8*hside^3/step_w;
            end

            if delta_n_w==1       % when delta_n_w == 1, we only compute it once
                if flag_delta_n_w==0
                    DOSarray(n_tmpt+1)=DOSarray(n_tmpt+1)+DOScontribution;
                    flag_delta_n_w=1;
                end
                continue;
            else
                if (n_tmpt>=0)&&(n_tmpt<=N_w)
                    DOSarray(n_tmpt+1)=DOSarray(n_tmpt+1)+DOScontribution;
                end
            end
        end
    end
end

% output DOS data into output.txt
if num_kpoints(3)==1    % the structure is 2 dimension
    DOSarray=DOSarray*num_kpoints(2);
end
if num_kpoints(1)*2==num_kpoints(2)     % the structure has time-reversal symmetry
    DOSarray=DOSarray*2;
end
    
file_output=fopen(file_DOSdata,'wt');
for nprint_w=1:(N_w+1)
    fprintf(file_output,'%.10f %.10f\n',w_min+step_w*(nprint_w-1),DOSarray(nprint_w));
end
fclose(file_output);
%% Band plot and DOS
% import band data
if draw_band==0
    return;
end

data_band = dlmread(file_bandline,' ',0,0); % the format corresponds to bash file

nbands = size(data_band,2)-3;

if nbands~=N_band       % numbers of band in two file are unequal!
    exit('error:numbers of band in two file are unequal!');
end

kindex = 1:size(data_band(:,1),1);
Ks = 0;   % record all
kidx=[];    % record the nodes of band plot

%scale each BZ section with the tight proportion
b1=reciprocalvector1;
b2=reciprocalvector2;
b3=reciprocalvector3;
bs = [b1;b2;b3];
imax = (length(kindex)-1)/(kinter+1); % how many sections
for i=1:imax
    k1=(i-1)*(kinter+1)+1; %starting section k
    k2=k1+(kinter+1); %ending section k
    A=kindex(k1:k2);
    
    % compute the length of a section (rr)
    coor_k1=data_band(k1,1:3)*bs;
    coor_k2=data_band(k2,1:3)*bs;
    r1=norm(coor_k1);
    r2=norm(coor_k2);
    if r1*r2==0
        rr=abs(r1-r2);
    else
        cost=sum(coor_k1.*coor_k2)/r1/r2;
        if r1^2+r2^2-2*r1*r2*cost<0
            error('Error:|r|<0\n');
        end
        rr=sqrt(r1^2+r2^2-2*r1*r2*cost);
    end
    
    A1=Ks(end) + (A-k1)*rr/(k2-k1);
    Ks(k1:k2) = A1;
    kidx = [kidx,k1];
end
kidx=[kidx,k2];

fs=10;
bandcolor=[34 34 120]/255;
doscolor=[122 122 174]/255;
bottomcolor='k';
bandlinewidth=1;
doslinewidth=1;
w_max_lim=w_max*w_max_dsp_coefficient;
figure
for i = 1:nbands 
    plot(Ks,data_band(:,3+i),'-','color',bandcolor,'LineWidth',bandlinewidth);
    hold on;
end

if maxDOS_custom==-1
    maxDOS=ceil(max(DOSarray));
else
    maxDOS=maxDOS_custom;
end

DOSarray(DOSarray>maxDOS)=maxDOS;
w_var=w_min+step_w*((1:(N_w+1))-1);   % frequency -- the variable of DOS
DOS_nrm=DOSratio*(Ks(end)-Ks(1))*DOSarray/maxDOS+Ks(end);
plot(DOS_nrm,w_var,'Color',bandcolor,'LineWidth',doslinewidth);
fill(DOS_nrm,w_var,doscolor);
plot(DOS_nrm(1)*ones(size(w_var,2),1),w_var,'color',bottomcolor);
sequence_points{end}=strcat(sequence_points{end},' | 0');
set(gca,'xTick', [Ks(kidx),Ks(end)*(1+DOSratio)],'XTickLabel',[sequence_points, num2str(maxDOS)],...
    'XGrid','on','GridLineStyle','-','layer','bottom');
xlim([Ks(1),(1+DOSratio)*Ks(end)]);
ylim([w_min,w_max_lim]);
ylabel('Frequency \omegaa/(2\pic)');
title('Band structure and Density of states D2\pic/a');
set(gca,'FontSize',fs,'FontName','Helvetica','Layer','top');
annotation(gcf,'textbox',...
    [0.71482722007722 0.0602135060573348 0.120428571428572 0.0476190476190477],...
    'String',{'per cell'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off');
hold off

saveas(gcf,'BandFigure.fig');
print('-depsc','-painters','BandFigure');
