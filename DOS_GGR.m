function DOS_GGR
%% DOS calculation using generalized GR method
% The program is for DOS calculation using generalized GR method, referring
% to article "Generalized Gilat-Raubenheimer Method for Density-of-States 
% Calculation in Photonic Crystals".

% This is edition 1, using for loop to calculate DOS just like C program.

%% Important notice for initial parameters!!!
% Necessary parameters: 
% 0. three files include band frequencies on high symmetry points,
% band frequencies in the whole Brillouin zone, group velocities in the
% whole Brillouin zone; (Notice: The band and velocity data could also be the half of 
% the Brillouin zone in latter two files if there is time-reversal symmetry, which will save half of
% the computing time for bands. The values of DOS only need to be multiplied 
% by 2 if you choose to use half of the Brillouin zone.)
% 1. the reciprocal vectors; 2. the number of k points; 3. number of bands.
% Optional parameters: 4. maximum and minimum of band frequency (w_max, w_min); 
% 5. resolution about the frequency  (N_w); 6. parameters about plot like color, fontsize, etc.

draw_band=1;    % draw the band structure with DOS
% draw_band=0;    % only calculate DOS and output the DOS data

file_bandline='band.txt';
file_bandmap='frequency_GGR.txt';
file_velocity='velocity_GGR.txt';

file_DOSdata='output.txt';    % save file for Density of states data

% notice the velocity should be compatible with bi basis length
reciprocalvector1=[0 1 1];
reciprocalvector2=[1 0 1];
reciprocalvector3=[1 1 0];

% num_kpoints(i) is the number of k points along the bi axis, 
num_kpoints=[12,12,12];

N_band=10;      % the total number of frequency bands

w_max_custom=-1;   % the range of frequency, '-1' denotes default settings
w_min_custom=-1;

N_w_custom=20000;       % denotes the resolution of frequency : dw = (w_max - w_min) / N_w

kinter = 30;       % the inter quantity of k points between two high symmetry points
maxDOS_custom=-1;        % the parameters about plot, '-1' denotes default settings
fs_custom=10;
bandcolor_custom='b';
bottomcolor_custom='k';
thelinewidth_custom=1;

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
hside=1/num_kpoints(1)/2;    % half of side length of one transfromed cube
DOS=zeros(N_w+1,1);       % initialze the density of states array

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
                    DOS(n_tmpt+1)=DOS(n_tmpt+1)+DOScontribution;
                    flag_delta_n_w=1;
                end
                continue;
            else
                if (n_tmpt>=0)&&(n_tmpt<=N_w)
                    DOS(n_tmpt+1)=DOS(n_tmpt+1)+DOScontribution;
                end
            end
        end
    end
end
% output DOS data into output.txt
file_output=fopen(file_DOSdata,'wt');
for nprint_w=1:N_w+1
    fprintf(file_output,'%.10f %.10f\n',w_min+step_w*(nprint_w-1),DOS(nprint_w));
end

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
b1=vectorb1;
b2=vectorb2;
b3=vectorb3;
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
        rr=sqrt(r1^2+r2^2-2*r1*r2*cost);
    end
    
    A1=Ks(end) + (A-k1)*rr/(k2-k1);
    Ks(k1:k2) = A1;
    kidx = [kidx,k1];
end
kidx=[kidx,k2];

fs=fs_custom;
bandcolor=bandcolor_custom;
bottomcolor=bottomcolor_custom;
thelinewidth=thelinewidth_custom;
figure
for i = 1:nbands 
    plot(Ks,data_band(:,3+i),'-','color',bandcolor,'LineWidth',thelinewidth);
    hold on;
end

if maxDOS_custom==-1
    maxDOS=ceil(max(DOS));
else
    maxDOS=maxDOS_custom;
end

DOS(DOS>maxDOS)=maxDOS;
w_var=w_min+step_w*((1:(N_w+1))-1);   % frequency -- the variable of DOS
DOS_nrm=(Ks(end)-Ks(1))*DOS/maxDOS+Ks(end);
plot(DOS_nrm,w_var,'Color',bandcolor);
fill(DOS_nrm,w_var,bandcolor);
plot(DOS_nrm(1)*ones(size(w_var,2),1),w_var,'color',bottomcolor);
set(gca,'FontSize',fs,'FontName','Helvetica','Layer','top');
set(gca,'xTick', [Ks(kidx),Ks(end)*2],'XTickLabel',{'H','\Gamma','N','P','\Gamma (0)',...
    num2str(maxDOS)},'XGrid','on','GridLineStyle','-','layer','bottom');
xlim([Ks(1),2*Ks(end)]);
ylim([w_min,w_max]);
ylabel('Normalized frequency \omega (a/\lambda_0)');
title('Band structure and its corresponding DOS');
hold off

saveas(gcf,'BandFigure.fig');
print('-depsc','-painters','BandFigure');
