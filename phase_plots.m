%% NEEDS fullpeaks.m, peaks.m, ThresholdingAlgo.m, ecatdir (1_Ecatalogue, files should have extension .e), momdir (Moments, files
%% should have extension .m), vmaxdir (1_vmax, files should have extension .v), popdir (Popr, files should have extension .p)
%% All files should follow the naming convention - Pmax_Pmaxdecimal____rP_rPdecimal.extension
%% Rename all the files (1_Ecatalogue, 1_vmax, Popr, Moments) into the naming convention above
%%

%% Variable listing
% a           |    split filename between pore pressure and rate
% adv         |    matrix containing Pmax,rP,advancement
% advancement |    advancement of first perturbed earthquake
% amom2       |    cumulative aseismic moment
% arec        |    average aseismic+seismic recurrence interval (preturbed)
% arecavg     |    matrix containing Pmax,rP,arec
% as_ratio    |    ratio of seismic to aseismic events (perturbed)
% ase         |    list of aseismic events (perturbed)
% asr         |    matrix containing Pmax,rP,as_ratio
% avg         |    average peak slip rate (perturbed)
% b           |    save filename of pore pressure file
% c           |    save filename of vmax file
% c1moments   |    only category 1 aseismic moment matrix
% cb          |    colorbar variable for labeling
% d           |    list of all files in ecatdir
% def_rec     |    recurrence interval without perturbaation
% e           |    1_Ecatalogue file
% em          |    starting times of all perturbed earthquakes
% et          |    start time of first perturbed earthquake
% et2         |    start time of last unperturbed earthquake
% f           |    save filename of moments file
% f1          |    waitbar handle
% f2          |    general figure handle
% g1          |    grid of thresholds
% g2          |    grid of smoms
% g3          |    grid of recavg
% g4          |    grid of asr
% g5          |    grid of naes
% g6          |    grid of nses
% g7          |    grid of moments
% g8          |    grid of arecavg
% g9          |    grid of adv
% g10         |    grid of c1moments
% h1          |    plot handle of threshold (p.p. at x,t of first eq)
% h2          |    plot handle of smoms (average peak slip rate)
% h3          |    plot handle of recavg (seismic rec. int.)
% h4          |    plot handle of asr (aseismic to seismic event ratio)
% h5          |    plot handle of naes (# of aseismic events)
% h6          |    plot handle of nses (# of seismic events)
% h7          |    plot handle of moments (cumulative aseismic moment)
% h8          |    plot handle of arecavg (aseismic rec. int.)
% h9          |    plot handle of adv (advancement)
% h10         |    plot handle of c1moments (category 1 aseismic moments)
% i           |    loop variable for iterating over files in ecatdir
% l           |    number of p,r combinations
% loc         |    location of first perturbed earthquake
% m           |    Moments file
% moments     |    matrix containing Pmax,rP,amom2
% nae         |    number of aseismic events
% naes        |    matrix containing Pmax,rP,nae
% nse         |    number of seismic events
% nses        |    matrix containing Pmax,rP,nse
% p           |    Pmax
% peak        |    output of peaks function
% perturbed   |    events during perturbation taken from peaks() output
% pop         |    pore pressure at the start of the first eq (perturbed)
% ppp         |    Popr file
% q           |    split filename at . to change extensions
% r           |    rP
% rec         |    seismic recurrence interval (perturbed)
% recavg      |    matrix containing Pmax,rP,rec
% slice       |    time series of moments during perturbation
% smoms       |    matrix containing Pmax,rP,avg
% sse         |    list of seismic events (perturbed)
% t1          |    start time of perturbation
% t2          |    end of rise of pore pressure
% t3          |    end of pore pressure plateau
% t4          |    end of pore pressure perturbation
% th          |    threshold pore pressure at the point of nucleation
% thresholds  |    matrix containing Pmax,rP,th
% v           |    Vmax file

warning('off','all');
%% Initialise empty arrays of size equal to the number of files, initialise array increment and defaults
[~,~]=mkdir('output'); % create output directory
[~,~]=mkdir('output/vmax_plots');
d=dir("copyecatdir/*__*.e"); % Using ecatdir as the directory to loop over because it has the smallest files
l=length(d);
% fprintf("Net size of parameter space = %d\n\n",l);
thresholds=zeros(l,3);
smoms=zeros(l,3);
recavg=zeros(l,3);
arecavg=zeros(l,3);
asr=zeros(l,3);
naes=zeros(l,3);
nses=zeros(l,3);
moments=zeros(l,3);
adv=zeros(l,3);
c1moments=[];

%f1 = waitbar(0,'0% (P_{max} = 0.000000, r_P = 0.000000)','Name','Progress');
parpool(4); % change this based on the number of logical cores that can be spared
parfor k=1:length(d) % loop over all the files in the ecatdir
    i=d(k)
    %% Get Pmax and rp from filename
    a=split(i.name,'___');
    q=split(i.name,'.');
    f=[q{1},'.m'];
    c=[q{1},'.v'];
    b=[q{1},'.p'];
    p=join(split(a{1},'_'),'.');
    p=str2double(p{1}); % Pmax
    r=join(split(a{2},'_'),'.');
    r=str2double(r{1}(1:end-2)); % rP
    % fprintf("%d)  (%d%%) Pmax = %f, rP = %f\n",k,int8(k*100/l),p,r);
    %waitbar(k/l,f1,sprintf("%d%% (P_{max} = %f, r_P = %f)",int8(k*100/l),p,r))
    if p>0 % only calculate anything if there is any perturbation
        %% Load Ecatalogue, Vmax, Moments and Popr files
        % These directories need to exist in the same folder as this program,
        % along with peaks.m
        e=load(['copyecatdir/',i.name]);
        v=load(['copyvmaxdir/',c]);
        ppp=load(['copypopdir/',b]);
        m=load(['copymomdir/',f]);

        %% Get starting and ending times of injection 
        % This following  section assumes 8th earthquake is perturbed with a 1
        % year interval and 4 year outflow period
        t1=1.325e10; % Perturbation start time
        t2=t1+(p*1e6/r);
        t3=t2+1.354e7;
        t4=t3+1.2626e8; % Perturbation end time
        et=e(8,4); % Timing of first earthquake after perturbation starts
        et2=e(7,4); % Timing of last earthquake before perturbation starts
        def_rec=56;
        advancement=def_rec*3.154e7-(et-et2); % advancement of first perturbed earthquake w.r.t an average unperturbed recurrence interval calculated from an unperturbed case (outermost else condition)
    
        %% Calculate threshold
        pop=ppp(ppp(:,1)>=et,3); % Pore pressure at the start of first earthquake (get all pore pressure values after the first earthquake,
        pop=pop(1);              % and then take the first one)
        loc=v(v(:,1)>=et,5);     % Location of the first earthquake in terms of cell number
        loc=0.05*loc(1)-50;      % Convert cell number to km (multiply by cell size and add value of starting cell)
        % Calculate pore pressure at the location of the earthquake by using the spatial function used in the simulation (appropriately converted to km)
        th=pop*(max(min((2e-2*(loc+10)+1),(2e-2*(loc+10)+1)),0)); % Pore pressure threshold at the point of nucleation
    
        %% Calculate number of events
        [perturbed,~,~]=fullpeak(q{1},v(:,4),v(:,1),t1,t4,1); % (see fullpeaks.m) Call peaks function which returns for all events between ts and te - peaks array which contains 1) peak time 2) peak slip rate 3) label (0=aseismic, 1=seismic)
        if ~isempty(perturbed) % If there are perturbed events only then calculate the rest, otherwise set some defaults
            ase=perturbed(perturbed(:,3)==0,:); % List of aseismic events during perturbation (label=0)
            sse=perturbed(perturbed(:,3)==1,:); % List of seismic events during perturbation (label=1)
            avg=mean(perturbed(:,2)); % Averave of Max V during perturbation (average peak slip rate)
            nae=numel(ase)/3; % Number of aseismic events (there are 3 columns, so number of elements/3 = number of events since each row is one event)
            nse=numel(sse)/3; % Number of seismic events
            as_ratio=nse/nae; % Ratio of seismic to aseismic events
        else
            avg=0;
            as_ratio=0;
            nae=0;
            nse=0;
        end
    
        %% Calculate average recurrence interval
        em=e(e(:,4)>=t1 & e(:,4)<=t4,4); % Starting times of all perturbed earthquakes (i.e. only seismic events)
        if isempty(em) % If there are no perturbed earthquakes, then recurrence interval is set to zero
            rec=0;
        else
            rec=(em(end)-em(1))/numel(em); % Average recurrence interval (s) = Starting time of first perturbed earthquake - starting time of last perturbed earthquake / number of perturbed earthquakes
        end
        if ~isempty(perturbed)  % If there are perturbed events only then calculate arec
            arec=(perturbed(end,1)-perturbed(1,1))/(nae+nse); % Average recurrence intervals including aseismic transients (s)
        else
            arec=0;
        end
    
        %% Calculate cumulative aseismic moment
        slice=m(m(:,1)>=t1 & m(:,1)<=t4,:); % Get slice of the recorded moments during perturbation
        amom2=log10(sum(slice(:,2))); % Cumulative aseismic moment (log10)
    
        %% Save all the variables
        adv(k,:)=[p,r,advancement/3.154e7];
        thresholds(k,:)=[p,r,th];
        smoms(k,:)=[p,r,avg];
        recavg(k,:)=[p,r,rec/3.154e7];
        arecavg(k,:)=[p,r,arec/3.154e7];
        asr(k,:)=[p,r,as_ratio];
        naes(k,:)=[p,r,nae];
        nses(k,:)=[p,r,nse];
        moments(k,:)=[p,r,amom2];


    else
        e=load(['ecatdir/',i.name]);
        def_rec=(e(end,4)-e(1,4))/(numel(e)*3.154e7/18);
    end
end

%% Save data to files
save 'output/nuc_thresh.out' thresholds -ascii;
save 'output/avg_max_v.out' smoms -ascii;
save 'output/avg_rec_int.out' recavg -ascii;
save 'output/avg_as_rec_int.out' arecavg -ascii;
save 'output/as_to_se_r' asr -ascii;
save 'output/n_as_e.out' naes -ascii;
save 'output/n_se_e.out' nses -ascii;
save 'output/log_cum_as_mom.out' moments -ascii;
save 'output/c1_log_cum_as_mom.out' c1moments -ascii;
save 'output/adv_first_eq.out' adv -ascii;

%% Plotting the data

f2=figure('visible','off');
g1=griddata(thresholds(:,1),thresholds(:,2),thresholds(:,3),thresholds(:,1),thresholds(:,2)');
h1=pcolor(thresholds(:,1),thresholds(:,2),g1);
h1.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Pore pressure threshold at nucleation point (MPa)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/threshold_plot.png");
saveas(f2, 'output/threshold_plot.fig');
close(f2);


f2=figure('visible','off');
g2=griddata(smoms(:,1),smoms(:,2),smoms(:,3),smoms(:,1),smoms(:,2)');
h2=pcolor(smoms(:,1),smoms(:,2),g2);
h2.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Average peak slip rate (log_{10}m/s)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/smom_plot.png");
saveas(f2,"output/smom_plot.fig");
close(f2);


f2=figure('visible','off');
g3=griddata(recavg(:,1),recavg(:,2),recavg(:,3),recavg(:,1),recavg(:,2)');
h3=pcolor(recavg(:,1),recavg(:,2),g3);
h3.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Average perturbed recurrence interval (years)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/recavg_plot.png");
saveas(f2,"output/recavg_plot.fig");
close(f2);


f2=figure('visible','off');
g4=griddata(asr(:,1),asr(:,2),asr(:,3),asr(:,1),asr(:,2)');
h4=pcolor(asr(:,1),asr(:,2),g4);
h4.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Ratio of seismic to aseismic events");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/asr_plot.png");
saveas(f2,"output/asr_plot.fig");
close(f2);


f2=figure('visible','off');
g5=griddata(naes(:,1),naes(:,2),naes(:,3),naes(:,1),naes(:,2)');
h5=pcolor(naes(:,1),naes(:,2),g5);
h5.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Number of aseismic transients");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/nase_plot.png");
saveas(f2,"output/nase_plot.fig");
close(f2);


f2=figure('visible','off');
g6=griddata(nses(:,1),nses(:,2),nses(:,3),nses(:,1),nses(:,2)');
h6=pcolor(nses(:,1),nses(:,2),g6);
h6.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Number of seismic events");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/nsse_plot.png");
saveas(f2,"output/nsse_plot.fig");
close(f2);


f2=figure('visible','off');
g7=griddata(moments(:,1),moments(:,2),moments(:,3),moments(:,1),moments(:,2)');
h7=pcolor(moments(:,1),moments(:,2),g7);
h7.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Aseismic Moment (log_{10} N/m^2)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/asmom_plot.png");
saveas(f2,"output/asmom_plot.fig");
close(f2);


f2=figure('visible','off');
g8=griddata(arecavg(:,1),arecavg(:,2),arecavg(:,3),arecavg(:,1),arecavg(:,2)');
h8=pcolor(arecavg(:,1),arecavg(:,2),g8);
h8.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Average perturbed recurrence interval (years)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/arecavg_plot.png");
saveas(f2,"output/arecavg_plot.fig");
close(f2);


f2=figure('visible','off');
g9=griddata(adv(:,1),adv(:,2),adv(:,3),adv(:,1),adv(:,2)');
h9=pcolor(adv(:,1),adv(:,2),g9);
h9.EdgeColor='none';
cb=colorbar;
xlabel("P_{max} (MPa)");
ylabel("r_P (Pa/s)");
ylabel(cb,"Advancement of first perturbed earthquake (years)");
colormap('plasma');
set(gca,'FontSize',30);
xlim([0.6,3]);
ylim([0.1,1]);
set(f2, 'visible', 'on');
set(f2, 'Position', get(0, 'Screensize'));
saveas(f2,"output/advancement_plot.png");
saveas(f2,"output/advancement_plot.fig")
close(f2);

%close(f1);

delete(gcp('nocreate'));