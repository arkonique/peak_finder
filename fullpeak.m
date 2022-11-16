function [fpeak,spans,signals] = fullpeak(fname,rates,times,ts,te,plotting,aevtime,sevtime,savedir,alimit,slimit,background,lag,threshold,influence)
    arguments
        fname
        rates
        times
        ts
        te
        plotting=0
        aevtime=864000 % 10 days
        sevtime=300
        savedir='output/vmax_plots/'
        alimit=1e-8
        slimit=1e-3
        background=1e-9
        lag=50
        threshold=2
        influence=0.05
    end
    [~,~]=mkdir(savedir);
    spans=[]; % spans are regions where the signal is a peak, extracted from the thresholding signal
    r1=rates(times>=ts & times<=te);
    t1=times(times>=ts & times<=te);
    [signals,~,~]=ThresholdingAlgo(r1,lag,threshold,influence); % run moving threshold algorithm
    signals(signals<0)=0;
    % find spans
    i=1;
    while i<=length(signals)
        if (signals(i)>0)
            ist=i;
            j=i+1;
            while j<=length(signals)
                if (signals(j)==0)
                    ied=j;
                    i=j;
                    break;
                end
                j=j+1;
            end
            spans=[spans;ist,ied];
        end
        i=i+1;
    end

    window=aevtime;

     % expand spans to be as long as a window
    for i=1:numel(spans)/2
        width=t1(spans(i,2))-t1(spans(i,1));
        if (width<window)
            midpoint=(t1(spans(i,2))+t1(spans(i,1)))/2;
            window_start=midpoint-window/2;
            window_end=midpoint+window/2;
            t1(spans(i,1));
            t1(spans(i,2));
            indexes=find(t1>=window_start & t1<=window_end);
            indexes(1);
            indexes(end);
            spans(i,1)=indexes(1);
            spans(i,2)=indexes(end);
        end
    end

     % join spans if they overlap
     for i=1:numel(spans)/2-1 % combine spans
        st1=spans(i,1);
        ed1=spans(i,2);
        for j=i+1:numel(spans)/2
            st2=spans(j,1);
            ed2=spans(j,2);
            if(st2>=st1 && st2<=ed1 && ed2>=ed1)
                st2=st1;
                ed1=ed2;
                spans(i,2)=ed1;
                spans(j,1)=st2;
            end
        end
    end
    
    spans=unique(spans);
    
    spans=reshape(spans,[2 length(spans)/2])';
    
    fpeak=[];

    % The max value within each span should be a peak
    for i=1:numel(spans)/2 % all aseismic peaks should be here
        [c,~]=max(r1(spans(i,1):spans(i,2)));
        tp=t1(r1==c);
        if (c>=slimit)
            l=1;
        else
            l=0;
        end
        if (c>=alimit && c<=slimit) % only save aseismic peaks because there might be multiple seismic peaks within one window
            fpeak=[fpeak;tp,c,l];
        end
    end


    % distance filtering for any clusters still left
    for i=1:numel(fpeak)/3-1
        d=sqrt((fpeak(i,1)/3.154e7-fpeak(i+1,1)/3.154e7)^2+(log10(fpeak(i,2))-log10(fpeak(i+1,2)))^2);
        if (d<0.05 && fpeak(i,2)<1e-3) % the limit for d is something that can be tweaked but I am keeping it out of the users' hands for now
            fpeak(i,3)=2;
        end
    end
    if (numel(fpeak)>0)
        fpeak(fpeak(:,3)==2,:)=[];
    end
    
    % seismic peaks
    [peak,~] = peaks(rates,times,ts,te,alimit,slimit,aevtime,sevtime,background,0);
    if (numel(peak)>0)
        peak=peak(peak(:,3)==1,:);
    end
    fpeak=[fpeak;peak];
    

    if (plotting==1)
        f=figure('Visible','off');
        f.Renderer='painter';
        plot(t1/3.154e7,log10(r1));
        hold on
        if (numel(fpeak)>0)
            h=plot(fpeak(:,1)/3.154e7,log10(fpeak(:,2)),'.');
            h.MarkerSize=40;
            ttl=['Seismic, aseismic = ' num2str(numel(fpeak(fpeak(:,3)==1))) ',' num2str(numel(fpeak(fpeak(:,3)==0)))];
            
        else
            ttl="Seismic, aseismic = 0,0";
        end
        title(ttl);
        xlabel('Time (years)')
        ylabel('log_{10} V_{max} (m/s)')
        hold off;
        legend('Original Data','Peaks')
        set(gca,'FontSize',30);
        set(f, 'Position', get(0, 'Screensize'));
        name=[savedir '/' num2str(fname) '.png'];
        saveas(f,name,'png');
        % close(f);
    end
end