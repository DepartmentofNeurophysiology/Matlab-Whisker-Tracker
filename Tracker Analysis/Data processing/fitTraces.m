function Traces_fit = fitTraces(Traces)
%%


nframes = size(Traces,  1);
Traces_fit = cell(1,nframes);

h =  waitbar(0, 'Fitting traces...');
for i = 1:nframes
    ntraces = size(Traces{i},2);
    Tsave = {};
    for j = 1:ntraces
        trace = Traces{i}{j};
           t1 = [];
        t2 = [];
        t3 = [];
        fitax = [1:length(trace)]';
        px1 = polyfit(fitax, trace(:,1), 1);
        py1 = polyfit(fitax, trace(:,2), 1);
        t1(:,1) = polyval(px1, fitax);
        t1(:,2) = polyval(py1, fitax);
        
        px2 = polyfit(fitax, trace(:,1), 2);
        py2 = polyfit(fitax, trace(:,2), 2);
        t2(:,1) = polyval(px2, fitax);
        t2(:,2) = polyval(py2, fitax);
        
        px3 = polyfit(fitax, trace(:,1), 3);
        py3 = polyfit(fitax, trace(:,2), 3);
        t3(:,1) = polyval(px3, fitax);
        t3(:,2) = polyval(py3, fitax);
        
        r(1) = getRSlocal(trace(:,1), t1(:,1), 1) + getRSlocal(trace(:,2), t1(:,2), 1);
        r(2) = getRSlocal(trace(:,1), t2(:,1), 2) + getRSlocal(trace(:,2), t2(:,2), 2);
        r(3) = getRSlocal(trace(:,1), t3(:,1), 3) + getRSlocal(trace(:,2), t3(:,2), 3);
        
        [~, keep_idx] = max(r);
        
        switch(keep_idx)
            case 1
                Tsave{j} = t1;
            case 2
                Tsave{j} = t2;
            case 3 
                Tsave{j} = t3;
        end
    end
    
    Traces_fit{i} = Tsave;
    waitbar(i/nframes)
end
close(h)

end



function rscor = getRSlocal(Y, Yfit, fitdeg)
SSres = sum( (Y-Yfit).^2);
SStot = sum( (Y-mean(Y)).^2);
Rs = 1 - SSres/SStot;
rscor = 1-(1-Rs)*(length(Y)-1)/(length(Y)-fitdeg-1);
end
