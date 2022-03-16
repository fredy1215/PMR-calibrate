% relexsation 10 sec senstivity map calculation
[Image,Path] = uigetfile('C:\Users\lvmen\Box\phantom\xie10_28\Xie_25881_2891\Biri_Research_Yibin - 91bd136439ec4f\*.dcm');
DCM_I = dicomread(fullfile(Path,Image));
dcm_info = dicominfo(fullfile(Path,Image));
DCM_I = double(DCM_I);
[N1,N2] = size(DCM_I);
for q = 1:N1
    for p = 1:N2
        if DCM_I(q,p)<800
            DCM_I_Masked(q,p) = 0;
        else
            DCM_I_Masked(q,p) = DCM_I(q,p);
        end
    end
end
I_median = median(nonzeros(DCM_I_Masked));
DCM_I_correction_value = DCM_I_Masked./I_median;
%%                  
for q = 1:10
    imagesc(DCM_I_Masked)
    colormap jet
    ROI(:,:,q) = roipoly;
end
close
%%
for p = 1:4
    Fpath =uigetdir('C:\Users\lvmen\Box\phantom\xie10_28\Xie_25881_2891\Biri_Research_Yibin - 91bd136439ec4f\normalized');
    cd(Fpath)
    imagelist = dir('*.dcm');
    findletter = isletter(imagelist(1).name);
    CATA = findletter.*imagelist(1).name;
    CATA = nonzeros(CATA);
    CATA = char(CATA(1:end-3));
    CATA = upper(sprintf('%c',CATA));
    fcount = numel(imagelist);
    Rawdata(p).catalog = CATA;
        for q = 1:fcount
            I = double(dicomread(imagelist(q).name));
            findletter = isletter(imagelist(q).name);
            for a = 1:numel(findletter)
                if findletter(a)~=0
                    notletter(a) = 0;
                else
                    notletter(a) = 1;
                end
            end
            condition = notletter.*imagelist(q).name;
            condition = char(nonzeros(condition));
            condition = sscanf(condition,'%d');
            Rawdata(p).condition(q) = condition;
            for qq = 1:10
                Tube(:,:,qq) = ROI(:,:,qq).*I;
                Tube_correctedNAN = ROI(:,:,qq).*(Tube(:,:,qq)./DCM_I_correction_value);
                Tube_correctedNAN(isnan(Tube_correctedNAN)) = 0;
                Tube_corrected(:,:,qq) = Tube_correctedNAN;
                Avgintensity(qq) = mean(nonzeros(Tube(:,:,qq)));
                Avgintensity_corrected(qq) =  mean(nonzeros(Tube_corrected(:,:,qq)));
            end
            Before_Senstivity_Correct = I;
            After_Senstivity_Correct = I./DCM_I_correction_value;
            After_Senstivity_Correct(isnan(After_Senstivity_Correct)) = 0;
            Rawdata(p).AvgTubeIntensity(q,:) = Avgintensity;
            Rawdata(p).AvgTubeIntensity_corrected(q,:) = Avgintensity_corrected;
            clear I findletter notletter condition
        end
    clear imagelist Fpath CATA Fcount
end
%% Scanner PMR
for p = 1:4
    ScannerPMR(p).catalog = Rawdata(p).catalog;
    [x,y] = size(Rawdata(p).AvgTubeIntensity_corrected);
    for q = 1:x
        for qq = 1:y
        ScannerPMR(p).value(q,qq) = Rawdata(p).AvgTubeIntensity_corrected(q,qq)/Rawdata(p).AvgTubeIntensity_corrected(q,10);
        end
    end
end
%% Predict T1
syms T1_pl
for p = 1:4
    Predict_T1(p).catalog = Rawdata(p).catalog;
    switch Predict_T1(p).catalog
        case 'FA'
[x,y] = size(ScannerPMR(p).value);
for q = 1:x
    for qq = 1:y
TI =600;
echos = 5;
n = 30;
HR = 65;
TR = 1/HR*60*1000;
flipA = Rawdata(p).condition(q);
T1_myo = 1143.9;
M0 = 1;
PMR = ScannerPMR(p).value(q,qq);
eq = PMR == abs(M0*(exp(-TI/T1_pl) - 1) - (exp(-TI/T1_pl)*(M0*(exp((echos*n - TR + TI)...
/T1_pl) - 1) + exp((echos*n - TR + TI)/T1_pl)*(M0*exp(-(echos*n)/T1_pl)*cos((pi*flipA)...
/180)^n*(exp(-TI/T1_pl) - 1) + (M0*(exp(-echos/T1_pl) - 1)*((exp(-echos/T1_pl)*cos((pi*flipA)/180))^n - 1))...
/(exp(-echos/T1_pl)*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_pl)*exp(-TI/T1_pl)*exp((echos*n - TR + TI)...
/T1_pl)*cos((pi*flipA)/180)^n + 1))/abs(M0*(exp(-TI/T1_myo) - 1) - (exp(-TI/T1_myo)*(M0*(exp((echos*n - TR + TI)...
/T1_myo) - 1) + exp((echos*n - TR + TI)/T1_myo)*(M0*exp(-(echos*n)/T1_myo)*cos((pi*flipA)/180)^n*(exp(-TI/T1_myo)...
- 1) + (M0*(exp(-echos/T1_myo) - 1)*((exp(-echos/T1_myo)*cos((pi*flipA)/180))^n - 1))/(exp(-echos/T1_myo)...
*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_myo)*exp(-TI/T1_myo)*exp((echos*n - TR + TI)/T1_myo)...
*cos((pi*flipA)/180)^n + 1));
solution = vpasolve(eq,T1_pl,[0 3000],'random',true);
Predict_T1(p).value(q,qq) = double(solution);
    end
end
    case 'HR'
[x,y] = size(ScannerPMR(p).value);
for q = 1:x
    for qq = 1:y
TI =600;
echos = 5;
n = 30;
HR = Rawdata(p).condition(q);
TR = 1/HR*60*1000;
flipA = 15;
T1_myo = 1143.9;
M0 = 1;
PMR = ScannerPMR(p).value(q,qq);
eq = PMR == abs(M0*(exp(-TI/T1_pl) - 1) - (exp(-TI/T1_pl)*(M0*(exp((echos*n - TR + TI)...
/T1_pl) - 1) + exp((echos*n - TR + TI)/T1_pl)*(M0*exp(-(echos*n)/T1_pl)*cos((pi*flipA)...
/180)^n*(exp(-TI/T1_pl) - 1) + (M0*(exp(-echos/T1_pl) - 1)*((exp(-echos/T1_pl)*cos((pi*flipA)/180))^n - 1))...
/(exp(-echos/T1_pl)*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_pl)*exp(-TI/T1_pl)*exp((echos*n - TR + TI)...
/T1_pl)*cos((pi*flipA)/180)^n + 1))/abs(M0*(exp(-TI/T1_myo) - 1) - (exp(-TI/T1_myo)*(M0*(exp((echos*n - TR + TI)...
/T1_myo) - 1) + exp((echos*n - TR + TI)/T1_myo)*(M0*exp(-(echos*n)/T1_myo)*cos((pi*flipA)/180)^n*(exp(-TI/T1_myo)...
- 1) + (M0*(exp(-echos/T1_myo) - 1)*((exp(-echos/T1_myo)*cos((pi*flipA)/180))^n - 1))/(exp(-echos/T1_myo)...
*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_myo)*exp(-TI/T1_myo)*exp((echos*n - TR + TI)/T1_myo)...
*cos((pi*flipA)/180)^n + 1));
solution = vpasolve(eq,T1_pl,[0 3000],'random',true);
Predict_T1(p).value(q,qq) = double(solution);
    end
end
    case 'TI'
[x,y] = size(ScannerPMR(p).value);
for q = 1:x
    for qq = 1:y
TI = Rawdata(p).condition(q);
echos = 5;
n = 30;
HR = 65;
TR = 1/HR*60*1000;
flipA = 15;
T1_myo = 1143.9;
M0 = 1;
PMR = ScannerPMR(p).value(q,qq);
eq = PMR == abs(M0*(exp(-TI/T1_pl) - 1) - (exp(-TI/T1_pl)*(M0*(exp((echos*n - TR + TI)...
/T1_pl) - 1) + exp((echos*n - TR + TI)/T1_pl)*(M0*exp(-(echos*n)/T1_pl)*cos((pi*flipA)...
/180)^n*(exp(-TI/T1_pl) - 1) + (M0*(exp(-echos/T1_pl) - 1)*((exp(-echos/T1_pl)*cos((pi*flipA)/180))^n - 1))...
/(exp(-echos/T1_pl)*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_pl)*exp(-TI/T1_pl)*exp((echos*n - TR + TI)...
/T1_pl)*cos((pi*flipA)/180)^n + 1))/abs(M0*(exp(-TI/T1_myo) - 1) - (exp(-TI/T1_myo)*(M0*(exp((echos*n - TR + TI)...
/T1_myo) - 1) + exp((echos*n - TR + TI)/T1_myo)*(M0*exp(-(echos*n)/T1_myo)*cos((pi*flipA)/180)^n*(exp(-TI/T1_myo)...
- 1) + (M0*(exp(-echos/T1_myo) - 1)*((exp(-echos/T1_myo)*cos((pi*flipA)/180))^n - 1))/(exp(-echos/T1_myo)...
*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_myo)*exp(-TI/T1_myo)*exp((echos*n - TR + TI)/T1_myo)...
*cos((pi*flipA)/180)^n + 1));
solution = vpasolve(eq,T1_pl,[0 3000],'random',true);
Predict_T1(p).value(q,qq) = double(solution);
    end
end
    case 'SEG'
[x,y] = size(ScannerPMR(p).value);
for q = 1:x
    for qq = 1:y
TI = 600;
echos = 5;
n = Rawdata(p).condition(q);
HR = 65;
TR = 1/HR*60*1000;
flipA = 15;
T1_myo = 1143.9;
M0 = 1;
PMR = ScannerPMR(p).value(q,qq);
eq = PMR == abs(M0*(exp(-TI/T1_pl) - 1) - (exp(-TI/T1_pl)*(M0*(exp((echos*n - TR + TI)...
/T1_pl) - 1) + exp((echos*n - TR + TI)/T1_pl)*(M0*exp(-(echos*n)/T1_pl)*cos((pi*flipA)...
/180)^n*(exp(-TI/T1_pl) - 1) + (M0*(exp(-echos/T1_pl) - 1)*((exp(-echos/T1_pl)*cos((pi*flipA)/180))^n - 1))...
/(exp(-echos/T1_pl)*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_pl)*exp(-TI/T1_pl)*exp((echos*n - TR + TI)...
/T1_pl)*cos((pi*flipA)/180)^n + 1))/abs(M0*(exp(-TI/T1_myo) - 1) - (exp(-TI/T1_myo)*(M0*(exp((echos*n - TR + TI)...
/T1_myo) - 1) + exp((echos*n - TR + TI)/T1_myo)*(M0*exp(-(echos*n)/T1_myo)*cos((pi*flipA)/180)^n*(exp(-TI/T1_myo)...
- 1) + (M0*(exp(-echos/T1_myo) - 1)*((exp(-echos/T1_myo)*cos((pi*flipA)/180))^n - 1))/(exp(-echos/T1_myo)...
*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_myo)*exp(-TI/T1_myo)*exp((echos*n - TR + TI)/T1_myo)...
*cos((pi*flipA)/180)^n + 1));
solution = vpasolve(eq,T1_pl,[0 3000],'random',true);
Predict_T1(p).value(q,qq) = double(solution);
    end
end
    end
end

%% Coverted PMR
for p = 1:4
    PMR_Covert(p).catalog = Rawdata(p).catalog;
    [x,y] = size(Predict_T1(p).value);
     for q = 1:x
        for qq = 1:y
            TI = 600;
            echos = 5;
            n = 30;
            HR = 65;
            TR = 1/HR*60*1000;
            flipA = 15;
            T1_myo = 1143.9;
            T1_pl = Predict_T1(p).value(q,qq);
            M0 = 1;
            PMR = abs(M0*(exp(-TI/T1_pl) - 1) - (exp(-TI/T1_pl)*(M0*(exp((echos*n - TR + TI)...
                /T1_pl) - 1) + exp((echos*n - TR + TI)/T1_pl)*(M0*exp(-(echos*n)/T1_pl)*cos((pi*flipA)...
                /180)^n*(exp(-TI/T1_pl) - 1) + (M0*(exp(-echos/T1_pl) - 1)*((exp(-echos/T1_pl)*cos((pi*flipA)/180))^n - 1))...
                /(exp(-echos/T1_pl)*cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_pl)*exp(-TI/T1_pl)*exp((echos*n - TR + TI)...
                /T1_pl)*cos((pi*flipA)/180)^n + 1))/abs(M0*(exp(-TI/T1_myo) - 1) - (exp(-TI/T1_myo)*(M0*(exp((echos*n - TR + TI)...
                /T1_myo) - 1) + exp((echos*n - TR + TI)/T1_myo)*(M0*exp(-(echos*n)/T1_myo)*cos((pi*flipA)/180)^n*(exp(-TI/T1_myo)...
                - 1) + (M0*(exp(-echos/T1_myo) - 1)*((exp(-echos/T1_myo)*cos((pi*flipA)/180))^n - 1))/(exp(-echos/T1_myo)...
                *cos((pi*flipA)/180) - 1))))/(exp(-(echos*n)/T1_myo)*exp(-TI/T1_myo)*exp((echos*n - TR + TI)/T1_myo)...
                *cos((pi*flipA)/180)^n + 1));
            PMR_Covert(p).value(q,qq) = PMR;
        end
    end
end

%% Bland Altman ([fa,hr,seg,TI])
PMR_True = [2.96227499069875,2.11056369259461,1.58334254096310,1.27837604174621,1.03886811324446,0.907521215097836,0.744889204157644,0.725237306775392,0.650223846059464,1];
T1_True = [367.46, 556.17, 744.76, 932.97, 1143.9, 1310.9, 1524.4, 1583.6, 1767.5,1143.9];
for p = 1:4
    Bland_PMR(p).catalog = Rawdata(p).catalog;
    Bland_T1(p).catalog = Rawdata(p).catalog;
    MSE(p).catalog = Rawdata(p).catalog;
    [x,y] = size(Predict_T1(p).value);
    for q = 1:x
         Bland_PMR(p).diffvalue(q,:) = PMR_Covert(p).value(q,:)- PMR_True;
         Bland_T1(p).diffvalue(q,:) = Predict_T1(p).value(q,:)- T1_True;
         Bland_PMR(p).meanvalue(q,:) = (PMR_Covert(p).value(q,:)+ PMR_True)/2;
         Bland_T1(p).meanvalue(q,:) = (Predict_T1(p).value(q,:)+ T1_True)/2;
         Bland_PMR(p).diffstd(q) = std(Bland_PMR(p).diffvalue(q,:));
         Bland_T1(p).diffstd(q) = std(Bland_T1(p).diffvalue(q,:));
         Bland_PMR(p).upperbound(q) = mean( Bland_PMR(p).diffvalue(q,:))+1.96*Bland_PMR(p).diffstd(q);
         Bland_T1(p).upperbound(q) = mean( Bland_T1(p).diffvalue(q,:))+1.96*Bland_T1(p).diffstd(q);
         Bland_PMR(p).lowerbound(q) = mean( Bland_PMR(p).diffvalue(q,:))-1.96*Bland_PMR(p).diffstd(q);
         Bland_T1(p).lowerbound(q) = mean( Bland_T1(p).diffvalue(q,:))-1.96*Bland_T1(p).diffstd(q);
         MSE(p).PMR(q) = sum(Bland_PMR(p).diffvalue(q,:).*Bland_PMR(p).diffvalue(q,:))/9;
         MSE(p).T1(q) = sum(Bland_T1(p).diffvalue(q,:).*Bland_T1(p).diffvalue(q,:))/9;
         PE(p).PMR(q,:) = abs(Bland_PMR(p).diffvalue(q,:))./PMR_True;
         PE(p).T1(q,:) = abs(Bland_T1(p).diffvalue(q,:))./T1_True;
    end
end
for p = 1:4
    [x,y] = size(Predict_T1(p).value);
    Bland_all_PMR(p).diffvalue = [];
    Bland_all_T1(p).diffvalue = [];
    Bland_all_PMR(p).meanvalue = [];
    Bland_all_T1(p).meanvalue = [];
    TruePMR_Mesh{p} = [];
    TrueT1_Mesh{p} = [];
    PE_PMR{p} = [];
    
    for q = 1:x
        Bland_all_PMR(p).diffvalue = [Bland_all_PMR(p).diffvalue,Bland_PMR(p).diffvalue(q,:)];
        Bland_all_T1(p).diffvalue = [Bland_all_T1(p).diffvalue,Bland_T1(p).diffvalue(q,:)];
        Bland_all_PMR(p).meanvalue = [Bland_all_PMR(p).meanvalue,Bland_PMR(p).meanvalue(q,:)];
        Bland_all_T1(p).meanvalue = [Bland_all_T1(p).meanvalue,Bland_T1(p).meanvalue(q,:)];
        TruePMR_Mesh{p} = [TruePMR_Mesh{p},PMR_True];
        TrueT1_Mesh{p} = [TrueT1_Mesh{p},T1_True];
        
    end
     Bland_all_PMR(p).diffvalue(abs(Bland_all_PMR(p).diffvalue)<0.0001) = 0;
    diff{p} = Bland_all_PMR(p).diffvalue;
    Bland_all_PMR(p).meandiff = mean(Bland_all_PMR(p).diffvalue);
    Bland_all_PMR(p).diffstd = std(Bland_all_PMR(p).diffvalue);
    Bland_all_T1(p).diffstd = std(Bland_all_T1(p).diffvalue);
    Bland_all_PMR(p).upperbound = mean( Bland_all_PMR(p).diffvalue)+1.96*Bland_all_PMR(p).diffstd;
    Bland_all_T1(p).upperbound = mean( Bland_all_T1(p).diffvalue)+1.96*Bland_all_T1(p).diffstd;
    Bland_all_PMR(p).lowerbound = mean( Bland_all_PMR(p).diffvalue)-1.96*Bland_all_PMR(p).diffstd;
    Bland_all_T1(p).lowerbound = mean( Bland_all_T1(p).diffvalue)-1.96*Bland_all_T1(p).diffstd;
    Bland_all_PMR(p).MPE = mean(nonzeros(abs(Bland_all_PMR(p).diffvalue./TruePMR_Mesh{p})));
    Bland_all_T1(p).MPE = mean(nonzeros(abs(Bland_all_T1(p).diffvalue./TrueT1_Mesh{p})));
    PE_PMR{p} = nonzeros(abs(diff{p})./TruePMR_Mesh{p});
end
%% Scanner PMR Bland altman
PMR_True = [2.96227499069875,2.11056369259461,1.58334254096310,1.27837604174621,1.03886811324446,0.907521215097836,0.744889204157644,0.725237306775392,0.650223846059464,1];
for p = 1:4
    Bland_PMRS(p).catalog = Rawdata(p).catalog;   
    [x,y] = size(Predict_T1(p).value);
    for q = 1:x
         Bland_PMRS(p).diffvalue(q,:) = ScannerPMR(p).value(q,:)- PMR_True;
         Bland_PMRS(p).meanvalue(q,:) = (ScannerPMR(p).value(q,:)+ PMR_True)/2;
    end
end
for p = 1:4
    [x,y] = size(ScannerPMR(p).value);
    Bland_all_PMRS(p).diffvalue = [];
    Bland_all_PMRS(p).meanvalue = [];
    TruePMR_Mesh{p} = [];
    PE_PMRS{p} = [];
    for q = 1:x
        Bland_all_PMRS(p).diffvalue = [Bland_all_PMRS(p).diffvalue,Bland_PMRS(p).diffvalue(q,:)];
        Bland_all_PMRS(p).meanvalue = [Bland_all_PMRS(p).meanvalue,Bland_PMRS(p).meanvalue(q,:)];
        TruePMR_Mesh{p} = [TruePMR_Mesh{p},PMR_True];
    end
    
    Bland_all_PMRS(p).diffvalue(abs(Bland_all_PMRS(p).diffvalue)<0.0001) = 0;
    diffs{p} = Bland_all_PMRS(p).diffvalue;
    Bland_all_PMRS(p).meandiff = mean(Bland_all_PMRS(p).diffvalue);
    Bland_all_PMRS(p).diffstd = std(Bland_all_PMRS(p).diffvalue);
    Bland_all_PMRS(p).upperbound = mean( Bland_all_PMRS(p).diffvalue)+1.96*Bland_all_PMRS(p).diffstd;
    Bland_all_PMRS(p).lowerbound = mean( Bland_all_PMRS(p).diffvalue)-1.96*Bland_all_PMRS(p).diffstd;
    PE_PMRS{p} = nonzeros(abs(diffs{p})./TruePMR_Mesh{p});
    Bland_all_PMRS(p).MPE = mean(nonzeros(abs(Bland_all_PMRS(p).diffvalue./TruePMR_Mesh{p})));
end
%% ICC
% for k = 1:4
%     M = Predict_T1(k).value';
% [r_T1(k), LB_T1(k), UB_T1(k), F_T1(k), df1_T1(k), df2_T1(k), p_T1(k)] = ICC(M,'1-1');
% end
% for k = 1:4
%     M = PMR_Covert(k).value';
% [r_PMR(k), LB_PMR(k), UB_PMR(k), F_PMR(k), df1_PMR(k), df2_PMR(k), p_PMR(k)] = ICC(M,'1-1');
% end
%% ttest
for q = 1:3
    [h_PMR(q), p_PMR(q)] = ttest(PE_PMRS{q},PE_PMR{q});
    Ma(q) = max(abs(PE_PMR{q}));
    Mi(q) = min(abs(PE_PMR{q}));
    Mas(q) = max(abs(PE_PMRS{q}));
    Mis(q) = min(abs(PE_PMRS{q}));
end
%% Pearson Correlation
for p = 1:4
    [x,~] = size(Predict_T1(p).value);
    T1_prediction{p} = [];
    True_T1{p} = [];
    for q = 1:x   
        T1_prediction{p} = [T1_prediction{p},Predict_T1(p).value(q,:)];
        True_T1{p} = [True_T1{p},T1_True]; 
    end
    [PearsonT1_c{p},PearsonT1_p{p}] = corr(transpose(T1_prediction{p}),transpose(True_T1{p}));
end
for p = 1:4
    [x,~] = size(PMR_Covert(p).value);
    CovertedPMR{p} = [];
    True_PMR{p} = [];
    for q = 1:x   
        CovertedPMR{p} = [CovertedPMR{p},PMR_Covert(p).value(q,:)];
        True_PMR{p} = [True_PMR{p},PMR_True]; 
    end
    [PearsonPMR_c{p},PearsonPMR_p{p}] = corr(transpose(CovertedPMR{p}),transpose( True_PMR{p}));
end
for p = 1:4
    [x,~] = size(ScannerPMR(p).value);
    CovertedPMRS{p} = [];
    True_PMRS{p} = [];
    for q = 1:x   
        CovertedPMRS{p} = [CovertedPMRS{p},ScannerPMR(p).value(q,:)];
        True_PMRS{p} = [True_PMRS{p},PMR_True]; 
    end
    [PearsonPMRS_c{p},PearsonPMRS_p{p}] = corr(transpose(CovertedPMRS{p}),transpose( True_PMRS{p}));
end
%%
for q = 1:3
figure;
plot(True_PMRS{q},CovertedPMRS{q},'b*')
xlabel('Standard PMR','FontSize',14)
ylabel('Uncalibrated PMR','FontSize',14)
hold on

plot(0:1900,0:1900,'r-')
hold off
xlim([0 5])
ylim([0 5])
text(1,4,sprintf('r^2 = %f p < 0.0001',PearsonPMRS_c{q}^2));
end
%%
for q = 1:3
figure;
plot(True_PMR{q},CovertedPMR{q},'b*')
xlabel('Standard PMR','FontSize',14)
ylabel('Calibrated PMR','FontSize',14)
hold on
plot(0:5,0:5,'r-')
hold off
xlim([0 5])
ylim([0 5])
text(1,4,sprintf('r^2 = %f p < 0.0001',PearsonPMR_c{q}^2));
end
%%
figure; imagesc(Before_Senstivity_Correct)
figure; imagesc(After_Senstivity_Correct)

%%
for q = 1:3
figure;
plot(Bland_all_PMR(q).meanvalue,Bland_all_PMR(q).diffvalue,'b*')
hold on
plot(0:4, ones(5)*Bland_all_PMR(q).upperbound,'r--')
hold on
plot(0:4, ones(5)*Bland_all_PMR(q).lowerbound,'r--')
hold on
plot(0:4, ones(5)*Bland_all_PMR(q).meandiff,'b-')
xlim([0 4])
ylim([-0.85 0.85])
xlabel('Mean between PMR__calibrated and PMR__standard','FontSize',18)
ylabel('Differences between PMR__calibrated and PMR__standard','FontSize',18)
hold off
end
%cd('C:\Users\lvmen\Box\Thesis\code')
%%
for q = 1:3
figure;
plot(Bland_all_PMRS(q).meanvalue,Bland_all_PMRS(q).diffvalue,'b*')
hold on
plot(0:4, ones(5)*Bland_all_PMRS(q).upperbound,'r--')
hold on
plot(0:4, ones(5)*Bland_all_PMRS(q).lowerbound,'r--')
hold on
plot(0:4, ones(5)*Bland_all_PMRS(q).meandiff,'b-')
xlim([0 4])
ylim([-0.85 0.85])
xlabel('Mean between PMR__uncalibrated and PMR__standard','FontSize',18)
ylabel('Differences between PMR__uncalibrated and PMR__standard','FontSize',18)
hold off
end
%cd('C:\Users\lvmen\Box\Thesis\code')
%% table data

%%
cd(Path)
save('Phantom_anaylsis_data_normal')

