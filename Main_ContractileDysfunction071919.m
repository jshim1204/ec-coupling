% arrhythmogenic risk & contractile dysfunction evalution 
% main function
% utility functions are found in the added path

% JHS
% 073018 

clear all;close all;clc


format compact

% Function Path
addpath('./Library/')
addpath('./UMI_data/')



% % Defining color table for intuitive display
igorct = zeros(256,3) ;
igorct(1,:) = 0 ;
igorct(256,:) = 0 ;
igorct(1:128,1) = 1 ;
igorct(1:128,2) = (1:128)/128 ;
igorct(1:128,3) = (1:128)/128 ;
igorct(129:256,1) = 1 - (1:128)/128 ;
igorct(129:256,2) = 1 - (1:128)/128 ;
igorct(129:256,3) = 1 ;

%% data loading

% load gene exp data
% plate={'plate6'};
% file='plate6.csv'; 
% plate6= importdata(file);
% genes=plate6.textdata(2:end);
% plate6=plate6.data;


% % ------------TMM method testing ---------------- % 
% verified to be the same as umi counts -- JHS 072518 saved output in
% output folder
% 
% % load gene exp data
% plate={'plate6-TMM'};
% file='normalized_plate6.csv';
% plate6=importdata(file);
% genes=plate6.textdata(2:end);
% plate6=plate6.data; 
% 
% % ------------------------------------------------ %
% % drugs & cell combo
% plate6_design='plate6_design.csv'; 
% plate6_design=readtable(plate6_design);
% plate6_design=table2cell(plate6_design);
% cell6=plate6_design(:,1);
% drug6=plate6_design(:,2);

% % NEW STAR aligned------------------------------------------------ %
plate={'plate6_STAR_UMIcounts'};
file='plate6_STAR_UMIcounts.csv'; 
plate6= importdata(file);
genes=plate6.textdata(2:end,1);
plate6=plate6.data;

% % ------------------------------------------------ %
% drugs & cell combo
plate6_design='plate6_design_STAR.csv'; 
plate6_design=readtable(plate6_design);
plate6_design_cell=strrep(plate6_design.iPSC_Line,'-','_');
plate6_design=[plate6_design_cell,plate6_design.State];
cell6=plate6_design(:,1);
drug6=plate6_design(:,2);

%rename the data -- these will be used throughout the script
data=plate6; % UMI count
CellLine_all=strrep(cell6,'-','_');
Drug_all=drug6;
genes=genes;







%% Data extraction

result=ExtractData(data,genes,CellLine_all,Drug_all,igorct,1);
[MSN02_4_CM_drugs,MSN09_4_CM_drugs,cellMSN02_4_CM_triplicates,cellMSN09_4_CM_triplicates]=result{:};

%% plot actinx myosin -- see function for other options
% plot_actinxmyosin('crossprodct_scalefactr_distr', cellMSN02_4_CM_triplicates,cellMSN09_4_CM_triplicates)

%% Execute Simulation

cell_var_pull=unique(CellLine_all);

for i=1:length(cell_var_pull)
    
    druglist=eval(strcat(cell_var_pull{i},'_drugs'));
    varname=strcat('cell',cell_var_pull{i},'_triplicates');
    triplicates_data=eval(varname);
    
   
    
    tic 
    
    for k=1:3
        switch k
            case 1
                stim='IKr_block';
                Kr_scalePercent_all=linspace(0.5,1,10); % from 50% block to 100%
                gradient_stim=Kr_scalePercent_all;
                
            case 2
                stim='hypokalemia';
                K0_scalePercent_all=linspace(0,0.7,10); % 0-70% reduction % increased from 0-50 bc most drugs had 0.5 as threshold
                gradient_stim=K0_scalePercent_all;
                
            case 3
                stim='extra_Ca';
                Ca_scalePercent_all=linspace(1,10,10); % 1-10X increase
                gradient_stim=Ca_scalePercent_all;
                
        end
        
        
        for kk=1:10 % 10 different stims
            
            switch k
                case 1
                    Kr_var=Kr_scalePercent_all(kk);
                    K0_var=0;
                    Ca_var=1;
                    
                case 2
                    Kr_var=0;
                    K0_var=K0_scalePercent_all(kk);
                    Ca_var=1;
                    
                case 3
                    Kr_var=0;
                    K0_var=0;
                    Ca_var=Ca_scalePercent_all(kk);
            end
            
            n_iterations=length(druglist);
            
            
            parfor ii=1 :n_iterations
                scalefactors=triplicates_data(1:13,ii);
                XB_scalefactors=(1+( (triplicates_data(14,ii) *triplicates_data(15,ii)) -1 )*.5);
                
                %% secondary insult
                
                Ko_scalePercent=K0_var; % percent 0.15,0.25, 0.35
                Ca_scalePercent=Ca_var; % X2, X3
                Kr_scalePercent=Kr_var; % 30 percent, 50 percent, 70 percent reduction
                
                constants=[Ko_scalePercent,Ca_scalePercent,Kr_scalePercent];
                
                ODE_results =RunODEfunctions(scalefactors,XB_scalefactors,constants);
                [time_data,statevalues]=ODE_results{:};
                
                %% save ouptut
                outputcell_all{ii}=statevalues;
                t_all{ii}=time_data;
            end
            
            %% save per diff stim level
            var_output=strcat(stim,'_',cell_var_pull{i},'_output_',num2str(kk));
            eval([ var_output '=outputcell_all;' ])
            
           
            var_time=strcat(stim,'_',cell_var_pull{i},'_time_',num2str(kk));
            eval([ var_time '=t_all;' ])
            
            %% save each output in cell array
            output_all_levels{kk}=outputcell_all;
            time_all_levels{kk}=t_all;
            
        end
        
        
        %% combine outputs in workspace and save it in .mat file
        
        var_combined=strcat(stim,'_',cell_var_pull{i},'_all');
        eval([var_combined '= output_all_levels;' ])
        
        var_combined=strcat(stim,'_',cell_var_pull{i},'_time');
        eval([var_combined '= time_all_levels;' ])
        
        
        saveODE_output(stim,cell_var_pull{i},output_all_levels,time_all_levels); % saved in the output folder
        
        
        %% compute APD 
        
        APD_CaTauc= Detect_APD(output_all_levels,time_all_levels,stim,gradient_stim,cell_var_pull{i},druglist,{'NIL','SUN','PON','TRS','LAP'},0,[]); %last input arg = plotting option for drug list , % CaT_plot_yes,CaT_plot_list
        
        [APD_stim,CaT_AUC]=APD_CaTauc{:};
       
        %% compute threshold 
        
        threshold_stim=computeThreshold(APD_stim,druglist,gradient_stim); % output is vertical
        threshold_all(:,k)=threshold_stim; % col 1-3: ikr, hypokalemia, ca aug
       
        %% compute SL shortening
        if k==2 % hypokalemia or extraCA bc both start from baseline
        
        SL_stim=Detect_SL(output_all_levels,time_all_levels,stim,cell_var_pull{i},druglist, []); % last space is for the list of drug for simulation plots
%         SL_stim_all(:,k)=SL_stim; % col 1-3: ikr, hypokalemia, ca aug
        end
        
        %% save APD & CaT AUC
        
        var_APD=strcat(cell_var_pull{i},'_',stim,'_APD');
        eval([ var_APD '=APD_stim;'])
        var_CaT=strcat(cell_var_pull{i},'_',stim,'_CaT');
        eval([ var_CaT '=CaT_AUC;'])
        
    end
        
        %% save threshold
        var_threshold=strcat(cell_var_pull{i},'_threshold');
        eval([ var_threshold '= threshold_all;' ]) % save threshold per stim per cell line (col-stim, rows-drugs)
        saveODE_output('threshold',cell_var_pull{i},threshold_all); % saved in the output folder
        
        %% save SL shortening
        var_SL=strcat(cell_var_pull{i},'_SL_shortening');
        eval([ var_SL '= SL_stim;' ])
        saveODE_output('SL',cell_var_pull{i},SL_stim);
        
                
    toc
    
        % compute arrhythmogenic index - EAD susceptibility ranking
        % generates plot
        
        arrhy_idx=computeEADsusceptibility(threshold_all,cell_var_pull{i},druglist,1,[1:10],1); % third arg : 1- plot, 0- no plot
        
        % save output
        saveRankings(cell_var_pull{i},arrhy_idx);
        
        
        [d3_all, sorted_log_d3, I_d3, arrhythmoScores_threeStim] = arrhy_idx{:};
        % d3_all : euclidean distance for each drug
        % sorted_log_d3= neg. log normalized(wrt ctrl) score
        % I_d3= sorted index for original druglist
        
        hypokal_score=arrhythmoScores_threeStim(:,1);
        extraCa_score=arrhythmoScores_threeStim(:,2);
        IKrblock_score=arrhythmoScores_threeStim(:,3);

 
        % contractile dysfunction- failure ranking
        
        triplicates=eval(strcat('cell',cell_var_pull{i},'_triplicates'));
        list_fail=generateCaT_actinXmyosin_analysis(.75,SL_stim,CaT_AUC(:,1),triplicates,cell_var_pull{i},druglist,1,{});
        plotContractileDysfunction(cell_var_pull{i},druglist,SL_stim,[1:10]); %[1:10] top 10 rankings 
    
        
end



