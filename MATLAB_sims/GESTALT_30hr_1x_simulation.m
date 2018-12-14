%    Copyright (C) 2017-2018 Irepan Salvador-Martínez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as
%    published by the Free Software Foundation, either version 3 of the
%    License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>
%
%    Contact: Irepan Salvador-Martinez <i.salvador@ucl.ac.uk>,
%    Department of Genetics, Evolution and Environment,
%    University College London,
%    Gower Street, London WC1E 6BT, England
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% This script is associated to the Figure 6B of the manuscript 
% "Is it possible to reconstruct an accurate cell lineage using CRISPR
% recorders?
% AUTHORS: Irepan Salvador-Martinez, Marco Grillo, Michalis Averof and
% Maximilian J Telford (preprint doi:https://doi.org/10.1101/37335)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% DESCRIPTION
% simulate the CRISPR mutations using parameters from the experimental data
% of GESTALT (McKenna et al 2016)
% for the explanation see the notebook :
% "GESTALT_30hr_1x_Mutation_rate_estimation.ipynb"
%-------------------------------------------------------------------------
% set the random seed to 123 so it can be reproducible
rng(123);

% Define the MAIN parameters
%-------------------------------------------------------------------------
% cells divisions
Ndiv = 16;        
% Number of targets
targets = 10;
% Number of simulations (repeats)
Nsim = 1000;

mut_rel = zeros(Nsim,targets);
mut_inter = zeros(Nsim,targets);
mut_intra = zeros(Nsim,targets);

Sat_total = zeros(Nsim,targets);
Sat_Norm = zeros(Nsim,targets);

Alleles = zeros(1,Nsim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MUTATION RATE FOR EACH TARGET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on Observed saturation
Mu = [0.135,0.044,0.026,0.117,0.087,0.105,0.232,0.022,0.152,0.007];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the number of states in the analysis. PAUP alows max 64 states.
states =60;
% create a matrix to store the number of mutations that occur in each cell division
muts = zeros(Ndiv,targets);
%-------------------------------------------------------------------------


for r = 1:Nsim                      % TO HAVE MULTIPLE SAMPLES
    
    %%%%%%% CREATE A TABLE WITH THE PROBABILITY OF MUTATION IN EACH TARGET,
    %%%%%%% WHICH WILL FOLLOW A GAMMA DISTRIBUTION with parameters 0.5 and 2
    MUT_TABLE_ALL = zeros(targets,states);
    for i = 1:targets
        %%%%%%%%%%%%%%%%%%  FOR EACH TARGET
        rvs = gamrnd(0.1,2, [1 states]);    rvs = rvs/sum(rvs);
        Nmer_prob = sort(rvs);
        MUT_TABLE_ALL(i,:) = Mu(i) * cumsum(Nmer_prob(1:states));
        MUT_TABLE_ALL(i,:) = sort(MUT_TABLE_ALL(i,:),'descend');
    end    
    
    %the array is a 2D array with the targets (within the cell) as columns and
    % each row is a daughter cell (matrix changes in each iteration)
    % ----------Initialize some variables ------------------------------------
    sequences = struct('Header',{},'Sequence',{});
    cell_name = ones(1);
    A1 = zeros(1,targets);
    
    %-----------------------------------------------------------------
    % -------------- MAIN LOOP ---------------------------------------
    %------------------------------------------------------------------
    for n_div = 1:Ndiv                           %for each division
        
        Ndaug_1 = 2^(n_div-1);
        Ndaugh  = 2^n_div;                            % calculate daughters number
        ncells_1 = size(A1,2);
                                      % to sum the mutations
        % create the matrix of daughters by copying the mother
        A2 = vertcat(A1,A1);
        A2_old = A2;
        %A3 = char(zeros(size(A2)));
        
        % keep track of the cell names
        cell_name = horzcat(cell_name,cell_name);
        for cn = 1:Ndaug_1
            cell_name(cn)   = (cell_name(cn) *2)-1;
            cell_name(cn+Ndaug_1)   = cell_name(cn+Ndaug_1) *2;
        end
        
        %%%%%%%%%%%%% FIRST ROUND OF CRISPR   %%%%%%%%%%%%%%%%%%%%%%%%%
        % for each copied target, mutate randomly with a fixed probability
        % SEPARATELY FOR EACH OF 10 TARGETS
        m1 = 0;
        %%%%%%%% FOR EACH TARGET 
        for i = 1:targets
            Aind = find(A2(:,i) == 0);
            AR = rand(size(Aind));
            
            AR = arrayfun(@(z)sum(z <= MUT_TABLE_ALL(i,:)), AR);
            
            A2(Aind,i) = AR;
            m1 = sum(AR > 0);
            % assign the numer of mutations to the matrix "muts"
            muts(n_div,i) = m1 /(Ndaugh);
            
        end
   
        %%%%%%%%% MAKE THE INTER-TARGETS HAPPEN  %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for j = 1:Ndaugh
            C_ne = find(A2(j,:)~=A2_old(j,:));
             if length(C_ne) > 1 
                C_ne = sort(datasample(C_ne,2,'Replace',false));
                % sum 100 to the targets at the ends (will modify to letters
                % later)
                A2(j,C_ne(1)) = A2(j,C_ne(1)) + 100;
                A2(j,C_ne(2)) = A2(j,C_ne(2)) + 100;
                
                % assign to the targets inbetween the state 200 
                A2(j,C_ne(1)+1:C_ne(2)-1) = 200;
                 
             end
        end      
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % make the daughter matrix the mother matrix for the next iteration
        A1 = A2;
        
        %end
        % ------------- IO --------------------------------------------
        % name of the simulation
        sim = ['GESTALT_',sprintf('%02d',n_div),'div_',...
            sprintf('%02d',targets),'_targets_'...
            sprintf('%02d',states),'_states_'];

        %define the output file name for the alignment (fasta)
        out = [sim, sprintf('%04d', r),'rep.fas'];

        %---------- produce the output for each division!!
        %%%% Converts to ASCII characters, including 0-9, A-Z, a-z
        %%%% Avoids troublesome characters as ?, \, etc..
        A2(A2>=36 & A2<100)=  A2(A2>=36 & A2<100) + 61;
        A2(A2>=10 & A2 <36) = A2(A2>=10 & A2 <36) + 55;
        A2(A2<10) = A2(A2<10) + 48;
        
        A2(A2>=136 & A2 <200) = A2(A2>=136 & A2 <200) -100 + 61;
        A2(A2>=110 & A2 <136) = A2(A2>=110 & A2 <136) -100 + 55;
        A2(A2>100 &  A2 <110) = A2(A2>100 &  A2 <110) -100 + 48;     
        
        %DELETED (MISSING) targets
        A2(A2==200)= 45;
        
                
        A3 = char(A2);
        
        % %INTRA-TARGET deletions
        % A3(A2>=1 & A2<=20)= char(A2(A2>=1 & A2<=20) + 64);
        % %INTRA-TARGET deletions
        % A3(A2>=21 & A2<=40)= char(A2(A2>=21 & A2<=40) + 76);
        % %DELETED (MISSING) targets
        % A3(A2==41)= char(45);
        % %UNMUTATED
        % A3(A2==0)= char(A2(A2==0) + 48);



        %%%%%% EXPORT THE SIMULATED FREQS TO A FILE %%%%%%%%%%%%

    %     unqA = unique(A3(1:numel(A3)));
    %     Nmers = cellstr(unqA');
    %     countNmers=histc(A3(1:numel(A3)),unqA);
    %     relFreq=countNmers/numel(A3);
    %     %transp(relFreq);
    %     T = table(Nmers, transp(relFreq));

        %clear A2 AR;

        % define the "true" tree with the number of final daughter cells
        Totdaugh = 2^n_div;
        Nbranch = Totdaugh -1;
        canonic = 1:(2*Nbranch);
        canonic = reshape(canonic,[2,Nbranch]);
        canonic = canonic';
        % create the name of the cells
        myfasta = fopen(['../simulations/',out],'w');
        fclose(myfasta);

        myfasta = fopen(['../simulations/',out],'a');
        for ii = 1:Totdaugh
            sequences(ii).Header = ['c_',sprintf('%07d',cell_name(ii))];
            sequences(ii).Sequence = strjoin(cellstr(A3(ii,:)));
        end
        %%%% SORT THE CELLS IN THE OUTPUT SO IT IS EASY TO READ
        [~,index] = sortrows({sequences.Header}.');
        sequences = sequences(index); clear index

        for ii = 1:Totdaugh
            % append to a file
            fprintf(myfasta,'%s\t',sequences(ii).Header);
            fprintf(myfasta,'%s\n',sequences(ii).Sequence);
        end
        fclose(myfasta);
        
    
    %define the output file name for the reference tree (newick)
    RefOut = [sim, 'REF.nw'];
    % construct a nw file to compare it with the simulated ones
    RefTree = phytree(canonic);
    d = get(RefTree,'distances');
    p = get(RefTree,'pointers');
    n = sort({sequences.Header});
    RefTree2 = phytree(p,d,n);
    %view(RefTree2)
    % write newick files to compare with ete software externally
    phytreewrite(['../simulations/',RefOut], RefTree2, ...
        'Distances', 'false','BranchNames','false');

    end
    %%%%%%%%%%% ANALYSIS OF THE DATA
    for i = 1:targets
        Sat_total(r,i) = sum(A1(:,i)>0)/length(A1(:,i));
        Sat_Norm(r,i)  = sum(A1(:,i)>0 & A1(:,i) <200) / (sum(A1(:,i)==0) +sum(A1(:,i)>0 & A1(:,i)<200));
        mut_rel(r,i)   = sum(A1(:,i)>0 & A1(:,i) <200)/length(A1(:,i));
        mut_intra(r,i) = sum(A1(:,i)>0 & A1(:,i) <100)/length(A1(:,i));
        mut_inter(r,i) = sum(A1(:,i)>100 & A1(:,i)<200)/length(A1(:,i)); 
    end
    
    %%%% determine the number of alleles with 100 samples
    allele_samples = zeros(1,100);
    for s = 1:100
        allele_samples(1,s) = length(unique(datasample(A3,10000,'Replace',false),'rows'));
    end
    %%%%%% Number of alleles
    Alleles(1,r) = median(allele_samples);
    
end

%%
%%%%%%%%%%% PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Title of the whole fig
f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'GESTALT 30hr 1x (1000 simulations)'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

    
% plot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1,'Parent',p);
plot(median(mut_rel)/sum(median(mut_rel)), 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'Marker','o','LineWidth',1.5,'Color',[0 0 0]);
hold on ;
% Titlr and labels
xlabel('Target'); ylabel('Mutated proportion');
title('Saturation (simulation)');

plot(median(mut_inter)/sum(median(mut_rel)), 'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 0.5 0],...
    'Marker','v','LineWidth',1,'Color',[1 0.5 0]);
plot(median(mut_intra)/sum(median(mut_rel)), 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
    'Marker','^','LineWidth',1,'Color',[0 0 1]);
hold off;
title('Inter + Intra Mutations');

% Total Saturation figure %%%%%%%%%%%%%%%%%%%%
%figure1 = figure('Name','Target Saturation');
subplot(1,3,2,'Parent',p) 
hold('on');
% Create plot
boxplot(Sat_Norm,'PlotStyle','compact','Symbol','','Color',[0 0 0])
plot(median(Sat_Norm,'omitnan'),'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],...
    'Marker','o','LineWidth',2,'Color',[1 0 1]);
% Titlr and labels
xlabel('Target'); ylabel('Mutated proportion');
title('Saturation (simulation)');
hold('off');


% Total Saturation figure %%%%%%%%%%%%%%%%%%%%
%figure1 = figure('Name','Target Saturation');
subplot(1,3,3,'Parent',p) 
hold('on');
% Create plot
boxplot(Sat_total,'PlotStyle','compact','Symbol','','Color',[0 0 0])
plot(median(Sat_total),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'Marker','o','LineWidth',2,'Color',[1 0 0]);
% Titlr and labels
xlabel('Target'); ylabel('Mutated proportion');
title('Saturation (simulation)');
hold('off');

%%
% PLOT THE FREQ OF EACH ALLELE
% [b,i,j] = unique(A3,'rows')
% bn = histc(j,1:length(b))
% semilogy(sort(bn))
x= mean(Alleles);
%%%% PLOT NUMBER OF ALLELES
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
% Create line
hist(Alleles)
line([x x],[0 250],'Parent',axes1,'DisplayName','Mean',...
    'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
% Create xlabel
xlabel('# Alleles');
% Create title
title({'15 div; 5 states (1K simulations)','(mean #Alleles from 100 samples of 10K cells)'});
% Create ylabel
ylabel('Freq');
box(axes1,'on');
% Create legend
legend(axes1,'show');
