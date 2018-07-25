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
% This script is associated to the Figure 2D of the manuscript 
% "Is it Is it possible to reconstruct an accurate cell lineage using CRISPR
% recorders?
% AUTHORS: Irepan Salvador-Martinez, Marco Grillo, Michalis Averof and
% Maximilian J Telford (preprint doi:https://doi.org/10.1101/37335)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% DESCRIPTION
% Explore the accuracy of tree reconstructions exploring different:
% mutation rate (from 0.01 to 0.3)
% number of character states (2, 4, 8, 16, 32)
% I will use 16 cell divisions, as this is the max number of divisions 
%-------------------------------------------------------------------------
% set the random seed to 123 so it can be reproducible
rng(123);

% Define the MAIN parameters
%-------------------------------------------------------------------------
% cells divisions
Ndiv = 16;

% Number of targets
targets = (10:10:150);

% Mutation rate for each number of divisions
Mu = (0.01:0.01:0.3);

% Number of simulations (repeats)
Nsim = 10;

% the number of states in the analysis. PAUP alows max 64 states.
s = 16;

% CREATE A TABLE WITH THE PROBABILITY OF MUTATION IN EACH TARGET,
% WHICH WILL BE THE SAME FOR EACH ONE
Nmer_prob = repmat(1/s,1,s);

%-------------------------------------------------------------------------
%
for t = targets

    for mu = Mu

        mutable = mu * cumsum(Nmer_prob);

        for r = 1:Nsim                      % TO HAVE MULTIPLE SAMPLES

            %the array is a 2D array with the targets (within the cell) as columns and
            % each row is a daughter cell (matrix changes in each iteration)
            % ----------Initialize some variables ------------------------------------
            sequences = struct('Header',{},'Sequence',{});
            cell_name = ones(1);
            A1 = zeros(1,t);
            
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
                A3 = char(zeros(size(A2)));
                
                % keep track of the cell names
                cell_name = horzcat(cell_name,cell_name);
                for cn = 1:Ndaug_1
                    cell_name(cn)   = (cell_name(cn) *2)-1;
                    cell_name(cn+Ndaug_1)   = cell_name(cn+Ndaug_1) *2;
                end
                
            % for each copied target, mutate randomly with a fixed probability
            % the site can mutate only once, into one of N possible states
            Aind = find(A2 == 0); AR = rand(size(Aind));
            AR = arrayfun(@(z)sum(z <= mutable), AR);
            A2(Aind) = AR;
            % make the daughter matrix the mother matrix for the next iteration
            A1 = A2;              
            end
            
            % ------------- IO --------------------------------------------
            % name of the simulation
            sim = [sprintf('%02d',t),'/',...
                sprintf('%02d',Ndiv),'div_',...
                sprintf('%02d',s),'states_',...
                sprintf('%02d',t),'_targets_'...
                sprintf('%0.2f',mu),'_Mu_'];
            
            %define the output file name for the alignment (fasta)
            out = [sim, sprintf('%04d', r),'rep.fas'];
            
            %---------- produce the output for each division!!
            %%%% Converts to ASCII characters, including 0-9, A-Z, a-z
            %%%% Avoids troublesome characters as ?, \, etc..

            %%% (a-x)
            %%% ----------------------------------------------------
            %A3(A2>=36 ) = char(A2(A2>=36) + 61);

            %   (A-Z)
            A3(A2>=10 & A2 <=35) = char(A2(A2>=10 & A2 <=35) + 55);

            %   (1-9)
            A3(A2< 10 & A2 >0 )= char(A2(A2 < 10 & A2 >0) + 48);

	    	% UNMUTATED- It can be changed to a mutated state (e.g., A) ALSO with perl 
            % inside the folder "perl -p -i -e 's/-/A/g' `find ./ -name *.fas`" 
            A3(A2==0)= char(A2(A2==0) + 48);
            %%% ----------------------------------------------------
            %clear A2 AR;
            
            % define the "true" tree with the number of final daughter cells
            Totdaugh = 2^n_div;

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
        end
    end
end
