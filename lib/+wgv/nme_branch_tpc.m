classdef nme_branch_tpc < mp.nme_branch & wgv.form_tpc
% wgv.nme_branch_tpc - Network model element for branch for TPC formulations.
%
% Implements building of the branch parameters :math:`\WQu`, math:`\M`, and
% :math:`\s`, and inherits from :class:`wgv.form_tpc`.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build_params(obj, nm, dm)% IN PROGRESS - WGV
            %
            build_params@mp.nme_branch(obj, nm, dm);    %% call parent

            dme = obj.data_model_element(dm);
            nl = obj.nk;                     %% number of branches
            if isfield(nm.node.idx.N,'bus')  %% total number of 1ph buses
                nb1p = nm.node.idx.N.bus;
            else
                nb1p =  0;
            end            
            np = obj.np;                     %% number of ports for individual branch 

            tm = ones(nl, 1);                %% default tap ratio = 1
            i = find(dme.tm);                %% indices of non-zero tap ratios
            tm(i) = dme.tm(i);               %% assign non-zero tap ratios
            T = tm .* exp(1j * dme.ta);      %% add phase shifters

            ys = 1 ./ (dme.r + 1j * dme.x);     %% series admittance
            yff = (ys + dme.g_fr + 1j * dme.b_fr) ./ (T .* conj(T));
            ytt = ys + dme.g_to + 1j * dme.b_to;
            yft = - ys ./ conj(T);
            ytf = - ys ./ T;
            
            %% 1) Qu parameter
            if nm.userdata.tpc.quad
                % Get row and column indices 
                [If, Jf, It, Jt] = get_indices_Qu(obj, nb1p, dme);
    
                % Qu matrices - 'from' ports 
                ww = [-1 1 1j 1j].';
                weights = [ww; -ww; -1j*ww; -1j*ww];
                weights = repmat(weights, 1, nl);
                Wvals_f = repmat(yft', (2*np)^2, 1);
                Wvals_f = Wvals_f.*weights;
                id_offset_f = (11:(2*np)^2:(nl*(2*np)^2)-5);
                offset_Wvals_f = spalloc(nl*(2*np)^2, 1, nl);                
                offset_Wvals_f(id_offset_f) = 4*yff';
                Wvals_f =  Wvals_f(:) + offset_Wvals_f;

                % Qu matrices - 'to' ports 
                ww = [-1 1 -1j -1j].';
                weights = [ww; -ww; 1j*ww; 1j*ww];
                weights = repmat(weights, 1, nl);
                Wvals_t = repmat(ytf', (2*np)^2, 1);
                Wvals_t = Wvals_t.*weights;
                id_offset_t = (16:(2*np)^2:(nl*(2*np)^2));                
                offset_Wvals_t = spalloc(nl*(2*np)^2, 1, nl);
                offset_Wvals_t(id_offset_t) = 4*ytt';
                Wvals_t =  Wvals_t(:) + offset_Wvals_t;               
    
                % Build Qu parameter
                Wvals = [Wvals_f; Wvals_t];                
                obj.Qu = mat2cell([[If Jf; It Jt] Wvals(:)], (2*np)^2 * ones(nl*np,1));
            end

            %% 2) M parameter
            % Get row and column indices
            [If, Jf, It, Jt] = get_indices_M(obj, nb1p, dme);

            % M matrix - 'from' ports
            ww = [1j -1j 1 1];
            weights = repmat(ww, nl, 1);
            weights = weights(:).';
            Mvals_f = repmat(yft', 1, 2*np);
            Mvals_f = weights.*Mvals_f;
            id_offset_f = (2*nl+1:3*nl);            
            offset_Wvals_f = spalloc(1, nl*2*np, nl);
            offset_Wvals_f(id_offset_f) = 2*yff';
            Mvals_f = (Mvals_f + offset_Wvals_f).';

            % M matrix - 'to' ports
            ww = conj(ww);
            weights = repmat(ww, nl, 1);
            weights = weights(:).';
            Mvals_t = repmat(ytf', 1, 2*np);
            Mvals_t = weights.*Mvals_t;
            id_offset_t = (3*nl+1:4*nl);            
            offset_Wvals_t = spalloc(1, nl*2*np, nl);
            offset_Wvals_t(id_offset_t) = 2*ytt';
            Mvals_t = (Mvals_t + offset_Wvals_t).';

            % Build M parameter            
            obj.M = sparse([If; It],[Jf; Jt],[Mvals_f; Mvals_t], nl*np, 2*nb1p);

            %% 3) s parameter
            obj.s = [(yff'+yft').'; (ytf'+ytt').'];
        end    
        
        function [If, Jf, It, Jt] = get_indices_Qu(obj, nb, dme)
            % Get the row and column indices for building Qu parameter
            % using sparse command
            %
            % Inputs:                        
            %   nb (integer) : number of system buses (nm.node.N)
            %   dme (mp.data_model) : data model element object
            %
            % Outputs:
            %   If (vector): row indices used to create quadratic matrices
            %                for ports related with from end of branches
            %   Jf (vector): column indices used to create quadratic matrices
            %                for ports related with from end of branches
            %   It (vector): row indices used to create quadratic matrices
            %                for ports related with to end of branches
            %   Jt (vector): column indices used to create quadratic matrices
            %                for ports related with to end of branches

            nl = obj.nk;                           %% Number of branches
            np = obj.np;                           %% Number of ports
            ft_row = [dme.fbus(dme.on) dme.tbus(dme.on)]';         %% Get the from-to indices as row vectors            

            %% Indices for from end of branches
            ii = [ft_row ; ft_row + nb];
            ii = repmat(ii(:), 1, 2*np);
            jj = ii';
            Jf = jj(:);

            ii = mat2cell(ii, 2*np*ones(nl,1));
            ii = ii';
            ii = cell2mat(ii);
            If = ii(:);

            %% Indices for to end of branches (the same as those for from end)
            It = If; 
            Jt = Jf;
        end

        function [If, Jf, It, Jt] = get_indices_M(obj, nb, dme)
            % Get the row and column indices for building M parameter
            % using sparse command
            %
            % Inputs:           
            %   dme (mp.data_model) : data model object element
            %   nb (integer) : number of system buses (nm.node.N)
            %
            % Outputs:
            %   If (vector): row indices used to create entries of M matrix
            %                for ports related with from end of branches
            %   Jf (vector): column indices used to create entries of M matrix
            %                for ports related with from end of branches
            %   It (vector): row indices used to create entries of M matrix
            %                for ports related with to end of branches
            %   Jt (vector): column indices used to create entries of M matrix
            %                for ports related with to end of branches

            nl = obj.nk;                           %% Number of branches            
            np = obj.np;                           %% Number of ports
            ft_col = [dme.fbus(dme.on) dme.tbus(dme.on)];          %% Get the from-to indices as column vectors (as in the mp case)
            
            %% Indices for from end of branches
            ii = repmat((1:nl)', 1, 2*np);
            If = ii(:);

            jj = [ft_col  (ft_col + nb)];
            Jf = jj(:);
            
            %% Indices for to end of branches (the same as those for from end)
            It = If + nl; 
            Jt = Jf;
        end

        %
    end     %% methods
end         %% classdef
