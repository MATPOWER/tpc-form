classdef nme_branch3p_tpc < mp.nm_element
% wgv.nme_branch3p_tpc - Network model element for 3-phase branches for TPC 
%                        formulation.
%
% Implements methods to support building the quadratic matrices stored in 
% parameter Qu for 3-phase branches such as lines and transformers.

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
        function id = p2i(obj, p1, p2)
            %
            % Define some scalars for indexing the information of the
            % series and shunt matrices of 3-phase lines
    
            if length(p1) ~= length(p2)
                error('wgv.nme_line3p_tpc:p2i: inputs p1 and p2 must be same size')
            end
    
            id_matrix = [1  2  3
                2  4  5
                3  5  6];
    
            if isscalar(p1) && isscalar(p2)
                id = id_matrix(p1, p2);
            else
                n = length(p1);
                id = zeros(1,n);
                for i = 1:n
                    id(i) = id_matrix(p1(i),p2(i));
                end
            end
        end
    
        function [ii, jj, vals] = make_triu(obj, vals_diag, id_var_diag, side)
            %
            % Inputs:
            %   ys  (matrix)  : 6 x nl matrix with the elements of the
            %                   series admitance matrix. Each column holds
            %                   the information of a 3-phase line
            %   ysh (matrix)  : 6 x nl matrix with the elements of the
            %                   shunt admitance matrix. Each column holds
            %                   the information of a 3-phase line
            %   phase (string): wheter 'a', 'b' or 'c' identifier for the
            %                   phase to be considered when building the
            %                   diagonal elements of the quadratic matrices
            %
            % Outputs:
            %
    
            %nl = size(vals_diag{1},2);
    
            % 1) Elments of Group 1: theta-theta (except pivot-pivot)
    
            if strcmp(side, 'from')
                %       phase         id_proj        types
                sce = {  'a',         (2:6)',        'all'
                    'b',      [1 3 4 5 6]',     'all'
                    'c',      [1 2 4 5 6]',     'all'
                    };
    
            elseif strcmp(side, 'to')
                %       phase         id_proj        types
                sce = {  'a',      [1 2 3 5 6]'      'all'
                    'b',      [1 2 3 4 6]',     'all'
                    'c',         (1:5)',        'all'
                    };
            else
                error('wgv.nme_line3p_tpc:make_triu: side %s not recognized',side);
            end
    
            vals_from = {};
            ii_from = {};
            jj_from = {};
    
            for k = 1:size(sce,2)
                [ii_aux, jj_aux, vals_aux] = obj.project(vals_diag{k}, id_var_diag, sce{k,1}, sce{k,2}, sce{k,3}, side);
                vals_from = horzcat(vals_from, vals_aux);
                ii_from = horzcat(ii_from, ii_aux);
                jj_from = horzcat(jj_from, jj_aux);
            end
    
            % 2) Elments of Group 2: lnvm-lnvm (except pivot-pivot)
    
            if strcmp(side, 'from')
                %       phase       id_proj             types
                sce = {  'a',       (8:12)',      {'self','nonself'}
                    'b',   [7 9 10 11 12]',  {'self','nonself'}
                    'c',   [7 8 10 11 12]',  {'self','nonself'}
                    };
            elseif strcmp(side, 'to')
                %       phase       id_proj             types
                sce = {  'a',   [7 8 9 11 12]'   {'self','nonself'}
                    'b',   [7 8 9 10 12]',  {'self','nonself'}
                    'c',       (7:11)',     {'self','nonself'}
                    };
            else
                error('wgv.nme_line3p_tpc:make_triu: side %s not recognized',side);
            end
    
            vals_to = {};
            ii_to = {};
            jj_to = {};
    
            for k = 1:size(sce,2)
                [ii_aux, jj_aux, vals_aux] = obj.project(vals_diag{k}, id_var_diag, sce{k,1}, sce{k,2}, sce{k,3}, side);
                vals_to = horzcat(vals_to, vals_aux);
                ii_to = horzcat(ii_to, ii_aux);
                jj_to = horzcat(jj_to, jj_aux);
            end
    
            % 3) Elements of pivot-pivot:
    
            if strcmp(side, 'from')
                vals_pp = {-1j*vals_diag{1}(1,:) -1j*vals_diag{2}(2,:) -1j*vals_diag{3}(3,:)};
                ii_pp = {id_var_diag(1,:) id_var_diag(2,:) id_var_diag(3,:)};
                jj_pp = {id_var_diag(7,:) id_var_diag(8,:) id_var_diag(9,:)};
            elseif strcmp(side, 'to')
                vals_pp = {-1j*vals_diag{1}(4,:) -1j*vals_diag{2}(5,:) -1j*vals_diag{3}(6,:)};
                ii_pp = {id_var_diag(4,:) id_var_diag(5,:) id_var_diag(6,:)};
                jj_pp = {id_var_diag(10,:) id_var_diag(11,:) id_var_diag(12,:)};
            else
                error('wgv.nme_line3p_tpc:make_triu: side %s not recognized',side);
            end
    
            % Finally, compute whole set of indices
            vals = vertcat(vals_from, vals_to, vals_pp);
            ii = vertcat(ii_from, ii_to, ii_pp);
            jj = vertcat(jj_from, jj_to, jj_pp);
        end
    
        function [ii, jj, vals] = project(obj, vals_diag, id_var_diag, phase, id_proj, types, side)
            %
            % Inputs:
            %   vals_diag (matrix) : values computed using the method
            %                        obj.make_diag for a given phase
            %       phase (string) : wheter 'a', 'b' or 'c' identifier for
            %                        the  phase to be considered when building
            %                        the upper-triangular values of the
            %                        quadratic matrices.
            %     id_proj (vector) : vector of indices in [1:6] or [7:12]
            %                        denoting the positions on the diagonal
            %                        to be projected in order to calculate
            %                        the upper-triangular values of the
            %                        quadratic matrices.
            %         dir (string) : whether 'up' or 'right' directions of
            %                        projection are supported. The
            %                        reference of the projection is the
            %                        pivot, that is, the element of the
            %                        diagonal indicated by the phase under
            %                        analysis
            %         types (cell) : cell array whose elements can be 'self',
            %                        'nonself' or 'peer' indicating projection
            %                        using the self pivot, the nonself pivot or
            %                        a projection to the peer of the element,
            %                        respectively. The pivot is determined
            %                        through the phase while the peer is the
            %                        element on the other group that is stored
            %                        in the analogous position. Group 1 goes
            %                        from column 1 to 6 and Group 2 goes from
            %                        colum 7 to 12. Using 'all' indicates
            %                        the user ask for projection though all
            %                        the three directions.
            %
            %
            % Outputs:
            %
    
            if any(id_proj > 12) || numel(id_proj) ~= numel(unique(id_proj))
                error('wgv.nme_line3p_tpc:poject: input id_proj must be formed as non repeated ids between 0 and 12')
            end
            if any(id_proj > 6) && any(id_proj <=6)
                error('wgv.nme_line3p_tpc:poject: id_proj must a vector of indices in [1:6] or [7:12], but not both')
            end
            if size(id_proj, 2) > 1
                id_proj = id_proj';
            end
    
            vals_to_proj = vals_diag(id_proj,:);
    
            if any(id_proj <= 6)
                group = 0;
            else
                group = 6;
            end
    
            if strcmp(types, 'all')
                types = {'self'; 'nonself'; 'peer'};
            end
    
            % 1) Project values
    
            vals = {};
    
            for t = 1:length(types)
                type = types{t};
                switch type
                    case {'self'}
                        if group == 0
                            vals_aux = -1 * vals_to_proj;
                        else
                            vals_aux =  vals_to_proj;
                        end
                    case {'nonself'}
                        vals_aux = 1j * vals_to_proj;
                    case {'peer'}
                        if group == 0
                            vals_aux = 1j * vals_to_proj;
                        else
                            vals_aux = -1j * vals_to_proj;
                        end
                end
                vals = vertcat(vals, {vals_aux});
            end
    
            % 2) Compute indices of projected values
    
            nl = size(vals_diag,2);
    
            if strcmp(side, 'from')
                switch phase
                    case {'a'}
                        piv = 1;
                    case {'b'}
                        piv = 2;
                    case {'c'}
                        piv = 3;
                end
            elseif strcmp(side, 'to')
                switch phase
                    case {'a'}
                        piv = 4;
                    case {'b'}
                        piv = 5;
                    case {'c'}
                        piv = 6;
                end
            else
                error('wgv.nme_line3p_tpc:project: side %s not recognized',side);
            end
    
            below_pivot = id_proj(id_proj > piv + group);
            above_pivot = id_proj(id_proj < piv + group);
    
            ii = {};
            jj = {};
    
            for t = 1:length(types)
                type = types{t};
                if group == 0     %% (theta-theta)
                    %ii_above = repmat(above_pivot, 1, nl);
                    ii_above = id_var_diag(above_pivot, :);
                    jj_above = [];
                    switch type
                        case {'self'}
                            %jj_above = repmat(piv, length(above_pivot), nl);
                            if ~isempty(ii_above)
                                jj_above = repmat(id_var_diag(piv,:), length(above_pivot), 1);
                            end
                            %ii_below = repmat(piv, length(below_pivot), nl);
                            ii_below = repmat(id_var_diag(piv,:), length(below_pivot), 1);
                            %jj_below = repmat(below_pivot, 1, nl);
                            jj_below = id_var_diag(below_pivot,:);
                        case {'nonself'}
                            %jj_above = repmat(piv+6, length(above_pivot), nl);
                            if ~isempty(ii_above)
                                jj_above = repmat(id_var_diag(piv+6,:), length(above_pivot), 1);
                            end
                            %ii_below = repmat(below_pivot, 1, nl);
                            ii_below = id_var_diag(below_pivot,:);
                            %jj_below = repmat(piv+6, length(below_pivot), nl);
                            jj_below = repmat(id_var_diag(piv+6,:), length(below_pivot), 1);
                        case {'peer'}
                            %jj_above = repmat(above_pivot+6, 1, nl);
                            if ~isempty(ii_above)
                                jj_above = id_var_diag(above_pivot+6,:);
                            end
                            %ii_below = repmat(below_pivot, 1, nl);
                            ii_below = id_var_diag(below_pivot,:);
                            %jj_below = repmat(below_pivot+6, 1, nl);
                            jj_below = id_var_diag(below_pivot+6,:);

                    end
                else            %% (lnvm-lnvm)
                    jj_below = id_var_diag(below_pivot, :);
                    ii_below = [];
                    switch type
                        case {'self'}
                            % ii_below = repmat(piv+6, length(below_pivot), nl);
                            if ~isempty(jj_below)
                                ii_below = repmat(id_var_diag(piv+6,:), length(below_pivot), 1);
                            end
                            %ii_above = repmat(above_pivot, 1, nl);
                            ii_above = id_var_diag(above_pivot,:);
                            %jj_above = repmat(piv+6, length(above_pivot), nl);
                            jj_above = repmat(id_var_diag(piv+6,:), length(above_pivot), 1);
                        case {'nonself'}
                            % ii_below = repmat(piv, length(below_pivot), nl);
                            if ~isempty(jj_below)
                                ii_below = repmat(id_var_diag(piv,:), length(below_pivot), 1);
                            end
                            % ii_above = repmat(piv, length(above_pivot), nl);
                            ii_above = repmat(id_var_diag(piv,:), length(above_pivot), 1);
                            %jj_above = repmat(above_pivot, 1, nl);
                            jj_above = id_var_diag(above_pivot,:);
                        case {'peer'}
                            %ii_below = repmat(above_pivot-6, 1, nl);
                            if ~isempty(jj_below)
                                ii_below = id_var_diag(above_pivot-6,:);
                            end
                            %ii_above = repmat(above_pivot-6, 1, nl);
                            ii_above = id_var_diag(above_pivot-6,:);
                            %jj_above = repmat(above_pivot, 1, nl);
                            jj_above = id_var_diag(above_pivot,:);
                    end
                end
    
                ii = vertcat(ii, {[ii_above; ii_below]});
                jj = vertcat(jj, {[jj_above; jj_below]});
            end
        end

        function vals_M = make_vals_M(obj, vals_diag, ys, yhs)
            %

            nl = obj.nk;

            vals_M = cell2mat(vals_diag);
            
            for p = 1:3       %% for each phase (a, b, and c)
                not_p_from = setdiff((1:6), p);
                not_p_to = setdiff((1:6), p+3);

                % Correct values on nonself elements of Group 1 (\theta - \theta)
                vals_M(not_p_from, (p-1)*nl+1 : p*nl) = 1j* vals_M(not_p_from, (p-1)*nl+1 : p*nl);      % from ports
                vals_M(not_p_to, (p+2)*nl+1 : (p+3)*nl) = 1j* vals_M(not_p_to, (p+2)*nl+1 : (p+3)*nl);  % to ports

                % Correct values on pivots of Group 1 (\theta - \theta)
                vals_M(p, (p-1)*nl+1 : p*nl) = -1j * vals_M(p, (p-1)*nl+1 : p*nl);                      % from ports
                vals_M(3+p, (p+2)*nl+1 : (p+3)*nl) = -1j * vals_M(3+p, (p+2)*nl+1 : (p+3)*nl);          % to ports

                % Correct values on pivots of Group 2 (lnvm-lnvm)
                correction = 2*conj(ys(obj.p2i(p,p),:) + yhs(obj.p2i(p,p),:));
                vals_M(6+p, (p-1)*nl+1 : p*nl) = vals_M(6+p, (p-1)*nl+1 : p*nl) - correction;           % from ports
                vals_M(9+p, (p+2)*nl+1 : (p+3)*nl) = vals_M(9+p, (p+2)*nl+1 : (p+3)*nl) - correction;   % to ports
            end
        end
    end     %% methods
end         %% classdef