classdef nme_line3p_tpc < wgv.nme_branch3p_tpc & wgv.form_tpc
% wgv.nme_line3p_tpc - Network model element for 3-phase line for TPC formulation
%
% Implements the network model element for 3-phase line elements, with
% 6 ports per 3-phase line.
%
% Implements building of the quadratic matrices stored in parameter Qu for 
% 3-phase lines and inherits from wgv.form_tpc.

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
        function name = name(obj)
            %
            name = 'line3p';
        end

        function np = np(obj)
            %
            np = 6;     %% this is a 6 port element
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);    %% call parent

            dme = obj.data_model_element(dm);
            bus_dme = dm.elements.(dme.cxn_type);
            nb3p = nm.node.idx.N.bus3p;         % number of 3-phase buses per phase
            nl = obj.nk;
            np = obj.np;           

            base_kv = bus_dme.tab.base_kv(bus_dme.on(dme.fbus(dme.on))) / sqrt(3);
            base_z = 1000 / dm.base_kva * base_kv .^ 2;

            sf = (base_z ./ dme.len) * ones(1, 6);      %% scale factor
            y_series = (dme.ys(dme.lc, :) .* sf).';     %% - of off-diagonal
            y_hshunt = (dme.yc(dme.lc, :)/2 .* sf).';   %% diagonal

            %% Compute the diagonal used both for parameters Wu and M
            
            phases = {'a', 'b', 'c'};
            rot_mat = [exp(1j*2*pi/3)    exp(-1j*2*pi/3)   exp(1j*2*pi/3)
                       exp(-1j*2*pi/3)   exp(1j*2*pi/3)    exp(-1j*2*pi/3)
                                       ];
            id_rot = [2  2  3
                      3  5  5];            

            vals_diag_from = {};
            ii_diag_from = {};
            jj_diag_from = {};
            for p = 1:length(phases)
                rot = ones(6,1);
                rot(id_rot(:,p)) = rot_mat(:,p);
                all_rot = repmat(rot,1,nl);
                [ii_aux, jj_aux, vals_aux] = ...
                    obj.make_diag(dme.fbus, dme.tbus, nb3p, all_rot.*y_series, all_rot.*y_hshunt, phases{p}, 'from');
                vals_diag_from = horzcat(vals_diag_from, {vals_aux});
                ii_diag_from = horzcat(ii_diag_from, {ii_aux});
                jj_diag_from = horzcat(jj_diag_from, {jj_aux});
            end

            vals_diag_to = {};
            ii_diag_to = {};
            jj_diag_to = {};
            for p = 1:length(phases)
                rot = ones(6,1);
                rot(id_rot(:,p)) = rot_mat(:,p);
                all_rot = repmat(rot,1,nl);
                [ii_aux, jj_aux, vals_aux] = ...
                    obj.make_diag(dme.fbus, dme.tbus, nb3p, all_rot.*y_series, all_rot.*y_hshunt, phases{p}, 'to');
                vals_diag_to = horzcat(vals_diag_to, {vals_aux});
                ii_diag_to = horzcat(ii_diag_to, ii_aux);
                jj_diag_to = horzcat(jj_diag_to, jj_aux);
            end
            
            %% 1) Wu parameter
            if nm.userdata.tpc.quad
                nval = 64;   %% number of nonzero values per matrix

                %% 1.1) Wu matrices - 'from' ports 
                
                % Compute the elements in the upper triangular
                [ii_triu, jj_triu, vals_triu_from] = obj.make_triu(vals_diag_from, ii_diag_from{1}, 'from');
    
                % Compute the elements on the lower triangular (symmetric matrices)
                vals_tril = vals_triu_from;
                ii_tril = jj_triu;
                jj_tril = ii_triu;
    
                % Compute final triplets (in the form of row-col-val)
                vals_aux = cell2mat(vertcat(vals_diag_from, vals_triu_from, vals_tril));
                %vals_aux = cell2sym(vertcat(vals_diag_from, vals_triu_from, vals_tril));
                ii_aux = cell2mat(vertcat(ii_diag_from, ii_triu, ii_tril));
                jj_aux = cell2mat(vertcat(jj_diag_from, jj_triu, jj_tril));
                
                vals_from = mat2cell(vals_aux(:), nval*ones(3*nl,1));
                % vals_from = cell(6,1);
                % for i = 1:6
                %     vals_from{i} = vals_aux(:,i);
                % end
                ii_from = mat2cell(ii_aux(:), nval*ones(3*nl,1));
                jj_from = mat2cell(jj_aux(:), nval*ones(3*nl,1));

                %% 1.2) Wu matrices - 'to' ports
                
                % Compute the elements in the upper triangular
                [ii_triu, jj_triu, vals_triu_to] = obj.make_triu(vals_diag_to, ii_diag_to{1}, 'to');
    
                % Compute the elements on the lower triangular (symmetric matrices)
                vals_tril = vals_triu_to;
                ii_tril = jj_triu;
                jj_tril = ii_triu;
    
                % Compute final triplets (in the form of row-col-val)
                vals_aux = cell2mat(vertcat(vals_diag_to, vals_triu_to, vals_tril));
                %vals_aux = cell2sym(vertcat(vals_diag_to, vals_triu_to, vals_tril));
                ii_aux = cell2mat(vertcat(ii_diag_to, ii_triu, ii_tril));
                jj_aux = cell2mat(vertcat(jj_diag_to, jj_triu, jj_tril));
                
                vals_to = mat2cell(vals_aux(:), nval*ones(3*nl,1));
                % vals_to = cell(6,1);
                % for i = 1:6
                %     vals_to{i} = vals_aux(:,i);
                % end
                ii_to = mat2cell(ii_aux(:), nval*ones(3*nl,1));
                jj_to = mat2cell(jj_aux(:), nval*ones(3*nl,1));
                
                %% 1.3) Build Wu parameter
                vals_Wu = [vals_from; vals_to];
                ii_Wu = [ii_from; ii_to];
                jj_Wu = [jj_from; jj_to];

                obj.Qu = mat2cell(cell2mat([ii_Wu jj_Wu vals_Wu]), nval*ones(6*nl,1));
            end
            
            %% 2) M parameter
            % Compute elements

            vals_diag = horzcat(vals_diag_from, vals_diag_to);
            vals_M = obj.make_vals_M(vals_diag, y_series, y_hshunt);

            % Compute indices 
            ii_M = repmat((1:nl*np), 2*np, 1);
            jj_M =  cell2mat(horzcat(ii_diag_from,  ii_diag_to));
            
            % Build M parameter
            obj.M = sparse(ii_M(:), jj_M(:), vals_M(:), nl*np, 2*sum(nb3p));
            
            %% 3) s parameter
            vals_s = [];

            rot_mat = [      1           exp(-1j*2*pi/3)   exp(1j*2*pi/3)
                       exp(1j*2*pi/3)         1            exp(-1j*2*pi/3)
                       exp(-1j*2*pi/3)   exp(1j*2*pi/3)          1];            
            
            for p = 1:length(phases)                
                aux_s = conj(sum(rot_mat(:,p).*y_hshunt(obj.p2i([1 2 3], [p p p]),:)));
                vals_s = [vals_s; aux_s(:)];
            end

            obj.s = repmat(vals_s, 2, 1);
        end

        function [ii, jj, vals] = make_diag(obj, fbus, tbus, nb3p, ys, ysh, phase, side)
            %
            % Inputs:
            %   ys  (matrix)      : 6 x nl matrix with the elements of the
            %                       series admitance matrix. Each column holds
            %                       the information of a 3-phase line
            %   ysh (matrix)      : 6 x nl matrix with the elements of the
            %                       shunt admitance matrix. Each column holds
            %                       the information of a 3-phase line
            %  phase (string)     : whether 'a', 'b' or 'c' identifier for the
            %                       phase to be considered when building the
            %                       diagonal elements of the quadratic matrices
            %   side (string)     : wheter 'from' or 'to' identifier for the
            %                       type of side being considered
            %
            % Outputs:
            %
    
            switch phase
                case {'a'}
                    p = 1;
                case {'b'}
                    p = 2;
                case {'c'}
                    p = 3;
            end
    
            id_phases = [1 2 3];
            not_p = setdiff(id_phases, p);
            nl = size(ys, 2);
            vals = zeros(12, nl);
            %vals = sym(zeros(12, nl));

            % 1) Compute values for quadratic matrices on the diagonal
    
            if strcmp(side, 'from')
                % theta_from - theta_from
                aux_id = obj.p2i(not_p, [p;p]);
                vals(p,:) = conj(ys(obj.p2i(p,p),:) - sum(ysh(aux_id,:)));
                vals(not_p,:) = -1 * conj(ys(aux_id,:) + ysh(aux_id,:));
    
                % theta_to - theta_to
                vals(id_phases+3,:) = conj(ys(obj.p2i([1 2 3], [p p p]),:));
    
                % lnvm_from - lnvm_from
                vals(p+6,:) = vals(p,:) + 2*conj(ys(obj.p2i(p,p),:) + 2*ysh(obj.p2i(p,p),:) + sum(ysh(aux_id,:)));
                vals(not_p+6,:) = -1  * vals(not_p,:);
    
                % lnvm_to - lnvm_to
                vals(id_phases+9,:) = -1 * vals(id_phases+3,:);
    
            elseif strcmp(side, 'to')
                % theta_from - theta_from
                vals(id_phases,:) = conj(ys(obj.p2i([1 2 3], [p p p]),:));
    
                % theta_to - theta_to
                aux_id = obj.p2i(not_p, [p;p]);
                vals(p+3,:) = conj(ys(obj.p2i(p,p),:) - sum(ysh(aux_id,:)));
                vals(not_p+3,:) = -1 * conj(ys(aux_id,:) + ysh(aux_id,:));
    
                % lnvm_from - lnvm_from
                vals(id_phases+6,:) = -1 * vals(id_phases,:);
    
                % lnvm_to - lnvm_to
                vals(p+9,:) = vals(p+3,:) + 2*conj(ys(obj.p2i(p,p),:) + sum(ysh(aux_id,:)) + 2*ysh(obj.p2i(p,p),:));
                vals(not_p+9,:) = -1  * vals(not_p+3,:);
            else
                error('wgv.nme_line3p_tpc:make_diag: side %s not recognized',side);
            end
    
    
            % 2) Compute row and column indices on the diagonal
            % 
            % The vector of variables is now: [\theta_1ph; lnvm_1ph; \theta_3ph; lnvm_3ph]
            
            % theta_from - theta_from
            Fbus = repmat(fbus', 3, 1);
            id_bias = repmat((id_phases'-1).*nb3p, 1, length(fbus));
            ii_tff = Fbus + id_bias;

            % theta_to - theta_to
            Tbus = repmat(tbus', 3, 1);
            ii_ttt = Tbus + id_bias;

            % lnvm_from - lnvm_from
            id_bias = repmat(3*nb3p + (id_phases'-1).*nb3p, 1, length(fbus));
            ii_lnvmff = Fbus + id_bias;

            % lnvm_to - lnvm_to
            ii_lnvmtt = Tbus + id_bias;
            
            % Compute row and column indices for one single phase
            ii = [ii_tff; ii_ttt; ii_lnvmff; ii_lnvmtt];
            jj = ii;
        end

        function vals_M = make_vals_M(obj, vals_diag, ys, yhs)
            %

            nl = obj.nk;

            vals_M = cell2mat(vals_diag);
            %vals_M = cell2sym(vals_diag);
            
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
