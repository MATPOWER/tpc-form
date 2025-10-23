classdef nme_xfmr3p_tpc < mp.nm_element & wgv.form_tpc
% wgv.nme_xfmr3p_tpc - Network model element for 3-phase transformer for
%                      TPC formulation.
%
% Implements the network model element for 3-phase transformer elements,
% with 6 ports per transformer.

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
            name = 'xfmr3p';
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
            nt = obj.nk;            

            % Change of basis
            base_kv = bus_dme.tab.base_kv(bus_dme.on(dme.fbus(dme.on)));
            z = dm.base_kva * (dme.r + 1j * dme.x) ./ dme.base_kva .* ...
                (dme.base_kv ./ ( base_kv / sqrt(3))) .^ 2;
            y = 1 ./ z;
            
            y_series = repmat(y.', 3, 1);  %% ys = [ya; yb; yc] (3 x nt)
            tm = dme.tab.tm;    %% tap ratios                

            %% Compute the diagonal used both for parameters Wu and M
            
            phases = {'a', 'b', 'c'};

            vals_diag_from = [];
            ii_diag_from = [];
            jj_diag_from = [];
            for p = 1:length(phases)
                [ii_aux, jj_aux, vals_aux] = ...
                    obj.make_diag(dme.fbus, dme.tbus, nb3p, y_series, tm', phases{p}, 'from');
                vals_diag_from = horzcat(vals_diag_from, vals_aux);
                ii_diag_from = horzcat(ii_diag_from, ii_aux);
                jj_diag_from = horzcat(jj_diag_from, jj_aux);
            end

            vals_diag_to = [];
            ii_diag_to = [];
            jj_diag_to = [];
            for p = 1:length(phases)
                [ii_aux, jj_aux, vals_aux] = ...
                    obj.make_diag(dme.fbus, dme.tbus, nb3p, y_series, tm', phases{p}, 'to');
                vals_diag_to = horzcat(vals_diag_to, vals_aux);
                ii_diag_to = horzcat(ii_diag_to, ii_aux);
                jj_diag_to = horzcat(jj_diag_to, jj_aux);
            end

            %% 1) Wu parameter
            if nm.userdata.tpc.quad
                nval = 16;   %% number of nonzero values per matrix

                %% 1.1) Qu matrices - 'from' ports 
                
                % Compute the elements in the upper triangular              
                [ii_triu, jj_triu, vals_triu] = obj.make_triu(dme.fbus, dme.tbus, y_series, nb3p, 'from');
    
                % Compute the elements on the lower triangular (symmetric matrices)
                vals_tril = vals_triu;
                ii_tril = jj_triu;
                jj_tril = ii_triu;
    
                % Compute final triplets (in the form of row-col-val)
                vals_aux = vertcat(vals_diag_from, vals_triu, vals_tril) ./ ...
                           repmat(repmat( tm', 1, 3) , nval , 1);                
                ii_aux = vertcat(ii_diag_from, ii_triu, ii_tril);
                jj_aux = vertcat(jj_diag_from, jj_triu, jj_tril);
                
                vals_from = mat2cell(vals_aux(:), nval*ones(3*nt,1));                
                ii_from = mat2cell(ii_aux(:), nval*ones(3*nt,1));
                jj_from = mat2cell(jj_aux(:), nval*ones(3*nt,1));

                %% 1.2) Qu matrices - 'to' ports
                
                % Compute the elements in the upper triangular
                [ii_triu, jj_triu, vals_triu] = obj.make_triu(dme.fbus, dme.tbus, y_series, nb3p, 'to');
    
                % Compute the elements on the lower triangular (symmetric matrices)
                vals_tril = vals_triu;
                ii_tril = jj_triu;
                jj_tril = ii_triu;
    
                % Compute final triplets (in the form of row-col-val)
                vals_aux = vertcat(vals_diag_to, vals_triu, vals_tril) ./ ...
                           repmat(repmat( tm', 1, 3) , nval , 1);                
                ii_aux = vertcat(ii_diag_to, ii_triu, ii_tril);
                jj_aux = vertcat(jj_diag_to, jj_triu, jj_tril);
                
                vals_to = mat2cell(vals_aux(:), nval*ones(3*nt,1));
                ii_to = mat2cell(ii_aux(:), nval*ones(3*nt,1));
                jj_to = mat2cell(jj_aux(:), nval*ones(3*nt,1));
                
                %% 1.3) Build Qu parameter
                vals_Qu = [vals_from; vals_to];
                ii_Wu = [ii_from; ii_to];
                jj_Wu = [jj_from; jj_to];

                obj.Qu = mat2cell(cell2mat([ii_Wu jj_Wu vals_Qu]), nval*ones(6*nt,1));
            end
            
            %% 2) M parameter            
            obj.M = obj.make_M(dme.fbus, dme.tbus, y_series, tm, nb3p);
            
            %% 3) Build s parameter
            s_blk1 = repmat((1-tm)./(tm.^2), 1, 3) .* y_series';
            s_blk2 = repmat((tm-1)./tm, 1, 3) .* y_series';
            obj.s = [s_blk1(:); s_blk2(:)];
        end
    
        function [ii, jj, vals] = make_triu(obj, fbus, tbus, ys, nb3p, side)
            
            nl = size(ys, 2);            
    
            % 1) Compute values in the upper triangular            
            Ypp = cell2mat(mat2cell(ys,ones(3,1))');
            if strcmp(side, 'from') 
                factors = repmat([-1 1j 1j -1j -1j -1].', 1, nl);
            elseif strcmp(side, 'to')
                factors = repmat([-1 -1j -1j 1j 1j -1].', 1, nl);
            else
                error('wgv.nme_xfmr3p_tpc:make_triu: side %s not recognized',side);
            end
            vals = conj(repmat(factors, 1, 3) .* repmat(Ypp, 6, 1));

            % 2) Compute row and column indices on the upper triangular
            % Block 1 of values (first row of three nonzero values)
            id_bias = cell2mat(mat2cell(repmat([0 1 2]'.*nb3p, 1, nl), ones(3,1))');
            ii_blk1 = repmat(repmat(fbus', 3, 1), 1, 3) + repmat(id_bias, 3, 1);
            jj_blk1 = repmat([tbus fbus+sum(nb3p) tbus+sum(nb3p)]', 1, 3) + repmat(id_bias, 3, 1);

            % Block 2 of values (second row of two nonzero values)            
            ii_blk2 = repmat(repmat(tbus', 2, 1), 1, 3) + repmat(id_bias, 2, 1);
            jj_blk2 = repmat([fbus+sum(nb3p) tbus+sum(nb3p)]', 1, 3) + repmat(id_bias, 2, 1);

            % Block 3 of values (third row of one nonzero value)            
            ii_blk3 = repmat(fbus'+sum(nb3p), 1, 3) + id_bias;
            jj_blk3 = repmat(tbus'+sum(nb3p), 1, 3) + id_bias;
            
            % Compute row and column indices for from/to end
            ii = [ii_blk1; ii_blk2; ii_blk3];
            jj = [jj_blk1; jj_blk2; jj_blk3];
        end

        function [ii, jj, vals] = make_diag(obj, fbus, tbus, nb3p, ys, tm, phase, side)
            %
            % Inputs:
            %   ys  (matrix)  : 6 x nl matrix with the elements of the
            %                   series admitance matrix. Each column holds
            %                   the information of a 3-phase line
            %    tm (vector)  : nl x 1 vector with the turns ratios of all
            %                   transformers
            %  phase (string) : wheter 'a', 'b' or 'c' identifier for the
            %                   phase to be considered when building the
            %                   diagonal elements of the quadratic matrices
            %   side (string) : wheter 'from' or 'to' identifier for the
            %                   type of side being considered
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
                        
            nl = size(ys, 2);
            vals = zeros(4, nl);        
    
            % 1) Compute values for quadratic matrices on the diagonal
        
            Ypp = conj(ys(p,:));     %% Values on the diagonal
            if strcmp(side, 'from')
                % theta_from - theta_from
                vals(1,:) = Ypp;
    
                % theta_to - theta_to
                vals(2,:) = Ypp;
    
                % lnvm_from - lnvm_from
                vals(3,:) = 4*Ypp./tm - Ypp;
    
                % lnvm_to - lnvm_to
                vals(4,:) = -1 * Ypp;
            elseif strcmp(side, 'to')
                % theta_from - theta_from
                vals(1,:) = Ypp;
    
                % theta_to - theta_to
                vals(2,:) = Ypp;
    
                % lnvm_from - lnvm_from
                vals(3,:) = -1 * Ypp;
    
                % lnvm_to - lnvm_to
                vals(4,:) = 4*Ypp.*tm - Ypp;
            else
                error('wgv.nme_xfmr3p_tpc:make_diag: side %s not recognized',side);
            end            
    
            % 2) Compute row and column indices on the diagonal
            % 
            % The vector of variables is now: [\theta_1ph; lnvm_1ph; \theta_3ph; lnvm_3ph]
            
            % theta_from - theta_from            
            id_bias = (p-1)*nb3p(p);
            ii_tff = fbus' + id_bias;

            % theta_to - theta_to            
            ii_ttt = tbus' + id_bias;

            % lnvm_from - lnvm_from
            id_bias = (p-1)*nb3p(p)+sum(nb3p);
            ii_lnvmff = fbus' + id_bias;

            % lnvm_to - lnvm_to
            ii_lnvmtt = tbus' + id_bias;
            
            % Compute row and column indices for one single phase
            ii = [ii_tff; ii_ttt; ii_lnvmff; ii_lnvmtt];
            jj = ii;
        end

        function C = make_M(obj, fbus, tbus, ys, tm, nb3p)
            %

            nt = obj.nk;
            
            %% 1) Compute values for C_\theta
            %              Diagonal                    Off-diagonal
            C_th = [1*repmat(ys.', 1, 2)    -1*repmat(ys.', 1, 2)];            

            %% 2) Compute values for C_\nu            
            C_nu_diag1 = 2*(ys.' ./ repmat(tm, 1, 3)) -1*ys.';
            C_nu_diag2 = 2*(ys.' .* repmat(tm, 1, 3)) -1*ys.';
            C_nu_offd  = -1*repmat(ys.', 1, 2);
            C_nu = [C_nu_diag1 C_nu_diag2 C_nu_offd];

            %% 3) Compute values for C
            C_vals = [-1j*conj(C_th) conj(C_nu)] ./ repmat(tm, 1, 24);

            %% 4) Compute indices
            % Indices for C_\theta
            ii_diag = reshape((1:6*nt)', nt, []);
            id_bias = repmat(repmat([0 1 2].*nb3p', nt, 1), 1, 2);
            jj1 = [repmat(fbus, 1, 3)   repmat(tbus, 1, 3)] + id_bias;
            jj2 = [repmat(tbus, 1, 3)   repmat(fbus, 1, 3)] + id_bias;

            ii_th = [ii_diag ii_diag];
            jj_th = [jj1 jj2];

            % Indices for C_\nu
            ii_nu = ii_th;
            jj_nu = jj_th + sum(nb3p);

            % Whole matrices of indices
            ii = [ii_th ii_nu];
            jj = [jj_th jj_nu];

            C = sparse(ii(:), jj(:), C_vals(:), nt*6, 2*sum(nb3p));
        end
    end     %% methods
end         %% classdef
