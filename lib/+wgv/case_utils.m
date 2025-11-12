classdef case_utils < mp.case_utils


    methods (Static)
        function mpc_out = unbalance_case(mpc_in,unbalancing,params)
            if nargin < 3
                params = struct('max_ldpf', 1, ...
                                'min_ldpf', 0.5, ...
                                'max_kW',   1.25, ...
                                'min_kW',   0.75, ...
                                'rnd_ldpf', 1);
            end

            mpc_out = mpc_in;
            if any(ismember(upper(unbalancing),{'LOAD'}))
                nl = size(mpc_in.load3p,1);
                id_phase_ldpf = randi(3,1,nl);
                id_phase_kW = randi(3,1,nl);
                old_ldpf = mpc_out.load3p(:,7:9);
                old_kW = mpc_out.load3p(:,4:6);
                val_max_ldpf = params.max_ldpf;
                val_min_ldpf = params.min_ldpf;
                new_ldpf = ((val_max_ldpf - val_min_ldpf)*rand(nl,1) + val_min_ldpf);
                if params.rnd_ldpf
                    new_ldpf = params.rnd_ldpf.*sign(randn(nl,1));
                end
                val_max_kW = params.max_kW; % percentage
                val_min_kW = params.min_kW; % percentage
                new_kW = (val_max_kW - val_min_kW)*rand(1,nl) + val_min_kW;

                old_ldpf(sub2ind(size(old_ldpf), 1:nl, id_phase_ldpf)) = new_ldpf;
                old_kW(sub2ind(size(old_kW), 1:nl, id_phase_kW)) = new_kW .* old_kW(sub2ind(size(old_kW), 1:nl, id_phase_kW));
                mpc_out.load3p(:,4:6) = old_kW;
                mpc_out.load3p(:,7:9) = old_ldpf;
            end
            if any(ismember(upper(unbalancing),{'LINE_PARAMS'}))
                nl = size(mpc_in.line3p,1);
                % Diagonal values
                min_factor_diagR = 0.98;
                max_factor_diagR = 1.1;
                min_factor_diagX = 1.18;
                max_factor_diagX = 1.3;
                min_factor_diagC = 0.14;
                max_factor_diagC = 0.19;
                vals_diagR = (max_factor_diagR - min_factor_diagR)*rand(nl,3) + min_factor_diagR;
                vals_diagX = (max_factor_diagX - min_factor_diagX)*rand(nl,3) + min_factor_diagX;
                vals_diagC = (max_factor_diagC - min_factor_diagC)*rand(nl,3) + min_factor_diagC;

                % Off-diagonal values
                %  Resistance
                min_factor_Rab = 0.89;
                max_factor_Rab = 1.2;
                min_factor_Rac = 0.86;
                max_factor_Rac = 1.13;
                min_factor_Rbc = 0.92;
                max_factor_Rbc = 1.27;
                vals_Rab = (max_factor_Rab - min_factor_Rab)*rand(nl,1)+min_factor_Rab;
                vals_Rac = (max_factor_Rac - min_factor_Rac)*rand(nl,1)+min_factor_Rac;
                vals_Rbc = (max_factor_Rbc - min_factor_Rbc)*rand(nl,1)+min_factor_Rbc;

                %  Reactance
                min_factor_Xab = 0.85;
                max_factor_Xab = 1.1;
                min_factor_Xac = 0.8;
                max_factor_Xac = 1.0;
                min_factor_Xbc = 0.9;
                max_factor_Xbc= 1.25;
                vals_Xab = (max_factor_Xab - min_factor_Xab)*rand(nl,1)+min_factor_Xab;
                vals_Xac = (max_factor_Xac - min_factor_Xac)*rand(nl,1)+min_factor_Xac;
                vals_Xbc = (max_factor_Xbc - min_factor_Xbc)*rand(nl,1)+min_factor_Xbc;

                %  Capacitance
                min_factor_Cab = 0.10;
                max_factor_Cab = 0.20;
                min_factor_Cac = 0.05;
                max_factor_Cac = 0.13;
                min_factor_Cbc = 0.11;
                max_factor_Cbc = 0.4;
                vals_Cab = (max_factor_Cab - min_factor_Cab)*rand(nl,1)+min_factor_Cab;
                vals_Cac = (max_factor_Cac - min_factor_Cac)*rand(nl,1)+min_factor_Cac;
                vals_Cbc = (max_factor_Cbc - min_factor_Cbc)*rand(nl,1)+min_factor_Cbc;

                % Final values
                valsR = [vals_diagR(:,1) vals_Rab vals_Rac vals_diagR(:,2) vals_Rbc vals_diagR(:,3)];
                valsX = [vals_diagX(:,1) vals_Xab vals_Xac vals_diagX(:,2) vals_Xbc vals_diagX(:,3)];
                valsC = [vals_diagC(:,1) vals_Cab vals_Cac vals_diagC(:,2) vals_Cbc vals_diagC(:,3)];

                mpc_out.lc(:,2:end) = mpc_out.lc(:,2:end) .* [valsR valsX valsC];
            end
        end
        
        function mpc = init_conf_eq_mpc(baseMVA, basekV)
            % Creates a small single-phase MATPOWER case with one PV bus and one genetaror
            
            if nargin < 2
                basekV = 230;
            end
            if nargin < 2
                baseMVA = 10;
            end

            mpc = loadcase('case4gs');   % Use the smallest available test case

            % --- 1) base power
            mpc.baseMVA = baseMVA;
            
            % --- 2) bus matrix
            mpc.bus(2:end,:) = [];            
            mpc.bus(2) = 2;           % Set the bus type as PV
            mpc.bus([3 4 5]) = 0;     % Zero out demand, susceptance, and conductance
            mpc.bus(10) = basekV;

            % --- 3) gen matrix
            mpc.gen(1,:) = [];        % This ensures generator is connected at single bus

            % --- 4) branch matrix
            mpc.branch = [];          % No branches in this case
        end

        function new_mpc = mpc_base_change(mpc,new_baseMVA,new_basekV)
            %

            % Ensure mpc is a struct
            mpc = loadcase(mpc);
            
            % define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN] = idx_bus;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX] = idx_brch;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            % Check inputs
            if nargin < 2
                error('wgv.case_utils.mpc_base_change: method require at least two inputs')
            end
            if nargin < 3
                new_basekV = mpc.bus(:,BASE_KV);
            else
                if numel(new_basekV) ~= size(mpc.bus,1) 
                    error('wgv.case_utils.mpc_base_change: new base voltages must be a vector of the same size as the number of system buses')
                end
                if any(new_basekV < 0)
                    error('wgv.case_utils.mpc_base_change: new base voltages must be positive')
                end
                if new_baseMVA < 0 || ~isscalar(new_baseMVA)
                    error('wgv.case_utils.mpc_base_change: new base power must be positive scalar in MVA')
                end
            end
            new_mpc = mpc;

            % Extract old bases and set new bases
            old_baseMVA = mpc.baseMVA;
            old_basekV = mpc.bus(:,BASE_KV);
            new_mpc.baseMVA = new_baseMVA;
            new_mpc.bus(:,BASE_KV) = new_basekV;

            % ---- 1) Change of bases in bus matrix
            zGs = mp.case_utils.z_base_change(1./mpc.bus(:,GS), ...
                        old_basekV,old_baseMVA,new_basekV,new_baseMVA);
            zBs = mp.case_utils.z_base_change(1./mpc.bus(:,BS), ...
                        old_basekV,old_baseMVA,new_basekV,new_baseMVA);
            new_mpc.bus(:,GS) = 1./zGs;
            new_mpc.bus(:,BS) = 1./zBs;
            new_mpc.bus(:,VM) = mpc.bus(:,VM) .* (old_basekV./new_basekV);
            new_mpc.bus(:,VMAX) = mpc.bus(:,VMAX) .* (old_basekV./new_basekV);
            new_mpc.bus(:,VMIN) = mpc.bus(:,VMIN) .* (old_basekV./new_basekV);

            % ---- 2) Change of bases in gen matrix
            [~,id_bus_gen] = ismember(mpc.gen(:,GEN_BUS),mpc.bus(:,BUS_I));
            new_mpc.gen(:,VG) = mpc.gen(:,VG) .* (old_basekV(id_bus_gen)./new_basekV(id_bus_gen));

            % ---- 3) Change of bases in branch matrix (we take voltage at the from sied of the branch)
            [~,id_branch_fbus] = ismember(mpc.branch(:,F_BUS),mpc.bus(:,BUS_I));
            new_mpc.branch(:,BR_R) = mp.case_utils.z_base_change(mpc.branch(:,BR_R), ...
                    old_basekV(id_branch_fbus),old_baseMVA, ...
                    new_basekV(id_branch_fbus),new_baseMVA);
            new_mpc.branch(:,BR_X) = mp.case_utils.z_base_change(mpc.branch(:,BR_X), ...
                    old_basekV(id_branch_fbus),old_baseMVA, ...
                    new_basekV(id_branch_fbus),new_baseMVA);
            zBR_B = mp.case_utils.z_base_change(1./mpc.branch(:,BR_B), ...
                    old_basekV(id_branch_fbus),old_baseMVA, ...
                    new_basekV(id_branch_fbus),new_baseMVA);
            new_mpc.branch(:,BR_B) = 1./zBR_B;
        end

        function [mpc_out, connect_info] = connect_cases(mpcT,mpcD,idbusT,idbusD,ishybrid,branch_links)
            %% Initial checks
            mpcT = loadcase(mpcT);  % ensures mpcT is a struct
            mpcD = loadcase(mpcD);  % ensures mpcD is a struct

            % define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN] = idx_bus;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX] = idx_brch;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            % check meaningful indices
            if ~all(ismember(idbusT,mpcT.bus(:,BUS_I)))
                error('connect_cases: all indices provided in input ''idbusT'' must belong to the set of bus indices of mpc1');
            end
            if ~all(ismember(idbusD,mpcD.bus(:,BUS_I)))
                error('connect_cases: all indices provided in input ''idbusD'' must belong to the set of bus indices of mpc2');
            end
            if length(idbusT) ~= length(idbusD)
                error('connect_cases: connection indices dimension mismatch')
            end

            if nargin < 6
                branch_links = [];
            end
            if nargin < 5
                ishybrid = 0;
            end
            if ~isempty(branch_links)
                if size(branch_links,1) ~= numel(idbusT) || size(branch_links,2) ~= size(mpcT.branch,2)
                    error('connect_cases: inconsistent information of branch links')
                end
            end

            % check for the presence of distribution network equivalents in mpcD
            if size(mpcD.bus,1) == numel(idbusD) && size(mpcD.gen,1) == numel(idbusD)
                eqDflag = 1;
            else
                eqDflag = 0;
            end

            % check for the presence of generators in the connection buses
            % of mpcD when this network is not an equivalent
            if any(ismember(idbusD,mpcD.gen(:,GEN_BUS))) && ~eqDflag
                error('connect_cases: generators connected in the connection buses of mpcD are not allowed')
            end

            %% Connect cases based on provided inputs
            nc = length(idbusT);                %% number of connections between cases
            max_id1 = max(mpcT.bus(:,BUS_I));   %% max id among buses of mpcT
            if ishybrid && ~eqDflag
                if isempty(branch_links)
                    % check consistency of base voltages
                    id_busT = ismember(mpcT.bus(:,BUS_I),idbusT);
                    lgcl_busD = ismember(mpcD.bus(:,BUS_I),idbusD);
                    base_kV1 = mpcT.bus(id_busT,BASE_KV);
                    base_kV2 = mpcD.bus(lgcl_busD,BASE_KV);
                    if ~isequal(base_kV1,base_kV2)
                        error('connect_cases: buses being connected must have the same BASE_KV values')
                    end

                    % create hybrid case, mpcT is single-phase and mpcD is three-phase
                    mpc_out = mp.case_utils.convert_1p_to_3p(mpcD);
                    mpc_out.buslink = [(1:nc)' idbusT idbusD ones(nc,1)];
                    mpc_out.bus = mpcT.bus;
                    mpc_out.gen = mpcT.gen;
                    mpc_out.branch = mpcT.branch;
                    mpc_out.gencost = mpcT.gencost;
                    mpc_out.baseMVA = mpcT.baseMVA;
                else

                end
                if nargout > 1
                    lgcl_branch_fbus2 = ismember(mpcD.branch(:,F_BUS),idbusD);
                    lgcl_branch_tbus2 = ismember(mpcD.branch(:,T_BUS),idbusD);
                    connect_info.interface.fname = cell(nc,1);
                    connect_info.interface.idx = cell(nc,1);
                    idx_branch_mpc2 = (1:size(mpcD.branch,1))';
                    idx_branch_mpc2 = idx_branch_mpc2(or(lgcl_branch_fbus2,lgcl_branch_tbus2));
                    for i = 1:numel(idx_branch_mpc2)
                        id_br = idx_branch_mpc2(i);
                        if mpcD.branch(id_br,TAP) ~= 0
                            fname = 'xfmr3p';
                        else
                            fname = 'line3p';
                        end
                        idx = find(mpc_out.(fname)(:,2)==mpcD.branch(id_br,F_BUS));
                        connect_info.interface.fname{i} = fname;
                        connect_info.interface.idx{i} = idx;
                    end
                end
            elseif ~eqDflag
                if isempty(branch_links)
                    % check consistency of base voltages
                    id_busT = ismember(mpcT.bus(:,BUS_I),idbusT);
                    lgcl_busD = ismember(mpcD.bus(:,BUS_I),idbusD);
                    base_kV1 = mpcT.bus(id_busT,BASE_KV);
                    base_kV2 = mpcD.bus(lgcl_busD,BASE_KV);
                    if ~isequal(base_kV1,base_kV2)
                        error('connect_cases: buses being connected must have the same BASE_KV values')
                    end

                    % merge bus connection info of mpcD into mpcT
                    % 1. bus matrix
                    mpcT.bus(id_busT,PD) = mpcT.bus(id_busT,PD) + mpcD.bus(lgcl_busD,PD);  % active power
                    mpcT.bus(id_busT,QD) = mpcT.bus(id_busT,QD) + mpcD.bus(lgcl_busD,QD);  % reactive power
                    zGs = mp.case_utils.z_base_change(1./mpcD.bus(lgcl_busD,GS), ...
                        base_kV2,mpcD.baseMVA,base_kV2,mpcT.baseMVA);
                    zBs = mp.case_utils.z_base_change(1./mpcD.bus(lgcl_busD,BS), ...
                        base_kV2,mpcD.baseMVA,base_kV2,mpcT.baseMVA);
                    mpcT.bus(id_busT,GS) = mpcT.bus(id_busT,GS) + 1./zGs;   % nodal conductances
                    mpcT.bus(id_busT,BS) = mpcT.bus(id_busT,BS) + 1./zBs;   % nodal susceptances
                    % 2. branch matrix
                    [~,id_fbus_2_all] = ismember(mpcD.branch(:,F_BUS),mpcD.bus(:,BUS_I));
                    base_kV2_all = mpcD.bus(id_fbus_2_all,BASE_KV);
                    mpcD.branch(:,BR_R) = mp.case_utils.z_base_change(mpcD.branch(:,BR_R), ...
                        base_kV2_all,mpcD.baseMVA,base_kV2_all,mpcT.baseMVA);
                    mpcD.branch(:,BR_X) = mp.case_utils.z_base_change(mpcD.branch(:,BR_X), ...
                        base_kV2_all,mpcD.baseMVA,base_kV2_all,mpcT.baseMVA);
                    zBR_B = mp.case_utils.z_base_change(1./mpcD.branch(:,BR_B), ...
                        base_kV2_all,mpcD.baseMVA,base_kV2_all,mpcT.baseMVA);
                    mpcD.branch(:,BR_B) = 1./zBR_B;
                    [~,id_fbus2] = ismember(idbusD,mpcD.branch(:,F_BUS));
                    [~,id_tbus2] = ismember(idbusD,mpcD.branch(:,T_BUS));
                    id_fbus2 = id_fbus2(id_fbus2~=0);
                    id_tbus2 = id_tbus2(id_tbus2~=0);
                    mpcD.branch(id_fbus2,F_BUS) = idbusT;          % update sending nodes in mpcD
                    mpcD.branch(id_tbus2,T_BUS) = idbusT;          % update receiving nodes in mpcD

                    mpcD.bus(lgcl_busD,:) = [];                         % after merging info, remove buses of mpcD
                    mpcD.bus(:,BUS_I) = mpcD.bus(:,BUS_I) + max_id1;    % update bus numbering of mpcD

                    id_branch_fbus2 = setdiff(1:size(mpcD.branch,1),id_fbus2);
                    id_branch_tbus2 = setdiff(1:size(mpcD.branch,1),id_tbus2);
                    mpcD.branch(id_branch_fbus2,F_BUS) = mpcD.branch(id_branch_fbus2,F_BUS) + max_id1;
                    mpcD.branch(id_branch_tbus2,T_BUS) = mpcD.branch(id_branch_tbus2,T_BUS) + max_id1;

                    if nargout > 1
                        nb1 = size(mpcT.bus,1);
                        nb2 = size(mpcD.bus,1);
                        nbr1 = size(mpcT.branch,1);
                        nbr2 = size(mpcD.branch,1);
                        connect_info.bus.mpc1.idx = [1 nb1];
                        connect_info.bus.mpc2.idx = [1 nb2]+nb1;
                        connect_info.branch.mpc1.fname = 'branch';
                        connect_info.branch.mpc1.idx = [1 nbr1];
                        connect_info.branch.mpc2.fname = 'branch';
                        connect_info.branch.mpc2.idx = [1 nbr2]+nbr1;
                        connect_info.interface.fname = 'branch';
                        connect_info.interface.idx = union(id_fbus2,id_tbus2)+nbr1;
                    end
                else

                end
                % combine cases into one single-phase case
                mpc_out = mpcT;
                mpc_out.bus = [mpc_out.bus; mpcD.bus];
                mpc_out.branch = [mpc_out.branch(:,F_BUS:ANGMAX); mpcD.branch];
            else % Distribution networks are equivalents
                mpc_out = mpcT;
                [~,id_busT] = ismember(idbusT,mpcT.bus(:,BUS_I));
                mpc_out.bus(id_busT,BS) = mpc_out.bus(id_busT,BS) + mpcD.bus(idbusD,BS);
                mpc_out.bus(id_busT,BUS_TYPE) = 2;
                mpcD.gen(:,GEN_BUS) = idbusT;
                mpc_out.gen = [mpc_out.gen; mpcD.gen];
                mpc_out.gencost = [mpc_out.gencost; repmat([2 zeros(1,6)],numel(idbusD),1)];
            end
        end

        function connect_info_3p = from_1p_to_3p_connect_info(connect_info_1p,case1p,case3p)
            idx_all_br = connect_info_1p.interface.idx;
            is_xfmr = case1p.branch(idx_all_br,9) ~= 0;
            idx_xfmr = idx_all_br(is_xfmr);
            idx_line = idx_all_br(~is_xfmr);

            connect_info_3p.interface.fname = cell(numel(idx_all_br),1);
            connect_info_3p.interface.idx = cell(numel(idx_all_br),1);
            for t = 1:length(idx_xfmr)
                idx = idx_xfmr(t);
                connect_info_3p.interface.fname{t} = 'xfmr3p';
                idx_fbus = find(case3p.xfmr3p(:,2) == case1p.branch(idx,1));
                connect_info_3p.interface.idx{t} = idx_fbus;
            end
            for l = 1:length(idx_line)
                idx = idx_line(l);
                connect_info_3p.interface.fname{l} = 'line3p';
                idx_fbus = find(case3p.line3p(:,2) == case1p.branch(idx,1));
                connect_info_3p.interface.idx{l} = idx_fbus;
            end
        end

        function [mmpc, merge_info] = merge_cases(cell_of_cases,tags)
            
            %
            
            define_constants;
            if nargout > 1
                merge_info = struct;
            end

            % Check inputs
            is_str = cellfun(@ischar,cell_of_cases);
            is_struct = cellfun(@isstruct,cell_of_cases);
            if ~all(or(is_str,is_struct))
                error('wgv.case_utils.merge_cases: input must be a cell array of valid case names and/or structs')
            end

            % Build by-default consecutive tags of the form 'mpc_#'
            n_mpc = numel(cell_of_cases);
            case_numbers = cellfun(@strtrim,cellstr(num2str((1:n_mpc)')),'UniformOutput',false);
            
            % Assign tags to all cases
            if nargin < 2                
                if all(is_struct)
                    tags = cellfun(@(x)(['mpc_' x]),case_numbers,'UniformOutput',false);
                else
                    tags = cell(n_mpc,1);
                    tags(is_str) = cell_of_cases(is_str);
                    tags(is_struct) = cellfun(@(x)(['mpc_' x]),case_numbers(is_struct),'UniformOutput',false);
                end                
            else
                if ~all(cellfun(@ischar,tags))
                    error('wgv.case_utils.mege_cases: input ''tags'' must be a cell array of char values')
                end
            end
            
            % Ensure each case is loaded as a struct
            cell_of_cases = cellfun(@(x)(loadcase(x)),cell_of_cases,'UniformOutput',false);
            
            % Find the number of buses and the maximmum bus ID per case            
            max_bus_id = cellfun(@(x)(max(x.bus(:,BUS_I))),cell_of_cases);
            [max_bus_id,sorted_max_bus_id] = sort(max_bus_id);
            nbus = cellfun(@(x)(size(x.bus,1)),cell_of_cases(sorted_max_bus_id));
            nbrch = cellfun(@(x)(size(x.branch,1)),cell_of_cases(sorted_max_bus_id));
            ngen = cellfun(@(x)(size(x.gen,1)),cell_of_cases(sorted_max_bus_id));
            
            % Move indices of each case according to max bus ID and stack mpc matrices
            mmpc = cell_of_cases{sorted_max_bus_id(1)};
            offset_id = cumsum(max_bus_id);
            cum_nbus = cumsum(nbus);
            cum_nbrch = cumsum(nbrch);
            cum_ngen = cumsum(ngen);
            for i = 2:n_mpc
                % Shift bus IDs
                mpc = cell_of_cases{sorted_max_bus_id(i)};
                mpc.bus(:,BUS_I) = mpc.bus(:,BUS_I) + offset_id(i-1);
                if ~isempty(mpc.branch)
                    mpc.branch(:,F_BUS) = mpc.branch(:,F_BUS) + offset_id(i-1);
                    mpc.branch(:,T_BUS) = mpc.branch(:,T_BUS) + offset_id(i-1);
                end                
                mpc.gen(:,GEN_BUS) = mpc.gen(:,GEN_BUS) + + offset_id(i-1);                

                % Stack info                
                mmpc.bus = [mmpc.bus; mpc.bus];             % bus matrix
                mmpc.branch = [mmpc.branch; mpc.branch];    % branch matrix
                mmpc.gen = [mmpc.gen; mpc.gen];             % gen matrix
                % mmpc.gencost = pending ...                % gencost matrix

                % Build merge info if requested
                if nargout > 1
                    fname = tags{sorted_max_bus_id(i)};
                    merge_info.(fname).busoffset = offset_id(i-1);
                    merge_info.(fname).bus.id1 = cum_nbus(i-1) + 1;
                    merge_info.(fname).bus.idN = cum_nbus(i);
                    merge_info.(fname).branch.id1 = cum_nbrch(i-1) + 1;
                    merge_info.(fname).branch.idN = cum_nbrch(i);
                    merge_info.(fname).gen.id1 = cum_ngen(i-1) + 1;
                    merge_info.(fname).gen.idN = cum_ngen(i);
                end
            end

            % Fill in merge info for first mpc
            if nargout > 1
                fname = tags{sorted_max_bus_id(1)};
                merge_info.(fname).busoffset = 0;
                merge_info.(fname).bus.id1 = 1;
                merge_info.(fname).bus.idN = cum_nbus(1);
                merge_info.(fname).branch.id1 = 1;
                merge_info.(fname).branch.idN = cum_nbrch(1);
                merge_info.(fname).gen.id1 = 1;
                merge_info.(fname).gen.idN = cum_ngen(1);
            end

            % Remove bus_name field and (temporally) gencost
            if isfield(mmpc,'bus_name')
                mmpc = rmfield(mmpc,'bus_name');
            end
            if isfield(mmpc,'gencost')
                mmpc = rmfield(mmpc,'gencost');
            end

            % Return the order in which input cases where merged
            merge_info.order.idx = sorted_max_bus_id;
            merge_info.order.str = tags(sorted_max_bus_id);

        end

        function merge_info_3p = from_1p_to_3p_merge_info(merge_info_1p,case1p,case3p)
            %
            
            define_constants;

            cases_D = fieldnames(merge_info_1p);
            id_order = find(strcmp('order',cases_D));
            cases_D(id_order) = [];

            merge_info_3p = struct();

            for c = 1:numel(cases_D)
                cname = cases_D{c};                
                    
                % ---- 1) buses
                id1 = merge_info_1p.(cname).bus.id1;
                idN = merge_info_1p.(cname).bus.idN;
                id_bus1p = case1p.bus(id1:idN,BUS_I);
                id_bus3p = ismember(case3p.bus3p(:,1),id_bus1p);

                % ---- 2) lines
                id1 = merge_info_1p.(cname).branch.id1;
                idN = merge_info_1p.(cname).branch.idN;                
                id_branch1p = (id1:idN)';
                id_line1p = case1p.branch(id_branch1p,TAP) == 0;
                id_line1p = id_branch1p(id_line1p);
                fbus_line1p = case1p.branch(id_line1p,F_BUS);
                tbus_line1p = case1p.branch(id_line1p,T_BUS);
                ft_str_1p = cellfun(@(x,y)([x '-' y]), ...
                    cellfun(@strtrim, cellstr(num2str(fbus_line1p)), 'UniformOutput', false), ...
                    cellfun(@strtrim, cellstr(num2str(tbus_line1p)), 'UniformOutput', false), ...
                    'UniformOutput',false);
                fbus_line3p = case3p.line3p(:,2);
                tbus_line3p = case3p.line3p(:,3);
                ft_str_3p = cellfun(@(x,y)([x '-' y]), ...
                    cellfun(@strtrim, cellstr(num2str(fbus_line3p)), 'UniformOutput', false), ...
                    cellfun(@strtrim, cellstr(num2str(tbus_line3p)), 'UniformOutput', false), ...
                    'UniformOutput',false);
                id_line3p = ismember(ft_str_3p,ft_str_1p);
                
                % ---- 3) transformers
                id_trafo1p = case1p.branch(id_branch1p,TAP) ~= 0;
                id_trafo1p = id_branch1p(id_trafo1p);                
                fbus_trafo1p = case1p.branch(id_trafo1p,F_BUS);
                tbus_trafo1p = case1p.branch(id_trafo1p,T_BUS);
                ft_str_1p = cellfun(@(x,y)([x '-' y]), ...
                    cellfun(@strtrim, cellstr(num2str(fbus_trafo1p)), 'UniformOutput', false), ...
                    cellfun(@strtrim, cellstr(num2str(tbus_trafo1p)), 'UniformOutput', false), ...
                    'UniformOutput',false);
                fbus_trafo3p = case3p.xfmr3p(:,2);
                tbus_trafo3p = case3p.xfmr3p(:,3);
                ft_str_3p = cellfun(@(x,y)([x '-' y]), ...
                    cellfun(@strtrim, cellstr(num2str(fbus_trafo3p)), 'UniformOutput', false), ...
                    cellfun(@strtrim, cellstr(num2str(tbus_trafo3p)), 'UniformOutput', false), ...
                    'UniformOutput',false);
                id_xfmr3p = ismember(ft_str_3p,ft_str_1p);

                % ---- 4) load
                id_load3p = ismember(case3p.load3p(:,2),id_bus1p);

                % store ouput info
                merge_info_3p.(cname).id_bus3p = id_bus3p;
                merge_info_3p.(cname).id_line3p = id_line3p;
                merge_info_3p.(cname).id_xfmr3p = id_xfmr3p;
                merge_info_3p.(cname).id_load3p = id_load3p;
            end
        end

        function plot_TD_sld(caseT,caseD,busesT,busesD,is_tiled)
            %
            if nargin < 5
                is_tiled = 0;
            end
            
            %% Create graph plots
            define_constants;
            
            mpc_T = loadcase(caseT);
            mpc_D = loadcase(caseD);
            
            f_T = wgv.case_utils.vec2char(mpc_T.branch(:,F_BUS));
            t_T = wgv.case_utils.vec2char(mpc_T.branch(:,T_BUS));
            f_D = wgv.case_utils.vec2char(mpc_D.branch(:,F_BUS));
            t_D = wgv.case_utils.vec2char(mpc_D.branch(:,T_BUS));
            busesT_cell_str = wgv.case_utils.vec2char(busesT);
            busesD_cell_str = wgv.case_utils.vec2char(busesD);

            graph_T = graph(f_T,t_T);
            graph_D = graph(f_D,t_D);

            if is_tiled
                figure;
                tl = tiledlayout(1,2);
                nexttile

                hT = plot(graph_T,'k');
                hT.Interpreter = "latex";
                
                %layout(hT,'layered','Direction','right','Sinks',busesT_cell_str,'AssignLayers','alap'),
                title('Transmission network',Interpreter='latex')
                
                nexttile
                hD = plot(graph_D,'k');
                hD.Interpreter = "latex";
                layout(hD,'layered','Direction','right','Sources',busesD_cell_str)
                title('Distribution network(s)',Interpreter='latex')

                tl.TileSpacing = "tight";
            else
                figure;                
                hT = plot(graph_T,'k');
                hT.Interpreter = "latex";
                %layout(hT,'layered','Direction','right','Sinks',busesT_cell_str,'AssignLayers','alap'),
                title('Transmission network',Interpreter='latex')

                figure;                
                hD = plot(graph_D,'k');
                hD.Interpreter = "latex";
                layout(hD,'layered','Direction','right','Sources',busesD_cell_str)
                title('Distribution network(s)',Interpreter='latex')
            end
            
            %% Customize graph plots
            % Classify nodes according to the presence of load, generation, or both
            nbT = size(mpc_T.bus,1);
            id_load_T = or(mpc_T.bus(:,PD)~=0, mpc_T.bus(:,QD)~=0);
            [~,id_bus_gen_T] = ismember(mpc_T.gen(:,GEN_BUS),mpc_T.bus(:,BUS_I));
            id_gen_T = zeros(nbT,1); id_gen_T(id_bus_gen_T) = 1; id_gen_T = logical(id_gen_T);
            id_mix_T = and(id_load_T,id_gen_T);
            id_load_T(id_mix_T) = 0;
            id_gen_T(id_mix_T) = 0;

            busesT_load = wgv.case_utils.vec2char(mpc_T.bus(id_load_T,BUS_I));
            busesT_gen = wgv.case_utils.vec2char(mpc_T.bus(id_gen_T,BUS_I));
            busesT_mix = wgv.case_utils.vec2char(mpc_T.bus(id_mix_T,BUS_I));

            nbD = size(mpc_D.bus,1);
            id_load_D = or(mpc_D.bus(:,PD)~=0, mpc_D.bus(:,QD)~=0);
            [~,id_bus_gen_D] = ismember(mpc_D.gen(:,GEN_BUS),mpc_D.bus(:,BUS_I));
            id_gen_D = zeros(nbD,1); id_gen_D(id_bus_gen_D) = 1; id_gen_D = logical(id_gen_D);
            id_mix_D = and(id_load_D,id_gen_D);
            id_load_D(id_mix_D) = 0;
            id_gen_D(id_mix_D) = 0;

            busesD_load = wgv.case_utils.vec2char(mpc_D.bus(id_load_D,BUS_I));
            busesD_gen = wgv.case_utils.vec2char(mpc_D.bus(id_gen_D,BUS_I));
            busesD_mix = wgv.case_utils.vec2char(mpc_D.bus(id_mix_D,BUS_I));

            % Customize load (red), generation (green) and mixed (yellow) color nodes
            highlight(hT,busesT_load,'NodeColor',[0.83 0.14 0.14]);
            highlight(hT,busesT_gen,'NodeColor',[0.25 0.80 0.54]);
            highlight(hT,busesT_mix,'NodeColor',[0.78 0.59 0.06]);

            highlight(hD,busesD_load,'NodeColor',[0.83 0.14 0.14]);
            highlight(hD,busesD_gen,'NodeColor',[0.25 0.80 0.54]);
            highlight(hD,busesD_mix,'NodeColor',[0.78 0.59 0.06]);

            % Identify branches with transformers
            id_xfmr_T = find(mpc_T.branch(:,TAP) ~= 0);
            id_xfmr_D = find(mpc_D.branch(:,TAP) ~= 0);

            % Customize colors for transformers
            highlight(hT,f_T(id_xfmr_T),t_T(id_xfmr_T),'EdgeColor',[0 0.44 0.74],'LineWidth',1.5)
            highlight(hD,f_D(id_xfmr_D),t_D(id_xfmr_D),'EdgeColor',[0 0.44 0.74],'LineWidth',1.5)

            % Highlight buses to be connected to form the T-D Interface
            highlight(hT,busesT_cell_str,'Marker','h','MarkerSize',8,'NodeColor',[0.54 0 0.54])
            highlight(hD,busesD_cell_str,'Marker','h','MarkerSize',8,'NodeColor',[0.54 0 0.54])
            %labelnode(hD,busesD_cell_str,cellfun(@(x,y)([x '(' y ')']),busesD_cell_str,busesT_cell_str,'UniformOutput',false))
        end

        function id_neigh = find_neighbors(mpc,id_bus)
            %

            define_constants;
            
            mpc = loadcase(mpc);

            f = mpc.branch(:,F_BUS);
            t = mpc.branch(:,T_BUS);

            G = graph(f,t);            
            id_neigh = neighbors(G,id_bus);
        end

        function charnum = vec2char(x)
            %
            if any(x)
                charnum = cellstr(num2str(x));
                charnum = cellfun(@strtrim, charnum, 'UniformOutput', false);
            else
                charnum = {};
            end
        end
    end     %% methods (Static)
end