classdef xt_tpc_3p_orpd < mp.extension
% wgv.xt_tpc_3p_orpd - |MATPOWER| extension to add unbalanced three-phase
%                      elements modeled through approximate formulations in
%                      Transformed Polar Coordinates (TPC) for Optimal
%                      Reactive Power Dispatch (ORPD) tasks
%
% For optimal reactive power dispatch tasks, overrides the mp.task_pf class 
% by the wgv.task_orpd class included in the +wgv package, which is part of 
% the tpc-form repository. Moreover, adds five data model an data model 
% converter elements:
%
%   - 'dmce_bus3p_orpd_mpc2' & 'dme_bus3p_orpd' - 3-phase bus
%   - 'dmce_gen3p_orpd_mpc2' & 'dme_gen3p_orpd' - 3-phase generator
%   - 'load3p_tpc' - 3-phase load
%   - 'dmce_line3p_orpd_mpc2' & 'dme_line3p_orpd' - 3-phase distribution line
%   - 'xfmr3p_tpc' - 3-phase transformer
%   - 'dmce_shunt3p_orpd_mpc2' & 'dme_shunt3p_orpd' - 3-phase shunt
%   - 'buslink_tpc' - 3-phase to single phase linking element

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%       properties
%       end

    methods
        function task_class = task_class(obj, task_class, mpopt)
            % Adds the task class for ORPD under TPC formulations

            task_class = @wgv.task_orpd;
        end
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Add six classes to data model converter elements.
            %
            % For ``'mpc2`` data formats, adds the classes:
            %
            %   - wgv.dmce_bus3p_orpd_mpc2
            %   - mp.dmce_gen3p_mpc2
            %   - mp.dmce_load3p_mpc2
            %   - mp.dmce_line3p_mpc2
            %   - mp.dmce_xfmr3p_mpc2
            %   - mp.dmce_buslink_mpc2

            switch fmt
                case 'mpc2'
                    dmc_elements = { ...
                        @wgv.dmce_bus3p_orpd_mpc2, @wgv.dmce_gen3p_orpd_mpc2, ...
                        @mp.dmce_load3p_mpc2, @wgv.dmce_line3p_orpd_mpc2, ...
                        @wgv.dmce_xfmr3p_orpd_mpc2, @wgv.dmce_shunt3p_orpd_mpc2, ...
                        @wgv.dmce_tdinter_orpd_mpc2 ...
                        };
                otherwise
                    dmc_elements = {};
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            % Add six classes to data model elements.
            %
            % For ``'PF'`` and ``'CPF'`` tasks, adds the classes:
            %
            %   - mp.dme_bus3p
            %   - mp.dme_gen3p
            %   - mp.dme_load3p
            %   - mp.dme_line3p
            %   - mp.dme_xfmr3p
            %   - mp.dme_buslink
            %
            % For ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.dme_bus3p_opf
            %   - mp.dme_gen3p_opf
            %   - mp.dme_load3p_opf
            %   - mp.dme_line3p_opf
            %   - mp.dme_xfmr3p_opf
            %   - mp.dme_buslink_opf

            switch task_tag                
                case {'ORPD'}
                    dm_elements = { ...
                        @wgv.dme_bus3p_orpd, @wgv.dme_gen3p_orpd, @mp.dme_load3p, ...
                        @wgv.dme_line3p_orpd, @wgv.dme_xfmr3p_orpd, ...
                        @wgv.dme_shunt3p_orpd, @wgv.dme_tdinter_orpd ...
                        };
                otherwise
                    error('xt_tpc_3p_orpd.nm_element_classes: curret implementation of ORPD tasks only works for TPC-based formulations.')
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            % Add six classes to network model elements.            
            %
            % For TPC formulations, add the classes:
            %
            %   - wgv.nme_bus3p_tpc
            %   - wgv.nme_gen3p_tpc
            %   - wgv.nme_load3p_tpc
            %   - wgv.nme_line3p_tpc
            %   - wgv.nme_xfmr3p_tpc
            %   - wgv.nme_buslink_tpc

            switch upper(mpopt.model)                
                case 'TPC'
                    changed_nm_elements = { ...
                        {@wgv.nme_branch_orpd_tpc, 'wgv.nme_branch_tpc'}, ...
                        {@wgv.nme_shunt_orpd_tpc, 'wgv.nme_shunt_tpc'}
                        };

                    new_nm_elements = { ...
                        @wgv.nme_bus3p_orpd_tpc, @wgv.nme_gen3p_orpd_tpc, ...
                        @wgv.nme_load3p_tpc, @wgv.nme_line3p_orpd_tpc, ...
                        @wgv.nme_xfmr3p_orpd_tpc, @wgv.nme_shunt3p_orpd_tpc, ...
                        @wgv.nme_tdinter_tpc ...
                        };

                    nm_elements = horzcat(changed_nm_elements,new_nm_elements);
                otherwise
                    error('xt_tpc_3p_orpd.nm_element_classes: curret implementation of ORPD tasks only works for TPC-based formulations.')
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Add five classes to mathematical model elements.            
            %
            % For TPC formulations, add the classes:
            %
            %   - wgv.mme_bus3p_orpd_tpc
            %   - wgv.mme_gen3p_tpc
            %   - wgv.mme_load3p_tpc
            %   - wgv.mme_line3p_tpc
            %   - wgv.mme_xfmr3p_tpc
            %   - wgv.mme_buslink_tpc

            switch task_tag                
                case {'ORPD'}
                    mm_elements = { ...
                        @wgv.mme_bus3p_orpd_tpc, @wgv.mme_gen3p_orpd_tpc, ...
                        @wgv.mme_line3p_orpd_tpc, @wgv.mme_xfmr3p_orpd_tpc, ...
                        @wgv.mme_shunt3p_orpd_tpc, @wgv.mme_tdinter_orpd_tpc ...
                        };
                otherwise
                    error('xt_tpc_3p_orpd.nm_element_classes: curret implementation of ORPD tasks only works for TPC-based formulations.')
            end
        end
    end     %% methods
end