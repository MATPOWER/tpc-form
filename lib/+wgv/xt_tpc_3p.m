classdef xt_tpc_3p < mp.extension
% wgv.xt_tpc_3p - |MATPOWER| extension to add unbalanced three-phase
%                 elements modeled through approximate formulations in
%                 Transformed Polar Coordinates (TPC)
%
% For power flow tasks, overrides the mp.task_pf class by the wgv.task_pf
% class included in the +wgv package, which is part of the tpc-form
% repository. Moreover, adds six new element types:
%
%   - 'bus3p_tpc' - 3-phase bus
%   - 'gen3p_tpc' - 3-phase generator
%   - 'load3p_tpc' - 3-phase load
%   - 'line3p_tpc' - 3-phase distribution line
%   - 'xfmr3p_tpc' - 3-phase transformer
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
            % Adds the task class for PF under TPC formulations

            task_class = @wgv.task_pf;
        end
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Add six classes to data model converter elements.
            %
            % For ``'mpc2`` data formats, adds the classes:
            %
            %   - mp.dmce_bus3p_mpc2
            %   - mp.dmce_gen3p_mpc2
            %   - mp.dmce_load3p_mpc2
            %   - mp.dmce_line3p_mpc2
            %   - mp.dmce_xfmr3p_mpc2
            %   - mp.dmce_buslink_mpc2

            switch fmt
                case 'mpc2'
                    dmc_elements = { ...
                        @mp.dmce_bus3p_mpc2, @mp.dmce_gen3p_mpc2, ...
                        @mp.dmce_load3p_mpc2, @mp.dmce_line3p_mpc2, ...
                        @mp.dmce_xfmr3p_mpc2, @mp.dmce_shunt3p_mpc2, ...
                        @wgv.dmce_tdinter_mpc2 ...
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
                case {'PF', 'CPF'}
                    dm_elements = { ...
                        @mp.dme_bus3p, @mp.dme_gen3p, @mp.dme_load3p, ...
                        @mp.dme_line3p, @mp.dme_xfmr3p, @mp.dme_shunt3p, ...
                        @wgv.dme_tdinter ...
                        };
                case 'OPF'
                    dm_elements = { ...
                        @mp.dme_bus3p_opf, @mp.dme_gen3p_opf, ...
                        @mp.dme_load3p_opf, @mp.dme_line3p_opf, ...
                        @mp.dme_xfmr3p_opf, @mp.dme_buslink_opf ...
                        };
                otherwise
                    dm_elements = {};
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            % Add six classes to network model elements.
            %
            % For AC *cartesian* voltage formulations, adds the classes:
            %
            %   - mp.nme_bus3p_acc
            %   - mp.nme_gen3p_acc
            %   - mp.nme_load3p
            %   - mp.nme_line3p
            %   - mp.nme_xfmr3p
            %   - mp.nme_buslink_acc
            %
            % For AC *polar* voltage formulations, adds the classes:
            %
            %   - mp.nme_bus3p_acp
            %   - mp.nme_gen3p_acp
            %   - mp.nme_load3p
            %   - mp.nme_line3p
            %   - mp.nme_xfmr3p
            %   - mp.nme_buslink_acp
            %
            % For TPC formulations, add the classes:
            %
            %   - wgv.nme_bus3p_tpc
            %   - wgv.nme_gen3p_tpc
            %   - wgv.nme_load3p_tpc
            %   - wgv.nme_line3p_tpc
            %   - wgv.nme_xfmr3p_tpc
            %   - wgv.nme_buslink_tpc


            switch task_tag
                case {'PF', 'CPF'}
                    v_cartesian = mpopt.pf.v_cartesian;
                case {'OPF'}
                    v_cartesian = mpopt.opf.v_cartesian;
            end
            switch upper(mpopt.model)
                case 'AC'
                    if v_cartesian
                        nm_elements = { ...
                            @mp.nme_bus3p_acc, @mp.nme_gen3p_acc, ...
                            @mp.nme_load3p, @mp.nme_line3p, ...
                            @mp.nme_xfmr3p, @mp.nme_buslink_acc ...
                            };
                    else
                        nm_elements = { ...
                            @mp.nme_bus3p_acp, @mp.nme_gen3p_acp, ...
                            @mp.nme_load3p, @mp.nme_line3p, ...
                            @mp.nme_xfmr3p, @mp.nme_buslink_acp ...
                            };
                    end
                case 'DC'
                    nm_elements = {};       %% no modifications
                case 'TPC'
                    nm_elements = { ...
                        @wgv.nme_bus3p_tpc, @wgv.nme_gen3p_tpc, ...
                        @wgv.nme_load3p_tpc, @wgv.nme_line3p_tpc, ...
                        @wgv.nme_xfmr3p_tpc, @wgv.nme_shunt3p_tpc, ...
                        @wgv.nme_tdinter_tpc ...
                        };
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Add five classes to mathematical model elements.
            %
            % For AC ``'PF'`` and ``'CPF'`` tasks, adds the classes:
            %
            %   - mp.mme_bus3p
            %   - mp.mme_gen3p
            %   - mp.mme_line3p
            %   - mp.mme_xfmr3p
            %   - mp.mme_buslink_pf_acc *(cartesian)* or
            %     mp.mme_buslink_pf_acp *(polar)*
            %
            % For AC ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.mme_bus3p_opf_acc *(cartesian)* or
            %     mp.mme_bus3p_opf_acp *(polar)*
            %   - mp.mme_gen3p_opf
            %   - mp.mme_line3p_opf
            %   - mp.mme_xfmr3p_opf
            %   - mp.mme_buslink_opf_acc *(cartesian)* or
            %     mp.mme_buslink_opf_acp *(polar)*
            %
            % For TPC formulations, add the classes:
            %
            %   - wgv.mme_bus3p_tpc
            %   - wgv.mme_gen3p_tpc
            %   - wgv.mme_load3p_tpc
            %   - wgv.mme_line3p_tpc
            %   - wgv.mme_xfmr3p_tpc
            %   - wgv.mme_buslink_tpc

            switch task_tag
                case {'PF', 'CPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.pf.v_cartesian
                                mm_elements = { ...
                                    @mp.mme_bus3p, @mp.mme_gen3p, ...
                                    @mp.mme_line3p, @mp.mme_xfmr3p, ...
                                    @mp.mme_buslink_pf_acc ...
                                    };
                            else
                                mm_elements = { ...
                                    @mp.mme_bus3p, @mp.mme_gen3p, ...
                                    @mp.mme_line3p, @mp.mme_xfmr3p, ...
                                    @mp.mme_buslink_pf_acp ...
                                    };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                        case 'TPC'
                            mm_elements = { ...
                                @wgv.mme_bus3p_pf_tpc, @wgv.mme_gen3p_pf_tpc, ...
                                @wgv.mme_line3p_pf_tpc, @wgv.mme_xfmr3p_pf_tpc ...
                                @wgv.mme_shunt3p_pf_tpc, @wgv.mme_tdinter_pf_tpc ...
                                };
                    end
                case {'OPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.opf.v_cartesian
                                mm_elements = { ...
                                    @mp.mme_bus3p_opf_acc, ...
                                    @mp.mme_gen3p_opf, ...
                                    @mp.mme_line3p_opf, ...
                                    @mp.mme_xfmr3p_opf, ...
                                    @mp.mme_buslink_opf_acc ...
                                    };
                            else
                                mm_elements = { ...
                                    @mp.mme_bus3p_opf_acp, ...
                                    @mp.mme_gen3p_opf, ...
                                    @mp.mme_line3p_opf, ...
                                    @mp.mme_xfmr3p_opf, ...
                                    @mp.mme_buslink_opf_acp ...
                                    };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                    end
            end
        end
    end     %% methods
end