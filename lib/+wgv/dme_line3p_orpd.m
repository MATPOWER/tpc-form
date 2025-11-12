classdef dme_line3p_orpd < mp.dme_line3p
% wgv.dme_line3p_orpd - Data model element for 3-phase line for Optimal Reactive 
%                       Power Dispatch (ORPD) tasks
%
% Implements the data element model for 3-phase distribution line elements
% for ORPD.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  =============================================
%   Name         Type       Description
%   ===========  =========  =============================================
%   ``smax``     *double*   maximum per-phase apparent power *(kVA)* for
%                           three-phase line
%   ===========  =========  =============================================

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties        
        smax      % maximum per-phase apparent power (kVA) for three-phase lines that are on
    end     

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus_fr', 'bus_to', 'lc', 'len', ...
                 'pl1_fr', 'ql1_fr', 'pl2_fr', 'ql2_fr', 'pl3_fr', 'ql3_fr', ...
                 'pl1_to', 'ql1_to', 'pl2_to', 'ql2_to', 'pl3_to', 'ql3_to', ...
                 'smax'});
        end
        function obj = build_params(obj, dm)
            %
            nlc = size(obj.lc, 1);
            obj.ys = zeros(nlc, 6);
            obj.yc = zeros(nlc, 6);
            obj.lc = obj.tab.lc(obj.on);
            obj.len = obj.tab.len(obj.on);
            obj.smax = obj.tab.smax(obj.on) / dm.base_kva;

            %% build Ys and Yc for relevant lines
            idx = unique(obj.lc);   %% unique codes used by lines that are on
            tmpmap = zeros(max(idx), 1);
            tmpmap(idx) = (1:length(idx))';
            obj.lc_y_idx = tmpmap(obj.lc);  %% index into idx & hence into rows
                                            %% of ys/yc for lines that are on
            rr = obj.lc_tab.r(idx, :);
            xx = obj.lc_tab.x(idx, :);
            cc = obj.lc_tab.c(idx, :);
            for k = 1:length(idx)
                R = obj.vec2symmat(rr(k, :));
                X = obj.vec2symmat(xx(k, :));
                C = obj.vec2symmat(cc(k, :));
                Ys = inv(R + 1j * X);
                Yc = 1j * 2*pi * obj.freq * 1e-9 * C;
                obj.ys(k, :) = obj.symmat2vec(Ys);
                obj.yc(k, :) = obj.symmat2vec(Yc);
            end
        end
    end     %% methods
end         %% classdef
