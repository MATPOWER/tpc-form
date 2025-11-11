classdef nme_tdinter_tpc < mp.nm_element & wgv.form_tpc
% wgv.nme_tdinter_tpc - Network model element abstract base class for 1-to-3-phase
%                       interface for TPC formulation
%
% Implements the network model element for 1-to-3-phase interface elements,
% with 4 ports and 3 non-voltage states per interface.
%
% Adds non-voltage state variables ``Pinter`` and ``Qinter`` to the network
% model, builds the parameter :math:`\NN`, and constructs voltage constraints.

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
            name = 'buslink';
        end

        function np = np(obj)
            %
            np = 4;     %% this is a 4 port element
        end

        function nz = nz(obj)
            %
            nz = 3;     %% 3 complex non-voltage states per element
        end

        function obj = add_zvars(obj, nm, dm, idx)
            %
            p = idx{1};
            ng = obj.nk;

            if p == 1
                nm.init_indexed_name('zr', 'Pinter', {obj.nz});
                nm.init_indexed_name('zi', 'Qinter', {obj.nz});
            end
            nm.add_var('zr', 'Pinter', {p}, obj.nk, 0, -Inf, Inf);
            nm.add_var('zi', 'Qinter', {p}, obj.nk, 0, -Inf, Inf);
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);      %% call parent
            
            nbuslink = obj.nk;           

            cob = (dm.base_kva / 1000) / dm.base_mva;     %% coefficient for change of base
            I1p = speye(nbuslink);
            I3p = speye(3*nbuslink);
            N_buslink = [ repmat(cob*I1p,1,3)  repmat(1j*cob*I1p,1,3);
                                -1*I3p                -1j*I3p ];

            obj.N = N_buslink;
        end

        function [A, b_va, b_vm, Istack_] = voltage_constraints(obj)
            %

            %% form constraint matrices for matching voltages
            nk = obj.nk;

            %% basic constraint matrix for voltage equality (all nodes)
            Istack = repmat(speye(nk), obj.nz, 1);  %% stacked identities
            A = [ Istack -speye(nk * obj.nz) ] * obj.C';            

            %% RHS
            b_va = zeros(nk*obj.nz, 1);             %% angle constraints
            b_vm = zeros(nk*obj.nz, 1);             %% magnitude constraints

            if nargout > 3
                Istack_ = Istack;
            end
        end
    end     %% methods
end         %% classdef
