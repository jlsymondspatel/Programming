%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 01702979 Jaimin Luke Symonds Patel (10/03/2021)                %%%
%%% Computational Methods Spring Coursework                        %%%
%%% Collapse of a redundant truss by incremental elastic analysis  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                %%%
%%%  Script file that utilises the Truss class to perform iL-MNA   %%%
%%%                         and iL-GMNA                            %%%
%%%                                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clearing any old data
clc;
clear all;
close all;

%% Assembly of properties which define a truss object

% define properties relevant to the bar elements in an ELEMENT struct.
ELEMENTS = struct('E',[],'R',[],'sigma_y',[],'L',[],'element_connectivity',[]);
% define the elastic Young's modulus [N/m2]:
ELEMENTS.E = 200e9;
% define the radius of each bar element [m]:
ELEMENTS.R = 20e-3;
% define the yield stress [N/m2]:
ELEMENTS.sigma_y = 250e6;
% define the horizontal length unit of each truss brace [m]:
ELEMENTS.L = 1000e-3;
% define element nodal connectivity
ELEMENTS.element_connectivity = ...
    [1 3;   % element 1 with element dofs u1,v1,u3,v3 (1,2,5,6)
    1 4;    % element 2 with element dofs u1,v1,u4,v4 (1,2,7,8)
    2 3;    % element 3 with element dofs u2,v2,u3,v3 (3,4,5,6)
    2 4;    % element 4 with element dofs u2,v2,u4,v4 (3,4,7,8)
    3 5;    % element 5 with element dofs u3,v3,u5,v5 (5,6,9,10)
    3 6;    % element 6 with element dofs u3,v3,u6,v6 (5,6,11,12)
    3 4;    % element 7 with element dofs u3,v3,u4,v4 (5,6,7,8)
    4 5;    % element 8 with element dofs u4,v4,u5,v5 (7,8,9,10)
    4 6;    % element 9 with element dofs u4,v4,u6,v6 (7,8,11,12)
    5 7;    % element 10 with element dofs u5,v5,u7,v7 (9,10,13,14)
    5 8;    % element 11 with element dofs u5,v5,u8,v8 (9,10,15,16)
    5 6;    % element 12 with element dofs u5,v5,u6,v6 (9,10,11,12)
    6 7;    % element 13 with element dofs u6,v6,u7,v7 (11,12,13,14)
    6 8;    % element 14 with element dofs u6,v6,u8,v8 (11,12,15,16)
    7 9;    % element 15 with element dofs u7,v7,u9,v9 (13,14,17,18)
    7 10;   % element 16 with element dofs u7,v7,u10,v10 (13,14,19,20)
    7 8;    % element 17 with element dofs u7,v7,u8,v8 (13,14,15,16)
    8 9;    % element 18 with element dofs u8,v8,u9,v9 (15,16,17,18)
    8 10];  % element 19 with element dofs u8,v8,u10,v10 (15,16,19,20)

% define properties relevant to the nodes in a NODES struct.
NODES = struct('coords',[],'dofs',[],'dofs_free',[],'dofs_restrained',[]);
% define the coordinates of each node:
NODES.coords = ...
    [0 ELEMENTS.L;          % node 1 x y coordinate
    0 0;                    % node 2 x y coordinate
    ELEMENTS.L ELEMENTS.L;  % node 3 x y coordinate
    ELEMENTS.L 0;           % node 4 x y coordinate
    2*ELEMENTS.L ELEMENTS.L;% node 5 x y coordinate
    2*ELEMENTS.L 0;         % node 6 x y coordinate
    3*ELEMENTS.L ELEMENTS.L;% node 7 x y coordinate
    3*ELEMENTS.L 0;         % node 8 x y coordinate
    4*ELEMENTS.L ELEMENTS.L;% node 9 x y coordinate
    4*ELEMENTS.L 0];        % node 10 x y coordinate
% define degrees of freedom for each node:
NODES.dofs = ...
    [1 2;   % node 1 - dofs u1,v1 - FIXED
    3 4;    % node 2 - dofs u2,v2 - FIXED
    5 6;    % node 3 - dofs u3,v3 - FREE
    7 8;    % node 4 - dofs u4,v4 - FREE
    9 10;   % node 5 - dofs u5,v5 - FREE
    11 12;  % node 6 - dofs u6,v6 - FREE
    13 14;  % node 7 - dofs u7,v7 - FREE
    15 16;  % node 8 - dofs u8,v8 - FREE
    17 18;  % node 9 - dofs u9,v9 - FIXED
    19 20]; % node 10 - dofs u10,v10 - FIXED
% define the free dofs (dofs 5 to 16):
NODES.dofs_free = [5 6 7 8 9 10 11 12 13 14 15 16];
% define the free dofs (dofs 1 to 4 and 17 to 20):
NODES.dofs_restrained = [1 2 3 4 17 18 19 20];

% define properties releveant to the loading in a LOADING struct.
LOADING = struct('P',[]);
% define unit point load [N]:
LOADING.P = 100e3;
% define the number of point loads and their distribution not taking into
% account the load proportionality factors.
LOADING.applied_load_distribution = zeros(3,2);
% define the left most point load acting on node 3, dof 6:
LOADING.applied_load_distribution(1,1) = -3*LOADING.P;
LOADING.applied_load_distribution(1,2) = 6;
% define the middle point load acting on node 5, dof 10:
LOADING.applied_load_distribution(2,1) = -1.5*LOADING.P;
LOADING.applied_load_distribution(2,2) = 10;
% define the right most point load acting on node 7, dof 14:
LOADING.applied_load_distribution(3,1) = -LOADING.P;
LOADING.applied_load_distribution(3,2) = 14;

%% Question 2 - MNA analysis

% constructing a Truss object for use with MNA analysis (so as to keep
% data):
TRUSS = Truss(NODES,ELEMENTS,LOADING);

% initialise the vector of element total axial forces as zero to begin
% with:
TRUSS.total_element_axial_forces = zeros(TRUSS.num_of_elements,1);

% initialise a vector of the total nodal displacements with zeros:
TRUSS.total_nodal_displacements = zeros(TRUSS.num_of_nodes*2,1);

% initialise the total load proportionality factor as zero to begin with
TRUSS.total_load_prop_factor_MNA = 0;

% define a vector which holds a list of all the intact (not failed)
% members, so that the element numbering may be preserved even when
% elements fail:
TRUSS.intact_elements_MNA = [1:TRUSS.num_of_elements]';

% assign an initial value of 1 to the increment of the load proportionality
% factor so that it is unscaled:
d_lambda_unscaled = 1;

% create a vector to hold a list of all the increments in the load
% proportionality factor:
TRUSS.d_lambda_list_MNA = [];

% define the plastic squash load:
yield_load = TRUSS.A * TRUSS.sigma_y;

% define an empty vector that holds the redundancies at each iteration
TRUSS.redundancies_list_MNA = [];

% at each iteration find the critical member, where the maximum possible
% iterations is the number of members which can yield
for I = 1:size(ELEMENTS.element_connectivity,1)
    % construct global stiffness matrix:
    TRUSS = TRUSS.constructGlobalStiffnessMatrix();
    % check if the global stiffness matrix is invertable, and so check if
    % the FE matrix system is solvable:
    TRUSS = TRUSS.checkSolvability();
    % if solvable, proceed to solve matrix system.
    if TRUSS.solvable
        % do nothing, and continue
    else
        % if matrix system is unsolvable, break out of the MNA iterative
        % process (for loop):
        fprintf('##Matrix system is now unsolvable##\n##(truss has too many failed members and has become a mechanism)##\n');
        break;
    end
    % construct global nodal force vector from unscaled loads
    TRUSS = TRUSS.constructGlobalNodalForceVector(d_lambda_unscaled);
    % solve matrix system to find nodal displacements and reaction forces:
    TRUSS = TRUSS.solveFiniteElementMatrixSystem();
    % find the new deformed nodal coordinates arising from the nodal
    % displacements:
    TRUSS = TRUSS.obtainNewNodalCoordinates();
    % Get increment in axial forces within each element from an unscaled
    % load proportionality factor on top of the previous total load
    % proportionality factor:
    TRUSS = TRUSS.obtainElementAxialForce();
    
    % initialise a vector of each element's increment in the load
    % proportionality factor in order to find the smallest one:
    element_d_lambdas = zeros(TRUSS.num_of_elements,1);
    % initialise a vector of the failure modes of all elements if they were
    % to fail:
    failure_modes = zeros(TRUSS.num_of_elements,1);
    % looping through all elements to find the lowest increment in the load
    % proportionality factor:
    for el = 1:TRUSS.num_of_elements
        % if the element is going to fail in tension.
        % find the elements increment in the load proportionality
        % factor:
        element_d_lambdas(el) = ...
            (yield_load - TRUSS.total_element_axial_forces(el))/TRUSS.element_axial_forces(el);
        % take note note of failure mode (compression):
        failure_modes(el) = 'T';
        if element_d_lambdas(el) < 0
            % if the element is going to fail in compression.
            % find the elements increment in the load proportionality
            % factor:
            element_d_lambdas(el) = ...
                (-yield_load - TRUSS.total_element_axial_forces(el))/TRUSS.element_axial_forces(el);
            % take note note of failure mode (compression):
            failure_modes(el) = 'C';
        end
    end
    
    % find the minimum increment in the load proportionality factor, and
    % which critical element it belongs to:
    [d_lambda, critical_element] = min(element_d_lambdas);
    % make not of the load proportionality factor increment:
    TRUSS.d_lambda_list_MNA(end+1,1) = d_lambda;
    % making note of the failure mode of the critical element:
    failure_mode = failure_modes(critical_element);
    % add the critial element to a list of critical elements that have
    % failed in each iteration as well as their failure mode:
    TRUSS.failed_elements_MNA(end+1,1) = TRUSS.intact_elements_MNA(critical_element);
    TRUSS.failure_modes_MNA(end+1,1) = char(failure_mode);
    % add the increment in the load proportionality factor to the total
    % proportionality factor:
    TRUSS.total_load_prop_factor_MNA = TRUSS.total_load_prop_factor_MNA + d_lambda;
    
    % apply the scaled load proportionality factor to find the axial
    % forces in all non failed elements.
    % construct global nodal force vector from scaled loads:
    TRUSS = TRUSS.constructGlobalNodalForceVector(d_lambda);
    % solve matrix system to find nodal displacements and reaction forces:
    TRUSS = TRUSS.solveFiniteElementMatrixSystem();
    % find the new deformed nodal coordinates arising from the nodal
    % displacements:
    TRUSS = TRUSS.obtainNewNodalCoordinates();
    % Get increment in axial forces within each element from the scaled
    % load proportionality factor:
    TRUSS = TRUSS.obtainElementAxialForce();
    % add the correctly scaled loads whereby the critical member is made
    % to fail:
    TRUSS.total_element_axial_forces = TRUSS.total_element_axial_forces + TRUSS.element_axial_forces;
    % add the increment in the nodal displacements to a vector of total
    % displacements:
    TRUSS.total_nodal_displacements = TRUSS.total_nodal_displacements + TRUSS.U;
    
    % replace the critical member, which has yielded, with a pair of nodal
    % forces for the next iteration, and so have four seperate forces on
    % four different dofs, depending on the orientation of the failed bar
    % member.
    % obtain the critical element's force as the axial force at failure:
    critical_element_force = TRUSS.total_element_axial_forces(critical_element);
    % replace the failed element with equivilent nodal forces in failure:
    TRUSS = TRUSS.replaceFailedElementWithForces(critical_element,critical_element_force);
    
    % get rid of the yielded element from the element connectivity
    % vector:
    TRUSS.element_connectivity(critical_element,:) = [];
    % subtract the failed element also from the list of intact elements:
    TRUSS.intact_elements_MNA(critical_element) = [];
    % update the number of elements remaining:
    TRUSS.num_of_elements = size(TRUSS.element_connectivity,1);
    % get rid of the failed element from the vector of axial element
    % forces:
    TRUSS.total_element_axial_forces(critical_element) = [];
    
    % work out number of redundancies and them to a list:
    redundancies = TRUSS.num_of_elements + 8 -(2*TRUSS.num_of_nodes);
    TRUSS.redundancies_list_MNA(end+1,1) = redundancies;
end

% assign the iteration counter to the truss object MNA counter, where the
% final iteration is the one after the last failed member has yielded, and
% so 1 must be taken from the final iteration counter to give the iteration
% which causes the whole truss to fail structurally:
TRUSS.MNA_counter = I-1;

% define a cumulative d_lambda for the printed table
d_lambda_cumulative = 0;
% print all the relevant information to screen with a table:
fprintf('Incremental Linear MNA Analysis results:\n');
fprintf('Event\t\tFailure mode\tRedundancies\td_lambda\t\t\tlambda_tot\n');
fprintf('--------------------------------------------------------------------------\n');
for i = 1:length(TRUSS.failed_elements_MNA)
    d_lambda_cumulative = d_lambda_cumulative + TRUSS.d_lambda_list_MNA(i);
    fprintf('Element %d \t\t %c \t\t\t %d \t\t\t\t %.3f \t\t\t\t %.3f \n',TRUSS.failed_elements_MNA(i),TRUSS.failure_modes_MNA(i),TRUSS.redundancies_list_MNA(i),TRUSS.d_lambda_list_MNA(i),d_lambda_cumulative);
end

%% Question 3 - GMNA analysis

% constructing a Truss object for GMNA analysis (to keep the previous
% data):
TRUSS2 = Truss(NODES,ELEMENTS,LOADING);

% initialise the vector of element total axial forces as zero to begin
% with:
TRUSS2.total_element_axial_forces = zeros(TRUSS2.num_of_elements,1);

% initialise a vector of the total nodal displacements with zeros:
TRUSS2.total_nodal_displacements = zeros(TRUSS2.num_of_nodes*2,1);

% initialise the total load proportionality factor as zero to begin with
TRUSS2.total_load_prop_factor_GMNA = 0;

% define a vector which holds a list of all the intact (not failed)
% members, so that the element numbering may be preserved even when
% elements fail:
TRUSS2.intact_elements_GMNA = [1:TRUSS2.num_of_elements]';

% assign an initial value of 1 to the increment of the load proportionality
% factor so that it is unscaled:
d_lambda_unscaled = 1;

% create a vector to hold a list of all the increments in the load
% proportionality factor:
TRUSS2.d_lambda_list_GMNA = [];

% define the plastic squash load:
yield_load = TRUSS2.A * TRUSS2.sigma_y;

% define an empty vector that holds the redundancies at each iteration
TRUSS2.redundancies_list_GMNA = [];

% at each iteration find the critical member, where the maximum possible
% iterations is the number of members which can yield
for I = 1:size(ELEMENTS.element_connectivity,1)
    % construct global stiffness matrix:
    TRUSS2 = TRUSS2.constructGlobalStiffnessMatrix();
    % check if the global stiffness matrix is invertable, and so check if
    % the FE matrix system is solvable:
    TRUSS2 = TRUSS2.checkSolvability();
    % if solvable, proceed to solve matrix system.
    if TRUSS2.solvable
        % do nothing, and continue
    else
        % if matrix system is unsolvable, break out of the MNA iterative
        % process (for loop):
        fprintf('__________________________________________________________________________\n');
        fprintf('##Matrix system is now unsolvable##\n##(truss has too many failed members and has become a mechanism)##\n');
        break;
    end
    % construct global nodal force vector from unscaled loads
    TRUSS2 = TRUSS2.constructGlobalNodalForceVector(d_lambda_unscaled);
    % solve matrix system to find nodal displacements and reaction forces:
    TRUSS2 = TRUSS2.solveFiniteElementMatrixSystem();
    % find the new deformed nodal coordinates arising from the nodal
    % displacements:
    TRUSS2 = TRUSS2.obtainNewNodalCoordinates();
    % Get increment in axial forces within each element from an unscaled
    % load proportionality factor on top of the previous total load
    % proportionality factor:
    TRUSS2 = TRUSS2.obtainElementAxialForce();
    
    % initialise a vector of each element's increment in the load
    % proportionality factor in order to find the smallest one:
    element_d_lambdas = zeros(TRUSS2.num_of_elements,1);
    % initialise a vector of the failure modes of all elements if they were
    % to fail:
    failure_modes = zeros(TRUSS2.num_of_elements,1);
    % looping through all elements to find the lowest increment in the load
    % proportionality factor:
    for el = 1:TRUSS2.num_of_elements
        % if the element is going to fail in tension.
        % find the elements increment in the load proportionality
        % factor:
        element_d_lambdas(el) = ...
            (yield_load - TRUSS2.total_element_axial_forces(el))/TRUSS2.element_axial_forces(el);
        % take note note of failure mode (compression):
        failure_modes(el) = 'T';
        if element_d_lambdas(el) < 0
            % define the length of the element.
            % obtain global node numers of local nodes:
            n1 = TRUSS2.element_connectivity(el,1);
            n2 = TRUSS2.element_connectivity(el,2);
            % obtaining original x and y coordinate of local node 1:
            x1_original = TRUSS2.coords(n1,1);
            y1_original = TRUSS2.coords(n1,2);
            % obtaining original x and y coordinate of local node 2:
            x2_original = TRUSS2.coords(n2,1);
            y2_original = TRUSS2.coords(n2,2);
            % obtaining original length of unstretched element:
            el_L = sqrt((x2_original-x1_original)^2 + (y2_original-y1_original)^2);
            % define the euler buckling load for each element:
            buckling_load = -(pi^2)*(TRUSS2.EI/(el_L^2));
            
            % if the element is going to fail in compression.
            % find the elements increment in the load proportionality
            % factor only taking into account plastic squash:
            element_d_lambdas(el) = ...
                (-yield_load - TRUSS2.total_element_axial_forces(el))/TRUSS2.element_axial_forces(el);
            % take note note of failure mode (compression):
            failure_modes(el) = 'C';
            
            % find the increment in load proportionality factor now taking
            % into account buckling:
            element_d_lambda_buckling = ...
                (buckling_load - TRUSS2.total_element_axial_forces(el))/TRUSS2.element_axial_forces(el);
            
            % if the increment in the load proportionality factor with
            % buckling is less than that with plastic squash, use the
            % increment obtained with buckling:
            if element_d_lambda_buckling < element_d_lambdas(el)
                element_d_lambdas(el) = element_d_lambda_buckling;
                failure_modes(el) = 'B';
            end
        end
    end
    
    % find the minimum increment in the load proportionality factor, and
    % which critical element it belongs to:
    [d_lambda, critical_element] = min(element_d_lambdas);
    % make not of the load proportionality factor increment:
    TRUSS2.d_lambda_list_GMNA(end+1,1) = d_lambda;
    % making note of the failure mode of the critical element:
    failure_mode = failure_modes(critical_element);
    % add the critial element to a list of critical elements that have
    % failed in each iteration as well as their failure mode:
    TRUSS2.failed_elements_GMNA(end+1,1) = TRUSS2.intact_elements_GMNA(critical_element);
    TRUSS2.failure_modes_GMNA(end+1,1) = char(failure_mode);
    % add the increment in the load proportionality factor to the total
    % proportionality factor:
    TRUSS2.total_load_prop_factor_GMNA = TRUSS2.total_load_prop_factor_GMNA + d_lambda;
    
    % apply the scaled load proportionality factor to find the axial
    % forces in all non failed elements.
    % construct global nodal force vector from scaled loads:
    TRUSS2 = TRUSS2.constructGlobalNodalForceVector(d_lambda);
    % solve matrix system to find nodal displacements and reaction forces:
    TRUSS2 = TRUSS2.solveFiniteElementMatrixSystem();
    % find the new deformed nodal coordinates arising from the nodal
    % displacements:
    TRUSS2 = TRUSS2.obtainNewNodalCoordinates();
    % Get increment in axial forces within each element from the scaled
    % load proportionality factor:
    TRUSS2 = TRUSS2.obtainElementAxialForce();
    % add the correctly scaled loads whereby the critical member is made
    % to fail:
    TRUSS2.total_element_axial_forces = TRUSS2.total_element_axial_forces + TRUSS2.element_axial_forces;
    % add the increment in the nodal displacements to a vector of total
    % displacements:
    TRUSS2.total_nodal_displacements = TRUSS2.total_nodal_displacements + TRUSS2.U;
    
    % replace the critical member, which has yielded, with a pair of nodal
    % forces for the next iteration, and so have four seperate forces on
    % four different dofs, depending on the orientation of the failed bar
    % member.
    % obtain the critical element's force as the axial force at failure:
    critical_element_force = TRUSS2.total_element_axial_forces(critical_element);
    % replace the failed element with equivilent nodal forces in failure:
    TRUSS2 = TRUSS2.replaceFailedElementWithForces(critical_element,critical_element_force);
    
    % get rid of the yielded element from the element connectivity
    % vector:
    TRUSS2.element_connectivity(critical_element,:) = [];
    % subtract the failed element also from the list of intact elements:
    TRUSS2.intact_elements_GMNA(critical_element) = [];
    % update the number of elements remaining:
    TRUSS2.num_of_elements = size(TRUSS2.element_connectivity,1);
    % get rid of the failed element from the vector of axial element
    % forces:
    TRUSS2.total_element_axial_forces(critical_element) = [];
    
    % work out number of redundancies and them to a list:
    redundancies = TRUSS2.num_of_elements + 8 -(2*TRUSS2.num_of_nodes);
    TRUSS2.redundancies_list_GMNA(end+1,1) = redundancies;
end

% assign the iteration counter to the truss object MNA counter, where the
% final iteration is the one after the last failed member has yielded, and
% so 1 must be taken from the final iteration counter to give the iteration
% which causes the whole truss to fail structurally:
TRUSS2.GMNA_counter = I-1;

% define a cumulative d_lambda for the printed table
d_lambda_cumulative = 0;
% print all the relevant information to screen with a table:
fprintf('Incremental Linear GMNA Analysis results:\n');
fprintf('Event\t\tFailure mode\tRedundancies\td_lambda\t\t\tlambda_tot\n');
fprintf('--------------------------------------------------------------------------\n');
for i = 1:length(TRUSS2.failed_elements_GMNA)
    d_lambda_cumulative = d_lambda_cumulative + TRUSS2.d_lambda_list_GMNA(i);
    fprintf('Element %d \t\t %c \t\t\t %d \t\t\t\t %.3f \t\t\t\t %.3f \n',TRUSS2.failed_elements_GMNA(i),TRUSS2.failure_modes_GMNA(i),TRUSS2.redundancies_list_GMNA(i),TRUSS2.d_lambda_list_GMNA(i),d_lambda_cumulative);
end

%% Question 5 - Varying radius
% in order to perform several analyses in a for loop, a vector for all the
% radii needs to be made:
radii = [1e-3:1e-3:40e-3]';
% a vector of all the total load proportionalities must also be made for
% MNA and GMNA analyses:
lambda_totals_MNA = zeros(length(radii),1);
lambda_totals_GMNA = zeros(length(radii),1);
% creating vectors to log the number of failed members within each
% analysis:
num_of_failed_elements_MNA = zeros(length(radii),1);
num_of_failed_elements_GMNA = zeros(length(radii),1);

% firstly for several MNA analyses:
for analysis = 1:length(radii)
    % replace the give radius with the varying radius:
    ELEMENTS.R = radii(analysis);
    % constructing a Truss object for use with MNA analysis (so as to keep
    % data):
    TRUSS3 = Truss(NODES,ELEMENTS,LOADING);
    
    % initialise the vector of element total axial forces as zero to begin
    % with:
    TRUSS3.total_element_axial_forces = zeros(TRUSS3.num_of_elements,1);
    
    % initialise a vector of the total nodal displacements with zeros:
    TRUSS3.total_nodal_displacements = zeros(TRUSS3.num_of_nodes*2,1);
    
    % initialise the total load proportionality factor as zero to begin with
    TRUSS3.total_load_prop_factor_MNA = 0;
    
    % define a vector which holds a list of all the intact (not failed)
    % members, so that the element numbering may be preserved even when
    % elements fail:
    TRUSS3.intact_elements_MNA = [1:TRUSS3.num_of_elements]';
    
    % assign an initial value of 1 to the increment of the load proportionality
    % factor so that it is unscaled:
    d_lambda_unscaled = 1;
    
    % create a vector to hold a list of all the increments in the load
    % proportionality factor:
    TRUSS3.d_lambda_list_MNA = [];
    
    % define the plastic squash load:
    yield_load = TRUSS3.A * TRUSS3.sigma_y;
    
    % define an empty vector that holds the redundancies at each iteration
    TRUSS3.redundancies_list_MNA = [];
    
    % at each iteration find the critical member, where the maximum possible
    % iterations is the number of members which can yield
    for I = 1:size(ELEMENTS.element_connectivity,1)
        % construct global stiffness matrix:
        TRUSS3 = TRUSS3.constructGlobalStiffnessMatrix();
        % check if the global stiffness matrix is invertable, and so check if
        % the FE matrix system is solvable:
        TRUSS3 = TRUSS3.checkSolvability();
        % if solvable, proceed to solve matrix system.
        if TRUSS3.solvable
            % do nothing, and continue
        else
            % if matrix system is unsolvable, break out of the MNA iterative
            % process (for loop):
            % (withhold from printing message due to several analyses):
            % fprintf('##Matrix system is now unsolvable##\n##(truss has too many failed members and has become a mechanism)##\n');
            break;
        end
        % construct global nodal force vector from unscaled loads
        TRUSS3 = TRUSS3.constructGlobalNodalForceVector(d_lambda_unscaled);
        % solve matrix system to find nodal displacements and reaction forces:
        TRUSS3 = TRUSS3.solveFiniteElementMatrixSystem();
        % find the new deformed nodal coordinates arising from the nodal
        % displacements:
        TRUSS3 = TRUSS3.obtainNewNodalCoordinates();
        % Get increment in axial forces within each element from an unscaled
        % load proportionality factor on top of the previous total load
        % proportionality factor:
        TRUSS3 = TRUSS3.obtainElementAxialForce();
        
        % initialise a vector of each element's increment in the load
        % proportionality factor in order to find the smallest one:
        element_d_lambdas = zeros(TRUSS3.num_of_elements,1);
        % initialise a vector of the failure modes of all elements if they were
        % to fail:
        failure_modes = zeros(TRUSS3.num_of_elements,1);
        % looping through all elements to find the lowest increment in the load
        % proportionality factor:
        for el = 1:TRUSS3.num_of_elements
            % if the element is going to fail in tension.
            % find the elements increment in the load proportionality
            % factor:
            element_d_lambdas(el) = ...
                (yield_load - TRUSS3.total_element_axial_forces(el))/TRUSS3.element_axial_forces(el);
            % take note note of failure mode (compression):
            failure_modes(el) = 'T';
            if element_d_lambdas(el) < 0
                % if the element is going to fail in compression.
                % find the elements increment in the load proportionality
                % factor:
                element_d_lambdas(el) = ...
                    (-yield_load - TRUSS3.total_element_axial_forces(el))/TRUSS3.element_axial_forces(el);
                % take note note of failure mode (compression):
                failure_modes(el) = 'C';
            end
        end
        
        % find the minimum increment in the load proportionality factor, and
        % which critical element it belongs to:
        [d_lambda, critical_element] = min(element_d_lambdas);
        % make not of the load proportionality factor increment:
        TRUSS3.d_lambda_list_MNA(end+1,1) = d_lambda;
        % making note of the failure mode of the critical element:
        failure_mode = failure_modes(critical_element);
        % add the critial element to a list of critical elements that have
        % failed in each iteration as well as their failure mode:
        TRUSS3.failed_elements_MNA(end+1,1) = TRUSS3.intact_elements_MNA(critical_element);
        TRUSS3.failure_modes_MNA(end+1,1) = char(failure_mode);
        % add the increment in the load proportionality factor to the total
        % proportionality factor:
        TRUSS3.total_load_prop_factor_MNA = TRUSS3.total_load_prop_factor_MNA + d_lambda;
        
        % apply the scaled load proportionality factor to find the axial
        % forces in all non failed elements.
        % construct global nodal force vector from scaled loads:
        TRUSS3 = TRUSS3.constructGlobalNodalForceVector(d_lambda);
        % solve matrix system to find nodal displacements and reaction forces:
        TRUSS3 = TRUSS3.solveFiniteElementMatrixSystem();
        % find the new deformed nodal coordinates arising from the nodal
        % displacements:
        TRUSS3 = TRUSS3.obtainNewNodalCoordinates();
        % Get increment in axial forces within each element from the scaled
        % load proportionality factor:
        TRUSS3 = TRUSS3.obtainElementAxialForce();
        % add the correctly scaled loads whereby the critical member is made
        % to fail:
        TRUSS3.total_element_axial_forces = TRUSS3.total_element_axial_forces + TRUSS3.element_axial_forces;
        % add the increment in the nodal displacements to a vector of total
        % displacements:
        TRUSS3.total_nodal_displacements = TRUSS3.total_nodal_displacements + TRUSS3.U;
        
        % replace the critical member, which has yielded, with a pair of nodal
        % forces for the next iteration, and so have four seperate forces on
        % four different dofs, depending on the orientation of the failed bar
        % member.
        % obtain the critical element's force as the axial force at failure:
        critical_element_force = TRUSS3.total_element_axial_forces(critical_element);
        % replace the failed element with equivilent nodal forces in failure:
        TRUSS3 = TRUSS3.replaceFailedElementWithForces(critical_element,critical_element_force);
        
        % get rid of the yielded element from the element connectivity
        % vector:
        TRUSS3.element_connectivity(critical_element,:) = [];
        % subtract the failed element also from the list of intact elements:
        TRUSS3.intact_elements_MNA(critical_element) = [];
        % update the number of elements remaining:
        TRUSS3.num_of_elements = size(TRUSS3.element_connectivity,1);
        % get rid of the failed element from the vector of axial element
        % forces:
        TRUSS3.total_element_axial_forces(critical_element) = [];
        
        % work out number of redundancies and them to a list:
        redundancies = TRUSS3.num_of_elements + 8 -(2*TRUSS3.num_of_nodes);
        TRUSS3.redundancies_list_MNA(end+1,1) = redundancies;
    end
    
    % log the total lambda:
    lambda_totals_MNA(analysis) = TRUSS3.total_load_prop_factor_MNA;
    % log the number of failed elements in this particular analysis:
    num_of_failed_elements_MNA(analysis) = length(TRUSS3.failed_elements_MNA);
end

% secondly for several GMNA analyses:
for analysis = 1:length(radii)
    % replace the give radius with the varying radius:
    ELEMENTS.R = radii(analysis);
    % constructing a Truss object for GMNA analysis (to keep the previous
    % data):
    TRUSS3 = Truss(NODES,ELEMENTS,LOADING);
    
    % initialise the vector of element total axial forces as zero to begin
    % with:
    TRUSS3.total_element_axial_forces = zeros(TRUSS3.num_of_elements,1);
    
    % initialise a vector of the total nodal displacements with zeros:
    TRUSS3.total_nodal_displacements = zeros(TRUSS3.num_of_nodes*2,1);
    
    % initialise the total load proportionality factor as zero to begin with
    TRUSS3.total_load_prop_factor_GMNA = 0;
    
    % define a vector which holds a list of all the intact (not failed)
    % members, so that the element numbering may be preserved even when
    % elements fail:
    TRUSS3.intact_elements_GMNA = [1:TRUSS3.num_of_elements]';
    
    % assign an initial value of 1 to the increment of the load proportionality
    % factor so that it is unscaled:
    d_lambda_unscaled = 1;
    
    % create a vector to hold a list of all the increments in the load
    % proportionality factor:
    TRUSS3.d_lambda_list_GMNA = [];
    
    % define the plastic squash load:
    yield_load = TRUSS3.A * TRUSS3.sigma_y;
    
    % define an empty vector that holds the redundancies at each iteration
    TRUSS3.redundancies_list_GMNA = [];
    
    % at each iteration find the critical member, where the maximum possible
    % iterations is the number of members which can yield
    for I = 1:size(ELEMENTS.element_connectivity,1)
        % construct global stiffness matrix:
        TRUSS3 = TRUSS3.constructGlobalStiffnessMatrix();
        % check if the global stiffness matrix is invertable, and so check if
        % the FE matrix system is solvable:
        TRUSS3 = TRUSS3.checkSolvability();
        % if solvable, proceed to solve matrix system.
        if TRUSS3.solvable
            % do nothing, and continue
        else
            % if matrix system is unsolvable, break out of the MNA iterative
            % process (for loop)
            % (withhold from printing message due to several analyses):
            % fprintf('__________________________________________________________________________\n');
            % fprintf('##Matrix system is now unsolvable##\n##(truss has too many failed members and has become a mechanism)##\n');
            break;
        end
        % construct global nodal force vector from unscaled loads
        TRUSS3 = TRUSS3.constructGlobalNodalForceVector(d_lambda_unscaled);
        % solve matrix system to find nodal displacements and reaction forces:
        TRUSS3 = TRUSS3.solveFiniteElementMatrixSystem();
        % find the new deformed nodal coordinates arising from the nodal
        % displacements:
        TRUSS3 = TRUSS3.obtainNewNodalCoordinates();
        % Get increment in axial forces within each element from an unscaled
        % load proportionality factor on top of the previous total load
        % proportionality factor:
        TRUSS3 = TRUSS3.obtainElementAxialForce();
        
        % initialise a vector of each element's increment in the load
        % proportionality factor in order to find the smallest one:
        element_d_lambdas = zeros(TRUSS3.num_of_elements,1);
        % initialise a vector of the failure modes of all elements if they were
        % to fail:
        failure_modes = zeros(TRUSS3.num_of_elements,1);
        % looping through all elements to find the lowest increment in the load
        % proportionality factor:
        for el = 1:TRUSS3.num_of_elements
            % if the element is going to fail in tension.
            % find the elements increment in the load proportionality
            % factor:
            element_d_lambdas(el) = ...
                (yield_load - TRUSS3.total_element_axial_forces(el))/TRUSS3.element_axial_forces(el);
            % take note note of failure mode (compression):
            failure_modes(el) = 'T';
            if element_d_lambdas(el) < 0
                % define the length of the element.
                % obtain global node numers of local nodes:
                n1 = TRUSS3.element_connectivity(el,1);
                n2 = TRUSS3.element_connectivity(el,2);
                % obtaining original x and y coordinate of local node 1:
                x1_original = TRUSS3.coords(n1,1);
                y1_original = TRUSS3.coords(n1,2);
                % obtaining original x and y coordinate of local node 2:
                x2_original = TRUSS3.coords(n2,1);
                y2_original = TRUSS3.coords(n2,2);
                % obtaining original length of unstretched element:
                el_L = sqrt((x2_original-x1_original)^2 + (y2_original-y1_original)^2);
                % define the euler buckling load for each element:
                buckling_load = -(pi^2)*(TRUSS3.EI/(el_L^2));
                
                % if the element is going to fail in compression.
                % find the elements increment in the load proportionality
                % factor only taking into account plastic squash:
                element_d_lambdas(el) = ...
                    (-yield_load - TRUSS3.total_element_axial_forces(el))/TRUSS3.element_axial_forces(el);
                % take note note of failure mode (compression):
                failure_modes(el) = 'C';
                
                % find the increment in load proportionality factor now taking
                % into account buckling:
                element_d_lambda_buckling = ...
                    (buckling_load - TRUSS3.total_element_axial_forces(el))/TRUSS3.element_axial_forces(el);
                
                % if the increment in the load proportionality factor with
                % buckling is less than that with plastic squash, use the
                % increment obtained with buckling:
                if element_d_lambda_buckling < element_d_lambdas(el)
                    element_d_lambdas(el) = element_d_lambda_buckling;
                    failure_modes(el) = 'B';
                end
            end
        end
        
        % find the minimum increment in the load proportionality factor, and
        % which critical element it belongs to:
        [d_lambda, critical_element] = min(element_d_lambdas);
        % make not of the load proportionality factor increment:
        TRUSS3.d_lambda_list_GMNA(end+1,1) = d_lambda;
        % making note of the failure mode of the critical element:
        failure_mode = failure_modes(critical_element);
        % add the critial element to a list of critical elements that have
        % failed in each iteration as well as their failure mode:
        TRUSS3.failed_elements_GMNA(end+1,1) = TRUSS3.intact_elements_GMNA(critical_element);
        TRUSS3.failure_modes_GMNA(end+1,1) = char(failure_mode);
        % add the increment in the load proportionality factor to the total
        % proportionality factor:
        TRUSS3.total_load_prop_factor_GMNA = TRUSS3.total_load_prop_factor_GMNA + d_lambda;
        
        % apply the scaled load proportionality factor to find the axial
        % forces in all non failed elements.
        % construct global nodal force vector from scaled loads:
        TRUSS3 = TRUSS3.constructGlobalNodalForceVector(d_lambda);
        % solve matrix system to find nodal displacements and reaction forces:
        TRUSS3 = TRUSS3.solveFiniteElementMatrixSystem();
        % find the new deformed nodal coordinates arising from the nodal
        % displacements:
        TRUSS3 = TRUSS3.obtainNewNodalCoordinates();
        % Get increment in axial forces within each element from the scaled
        % load proportionality factor:
        TRUSS3 = TRUSS3.obtainElementAxialForce();
        % add the correctly scaled loads whereby the critical member is made
        % to fail:
        TRUSS3.total_element_axial_forces = TRUSS3.total_element_axial_forces + TRUSS3.element_axial_forces;
        % add the increment in the nodal displacements to a vector of total
        % displacements:
        TRUSS3.total_nodal_displacements = TRUSS3.total_nodal_displacements + TRUSS3.U;
        
        % replace the critical member, which has yielded, with a pair of nodal
        % forces for the next iteration, and so have four seperate forces on
        % four different dofs, depending on the orientation of the failed bar
        % member.
        % obtain the critical element's force as the axial force at failure:
        critical_element_force = TRUSS3.total_element_axial_forces(critical_element);
        % replace the failed element with equivilent nodal forces in failure:
        TRUSS3 = TRUSS3.replaceFailedElementWithForces(critical_element,critical_element_force);
        
        % get rid of the yielded element from the element connectivity
        % vector:
        TRUSS3.element_connectivity(critical_element,:) = [];
        % subtract the failed element also from the list of intact elements:
        TRUSS3.intact_elements_GMNA(critical_element) = [];
        % update the number of elements remaining:
        TRUSS3.num_of_elements = size(TRUSS3.element_connectivity,1);
        % get rid of the failed element from the vector of axial element
        % forces:
        TRUSS3.total_element_axial_forces(critical_element) = [];
        
        % work out number of redundancies and them to a list:
        redundancies = TRUSS3.num_of_elements + 8 -(2*TRUSS3.num_of_nodes);
        TRUSS3.redundancies_list_GMNA(end+1,1) = redundancies;
    end
    
    % log the total lambda:
    lambda_totals_GMNA(analysis) = TRUSS3.total_load_prop_factor_GMNA;
    % log the number of failed elements in this particular analysis:
    num_of_failed_elements_GMNA(analysis) = length(TRUSS3.failed_elements_GMNA);
end

% plot the total load proportionality factors against the radius:
figure();
plot(radii,lambda_totals_MNA,'-k');
hold on;
plot(radii,lambda_totals_GMNA,'--k');
title('Total Load Proportionality Factor Dependent on Bar Element Radius');
legend('MNA','GMNA');
xlabel('Radius [m]');
ylabel('Total Load Proportionality Factor');

%% Question 6 - Range of R which is at risk

% plot the number of failed number of elements against the radius:
figure();
plot(radii,num_of_failed_elements_MNA,'-k');
hold on;
plot(radii,num_of_failed_elements_GMNA,'--k');
title('Number of Failed Elements Until Complete Failure for Varying Radii');
legend('MNA','GMNA');
xlabel('Radius [m]');
ylabel('Number of Failed Elements');
