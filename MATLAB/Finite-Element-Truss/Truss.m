%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collapse of a redundant truss by incremental elastic analysis  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                %%%
%%% Truss class which holds all the relevant FE-related algorithms %%%
%%%                                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INCLUDE DEPENDENCY ADVICE OF METHODS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

classdef Truss
    % define properties which a Truss object would hold.
    properties
        % define properties which are crucial to object construction.
        % the elastic Young's modulus [N/m2]:
        E
        % the second moment of area [m4]:
        I
        % the radius of each bar element [m]:
        R
        % the area of each bar element:
        A
        % the elastic young's modulus multiplied by second moment of area
        % (since both are constant, and can make things simpler):
        EI
        % the elastic youngs modulus multiplied by area (since both are
        % constant):
        EA
        % the yield stress [N/m2]:
        sigma_y
        % the horizontal length unit of each truss brace [m]:
        L
        % the element nodal connectivity
        element_connectivity
        % the number of elements:
        num_of_elements
        % the coordinates of each node:
        coords
        % the dofs of each node:
        dofs
        % the free degrees of freedom (dofs):
        dofs_free
        % the restrained degrees of freedom (dofs):
        dofs_restrained
        % the number of nodes:
        num_of_nodes
        % the unit point load [N]:
        P
        % global counter used for incremenetal MNA analysis:
        MNA_counter
        % global counter used for incremenetal GMNA analysis:
        GMNA_counter
        % global stiffness matrix
        K
        % applied point loads, and their positions in terms of the dofs:
        applied_load_distribution
        % global nodal force vector:
        F
        % flag for solvability of finite element matrix system:
        solvable
        % global nodal displacement vector
        U
        % new coordinates of nodes after an iteration of loading
        % deformation:
        new_coords
        % new amplified coordinates of nodes after an iteration of loading
        % deformation (for plotting purposes only):
        amp_coords
        % original nodal coordinates:
        original_coords
        % vector of axial forces in elements:
        element_axial_forces
        % the final load proportionality factor based on the MNA faiure
        % criterion:
        total_load_prop_factor_MNA
        % a vector of the total displacements of each node:
        total_nodal_displacements
        % the final load proportionality factor based on the GMNA faiure
        % criterion:
        total_load_prop_factor_GMNA
        % vector of total axial forces in each members after the Ith
        % iteration/increment:
        total_element_axial_forces
        % a vector of failed elements with MNA failure criterion:
        failed_elements_MNA
        % a vector of failed elements with GMNA failure criterion:
        failed_elements_GMNA
        % a vector that holds the failure modes of each failed element in
        % MNA:
        failure_modes_MNA
        % a vector that holds the failure modes of each failed element in
        % GMNA:
        failure_modes_GMNA
        % a vector of a list of all the intact (non failed) elements
        % throughout and after MNA analysis
        intact_elements_MNA
        % a vector of a list of all the intact (non failed) elements
        % throughout and after MNA analysis
        intact_elements_GMNA
        % a list of all load proportionality increments for MNA:
        d_lambda_list_MNA
        % a list of all load proportionality increments for GMNA:
        d_lambda_list_GMNA
        % a list of the redundancies in each iteration in MNA:
        redundancies_list_MNA
        % a list of the redundancies in each iteration in GMNA:
        redundancies_list_GMNA
    end
    methods
        function obj = Truss(NODES,ELEMENTS,LOADING)
            % set object elastic young's modulus to what is given:
            obj.E = ELEMENTS.E;
            % set object element bar radius to what is given:
            obj.R = ELEMENTS.R;
            % set object element bar cross sectional area (circular) to
            % what is given:
            obj.A = pi*obj.R^2;
            % define the second moment of area for each of the elements:
            obj.I = 0.25*pi*obj.R^4;
            % obtain the elastic modulus mutliplied by the second moment of
            % area for simplicity:
            obj.EI = obj.E*obj.I;
            % set object elastic young's modulus multiplied by area (since
            % both constant) to what is given:
            obj.EA = obj.E*obj.A;
            % set object yield stress to what is given:
            obj.sigma_y = ELEMENTS.sigma_y;
            % set object horizontal and vertical brace length to what is
            % given:
            obj.L = ELEMENTS.L;
            % set object element nodal connectivity to what is given:
            obj.element_connectivity = ELEMENTS.element_connectivity;
            % set object number of elements to what is given:
            obj.num_of_elements = size(obj.element_connectivity,1);
            
            % set object nodal coordinates to what are given:
            obj.coords = NODES.coords;
            % set the original nodal coordates of the unloaded structure:
            obj.original_coords = obj.coords;
            % set object to what is given:
            obj.dofs = NODES.dofs;
            % set object free dofs to what are given:
            obj.dofs_free = NODES.dofs_free;
            % set object restrained dofs to what are given:
            obj.dofs_restrained = NODES.dofs_restrained;
            % set object number of nodes to what is given:
            obj.num_of_nodes = size(obj.coords,1);
            
            % set object applied point load unit to what is given:
            obj.P = LOADING.P;
            % set object applied point loads, and their positions in terms
            % of the dofs to what is given:
            obj.applied_load_distribution = LOADING.applied_load_distribution;
            
            % initialise a vector of a list of the failed elements in MNA:
            obj.failed_elements_MNA = [];
            % initialise a vector of a list of the failed elements in GMNA:
            obj.failed_elements_GMNA = [];
        end
        function obj = constructGlobalStiffnessMatrix(obj)
            % this method constructs the global stiffness matrix K.
            % initialise a square global stiffness matrix K with zeros,
            % which row and column count is the same as number of dofs (2
            % times number of nodes):
            obj.K = zeros(2*obj.num_of_nodes);
            % loop though each element to build a global stiffness matrix:
            for el = 1:obj.num_of_elements
                % identify the element node numbers:
                n1 = obj.element_connectivity(el,1);
                n2 = obj.element_connectivity(el,2);
                % x and y coordinate of element node 1:
                x1 = obj.coords(n1,1);
                y1 = obj.coords(n1,2);
                % x and y coordinate of element node 2:
                x2 = obj.coords(n2,1);
                y2 = obj.coords(n2,2);
                % element node 1 dofs:
                n1_dof1 = obj.dofs(n1,1);
                n1_dof2 = obj.dofs(n1,2);
                % element node 2 dofs:
                n2_dof1 = obj.dofs(n2,1);
                n2_dof2 = obj.dofs(n2,2);
                % angle of inclination of element relative to the POSITIVE
                % direction of the x axis:
                alpha = atan2(y2-y1,x2-x1);
                % angle parameters:
                c = cos(alpha);
                c2 = c*c;
                s = sin(alpha);
                s2 = s*s;
                cs = c*s;
                % element length:
                el_L = sqrt((x2-x1)^2 + (y2-y1)^2);
                % element axial stiffness:
                el_stiffness = obj.EA/el_L;
                
                % updating global stiffness matrix [K] coefficients.
                
                % row 1 (element n1_dof1).
                % column 1 (element n1_dof1):
                obj.K(n1_dof1,n1_dof1) = obj.K(n1_dof1,n1_dof1) + el_stiffness*c2;
                % column 2 (element n1_dof2):
                obj.K(n1_dof1,n1_dof2) = obj.K(n1_dof1,n1_dof2) + el_stiffness*cs;
                % column 3 (element n2_dof1):
                obj.K(n1_dof1,n2_dof1) = obj.K(n1_dof1,n2_dof1) - el_stiffness*c2;
                % column 4 (element n2_dof2):
                obj.K(n1_dof1,n2_dof2) = obj.K(n1_dof1,n2_dof2) - el_stiffness*cs;
                
                % row 2 (element n1_dof2).
                % column 1 (element n1_dof1):
                obj.K(n1_dof2,n1_dof1) = obj.K(n1_dof2,n1_dof1) + el_stiffness*cs;
                % column 2 (element n1_dof2):
                obj.K(n1_dof2,n1_dof2) = obj.K(n1_dof2,n1_dof2) + el_stiffness*s2;
                % column 3 (element n2_dof1):
                obj.K(n1_dof2,n2_dof1) = obj.K(n1_dof2,n2_dof1) - el_stiffness*cs;
                % column 4 (element n2_dof2):
                obj.K(n1_dof2,n2_dof2) = obj.K(n1_dof2,n2_dof2) - el_stiffness*s2;
                
                % row 3 (element n2_dof1).
                % column 1 (element n1_dof1):
                obj.K(n2_dof1,n1_dof1) = obj.K(n2_dof1,n1_dof1) - el_stiffness*c2;
                % column 2 (element n1_dof2):
                obj.K(n2_dof1,n1_dof2) = obj.K(n2_dof1,n1_dof2) - el_stiffness*cs;
                % column 3 (element n2_dof1):
                obj.K(n2_dof1,n2_dof1) = obj.K(n2_dof1,n2_dof1) + el_stiffness*c2;
                % column 4 (element n2_dof2):
                obj.K(n2_dof1,n2_dof2) = obj.K(n2_dof1,n2_dof2) + el_stiffness*cs;
                
                % row 4 (element n2_dof2).
                % column 1 (element n1_dof1):
                obj.K(n2_dof2,n1_dof1) = obj.K(n2_dof2,n1_dof1) - el_stiffness*cs;
                % column 2 (element n1_dof2):
                obj.K(n2_dof2,n1_dof2) = obj.K(n2_dof2,n1_dof2) - el_stiffness*s2;
                % column 3 (element n2_dof1):
                obj.K(n2_dof2,n2_dof1) = obj.K(n2_dof2,n2_dof1) + el_stiffness*cs;
                % column 4 (element n2_dof2):
                obj.K(n2_dof2,n2_dof2) = obj.K(n2_dof2,n2_dof2) + el_stiffness*s2;
            end
        end
        function obj = constructGlobalNodalForceVector(obj,load_prop_factor)
            % initialise the global nodal force vector, with length same
            % as the number of dofs, with zeros:
            obj.F = zeros(2*obj.num_of_nodes,1);
            % inserting the point loads into the global nodal force vector:
            for i = 1:size(obj.applied_load_distribution,1)
                % insert the force into the correct dof index, and multiply
                % by the load proportionality factor to each applied load,
                % as well as adding what might already be applied to the
                % node:
                obj.F(obj.applied_load_distribution(i,2)) = ...
                    obj.F(obj.applied_load_distribution(i,2)) + load_prop_factor*obj.applied_load_distribution(i,1);
                
            end
        end
        function obj = checkSolvability(obj)
            % first, clarify whether the matrix of free dofs in K can be
            % solved using rcond whereby if rcond returns a value below
            % machine precision, the matrix cannot be solved.
            if rcond(obj.K(obj.dofs_free,obj.dofs_free)) > 1e-15
                % if rcond returns a value larger than machine precision,
                % the system can be solved:
                obj.solvable = true(1);
            else
                % if rcond returns a value that is smaller than machine
                % precision, then the system cannot be solved:
                obj.solvable = false(1);
            end
        end
        function obj = solveFiniteElementMatrixSystem(obj)
            % define a global nodal displacement vector with zeros:
            obj.U = zeros(2*obj.num_of_nodes,1);
            
            % solve for the unknown nodal dofs.
            % uF = KFF\(fF - KFR*uR):
            obj.U(obj.dofs_free) = ...
                obj.K(obj.dofs_free,obj.dofs_free)\(obj.F(obj.dofs_free) - (obj.K(obj.dofs_free,obj.dofs_restrained)*obj.U(obj.dofs_restrained)));
            
            % solve for the unknown reactions.
            % fR = KRF*uF + KRR*uR:
            obj.F(obj.dofs_restrained) = ...
                (obj.K(obj.dofs_restrained,obj.dofs_free)*obj.U(obj.dofs_free)) + (obj.K(obj.dofs_restrained,obj.dofs_restrained)*obj.U(obj.dofs_restrained));
        end
        function obj = obtainNewNodalCoordinates(obj)
            % vector of normal nodal coordinates after deformation:
            obj.new_coords = zeros(size(obj.coords));
            for I = 1:size(obj.coords,1)
                for J = 1:size(obj.coords,2)
                    % obtaining new non amplified nodal coordinates by
                    % adding the displacements:
                    obj.new_coords(I,J) = obj.coords(I,J) + obj.U(obj.dofs(I,J));
                end
            end
        end
        function obj = obtainElementAxialForce(obj)
            % define a vector of element axial forces:
            obj.element_axial_forces = zeros(obj.num_of_elements,1);
            
            % loop though each element to find its axial force:
            for el = 1:obj.num_of_elements
                % obtain global node numers of local nodes:
                n1 = obj.element_connectivity(el,1);
                n2 = obj.element_connectivity(el,2);
                
                % obtaining new x and y coordinate of local node 1:
                x1_new = obj.new_coords(n1,1);
                y1_new = obj.new_coords(n1,2);
                % obtaining new x and y coordinate of local node 2:
                x2_new = obj.new_coords(n2,1);
                y2_new = obj.new_coords(n2,2);
                
                % obtaining original x and y coordinate of local node 1:
                x1_original = obj.coords(n1,1);
                y1_original = obj.coords(n1,2);
                % obtaining original x and y coordinate of local node 2:
                x2_original = obj.coords(n2,1);
                y2_original = obj.coords(n2,2);
                
                % obtaining original length of unstretched element:
                el_L_original = sqrt((x2_original-x1_original)^2 + (y2_original-y1_original)^2);
                % define element axial stiffness:
                el_stiffness = obj.EA/el_L_original;
                
                % displacement of nodes in the direction of which the bar
                % element can do work against:
                alpha = atan2(y2_original-y1_original,x2_original-x1_original);
                u1 = x1_new - x1_original;
                v1 = y1_new - y1_original;
                u2 = x2_new - x2_original;
                v2 = y2_new - y2_original;
                up1 = cos(alpha)*u1 + sin(alpha)*v1;
                up2 = cos(alpha)*u2 + sin(alpha)*v2;
                
                % difference in displacements of nodes:
                dup = up2 - up1;
                
                % obtain axial force in element:
                obj.element_axial_forces(el) = el_stiffness*dup;
            end
        end
        function obj = replaceFailedElementWithForces(obj,critical_element,critical_element_force)
            % replace the critical member, which has yielded, with a pair of nodal
            % forces for the next iteration, and so have four seperate forces on
            % four different dofs, depending on the orientation of the failed bar
            % member.
            % obtain the nodes at each end of the bar element:
            n1 = obj.element_connectivity(critical_element,1);
            n2 = obj.element_connectivity(critical_element,2);
            % x and y coordinate of element node 1:
            x1 = obj.coords(n1,1);
            y1 = obj.coords(n1,2);
            % x and y coordinate of element node 2:
            x2 = obj.coords(n2,1);
            y2 = obj.coords(n2,2);
            % element node 1 dofs:
            n1_dof1 = obj.dofs(n1,1);
            n1_dof2 = obj.dofs(n1,2);
            % element node 2 dofs:
            n2_dof1 = obj.dofs(n2,1);
            n2_dof2 = obj.dofs(n2,2);
            % angle of inclination of element relative to the POSITIVE
            % direction of the x axis:
            alpha = atan2(y2-y1,x2-x1);
            if (0 <= alpha) && (alpha <= 0.5*pi)
                % if angle is in quadrant 1.
                if critical_element_force > 0
                    % if member is in tension
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                else
                    % if member is in compression
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                end
            elseif (0.5*pi < alpha) && (alpha <= pi)
                % if angle is in quadrant 2:
                if critical_element_force > 0
                    % if member is in tension
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                else
                    % if member is in compression
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                end
            elseif (-pi < alpha) && (alpha < -0.5*pi)
                % if angle is in quadrant 3:
                if critical_element_force > 0
                    % if member is in tension
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                else
                    % if member is in compression
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                end
            else
                % if angle is in quadrant 4:
                if critical_element_force > 0
                    % if member is in tension
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                else
                    % if member is in compression
                    % define the force of the element on the first node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof1;
                    % define the force of the element on the first node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n1_dof2;
                    % define the force of the element on the second node's horizontal
                    % dof:
                    obj.applied_load_distribution(end+1,1) = -abs(critical_element_force*cos(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof1;
                    % define the force of the element on the second node's vertical
                    % dof:
                    obj.applied_load_distribution(end+1,1) = abs(critical_element_force*sin(alpha));
                    obj.applied_load_distribution(end,2) = n2_dof2;
                end
            end
        end
        
    end
end