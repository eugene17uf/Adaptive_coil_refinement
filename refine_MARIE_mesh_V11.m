function refine_MARIE_mesh_V11(filename,filename2,iteration,faces_to_refine,faces_not_to_refine,vertices_STL,faces_STL)

%Inputs:
% filename: Name of starting mesh for adaptive mesh refinement (AMR)
% filename2: Name of output mesh/txt/smm file; the iteration number is appended to the file
% iteration: iteration number of adaptive refinement
% faces_to_refine: the faces (triangles) that will be refined in the mesh
% faces_not_to_refine: the faces not to be refined; note that some faces will eventually be refined from this list to eliminate any hanging nodes
% vertices_STL: vertices (nodes) from super refined mesh that represents STL model of object
% faces_STL: faces (triangles) from super refined mesh that represents STL model of object

% Outputs:
% .msh, .txt, .smm files of refined mesh for running in MARIE

% Import the coil with MARIE and have it preprocess to generate COIL
% structure, which the coordinates of the nodes and the triangles
% making up mesh
[COIL] = Import_COIL([filename,'.smm']);

% Get 3D coordinates of all the nodes from COIL.node
FV.vertices = COIL.node';

% Get all the triangles (faces of mesh) from COIL.elem; each triangle is
% given with the indices of the three nodes that make up that triangle
% (first three columns of FV.faces'
FV.faces = COIL.elem';

% Use the mesh_parse.m script from MARIE to obtain nodes and IDs that make
% up lines, which are used by MARIE to get port locations. Ports can be
% made up of multiple lines; these lines all have the same IDs (3rd row of
% "lines" variable)
[nodes,~,lines,~,~] = Mesh_Parse([filename,'.msh']);


% Organize the ports by the line(s) in the .msh file that define the port;
% The third row of the variable "lines" indicate groupings of the line(s)
% that make up a port; we want to first group all the nodes by the port
% they define; then find the endpoints of the port by finding the max
% distance within each group of nodes; we will use these two endpoints
% later to find if newly created nodes lie on that line
port_endpoints = [];
for ii = min(lines(3,:)):max(lines(3,:)) % index by 3rd row of "lines" variable
    
    port_indices_temp = lines(1:2,find(lines(3,:) == ii));
    
    port_indices = unique(port_indices_temp(:));
    
    port_nodes = nodes(:,port_indices)'; % nodes corresponding to the points on the line defining port
    
    % get distances of all points on line to one another; longest distances
    % indicate endpoints of line; we will use these to define the equation of
    % the line later to find which nodes fall on these lines
    [distances, index] = pdist2(port_nodes,port_nodes, 'euclidean', 'largest', 1);
    % Find max distance
    [maxDistance, indexOfMax] = max(distances);
    index2 = indexOfMax;
    index1 = index(indexOfMax);
    
    % All port endpoints, grouped every two rows
    port_endpoints = [port_endpoints; port_nodes(index2,:); port_nodes(index1,:)];
    
end


% Faces we want to refine and faces that are not refined to place back in
% our mesh after initial refinement
faces_mesh = faces_to_refine';
faces_remaining = faces_not_to_refine';

% Generate files needed to run triangulation_triangle_neighbors.m script to
% get the neighbors of all the faces (triangles)
fileID = fopen('SKYRA_nodes.txt','w');
% fprintf(fileID,'%6s\n','#Nodes');

    for i = 1:size(FV.vertices,1)
       
        fprintf(fileID,'%12.16f %12.16f %12.16f\n',FV.vertices(i,1),FV.vertices(i,2),FV.vertices(i,3));
        
    end
fclose(fileID);
fileID = fopen('SKYRA_elements.txt','w'); 
% fprintf(fileID,'%6s\n','#Elements');

    for i = 1:size(FV.faces,1)
       
        fprintf(fileID,'%d %d %d\n',FV.faces(i,1),FV.faces(i,2),FV.faces(i,3));
        
    end
fclose(fileID);

triangulation_triangle_neighbors('SKYRA');

% Get matrix of neighboring triangles to each triangle (face)
triangle_neighbors = importdata('SKYRA_element_neighbors.txt');

triangle_neighbors_to_refine = triangle_neighbors(faces_mesh,:); %neighbor elements of elements to be refined that need to be either split or refined via triangulation
triangle_neighbors_to_refine_reshaped = triangle_neighbors_to_refine(:); %reshape into vector
triangle_neighbors_to_refine_reshaped(triangle_neighbors_to_refine_reshaped == -1) = []; %remove any -1, which indicate no neighbor to triangle (out of 3 potential neighbors)
triangles_remaining = triangle_neighbors_to_refine_reshaped(~ismember(triangle_neighbors_to_refine_reshaped,faces_mesh)); %remove elements from neighbor list that were going to be refined anyways, leaving only elements that need to be split or possibly refined via triangulation
[unique_triangles, ia] = unique(triangles_remaining); % find all unique elements that need to be split/refined; the ones to be refined are the ones that appear more than once on list

% Get triangles that are repeated from triangles_remaining; these need to
% be refined by triangulation; need to be added to the faces_mesh and
% removed from faces_remaining

% indices of duplicate faces
duplicate_ind = setdiff(1:size(triangles_remaining, 1), ia);
% duplicate faces 
duplicate_value = unique(triangles_remaining(duplicate_ind));

faces_mesh_new = [faces_mesh duplicate_value']; % add the neighbor triangles that need to be refined by triangulation
faces_remaining_new = setdiff(faces_remaining,duplicate_value); % remove the neighbor triangles that need to be refined by triangulation
triangles_to_split = setdiff(unique_triangles,duplicate_value); % remove the neighbor triangles that need to be refined by triangulation, leaving only ones that need to be split

% Need to find the neighbors to split of just added triangles to
% faces_mesh_new
new_neighbors = triangle_neighbors(duplicate_value,:);
new_neighbors_2 = new_neighbors(:);
new_neighbors_2(new_neighbors_2 == -1) = [];
new_neighbors_3 = new_neighbors_2(~ismember(new_neighbors_2,faces_mesh_new));

% Need to find if new neighbors already appeared in triangles_to_split, if
% so, they need to be refined by triangulation. Should do this in while
% statement, exiting when there are no more repeating values in
% triangles_to_split

temp2 = [triangles_to_split;new_neighbors_3];
[temp3, ia] = unique(temp2);
% indices of duplicate faces
duplicate_ind = setdiff(1:size(temp2, 1), ia);
% duplicate faces
dup = unique(temp2(duplicate_ind));

while ~isempty(dup) 
    faces_mesh_new = [faces_mesh_new dup']; % add the neighbor triangles that need to be refined by triangulation
    faces_remaining_new = setdiff(faces_remaining_new,dup); % remove the neighbor triangles that need to be refined by triangulation
    triangles_to_split = setdiff(temp3,dup); % remove the neighbor triangles that need to be refined by triangulation, leaving only ones that need to be split
    
    new_neighbors = triangle_neighbors(dup,:);
    new_neighbors_2 = new_neighbors(:);
    new_neighbors_2(new_neighbors_2 == -1) = [];
    new_neighbors_3 = new_neighbors_2(~ismember(new_neighbors_2,faces_mesh_new));

    temp2 = [triangles_to_split;new_neighbors_3];
    [temp3, ia] = unique(temp2);
    % duplicate indices
    duplicate_ind = setdiff(1:size(temp2, 1), ia);
    % duplicate values
    dup = unique(temp2(duplicate_ind));
    triangles_to_split = unique([triangles_to_split;new_neighbors_3]);
end

triangles_to_split = temp3;
FV.faces = (COIL.elem(1:3,faces_mesh_new))';

% generate new .txt files for triangulation_refine_eugene_mod.m script
fileID = fopen('SKYRA_nodes.txt','w');
fprintf(fileID,'%6s\n','#Nodes');

    for i = 1:size(FV.vertices,1)
       
        fprintf(fileID,'%12.16f %12.16f %12.16f\n',FV.vertices(i,1),FV.vertices(i,2),FV.vertices(i,3));
        
    end
fclose(fileID);
fileID = fopen('SKYRA_elements.txt','w'); 
fprintf(fileID,'%6s\n','#Elements');

    for i = 1:size(FV.faces,1)
       
        fprintf(fileID,'%d %d %d\n',FV.faces(i,1),FV.faces(i,2),FV.faces(i,3));
        
    end
fclose(fileID);

% Get refined node and element (triangle) list
[node_xy2, triangle_node2] = triangulation_refine_eugene_mod('SKYRA');

FV.vertices = node_xy2';
FV.faces = triangle_node2';

% Put back elements not refined into our mesh
FV.faces((size(triangle_node2,2)+1):(size(triangle_node2,2))+size(faces_remaining_new,2),1:3) = (COIL.elem(1:3,faces_remaining_new))';

% Split triangles to eliminate any hanging nodes
if ~isempty(triangles_to_split)
    [element_node2] = split_triangles_V5(triangles_to_split,FV.vertices,FV.faces,COIL.elem(1:3,:)');
    FV.faces = element_node2';
else
    element_node2 = FV.faces';
end


% Snap back elements onto pre-refined mesh; use a really refined mesh as
% the surface to snap back to
XV.vertices = vertices_STL'; %vertices of pre-refined mesh
XV.faces = faces_STL(1:3,:)'; %faces of pre-refined mesh
[ distances, surface_points] = point2trimesh( XV,'QueryPoints',FV.vertices,'Algorithm','parallel' );

% Surface_points are snapped back points
FV.vertices = surface_points;

node_xy2 = surface_points';

% Need to find all the nodes that are along the lines corresponding to
% ports; this is why we got endpoints from each port previously
% (cc-aa)*null(bb-aa) should equal 0 if point cc is on the line defined by
% endpoints aa and bb; aa and bb were found earlier 

new_line_nodes_temp = [];

for ii = 1:size(node_xy2,2)
    
    for jj = 1:2:size(port_endpoints,1)
        
        if max(abs((node_xy2(:,ii)' - port_endpoints(jj,:)) * null(port_endpoints(jj+1,:) - port_endpoints(jj,:)))) <= 1e-9 && (pdist2(node_xy2(:,ii)',port_endpoints(jj,:)) <= pdist2(port_endpoints(jj+1,:),port_endpoints(jj,:))) && (pdist2(node_xy2(:,ii)',port_endpoints(jj+1,:)) <= pdist2(port_endpoints(jj+1,:),port_endpoints(jj,:)))
           
            temp(1,1) = jj; %create an index for grouping nodes for each port; will have same value here
            temp(2:4,1) = node_xy2(:,ii);
            temp(5,1) = ii; %store the node number to write out in .msh file
            new_line_nodes_temp = [new_line_nodes_temp temp];
            clear temp
            
        end
    end
end

% Sort the nodes together by port
new_line_nodes = sortrows(new_line_nodes_temp',1);

% Calculating the number of lines present in our refined mesh; need for
% fprintf statement that prints how many elements are present
% later when we are writing out new mesh file; number of lines per port is
% number of nodes - 1
num_lines = 0;

for ii = 1:2:size(port_endpoints,1)
    
    lines_per_port = size(find(new_line_nodes(:,1) == ii),1);
    
    num_lines = num_lines + (lines_per_port - 1);
    
end    
    
        
% Write new .msh file with refined triangles; format setup for reading in MARIE
fileID = fopen([filename2,'_',num2str(iteration),'.msh'],'w');
fprintf(fileID,'%6s\n','$MeshFormat');
fprintf(fileID,'%6s\n','2.2 0 8');
fprintf(fileID,'%6s\n','$EndMeshFormat');
fprintf(fileID,'%6s\n','$Nodes');
fprintf(fileID,'%d\n',size(FV.vertices,1));
for i = 1:size(FV.vertices,1)
    fprintf(fileID,'%d %12.16f %12.16f %12.16f\n',i,FV.vertices(i,1),FV.vertices(i,2),FV.vertices(i,3));
end
fprintf(fileID,'%6s\n','$EndNodes');
fprintf(fileID,'%6s\n','$Elements');
fprintf(fileID,'%d\n',size(FV.faces,1)+num_lines);

ll = 1; % counter for element (first column in .msh file in elements section)
port_num = 1; % counter to track of port number since all lines corresponding to a port need to have same ID (fourth column in elements section of .msh fle)

% Write out lines first
for ii = 1:2:size(port_endpoints,1)
    
     % indices of nodes for a particular port
    yy = new_line_nodes(find(new_line_nodes(:,1) == ii),2:4);
    
    % find distance between one endpoint of port and all other points
    % making up port
    [D] = pdist2(yy(1,:),yy,'euclidean', 'largest', 1);
    
    % Append distances to nodes matrix
    port_nodes_distances = [new_line_nodes(find(new_line_nodes(:,1) == ii),2:5) D'];
    
    % Sort port_nodes_distances by distance
    port_nodes_distances_ordered = sortrows(port_nodes_distances,5);
    
    % Print lines
    for jj = 1:(size(port_nodes_distances_ordered,1)-1)
        fprintf(fileID,'%d %d %d %d %d %d %d\n',ll,1,2,1000+port_num,1800+port_num,port_nodes_distances_ordered(jj,4),port_nodes_distances_ordered(jj+1,4)); %Lines corresponding to a port need to have same ID (1000+ii term)
        ll = ll+1;
    end
    
    port_num = port_num + 1; 
end

% Write out faces (triangles)
for i = 1:size(FV.faces,1)

    fprintf(fileID,'%d %d %d %d %d %d %d %d\n',i+ll-1,2,2,1000,1800,FV.faces(i,1),FV.faces(i,2),FV.faces(i,3));

end
fprintf(fileID,'%6s\n','$EndElements');
fclose(fileID);

% Copy .txt file with port information for MR_solver.m
copyfile([filename,'.txt'],[filename2,'_',num2str(iteration),'.txt']);

% Make new .smm file for loading next iteration of mesh
fileID = fopen([filename2,'_',num2str(iteration),'.smm'],'w');
fprintf(fileID,'%6s\n',['SKYRA with shield iteration ',num2str(iteration)]);
fprintf(fileID,'%12.16f %12.16f \n',1.2579e-8,0.01);
fprintf(fileID,'%6s\n',[filename2,'_',num2str(iteration),'.msh']);
fclose(fileID);

end

