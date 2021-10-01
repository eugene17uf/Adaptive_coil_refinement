function [faces_after_split] = split_triangles_V5(triangles_to_split,new_nodes,new_faces,old_faces_list)

% Inputs:

% triangles_to_split: list of neighboring triangles to be split into two
% triangles using hanging nodes; the values in the vector correspond to the
% faces listed in old_faces_list

% new_nodes: list of nodes after refinement by triangulation

% new_faces: list of faces after refinement by trianglulation

% old_faces_list: list of faces of original mesh that needed to be refined;
% will be used to get the nodal coordinates of the the triangles_to_split
% since the refinement by triangulation changed faces list


% Outputs:

% faces_after_split: list of faces with newly split triangles

%Split the neighboring triangles 
num_faces = size(new_faces,1); % Number of faces in mesh with only refinement by triangulation

% faces_after_split is the matrix containing the faces from the refinement by triangulation and will have new faces added after neighboring triangles are split
faces_after_split = new_faces'; 


for ii = 1:size(triangles_to_split,1)
    % Get nodal coordinates of the triangle to split; use the
    % old_faces_list to get the correct indices to use within new_nodes
     node1 = new_nodes(old_faces_list(triangles_to_split(ii),1),:);     
     node2 = new_nodes(old_faces_list(triangles_to_split(ii),2),:);  
     node3 = new_nodes(old_faces_list(triangles_to_split(ii),3),:); 
     
     % Find the midpoint of the edges of the triangle using the
     % aforementione nodes; one of these nodes has already been created by
     % the refinement by triangulation; the one that exists will serve as
     % the endpoint for bisecting the triangle
     middle_node1 = (node1 + node2)/2.0;
     middle_node2 = (node1 + node3)/2.0;
     middle_node3 = (node2 + node3)/2.0;
     
     % Find which of the midpoints (middle_node) already exists in
     % new_nodes; only one of the following three if statements is true
     
     if ismember(middle_node1,new_nodes,'rows')
         % use ismember function to get the index of that node (Loca);
        [Lia, Loca] = ismember(middle_node1,new_nodes,'rows');
        
        % Create first triangle; find the index of the triangle that needs to be
        % split (which column has the indices of the nodes of this
        % triangle) and place the triangle there since we don't want the
        % triangle to be split and still appear as a face in our final mesh
        % nodes to be oriented counter-clockwise to unify mesh normals
        faces_after_split(1:3,find(faces_after_split(1,:) == old_faces_list(triangles_to_split(ii),1) & faces_after_split(2,:) == old_faces_list(triangles_to_split(ii),2) & faces_after_split(3,:) == old_faces_list(triangles_to_split(ii),3))) = [  Loca; old_faces_list(triangles_to_split(ii),3); old_faces_list(triangles_to_split(ii),1) ];
        
        % Create second triangle; will need new index in the
        % faces_after_split list
        faces_after_split(1:3,num_faces+1) = [  Loca; old_faces_list(triangles_to_split(ii),2); old_faces_list(triangles_to_split(ii),3) ];
        
        % Keep counter of new indices
        num_faces = num_faces + 1;
        
     end
     
     if ismember(middle_node2,new_nodes,'rows') 
         
        [Lia, Loca] = ismember(middle_node2,new_nodes,'rows');
        faces_after_split(1:3,find(faces_after_split(1,:) == old_faces_list(triangles_to_split(ii),1) & faces_after_split(2,:) == old_faces_list(triangles_to_split(ii),2) & faces_after_split(3,:) == old_faces_list(triangles_to_split(ii),3))) = [ old_faces_list(triangles_to_split(ii),1); old_faces_list(triangles_to_split(ii),2); Loca ];
        faces_after_split(1:3,num_faces+1) = [  Loca; old_faces_list(triangles_to_split(ii),2); old_faces_list(triangles_to_split(ii),3) ];
        num_faces = num_faces + 1;
         
     end
     
     if ismember(middle_node3,new_nodes,'rows')
         
        [Lia, Loca] = ismember(middle_node3,new_nodes,'rows');
        faces_after_split(1:3,find(faces_after_split(1,:) == old_faces_list(triangles_to_split(ii),1) & faces_after_split(2,:) == old_faces_list(triangles_to_split(ii),2) & faces_after_split(3,:) == old_faces_list(triangles_to_split(ii),3))) = [ Loca; old_faces_list(triangles_to_split(ii),1); old_faces_list(triangles_to_split(ii),2)  ];
        faces_after_split(1:3,num_faces+1) = [ old_faces_list(triangles_to_split(ii),3); old_faces_list(triangles_to_split(ii),1); Loca  ];
        num_faces = num_faces + 1;
         
     end
          
          

end