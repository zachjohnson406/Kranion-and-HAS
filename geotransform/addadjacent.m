function tmsk = addadjacent(ttt)
%Input 
%ttt is a 3D logical or integer variable
%tmsk = the the output of adding all the neighboring points
tmsk = ttt + circshift(ttt,1,1)+ circshift(ttt,-1,1)+ circshift(ttt,1,2) ...
    + circshift(ttt,-1,2)+ circshift(ttt,1,3)+ circshift(ttt,-1,3);    %each point will have many neighbors
end