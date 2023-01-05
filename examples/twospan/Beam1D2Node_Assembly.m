% ---------------------------------------------------------------------
% K代表整体的刚度矩阵；
% k代表各个梁单元的刚度矩阵；
% i，j代表两个相连的梁元；
% 该函数的功能是实现两个相连的梁元之间的刚度矩阵的组合
% ---------------------------------------------------------------------
function K=Beam1D2Node_Assembly(K,k,i,j)
DOF(1)=2*i-1;
DOF(2)=2*i;
DOF(3)=2*j-1;
DOF(4)=2*j;
for n1=1:4;
    for n2=1:4;
        K(DOF(n1),DOF(n2))=K(DOF(n1),DOF(n2))+k(n1,n2);
    end
end
K=K;