% --------------------------------------------------------------------
% P���������ܶȣ�A�������ĺ�������L�������ĳ��ȣ�
% �ú����Ĺ��������ÿһ����Ԫ����������
% --------------------------------------------------------------------
function m=Beam1D2Node_Mass(P,A,L)
m=P*A*L/420*[156 22*L 54 -13*L;22*L 4*L*L 13*L -3*L*L;...
    54 13*L 156 -22*L;-13*L -3*L*L -22*L 4*L*L];
