# File "general_def.mpl";

d1:=g1-G1:d2:=g2-G2:

MatrixD:=Matrix([[-G1*y1, -z1-1,-G1+d1*z1,2*d1*y1],
                 [-g1*z1,y1,d1*y1,(G1-d1)-2*d1*z1],
                 [-G2*y2, -z2-1,-G2+d2*z2,2*d2*y2],
                 [-g2*z2,y2,d2*y2,(G2-d2)-2*d2*z2]]):

g1 := 1:
vars := [y1,y2,z1,z2]:
params := [G1,g2,G2]:

H := [1-y1^2-(z1+1)^2, 1-y2^2-(z2+1)^2]:

k := 4: n := 4: t := 3:
