# restart; read "water_cad.mpl";

read "water_def.mpl";

SDetD := [DetD,diff(DetD,y1),diff(DetD,y2),
          diff(DetD,z1),diff(DetD,z2)]:

Pols := [g2,
         G2,
         g2-2*G2,
         G2-1,
         3*G2-1-2*g2,
         3*G2^2-5*G2*g2+g2^2+2*G2-2*g2+1,
         2*G2^2-5*G2*g2+2*g2^2-2*G2+3*g2,
         (2*G2^3-3*G2^2*g2-3*G2*g2^2+2*g2^3
          +2*G2^2+9*G2*g2-11*g2^2-4*G2+6*g2),
         -2*g2+G2+1,
         2*G2-1-g2,
         G2-2+g2,
         2*G2^2-5*G2*g2+2*g2^2+1];

xi := product(Pols[i], i = 1 .. nops(Pols));

with(RegularChains); with(ChainTools); with(SemiAlgebraicSetTools);
R := PolynomialRing([G2,g2]);
cadfull := CylindricalAlgebraicDecompose(Pols,R,output=cadcell):
sols := [];
for j from 1 to nops(cadfull) do
    sols := [op(sols),
             subs(op(2,op(1,subs(op(2,op(1,cadfull[j])),SamplePoint))),
                  box_bwe)[1]
            ]:
od:

print(nops(sols));
print(sols[1]);

avg2 := proc(l):
    return (l[1]+l[2])/2;
end:



insideBB := x -> (evala(eval(y1^2+(z1+1)^2-1, x)) <= 0 and
                  evala(eval(y2^2+(z2+1)^2-1, x)) <= 0):

PolsRed := remove(x -> x = g2 or x = G2 or x = 2*G2-g2, Pols):

sAvg := map(x -> [g2 = avg2(eval(g2,x)),
                        G2 = avg2(eval(G2,x))],
                  sols):

sAvgValid := select(x ->(eval(xi,x) <> 0
                                  and 0 < eval(g2,x)
                                  and 0 < eval(G2,x)
                                  and eval(g2,x) < 2*eval(G2,x)),
                             sAvg):
print(nops(%)):
# 570




Sols := [seq(select(x -> x <> {y1 = 0, y2 = 0,
                               z1 = -1, z2 = -1},
                    [solve(Groebner[Basis](eval(SDetD,
                                                sAvgValid[i]),
                                           plex(y1, y2, z1, z2)))]),
             i = 1..nops(sAvgValid))]:
print(map(nops, Sols));
# [2..2]

SolsInBB := [seq(select(x -> (evalb(0 <= evala(eval(y1^2, x)))
                              and evalb(0 <= evala(eval(y2^2, x)))
                              and insideBB(x)), Sols[i]),
                 i = 1 .. nops(Sols))]:

print(nops(select(x -> x <> [], SolsInBB)));
# 187

SolsBBSym := [seq(select(x -> (evalb(evala(eval(y1, x)) = 0)
                               and evalb(evala(eval(y2, x)) = 0)
                               and insideBB(x)), Sols[i]),
                  i = 1 .. nops(Sols))]:

print(nops(select(x -> ( x <> []), SolsBBSym)));
# 156

IndBBSym := select(x -> ( SolsBBSym[x] <> []),
                   [seq(i, i = 1 .. nops(SolsBBSym))]):

SolsBBNoSym := [seq(select(x -> (evalb(0 < evala(eval(y1^2, x)))
                                 and evalb(0 < evala(eval(y2^2, x)))
                                 and insideBB(x)), Sols[i]),
                    i = 1 .. nops(Sols))]:

print(nops(select(x -> ( x <> []), SolsBBNoSym)));
# 31

IndBBNoSym := select(x -> ( SolsBBNoSym[x] <> []),
                     [seq(i, i = 1 .. nops(SolsBBNoSym))]):

g2G2_1 := sAvgValid[IndBBSym]:
g2G2_2 := sAvgValid[IndBBNoSym]:

filename1 := "pts_cad_1.txt";
fid1 := fopen(filename1,WRITE);
for s in g2G2_1 do
    fprintf(fid1,"%a %a\n",evalf(eval(g2,s)), evalf(eval(G2,s)));
od:
fclose(fid1);

filename2 := "pts_cad_2.txt";
fid2 := fopen(filename2,WRITE);
for s in g2G2_2 do
    fprintf(fid2,"%a %a\n",evalf(eval(g2,s)), evalf(eval(G2,s)));
od:
fclose(fid2);

f0,f1,f2,f3,f4,f5,f6,f7,f8,f9 := op(Pols[3..-1]):

test_all := x -> map(f -> evalb(eval(f,x) > 0),
                     [f1,f2,f3,f4,f5,f6,f7,f8,f9]);




crit_1 := x -> ((eval(f1,x) > 0
                 and eval(f4,x) > 0
                 and eval(f2,x) < 0)
                or (eval(f1,x) < 0
                    and eval(f2,x) > 0));
print({op(map(crit_1,g2G2_1))});

print(nops(g2G2_1));
# 156
print(nops(select(crit_1,sAvgValid)));
# 156



crit_2 := x -> ((eval(f1,x) < 0
                 and eval(f6,x) > 0
                 and eval(f3,x) < 0)
                or (eval(f1,x) > 0
                    and eval(f6,x) < 0
                    and eval(f5,x) > 0)):
print({op(map(crit_2,g2G2_2))});

print(nops(g2G2_2));
# 31
print(nops(select(crit_2,sAvgValid)));
# 31
