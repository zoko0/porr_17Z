// overloads '*' operator for arrays
proc multiplyMatixes(a:[], b:[]) {

    if (a.eltType != b.eltType) then
        writeln("type mismatch: ", a.eltType, " ", b.eltType);

    var ad = a.domain.dims();
    var bd = b.domain.dims();
    var (arows, acols) = ad;
    var (brows, bcols) = bd;
    if (arows != bcols) then
        writeln("dimension mismatch: ", ad, " ", bd);

    var c:[{arows, bcols}] a.eltType = 0;

    for i in arows do
        for j in bcols do
            for k in acols do
                c(i,j) += a(i,k) * b(k,j);

    return c;
}


var m1:[{1..2, 1..2}] int;
m1(1,1) = 1; m1(1,2) = 2;
m1(2,1) = 3; m1(2,2) = 4;
writeln("m1 ", m1);

var m2:[{1..2, 1..2}] int;
m2(1,1) = 2; m2(1,2) = 3;
m2(2,1) = 4; m2(2,2) = 5;
writeln("m2 ",m2);

var m3 = multiplyMatixes(m1, m2);
writeln("m3 ",m3);

var m4:[{1..2, 1..3}] int;
m4(1, 1) = 1; m4(1, 2) = 2; m4(1, 3) = 3;
m4(2, 1) = 4; m4(2, 2) = 5; m4(2, 3) = 6;
writeln("m4 ",m4);

var m5:[{1..3, 1..2}] int;
m5(1, 1) = 6; m5(1, 2) = -1;
m5(2, 1) = 3; m5(2, 2) =  2;
m5(3, 1) = 0; m5(3, 2) = -3;
writeln("m5 ",m5);

writeln("m4*m5 ",multiplyMatixes(m4,m5));
