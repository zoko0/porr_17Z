// function to multiply matrixes
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

// two dimensioned matrix
var m1:[{1..2, 1..2}] int;
m1(1,1) = 1; m1(1,2) = 2;
m1(2,1) = 3; m1(2,2) = 4;
// end of matrix

writeln("m1 ", m1);

var m2:[{1..2, 1..2}] int;
m2(1,1) = 2; m2(1,2) = 3;
m2(2,1) = 4; m2(2,2) = 5;
writeln("m2 ",m2);

var m3 = multiplyMatixes(m1, m2);
writeln("m3 ",m3);  
