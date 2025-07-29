set terminal png size 3000,1000
set output ofilename
set boxwidth 0.5
set style fill solid
set xtics rotate by 90 right noenhanced
plot filename using 2:xtic(1) with boxes
