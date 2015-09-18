# Used to create a visual representation of codon frequency
# Author: Jairo Navarro
# Date: July 21, 2015

# Create a PDF of the graph
set term jpeg
set output "Beta_glu.jpeg"
set boxwidth .3
set style fill solid


#Label axis and title
set title 'Haloarcula Hispanica β-glucosidase'
set ylabel "Codon Frequency"
set xlabel "Amino Acid postition"

set key top right

set yrange[.05:6]

set mxtics 10

set label 1 sprintf("Loop Region") at 190,5.5 right textcolor rgb 'gold'
set label 2 sprintf("β-Sheet") at 190,5.2 right textcolor rgb 'blue'
set label 3 sprintf("α-Helix") at 190,4.9 right textcolor rgb 'red'

#set label sprintf("W") at 185,5
#set label sprintf("W") at 199, 5

#Coil: 1 
#Helix: 2
#Sheet: 4
set palette model RGB defined ( 1 'gold', 2 'red', 4 'blue' )
unset colorbox

plot [1:200] 'GNU_Hispanica.txt' using 1:3:4 \
     title "" with boxes palette

unset label
set label 1 sprintf("Loop Region") at 387,5.5 right textcolor rgb 'gold'
set label 2 sprintf("β-Sheet") at 387,5.2 right textcolor rgb 'blue'
set label 3 sprintf("α-Helix") at 387,4.9 right textcolor rgb 'red'

#set label sprintf("W") at 295,5
#set label sprintf("W") at 315, 5

plot [200:400] 'GNU_Hispanica.txt' using 1:3:4\
     title "" w boxes palette
unset label 
set label 1 sprintf("Loop Region") at 590,5.5 right textcolor rgb 'gold'
set label 2 sprintf("β-Sheet") at 590,5.2 right textcolor rgb 'blue'
set label 3 sprintf("α-Helix") at 590,4.9 right textcolor rgb 'red'

#set label sprintf("W") at 465,5
#set label sprintf("W") at 469, 5

plot [400:600] 'GNU_Hispanica.txt' using 1:3:4\
     title "" w boxes palette
unset label
set label 1 sprintf("Loop Region") at 790,5.5 right textcolor rgb 'gold'
set label 2 sprintf("β-Sheet") at 790,5.2 right textcolor rgb 'blue'
set label 3 sprintf("α-Helix") at 790,4.9 right textcolor rgb 'red'

#set label sprintf("W") at 507,5
#set label sprintf("W") at 508, 5

plot [600:800] 'GNU_Hispanica.txt' using 1:3:4\
     title "" w boxes palette
unset label 
set label 1 sprintf("Loop Region") at 990,5.5 right textcolor rgb 'gold'
set label 2 sprintf("β-Sheet") at 990,5.2 right textcolor rgb 'blue'
set label 3 sprintf("α-Helix") at 990,4.9 right textcolor rgb 'red'

#set label sprintf("W") at 641,5
#set label sprintf("W") at 815,5

plot [800:1000] 'GNU_Hispanica.txt' using 1:3:4\
     title "" w boxes palette

