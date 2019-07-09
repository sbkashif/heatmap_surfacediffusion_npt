set encoding iso_8859_1
set terminal postscript solid color enhanced "Helvetica, 30"
set output "|ps2pdf - ./tmp.pdf"

name="plot.pdf"
margins="5"


unset key
#set key at graph -0.97,0.97 left samplen 1
set xrange[0:4.7]
set yrange[0:4.7]
p "cog_list_allComp.xvg" u 2:3, "cog_list_allComp.xvg" u 4:5 




set size ratio -1
set xrange[0:5]
set yrange[0:5]
set cbrange[0:1]
plot "freq_list_MBL.xvg" u 1:2:4 w p palette pt 5 ps 4.2

set xrange[0:5]
set yrange[0:5]
set cbrange[0:1]
plot "freq_list_pa.xvg" u 1:2:4 w p palette pt 5 ps 4.2

unset output
!pdfcrop ./tmp.pdf ./@name --margins @margins
!rm -f ./tmp.pdf


