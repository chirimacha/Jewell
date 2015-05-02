#!/bin/bash
#===============================================================================
#
#          FILE:  check_points_in_poly.sh
# 
#         USAGE:  ./check_points_in_poly.sh pointsfile.kmz blocksfile.kmz
# 
#   DESCRIPTION:  get kml of points and kml of polygons and return
#		list of points not in polygons
# 		use ogr2ogr (sudo apt-get install gdal-bin)
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  03/15/2012 05:36:02 PM EDT
#      REVISION:  ---
#===============================================================================
(
echo "5" ; 
poriginal=$(zenity --file-selection --title="Por favor selectiona kmz con casas")
pbasename=$(basename $poriginal .kmz)
echo "15" ; 
boriginal=$(zenity --file-selection --title="Por favor selectiona kmz con manzanas")
bbasename=$(basename $boriginal .kmz)
echo "25" ; 

echo casas en: $poriginal
echo manzanas en: $boriginal

# kmz to kml
unzip "$poriginal"
mv doc.kml $pbasename.kml
echo "30" ; 
unzip "$boriginal"
mv doc.kml "$bbasename.kml"
echo "35" ; 

echo "$pbasename.kml and $bbasename OK"
# kml to gmt for polygons
# ogr2ogr -f "GMT" $bbasename.gmt $bbasename.kml # for some reason import only a small fraction of the polygons
# kml to pseudo gmt for polygons
# grep "coordinates" "$bbasename.kml" | sed "s/^.*<coordinates>/>/" | sed "s#</coordinates># #" |sed "s/ \+/ /g" |tr ' ' '\n' |tr ',' ' ' > "$bbasename.pseudogmt"
awk '/coordinates/,/\/coordinates/' "$bbasename.kml" | sed "s/^.*<coordinates>/>/" | sed "s/<\/coordinates>//" |sed "s/[ \t]\+/ /g" | sed 's/^ //'|tr ' ' '\n' |tr ',' ' ' > "$bbasename.pseudogmt"
echo "$bbasename.pseudogmt OK"
echo "45" ; 

# kml to gmt for points
unified_kml.sh "$pbasename.kml"
echo "unified $pbasename.kml OK"
grep "[0-9]\+\.[0-9]\+\.[0-9]\+[A-Z]*\.[0-9]\|<coordinates>" "$pbasename-2.kml" | sed ':a;N;$!ba;s/<\/name>\n/ /g' |sed 's/<coordinates>/,/'| sed 's/<[^>]\+>//g' | sed 's/[\t ]*//g' | awk -F ',' '{print $3","$2","$4","$1}' > "$pbasename.csv"
echo "$pbasename.csv OK"
echo "65" ; 

## doesn't work well for some reason
# ogr2ogr -f "GMT" $pbasename.gmt $pbasename-2.kml
# # gmt to csv for points
# gmt2csv.sh $pbasename.gmt
echo -e "\nLaunch points_to_polygons"
points_to_polygons "$pbasename.csv" "$bbasename.pseudogmt" "$pbasename-blocks.csv" 2>points_to_polygons.err 1>points_to_polygons.out
echo "100" ; 
)|
zenity --progress \
	--title="Verificando correspondencia manzaneo/casas" \
	--text="Verificando..." \
	--percentage=0

error=$(grep ERROR points_to_polygons.out)
echo "$error"

if [ "$?" = -1 ] ; then
	zenity --error \
		--text="Verificacion anulada."
else
	if [[ $error != "" ]] ; then
		csvtokml.py points_to_polygons_errors.csv
		zenity --warning --text "$error. \nPara ver los problemas: abre el points_to_polygons_errors.kml en Googleearth."
	else
		zenity --info --text "Todo bien."
	fi
fi

# check houses out of manzanas
# grep " -1$" out.txt

