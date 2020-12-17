methratiobwfile=$1
prefix=$2
regionlen=$3
maxvalue=$4 ##0.2
bed6=$5

computeMatrix scale-regions -S ${methratiobwfile} -R ${bed6} --beforeRegionStartLength ${regionlen} --regionBodyLength ${regionlen} --afterRegionStartLength ${regionlen} -o ${prefix}.matrix.gz --skipZeros

plotProfile -m ${prefix}.matrix.gz -out ${prefix}.profile.pdf  --plotTitle "" --zMin 0

plotHeatmap -m ${prefix}.matrix.gz -out ${prefix}.heatmap.pdf --missingDataColor white  --whatToShow 'heatmap and colorbar' --zMin 0 --kmeans 1 --zMax ${maxvalue}
