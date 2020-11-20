BEGIN{ORS=""}
{
	if(NR!=1){
		print "\n"
	}
	for(i=1;i<ncol;i++){
		for(j=2;(j-1)*ncol<=NF;j++){
			$i+=$((j-1)*ncol+i)
		}
		print $i" "
	}
	for(j=2;(j-1)*ncol<=NF;j++){
		$ncol+=$((j-1)*ncol+ncol)
	}
	print $ncol
}

