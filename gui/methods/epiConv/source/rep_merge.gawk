BEGIN{ORS="";row=2;col=1;print "\n";if(infv=="") infv=0.0001}
{
	sum=0
	for(i=1;i<=NF;i++){
		if($i==0){
			$i=infv
		}
		sum+=log($i)/log(10)
	}
	print sum/NF
	col++
	if(col<row){
		print ","
	}else{
		print "\n"
		col=1
		row++
	}
}
