#!/usr/bin/awk -f

BEGIN{
	OFS = "\t"
}

{
	# set second chromosome
	
	if($7 == "="){
		$7 = $3
	}

	# get positions
	
	ci = $3
	cj = $7
	pi = int($4 / w) * w
	pj = int($8 / w) * w
	bi = ci OFS pi
	bj = cj OFS pj

	# print if new bin
	
	if(bi != bi_old){
		for(k in a){
			print bi_old, k, a[k]
		}
		split("", a) # release array
	}
		
	a[bj] += 1
	bi_old = bi
}

END{
	for(k in a){
		print bi_old, k, a[k]
	}
}
 
