#!/bin/bash
#for i in 1 2 5; do
for((i=1;i<=51;i++)); do
	U=$(awk -v U=$i 'BEGIN{printf "%.2f", (U-1)*0.10}')
	for((j=1;j<=51;j++)); do
		V=$(awk -v V=$j 'BEGIN{printf "%.2f", (V-1)*0.10}')
		for((k=1;k<=1;k++)); do
			mu=$(awk -v mu=$k 'BEGIN{printf "%.3f", (mu-1)*0.001}')
            dir=u${U}v${V}m${mu}
            if [ ! -d "$dir" ]; then
				mkdir ${dir}
	        fi
            cd ${dir}
	        if [ ! -e "flow.dat" ]; then
		        cp ../smfrg_seq.out .
	            echo ${U} > U.input
        		echo ${V} > Vnn.input
        		echo ${mu} > mu.input
        		./smfrg_seq.out > output
			fi

			if [ -e "flow.dat" ]; then
				#echo ${mu} >> ../mu.dat
	    		#sed -n '1,1p' cdwform.dat >> ../f0.dat
	    		tail -1 flow.dat >> ../phd.dat
				#echo $mu
				#tail -1 flow.dat
    		fi
	
			cd ..
		done
    done
done
