#!/bin/bash
#for i in 1 2 5; do
cd /cygdrive/f/workspace
rm *.dat
for((i=1;i<=1;i++)); do
U=$(awk -v U=$i 'BEGIN{printf "%.2f", 2.0+(U-1)*0.10}')
    for((j=1;j<=1;j++)); do
#V=$(awk -v V=$j -v V0=$U 'BEGIN{printf "%.2f", V0*0.5+(V-11)*0.01}')
        V=$(awk -v V=$j 'BEGIN{printf "%.2f", 0.0+(V-1)*0.10}')
		for((k=1;k<=101;k++)); do
			mu=$(awk -v mu=$k 'BEGIN{printf "%.3f", (mu-1)*0.001}')
            dir=u${U}v${V}m${mu}
            if [ ! -d "$dir" ]; then
				mkdir ${dir}
	        fi
            cd ${dir}
	        if [ ! -e "flow.dat" ]; then
		        cp ../../research/fortran/smfrg_seq_1d/smfrg_seq.out . 
	            echo ${U} > U.input
        		echo ${V} > Vnn.input
        		echo ${mu} > mu.input
        		./smfrg_seq.out > output
			fi

			if [ -e "flow.dat" ]; then
	    		sed -n '1,1p' cdwform.dat >> ../f0.dat
				tail -1 flow.dat >> ../la.dat
    		fi
	
			cd ..
		done
    done
done
