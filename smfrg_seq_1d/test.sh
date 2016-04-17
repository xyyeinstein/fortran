#!/bin/bash
#for i in 1 2 5; do
for((i=1;i<=26;i++)); do
	U=$(awk -v U=$i 'BEGIN{printf "%.2f", (U-1)*0.20}')
	for((j=1;j<=21;j++)); do
		V=$(awk -v V=$j -v V0=$U 'BEGIN{printf "%.2f", V0*0.5+(V-11)*0.01}')
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
#if 判定语句
#-e  文件存在
#-f  被测文件是一个regular文件（正常文件，非目录或设备）
#-s  文件长度不为0
#-d  被测对象是目录
#-b  被测对象是块设备
#-c  被测对象是字符设备
#-p  被测对象是管道
#-h  被测文件是符号连接
#-L  被测文件是符号连接
#-S(大写)  被测文件是一个socket
#-t  关联到一个终端设备的文件描述符。用来检测脚本的stdin[-t0]或[-t1]是一个终端
#-r                          文件具有读权限，针对运行脚本的用户
#-w                         文件具有写权限，针对运行脚本的用户
#-x                          文件具有执行权限，针对运行脚本的用户
#-u                          set-user-id(suid)标志到文件，即普通用户可以使用的root权限文件，通过chmod +s file实现
#-k                          设置粘贴位
#-O                         运行脚本的用户是文件的所有者
#-G                         文件的group-id和运行脚本的用户相同
#-N                         从文件最后被阅读到现在，是否被修改
#f1 -nt f2                文件f1是否比f2新
#f1 -ot f2                文件f1是否比f2旧
#f1 -ef f2                文件f1和f2是否硬连接到同一个文件

#二元比较操作符，比较变量或比较数字
#整数比较：
#-eq                       等于            if [ "$a" -eq "$b" ]
#-ne                       不等于         if [ "$a" -ne "$b" ]
#-gt                        大于            if [ "$a" -gt "$b" ]
#-ge                       大于等于      if [ "$a" -ge "$b" ]
#-lt                         小于            if [ "$a" -lt "$b" ]
#-le                        小于等于      if [ "$a" -le "$b" ]
#<                          小于（需要双括号）       (( "$a" < "$b" ))
#<=                        小于等于(...)                (( "$a" <= "$b" ))
#>                          大于(...)                      (( "$a" > "$b" ))
#>=                        大于等于(...)                (( "$a" >= "$b" ))

#字符串比较：
#=                          等于           if [ "$a" = "$b" ]
#==                        与=等价
#!=                         不等于        if [ "$a" = "$b" ]
#<                          小于，在ASCII字母中的顺序：
#                            if [[ "$a" < "$b" ]]
#                            if [ "$a" \< "$b" ]         #需要对<进行转义
#						>                          大于
#-z                         字符串为null，即长度为0
#-n                         字符串不为null，即长度不为0
