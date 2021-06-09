for i in [Bl]*[0-9]* ; do 
    echo $i
    ln -s $i `echo $i \
	      | sed \
	      -e 's/lnc0/NI/g' \
	      -e 's/lncM/mock/g' \
	      -e 's/lnc1/100pfu/g' \
	      -e 's/lnc2/010pfu/g' \
	      -e 's/B3/001pfu/g' \
	      -e 's/B4/000.1pfu/g'` ; 
done
 
