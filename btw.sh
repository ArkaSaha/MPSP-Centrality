
DIR="data/Real/Autism_int/"

files="${DIR}*.csv.txt"

cores=2
C=1
for i in ${files}; do
    filename=${i##*/};
    outputfile="${DIR}BTW/${filename}.btw";
    logfile="${DIR}log/${filename}.log";
    echo ${i};
    ./mpsp ${i} ${outputfile} > ${logfile} &
    ((C++==0)) && wait;
    ((C=C%cores));
done
