
DIR="data/Real/Autism_int/"

files="${DIR}*.csv.txt"

for i in ${files}; do
    filename=${i##*/};
    outputfile="${DIR}BTW/${filename}";
    echo ${i};
    ./mpsp ${i} ${outputfile};
done
