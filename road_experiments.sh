#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "You must enter the city name (Beijing, Porto or San_Francisco)";
    exit 1;
fi

city=$1

#./experiments data/Real/Road/${city}.graph data/Real/Road/${city}.queries output/Road/${city}_WISE_1.output 1 &
#./experiments data/Real/Road/${city}.graph data/Real/Road/${city}.queries output/Road/${city}_WISE_60.output 60 &
./experiments data/Real/Road/${city}.graph data/Real/Road/${city}.queries output/Road/${city}_WISE_10.output 10 &
#./mpsp data/Real/Road/${city}.graph data/Real/Road/${city}.queries output/Road/${city}_MPSP_1.output &

