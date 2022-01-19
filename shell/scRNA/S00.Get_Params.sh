#!/bin/sh

#####Get parameters#####
ARGV=`getopt -o hmn:p:u: -l human,mouse,nMaxMultimap,help,partition,user -- "$@"`
#echo $ARGV
eval set -- "$ARGV"
while true
do
case "$1" in
    -h|--human)
        species="human"
        Nmax=1
        shift
        ;;
    -m|--mouse)
        species="mouse"
        Nmax=5
        shift
        ;;
    -n|--nMaxMultimap)
        Nmax=$2
        shift 2
        ;;
    -p|--partition)
        qos=$2
        if [ $qos == "cnl" ];then
                part="cn-long"
        elif [ $qos == "cns" ];then
                part="cn-short"
        elif [ $qos == "cnnl" ];then
                part="cn_nl"
        else
                qos="cnl"
                part="cn-long"
                echo "Couldn't detect available PARTITION name, will run pipeline in cn-long"
        fi
        shift 2
        ;;
    -u|--user)
        User=$2
        shift 2
        ;;
    --help)
        cat ${RNA_PIPELINE}/Documentation
        exit 1
        ;;
    --)
        break
        ;;
    *)
        echo " Unkonwn parameters:"$1""
        ;;
esac
done

if [ ! $species ];then
        echo "Species is required !"
        exit 1
fi

if [ ! $User ];then
        echo "User name is reqiured !"
        exit 1
fi

echo "User = $User"
echo "Species = $species"
echo "Max multimap number = $Nmax"

if [ -s sample.xls ];then
        echo "Samples to be processed:"
        cat sample.xls
else
        echo "sample.xls is not found or it's empty!"
        exit 1
fi

if [ ! $qos ];then
        qos="cnl"
        part="cn-long"
fi

echo "Pipeline will be run in $part"
