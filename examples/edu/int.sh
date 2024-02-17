
touch $1.stp
echo "Stop request sent to $1. Please wait..."
sleep 3
while [ -e $1.loc ] ; do sleep 2 ; done
rm $1.stp
