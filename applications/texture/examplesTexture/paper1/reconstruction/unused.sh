users=`who | grep ":0 " | wc -l`
day=`date +%a`
hour=`date +%H`
echo $users
echo $day
echo $hour
[ $users -eq 0 ] || [ "$day" == "Sa" ] || [ "$day" == "So" ] || [ "$hour" \> "18" ] || [ "$hour" \< "08" ]
exit $?
