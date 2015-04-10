# save filenames to a array in linux system 
i=0
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -ls)

echo ${array[1]}
# write array into txt file 
for j in "${array[@]}"
do
      echo $j 
done >tmp.txt