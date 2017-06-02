#!/bin/bash
fun_names=(*.R)

# first function has not been added yet
first=true

# loop through all functions' file names
for f in ${fun_names[@]}
	do
	# store output of a grep on each of those filenames
	# to determine which of them do not belong (start witn "inner")
	temp=$(printf '%s' "$f" | grep ^internal)
	
	# if said output is of length 0, append function name to file
	if [ ${#temp} \< 1 ]
	then
		# if it's the first function, overwrite the file
		# with a new string and first function name
		if $first 
		then
			printf 'export("%s"\n' ${f::-2} > '../NAMESPACE'
			first=false
		else
			# else simply append to existing file
			printf '\t\t,"%s"\n' ${f::-2} >> '../NAMESPACE'
		fi
	fi
done

# append a closing bracket
printf '\t\t)\n' >> '../NAMESPACE'

