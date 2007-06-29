# Parmeters:

omega_file_name=omega1
output_file_name="omega1.2"

# Arguments:

N=$1
N1=$2
N2=$3
solver_algo=$4
first_case=$5
test_cases=$6
let last_case="$first_case + $test_cases - 1"

name=todo.$output_file_name.N1_$N1.N2_$N2.N_$N.algo_$solver_algo.$first_case-$last_case

echo $test_cases > $name
echo $first_case >> $name
echo $N >> $name
echo $omega_file_name >> $name
echo h$N1 >> $name
echo r$N2 >> $name
echo $output_file_name >> $name
echo $solver_algo >> $name
echo 2 stop_$name >> $name
echo >> $name
