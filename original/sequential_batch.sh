export TIMEFORMAT=%R
INPUT_SIZES=(100 1000 2000 4000 5000 8000 10000 12000 15000 20000 25000 30000)
TIMESTEPS=(10)
ITERATIONS=30

touch temp_time.txt

mkdir -p outputs
mkdir -p data

printf "particle_count,timesteps,sample,time\n" > data/times_results.txt

for INPUT in ${INPUT_SIZES[*]}; do
    for TIMESTEP in ${TIMESTEPS[*]}; do
        for ITERATION in `seq $ITERATIONS`; do
            echo "Running ($INPUT, $TIMESTEP) iteration $ITERATION"
            touch "outputs/${INPUT}_${TIMESTEP}.txt"

            printf "${INPUT}\n${TIMESTEP}\n" > sequential_batch_input.in

            if [ $ITERATION -eq 1 ]; then
                RUN_TIME=$({ time ./nbody < sequential_batch_input.in >"outputs/${INPUT}_${TIMESTEP}.txt"; } 2>&1)
            else
                RUN_TIME=$({ time ./nbody < sequential_batch_input.in >/dev/null; } 2>&1)
            fi

            printf "${INPUT},${TIMESTEP},${ITERATION},${RUN_TIME}\n" >> data/times_results.txt
        done
    done
done

rm temp_time.txt
rm sequential_batch_input.in
