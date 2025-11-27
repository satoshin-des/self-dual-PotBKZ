for dim in 100
do
    for seed in `seq 0 9`
    do
        printf "$dim $seed 0" | python3 Main.py;
    done
    printf "$dim 0 1" | python3 Main.py;
done
