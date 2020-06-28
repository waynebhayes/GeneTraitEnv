for seed in {1..5}; do
	for K in 2 4 5 10; do
		for traits in 1 3 5 10; do
		dir_name=Test-1/Seed_${seed}Traits_${traits}K_${K}
	        module load R/3.5.1;mkdir -p $dir_name; cd $dir_name; ls; Rscript ../../Run_shell_synthesis.R "../../Data_Matrix.csv" $seed $traits $K; cd ../..
		done
	done
done 
