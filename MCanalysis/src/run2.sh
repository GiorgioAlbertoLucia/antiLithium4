start_time=$(date +%s)

# Run the analysis 
root -l <<-EOF 
.x /Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/src/efficiency_ptresolution.cpp
.q
EOF
echo "All done!"

end_time=$(date +%s)
runtime=$((end_time - start_time))

hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))

printf "Runtime: %02d:%02d:%02d\n" $hours $minutes $seconds