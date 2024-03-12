inFilePath="/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowLauncher/AO2D_lit_mc.root"
nKeysFilePath="/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/output/nKeys.txt"

root -l <<-EOF
const char * inFilePath = "$inFilePath";
TFile * inFile = TFile::Open(inFilePath);
const int nKeys = inFile->GetNkeys();
inFile->Close();
std::ofstream outputFile("$nKeysFilePath");
outputFile << nKeys << endl;
outputFile.close();
.q
EOF

nKeys=$(cat $nKeysFilePath)-1
rm $nKeysFilePath
echo "nKeys = $nKeys"

outFilePath="/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/output/efficiency_ptresolution.root"
root -l <<-EOF
const char * outFilePath = "$outFilePath";
TFile * outFile = TFile::Open(outFilePath, "RECREATE");
outFile->Close();
EOF

start_time=$(date +%s)

# Run the analysis 
MAX_CURRENT_PROCESSES=8

function check_running_processes {
    local count=$(jobs -p | wc -l)
    echo $count
}


for ((i=0; i<$nKeys; i++))
do
    while [ $(check_running_processes) -ge $MAX_CURRENT_PROCESSES ]; do
        sleep 1
    done

    inputVar=$i
    root -l -b -q "/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/src/efficiency_ptresolution.cpp($inputVar)" &
done

wait
echo "All done!"

end_time=$(date +%s)
runtime=$((end_time - start_time))

hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))

printf "Runtime: %02d:%02d:%02d\n" $hours $minutes $seconds