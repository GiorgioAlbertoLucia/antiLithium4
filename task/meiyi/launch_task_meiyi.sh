LOGFILE="output.log"
CONF="-b --configuration json://configuration.json"
OUTPUT_DIR="OutputDirector.json"


o2-analysis-pid-tof-full $CONF | 
o2-analysis-pid-tof-base $CONF | 
o2-analysis-pid-tpc $CONF | 
o2-analysis-ft0-corrected-table $CONF | 
o2-analysis-timestamp $CONF | 
o2-analysis-tracks-extra-v002-converter $CONF | 
#o2-analysis-mccollision-converter $CONF |
o2-analysis-event-selection $CONF | 
o2-analysis-lf-he3hadronfemto $CONF | 
o2-analysis-pid-tpc-base $CONF | 
o2-analysis-multiplicity-table $CONF | 
o2-analysis-centrality-table $CONF | 
o2-analysis-track-propagation $CONF | 
o2-analysis-trackselection $CONF --aod-file @input_data.txt --aod-writer-json $OUTPUT_DIR > $LOGFILE


# report the status of the workflow
rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
    exit $rc
fi

    