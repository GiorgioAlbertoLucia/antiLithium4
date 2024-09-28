# inputs
# data: /data/galucia/lithium4/LHC24_pass1_skimmed/LHC24ag_pass1_skimmed_slice.root
# mc: /data/galucia/lithium4/MC/AO2D_injectedli4.root

LOGFILE="output.log"
# CONF="-b --configuration json://configuration.json"
# CONF="-b --configuration json://configuration_new.json"
CONF="-b --configuration json://configuration_newest.json"
# CONF="-b --configuration json://configuration_mc.json"
# CONF="-b --configuration json://configuration_mc_newest.json"
OUTPUT_DIR="OutputDirector.json"
# OUTPUT_DIR="OutputDirector_mc.json"

echo "********************"
echo ${CONF}
echo "********************"


o2-analysis-lf-lithium4analysis ${CONF}|\

# converters
#o2-analysis-tracks-extra-converter ${CONF}|
#o2-analysis-bc-converter ${CONF}|
#o2-analysis-centrality-table ${CONF}|
#o2-analysis-mccollision-converter ${CONF}|

# standard wagons
o2-analysis-timestamp ${CONF}|\
o2-analysis-multiplicity-table ${CONF}|\
o2-analysis-event-selection ${CONF}|\
o2-analysis-track-propagation ${CONF}|\
o2-analysis-trackselection ${CONF}|\
o2-analysis-pid-tof-full ${CONF}|\
o2-analysis-pid-tof-base ${CONF}|\
o2-analysis-pid-tpc-base ${CONF}|\
o2-analysis-ft0-corrected-table ${CONF}|\
o2-analysis-pid-tpc ${CONF} --aod-writer-json $OUTPUT_DIR --aod-file @input_data.txt > $LOGFILE


# report the status of the workflow
rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
    exit $rc
fi

    