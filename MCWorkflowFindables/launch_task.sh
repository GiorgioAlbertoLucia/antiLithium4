LOGFILE="output.log"
CONF="-b --configuration json://configuration_mc.json"


o2-analysis-lf-lithium4findables $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|

    # converters
    # o2-analysis-tracks-extra-converter $CONF|
    # o2-analysis-mc-converter $CONF|
    # o2-analysis-bc-converter $CONF|

    # standard wagons
    o2-analysis-timestamp $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-event-selection $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|      
    o2-analysis-track-propagation $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 --aod-writer-json OutputDirector_mc.json > $LOGFILE
    # o2-analysis-trackselection $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    # o2-analysis-multiplicity-table $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    # o2-analysis-pid-tpc-base $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    # o2-analysis-pid-tpc-full $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    # o2-analysis-pid-tof-base $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    # o2-analysis-pid-tof-full $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 


# report the status of the workflow
rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
    exit $rc
fi

    